"""Explicit sequence, confidence and geometric descriptors; no clinical claims."""
import re
import numpy as np
from Bio.PDB import NeighborSearch
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import ProtParamData
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from binderflow_utils import chain_data
from run_colabfold import validate_prediction


def sequence_features(sequence):
    if not sequence or set(sequence) - set('ACDEFGHIKLMNPQRSTVWY'):
        raise ValueError('Features require a nonempty canonical protein sequence.')
    analysis = ProteinAnalysis(sequence)
    n = len(sequence)
    width = min(7, n)
    windows = [sum(ProtParamData.kd[a] for a in sequence[i:i + width]) / width
               for i in range(n - width + 1)]
    return dict(length=n, molecular_weight_da=analysis.molecular_weight(),
                isoelectric_point=analysis.isoelectric_point(), charge_ph7=analysis.charge_at_pH(7),
                abs_charge_per_residue_ph7=abs(analysis.charge_at_pH(7)) / n,
                gravy=analysis.gravy(), aromatic_fraction=analysis.aromaticity(),
                instability_index=analysis.instability_index(), cysteine_count=sequence.count('C'),
                gly_pro_fraction=sum(sequence.count(a) for a in 'GP') / n,
                oxidation_prone_mw_count=sum(sequence.count(a) for a in 'MW'),
                deamidation_ng_ns_count=len(re.findall(r'(?=N[GS])', sequence)),
                n_glycosylation_motif_count=len(re.findall(r'(?=N[^P][ST])', sequence)),
                aggregation_window_hydropathy_max=max(windows), aggregation_window_length=width,
                aggregation_method='hydropathy_and_exposed_hydrophobic_patch_proxies',
                solubility_method='gravy_charge_and_exposure_descriptors_only')


def fit(mobile, reference):
    """Least-squares rotation for row-vector coordinates, with no reflection."""
    mobile, reference = np.asarray(mobile), np.asarray(reference)
    if mobile.shape != reference.shape or len(mobile) < 3:
        raise ValueError('RMSD requires matched sets of at least three CA atoms.')
    mc, rc = mobile.mean(0), reference.mean(0)
    u, _, vt = np.linalg.svd((mobile - mc).T @ (reference - rc))
    if np.linalg.det(u @ vt) < 0:
        u[:, -1] *= -1
    return u @ vt, mc, rc


def rmsd(a, b):
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def residue_sasa(chain):
    copy = chain.copy()
    ShrakeRupley(probe_radius=1.4, n_points=100).compute(copy, level='R')
    return [float(r.sasa) for r in copy if r.id[0] == ' ']


def structure_features(pdb, scores_path, reference, candidate):
    scores = validate_prediction(pdb, scores_path, candidate)
    model, chains = chain_data(pdb)
    _, reference_chains = chain_data(reference)
    bseq, br = chains['A']
    tseq, tr = chains['B']
    rb = reference_chains[candidate['reference_binder_chain']][1]
    rt_seq, rt = reference_chains[candidate['reference_target_chain']]
    if len(rb) != len(br) or rt_seq != tseq:
        raise ValueError('Reference and prediction do not have matching chains.')
    bc, tc, rbc, rtc = [np.array([r['CA'].coord for r in residues], dtype=float)
                         for residues in (br, tr, rb, rt)]
    rotation, center, refcenter = fit(tc, rtc)
    pose_rmsd = rmsd((bc - center) @ rotation + refcenter, rbc)
    brot, bcenter, rbcenter = fit(bc, rbc)
    fold_rmsd = rmsd((bc - bcenter) @ brot + rbcenter, rbc)

    # Exclude hydrogens from contacts/clashes and SASA consistently.
    for residue in model.get_residues():
        for atom in list(residue):
            if atom.element in ('H', 'D'):
                residue.detach_child(atom.id)
    atoms = list(model.get_atoms())
    contacts, interface_binder = set(), set()
    clashes = 0
    for first, second in NeighborSearch(atoms).search_all(5.0, level='A'):
        r1, r2 = first.get_parent(), second.get_parent()
        c1, c2 = r1.get_parent().id, r2.get_parent().id
        if c1 == c2:
            continue
        binder_res, target_res = (r1, r2) if c1 == 'A' else (r2, r1)
        contacts.add((binder_res.id, target_res.id))
        interface_binder.add(binder_res.id)
        if float(np.linalg.norm(first.coord - second.coord)) < 2.0:
            clashes += 1
    ShrakeRupley(probe_radius=1.4, n_points=100).compute(model, level='R')
    bound_binder = [float(r.sasa) for r in br]
    bound_target = [float(r.sasa) for r in tr]
    free_binder, free_target = residue_sasa(model['A']), residue_sasa(model['B'])
    binder_area, target_area = sum(free_binder), sum(free_target)
    complex_area = sum(bound_binder) + sum(bound_target)
    buried = binder_area + target_area - complex_area
    hydrophobic = set('AVILMFWY')
    hydrophobic_area = sum(area for aa, area in zip(bseq, free_binder) if aa in hydrophobic)
    # Connected exposed hydrophobic residues, CA distance <= 8 A, isolated binder.
    unvisited = {i for i, (aa, area) in enumerate(zip(bseq, free_binder)) if aa in hydrophobic and area >= 20}
    largest_patch = 0.0
    while unvisited:
        todo = [unvisited.pop()]
        patch = 0.0
        while todo:
            i = todo.pop()
            patch += free_binder[i]
            neighbours = {j for j in unvisited if np.linalg.norm(bc[i] - bc[j]) <= 8}
            unvisited -= neighbours
            todo.extend(neighbours)
        largest_patch = max(largest_patch, patch)
    n = len(br)
    plddt = np.asarray(scores['plddt'])
    pae = np.asarray(scores['pae'])
    interface_indices = [i for i, r in enumerate(br) if r.id in interface_binder]
    return dict(binder_plddt=float(plddt[:n].mean()), target_plddt=float(plddt[n:].mean()),
                binder_pae=float(pae[:n, :n].mean()),
                pae_binder_rows_target_cols=float(pae[:n, n:].mean()),
                pae_target_rows_binder_cols=float(pae[n:, :n].mean()),
                interchain_pae=float((pae[:n, n:].mean() + pae[n:, :n].mean()) / 2),
                ptm=float(scores['ptm']), iptm=float(scores['iptm']),
                interface_binder_plddt=float(plddt[interface_indices].mean()) if interface_indices else None,
                binder_sasa_a2=binder_area, target_sasa_a2=target_area, complex_sasa_a2=complex_area,
                buried_sasa_total_a2=max(0.0, buried), interface_area_a2=max(0.0, buried / 2),
                binder_buried_sasa_a2=max(0.0, binder_area - sum(bound_binder)),
                exposed_hydrophobic_sasa_a2=hydrophobic_area,
                exposed_hydrophobic_sasa_fraction=hydrophobic_area / binder_area if binder_area else 0.0,
                largest_exposed_hydrophobic_patch_a2=largest_patch,
                interface_residue_pairs=len(contacts), interface_binder_residues=len(interface_binder),
                interchain_clashing_atom_pairs=clashes,
                binder_ca_radius_gyration_a=float(np.sqrt(np.mean(np.sum((bc - bc.mean(0)) ** 2, axis=1)))),
                binder_fold_ca_rmsd_a=fold_rmsd, binder_target_aligned_ca_rmsd_a=pose_rmsd,
                target_ca_rmsd_a=rmsd((tc - center) @ rotation + refcenter, rtc))
