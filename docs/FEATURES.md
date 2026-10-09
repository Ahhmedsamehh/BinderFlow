# Feature definitions and ranking

Features describe the designed binder sequence and the rank-one predicted
binder–target complex. Confidence measures (pLDDT, pTM, ipTM and PAE) are model
outputs, not geometrical measurements or experimental binding measurements.
Immunogenicity is deliberately excluded.

## Confidence

- `binder_plddt`, `target_plddt`: mean residue confidence, 0–100, sliced using
  validated chain lengths and the input binder:target order.
- `binder_pae`: mean binder–binder PAE block, including its diagonal, in Å.
- `pae_binder_rows_target_cols` and its reverse: separate directional blocks.
  `interchain_pae` averages these two block means in Å; it is not restricted to
  contacting residues.
- `ptm`, `iptm`: complex pTM and interface pTM reported by ColabFold, 0–1.
- `interface_binder_plddt`: mean confidence over binder residues with an
  interchain heavy-atom contact within 5 Å. Missing if no interface exists.

PDB and score JSON must come from the same rank-one model. Mismatched sequences,
array sizes, nonfinite scores and out-of-range confidence cause a task failure.

## Geometry and solvent exposure

SASA uses [Biopython Shrake–Rupley](https://biopython.org/docs/latest/api/Bio.PDB.SASA.html),
a 1.4 Å water probe and 100 surface points per atom. Hydrogen atoms are excluded.
Areas are in Å² and are approximate; use more surface points for convergence
studies. Calculations use the unrelaxed predicted complex without solvent,
glycans or other binding partners.

| Fields | Definition |
| --- | --- |
| `binder_sasa_a2`, `target_sasa_a2` | Isolated chains in their predicted bound conformations |
| `complex_sasa_a2` | SASA of the two-chain complex |
| `buried_sasa_total_a2` | SASA(binder) + SASA(target) − SASA(complex) |
| `interface_area_a2` | Half the total buried SASA; avoids a factor-of-two ambiguity |
| `binder_buried_sasa_a2` | Isolated binder SASA minus binder SASA in the complex |
| `interface_residue_pairs` | Distinct interchain residue pairs with any heavy-atom separation ≤5 Å |
| `interface_binder_residues` | Number of binder residues in those contacts |
| `interchain_clashing_atom_pairs` | Interchain heavy-atom pairs separated by <2 Å; a crude screen, not a MolProbity clashscore |
| `binder_ca_radius_gyration_a` | Root mean squared CA distance from the binder CA centroid |
| `binder_fold_ca_rmsd_a` | Binder CA RMSD after independently aligning the binder to the designed backbone |
| `binder_target_aligned_ca_rmsd_a` | Binder CA RMSD after aligning only the target; captures changes in docking pose and binder conformation |
| `target_ca_rmsd_a` | Target CA RMSD after that target alignment |

RMSDs use ordered, matching chains and a least-squares rotation without
reflection. Target identity and binder length are checked against the reference.
These are whole-chain comparisons; no flexible-region trimming is applied.

## Sequence developability and aggregation/solubility proxies

[Biopython ProtParam](https://biopython.org/docs/latest/api/Bio.SeqUtils.ProtParam.html)
provides molecular weight (Da), estimated pI, estimated charge at pH 7, GRAVY,
aromatic fraction and the sequence instability index. Absolute charge divided
by binder length is also reported. These are context-free sequence descriptors;
the instability index is not a measured stability or melting temperature.

- `aggregation_window_hydropathy_max`: maximum mean Kyte–Doolittle hydropathy
  over overlapping seven-residue windows (the whole sequence if shorter).
  A high value flags a hydrophobic stretch; it is **not a validated aggregation
  probability, TANGO score or AGGRESCAN score**.
- `exposed_hydrophobic_sasa_a2` and its fraction: isolated-binder SASA assigned
  to A/V/I/L/M/F/W/Y residues. All atoms of these residues contribute, including
  backbone atoms. This is an operational descriptor, not a polar/nonpolar atom
  decomposition. The binding face counts because the free binder exposes it.
- `largest_exposed_hydrophobic_patch_a2`: largest connected set of those residues
  having isolated SASA ≥20 Å² and CA separations ≤8 Å. Patch area is their summed
  SASA. Thresholds are heuristic and do not define a validated aggregation model.
- **Solubility** is represented by GRAVY, charge density and exposed hydrophobic
  surface descriptors. No calibrated solubility prediction is emitted. Solubility
  and aggregation depend on pH, ionic strength, concentration and formulation,
  which are absent here. Dedicated predictors would require separate validation.
- Cysteine count, G/P fraction, M/W oxidation-prone residue count, N-G/N-S motif
  count and N-X-S/T motifs (X ≠ P) flag sequence liabilities. Counts do not prove
  oxidation, deamidation or glycosylation; the expression system and local
  structure matter.

`aggregation_method` and `solubility_method` in every row identify this scope.

## Default selection rules

The editable [ranking configuration](../config/ranking.json) first checks:
binder pLDDT ≥70, mean interchain PAE ≤15 Å, at least one interface residue pair,
at most ten interchain atom clashes, and target-aligned binder RMSD ≤5 Å.
These are initial engineering filters, not validated success thresholds.

Within each target, each ranking feature is converted to a percentile; lower-is-
better features are reversed, ties receive their average percentile, and a
single-candidate group receives 0.5. The weighted mean uses:

| Feature | Preferred direction | Weight |
| --- | --- | --- |
| ipTM | Higher | 0.25 |
| Interchain PAE | Lower | 0.20 |
| Binder pLDDT | Higher | 0.15 |
| Target-aligned binder RMSD | Lower | 0.15 |
| Binder buried SASA | Higher | 0.10 |
| Exposed hydrophobic SASA fraction | Lower | 0.10 |
| Maximum window hydropathy | Lower | 0.05 |

Filter-passing candidates precede failures; remaining ties sort by candidate ID.
Nothing is silently dropped. Missing required ranking metrics cause failure,
not an invented score. Correlated features can count related evidence twice,
and area measures depend on binder size. Weights are transparent defaults that
need calibration to your data. Do not compare scores across targets or runs:
percentiles change with the candidate set. No score is a probability of binding,
expression or therapeutic suitability.
