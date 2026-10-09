"""Scientific descriptor checks and end-to-end wiring tests (no model inference)."""
import copy
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
from Bio.PDB import PDBIO, PDBParser

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
from binderflow_utils import mock_pdb, write_json
from developability import fit, rmsd, sequence_features, structure_features
from run_mpnn import parse_sequences
from run_colabfold import validate_prediction
from rank_binders import rank


def fixture(folder):
    candidate = dict(candidate_id='test_0__s001', sample_id='test', backbone_id='test_0',
                     binder_sequence='AEKLAEKLAEKL', target_sequence='VTISCTGSSSNIGAGYDVH',
                     reference_binder_chain='A', reference_target_chain='B', stub=True)
    n = len(candidate['binder_sequence']) + len(candidate['target_sequence'])
    for name in ('reference.pdb', 'prediction.pdb'):
        mock_pdb(folder / name, candidate['binder_sequence'], candidate['target_sequence'])
    write_json(folder / 'scores.json', dict(plddt=[80] * n, pae=np.full((n, n), 4.0).tolist(), ptm=0.8, iptm=0.7))
    return candidate


class DescriptorTests(unittest.TestCase):
    def test_sequence_features_have_interpretable_scales(self):
        values = sequence_features('VVVVVVVV')
        self.assertAlmostEqual(values['gravy'], 4.2)
        self.assertAlmostEqual(values['aggregation_window_hydropathy_max'], 4.2)
        self.assertEqual(values['length'], 8)
        self.assertEqual(sequence_features('NATNGSCMW')['n_glycosylation_motif_count'], 2)
        with self.assertRaises(ValueError):
            sequence_features('AXZ')

    def test_alignment_removes_translation_and_rotation(self):
        ref = np.array([[0, 0, 0], [1, 0, 0], [0, 2, 0], [0, 0, 3]], dtype=float)
        rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        mobile = ref @ rotation + [7, 4, -8]
        r, m, c = fit(mobile, ref)
        self.assertLess(rmsd((mobile - m) @ r + c, ref), 1e-10)
        self.assertGreater(np.linalg.det(r), 0.99)

    def test_separated_chains_lose_interface_and_pose_agreement(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            candidate = fixture(folder)
            before = structure_features(folder / 'prediction.pdb', folder / 'scores.json', folder / 'reference.pdb', candidate)
            self.assertGreater(before['interface_residue_pairs'], 0)
            self.assertGreater(before['buried_sasa_total_a2'], 0)
            self.assertLess(before['binder_target_aligned_ca_rmsd_a'], 1e-5)
            model = PDBParser(QUIET=True).get_structure('x', folder / 'prediction.pdb')
            for atom in model[0]['B'].get_atoms():
                atom.coord += [100, 0, 0]
            writer = PDBIO()
            writer.set_structure(model)
            writer.save(str(folder / 'prediction.pdb'))
            after = structure_features(folder / 'prediction.pdb', folder / 'scores.json', folder / 'reference.pdb', candidate)
            self.assertEqual(after['interface_residue_pairs'], 0)
            self.assertIsNone(after['interface_binder_plddt'])
            self.assertAlmostEqual(after['buried_sasa_total_a2'], 0, places=3)
            self.assertGreater(after['binder_target_aligned_ca_rmsd_a'], 99)
            self.assertLess(after['binder_fold_ca_rmsd_a'], 1e-5)

    def test_confidence_dimensions_and_chain_sequences_are_checked(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            candidate = fixture(folder)
            validate_prediction(folder / 'prediction.pdb', folder / 'scores.json', candidate)
            wrong = dict(candidate, binder_sequence='A' * len(candidate['binder_sequence']))
            with self.assertRaisesRegex(ValueError, 'sequences'):
                validate_prediction(folder / 'prediction.pdb', folder / 'scores.json', wrong)
            scores = json.loads((folder / 'scores.json').read_text())
            scores['pae'] = [[1]]
            write_json(folder / 'scores.json', scores)
            with self.assertRaisesRegex(ValueError, 'pae'):
                validate_prediction(folder / 'prediction.pdb', folder / 'scores.json', candidate)

    def test_mpnn_fasta_skips_native_and_rejects_wrong_length(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'sequences.fa'
            path.write_text('>native, score=0.1\nGGGG\n>T=0.1, sample=1, score=0.5, global_score=1\nAEKL\n')
            self.assertEqual(parse_sequences(path, 1, 4), [('AEKL', 0.5)])
            with self.assertRaises(ValueError):
                parse_sequences(path, 2, 4)
            with self.assertRaises(ValueError):
                parse_sequences(path, 1, 5)

    def test_ranking_is_per_target_ties_are_equal_and_gates_come_first(self):
        config = {'gates': {'clashes': {'max': 0}}, 'metrics': {'quality': {'weight': 1, 'direction': 'higher'}}}
        rows = [dict(candidate_id=name, sample_id=sample, target_sequence=sample, stub=False,
                     quality=quality, clashes=clashes) for name, sample, quality, clashes in
                [('good', 'a', 5, 0), ('bad', 'a', 10, 1), ('tie', 'a', 5, 0), ('other', 'b', -100, 0)]]
        ranked = rank(rows, config)
        self.assertEqual([r['candidate_id'] for r in ranked], ['good', 'tie', 'bad', 'other'])
        self.assertEqual(ranked[0]['ranking_score'], ranked[1]['ranking_score'])
        self.assertEqual(ranked[-1]['rank_within_target'], 1)
        self.assertEqual(ranked[-1]['ranking_score'], 0.5)
        rows[0]['quality'] = None
        with self.assertRaisesRegex(ValueError, 'ranking metric'):
            rank(rows, config)
        rows[0]['quality'] = 5
        rows[0]['stub'] = True
        with self.assertRaisesRegex(ValueError, 'mix mock'):
            rank(rows, config)


class AdapterTests(unittest.TestCase):
    def test_nonstub_mpnn_and_colabfold_commands_with_fake_engines(self):
        with tempfile.TemporaryDirectory(prefix="binderflow's contract ") as tmp:
            folder = Path(tmp)
            candidate = fixture(folder)
            # The retained target is A, so the binder must be identified as B.
            original = (folder / 'reference.pdb').read_text()
            swapped = '\n'.join(line[:21] + ('B' if line[21] == 'A' else 'A') + line[22:]
                                if line.startswith('ATOM  ') else line for line in original.splitlines())
            (folder / 'reference.pdb').write_text(swapped + '\n')
            engine = folder / 'fake mpnn.py'
            engine.write_text('import sys\nfrom pathlib import Path\n'
                'a=dict(zip(sys.argv[1::2],sys.argv[2::2]))\n'
                'assert a["--pdb_path_chains"] == "B"\n'
                'out=Path(a["--out_folder"])/"seqs"\nout.mkdir(parents=True)\n'
                '(out/(Path(a["--pdb_path"]).stem+".fa")).write_text(">native\\nGGGG\\n>T=0.1, sample=1, score=0.5\\nAEKLAEKLAEKL\\n")\n')
            result = subprocess.run([sys.executable, str(ROOT / 'bin/run_mpnn.py'), '--backbone', 'reference.pdb',
                '--sample-id', 'sample', '--backbone-id', 'sample_0', '--contigs', 'A1-19/0 12',
                '--script', str(engine), '--count', '1'], cwd=folder, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            records = json.loads((folder / 'sequences/candidates.json').read_text())
            self.assertEqual(records[0]['target_sequence'], candidate['target_sequence'])
            self.assertEqual(records[0]['reference_binder_chain'], 'B')
            self.assertFalse(records[0]['stub'])
            fake_af = folder / 'fake colabfold'
            fake_af.write_text(f'#!{sys.executable}\nimport sys,json\nfrom pathlib import Path\n'
                f'sys.path.insert(0,{str(ROOT / "bin")!r})\nfrom binderflow_utils import mock_pdb,write_json\n'
                'a=dict(zip(sys.argv[3::2],sys.argv[4::2]))\nassert a["--msa-mode"] == "single_sequence"\n'
                'out=Path(sys.argv[2]);out.mkdir()\n'
                'for block in Path(sys.argv[1]).read_text().split(">")[1:]:\n'
                ' name,seq=block.strip().splitlines();b,t=seq.split(":");n=len(b)+len(t)\n'
                ' tag="rank_001_alphafold2_multimer_v3_model_1_seed_001"\n'
                ' mock_pdb(out/f"{name}_unrelaxed_{tag}.pdb",b,t)\n'
                ' write_json(out/f"{name}_scores_{tag}.json",dict(plddt=[80]*n,pae=[[5]*n for _ in range(n)],ptm=0.8,iptm=0.7))\n')
            fake_af.chmod(0o755)
            (folder / 'model weights').mkdir()
            result = subprocess.run([sys.executable, str(ROOT / 'bin/run_colabfold.py'), '--sequences', 'sequences',
                '--models', 'model weights', '--executable', str(fake_af)], cwd=folder, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((folder / 'predictions/sample_0__s001.pdb').is_file())
            self.assertFalse(json.loads((folder / 'predictions/candidates.json').read_text())[0]['stub'])


@unittest.skipUnless(os.environ.get('NEXTFLOW_BIN'), 'Set NEXTFLOW_BIN for full workflow integration.')
class FullWorkflowTests(unittest.TestCase):
    def test_multiple_targets_backbones_sequences_and_resume(self):
        with tempfile.TemporaryDirectory(prefix="binderflow's full workflow ") as tmp:
            folder = Path(tmp)
            (folder / 'target.pdb').write_bytes((ROOT / 'sample-data/ctla4.pdb').read_bytes())
            (folder / 'samples.csv').write_text('id,target_pdb,contigs,hotspots\n'
                    'one,target.pdb,C3-23/0 70-100,C3\n'
                    'two_long,target.pdb,C3-23/0 70-100,\n')
            command = [os.environ['NEXTFLOW_BIN'], 'run', str(ROOT / 'main.nf'), '-profile', 'stub', '-stub-run',
                       '--input', str(folder / 'samples.csv'), '--num_designs', '2', '--sequences_per_backbone', '2',
                       '--outdir', str(folder / 'published results'), '-ansi-log', 'false']
            env = dict(os.environ, PATH=str(Path(sys.executable).parent) + os.pathsep + os.environ['PATH'])
            run = subprocess.run(command, cwd=folder, env=env, text=True, capture_output=True)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            out = folder / 'published results'
            rows = json.loads((out / 'ranked_candidates.json').read_text())
            self.assertEqual(len(rows), 8)
            self.assertEqual(len({r['candidate_id'] for r in rows}), 8)
            for sample in ('one', 'two_long'):
                self.assertEqual(sorted(r['rank_within_target'] for r in rows if r['sample_id'] == sample), [1, 2, 3, 4])
            for row in rows:
                self.assertTrue(row['stub'])
                self.assertTrue((out / row['pdb']).is_file())
                self.assertTrue((out / row['scores']).is_file())
                self.assertFalse(any('immuno' in key.lower() for key in row))
            resumed = subprocess.run(command + ['-resume'], cwd=folder, env=env, text=True, capture_output=True)
            self.assertEqual(resumed.returncode, 0, resumed.stdout + resumed.stderr)
            self.assertEqual(resumed.stdout.count('Cached process'), 15)


if __name__ == '__main__':
    unittest.main()
