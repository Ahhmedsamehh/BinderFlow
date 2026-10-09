"""Contract tests; fake inference is not a scientific validation of RFdiffusion."""
import csv
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('runner', ROOT / 'bin/run_rfdiffusion.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)
TARGET = ROOT / 'sample-data/ctla4.pdb'


class RunnerTests(unittest.TestCase):
    def test_real_sample_target_and_hotspots(self):
        self.assertEqual(runner.validate_target(TARGET, 'C3-23/0 70-100', 'C3,C10,C20,C15'),
                         ['C3', 'C10', 'C20', 'C15'])

    def test_invalid_contigs_and_missing_residues(self):
        for contig in ('C23-3/0 70-100', 'C3-23/0 100-70', 'C3-999/0 70-100', 'garbage'):
            with self.subTest(contig=contig), self.assertRaises(ValueError):
                runner.validate_target(TARGET, contig, 'C3')

    def test_hotspots_must_be_retained_and_unique(self):
        for hotspots in ('C24', 'D3', 'C3,C3', 'C3,,C10', 'C3;echo bad'):
            with self.subTest(hotspots=hotspots), self.assertRaises(ValueError):
                runner.validate_target(TARGET, 'C3-23/0 70-100', hotspots)

    def test_output_pairs_must_exist_and_be_nonempty(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / 'test_0.pdb').write_text('ATOM  test\n')
            with self.assertRaises(ValueError):
                runner.collect_outputs(folder, 'test', 0, 1)
            (folder / 'test_0.trb').touch()
            with self.assertRaisesRegex(ValueError, 'Empty'):
                runner.collect_outputs(folder, 'test', 0, 1)
            (folder / 'test_0.trb').write_bytes(b'test metadata')
            runner.collect_outputs(folder, 'test', 0, 1)
            with (folder / 'manifest.csv').open() as handle:
                manifest = list(csv.DictReader(handle))
            self.assertEqual(manifest[0]['pdb'], '../backbones/test_0.pdb')

    def run_wrapper(self, tmp, *extra):
        return subprocess.run([
            sys.executable, str(ROOT / 'bin/run_rfdiffusion.py'),
            '--sample-id', 'test', '--target', str(TARGET), '--contigs', 'C3-23/0 70-100',
            '--hotspots', 'C3,C10', '--models', str(Path(tmp) / 'model weights'),
            '--script', str(Path(tmp) / 'fake inference.py'), '--num-designs', '2',
            '--design-startnum', '7', *extra
        ], cwd=tmp, text=True, capture_output=True)

    def test_stub_count_start_index_and_markers(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = self.run_wrapper(tmp, '--stub')
            self.assertEqual(result.returncode, 0, result.stderr)
            out = Path(tmp) / 'outputs'
            self.assertEqual(sorted(p.name for p in out.glob('*.pdb')), ['test_7.pdb', 'test_8.pdb'])
            self.assertTrue(json.loads((out / 'run.json').read_text())['stub'])
            self.assertIn('MOCK', (out / 'test_7.pdb').read_text())

    def test_real_command_contract_with_fake_engine_and_spaced_paths(self):
        # Exercise the real subprocess branch, but replace the model engine.
        with tempfile.TemporaryDirectory(prefix='binderflow test ') as tmp:
            folder = Path(tmp)
            (folder / 'model weights').mkdir()
            (folder / 'model weights/Complex_base_ckpt.pt').write_bytes(b'FAKE CHECKPOINT')
            (folder / 'fake inference.py').write_text(
                'import json, sys\nfrom pathlib import Path\n'
                'args = dict(arg.split("=", 1) for arg in sys.argv[1:])\n'
                'prefix = json.loads(args["inference.output_prefix"])\n'
                'assert Path(json.loads(args["inference.ckpt_override_path"])).is_file()\n'
                'assert json.loads(args["contigmap.contigs"]) == ["C3-23/0 70-100"]\n'
                'start = int(args["inference.design_startnum"])\n'
                'for i in range(start, start + int(args["inference.num_designs"])):\n'
                '    Path(f"{prefix}_{i}.pdb").write_text("ATOM  FAKE STRUCTURE\\n")\n'
                '    Path(f"{prefix}_{i}.trb").write_text("FAKE METADATA\\n")\n'
            )
            result = self.run_wrapper(tmp)
            self.assertEqual(result.returncode, 0, result.stderr)
            report = json.loads((folder / 'outputs/run.json').read_text())
            self.assertFalse(report['stub'])
            self.assertEqual(len(report['checkpoint_sha256']), 64)
            (folder / 'fake inference.py').write_text('import sys\nprint("engine failed")\nsys.exit(3)\n')
            result = self.run_wrapper(tmp)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('status 3', result.stderr)

    def test_missing_model_and_invalid_parameters(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / 'fake inference.py').write_text('raise RuntimeError("must not run")')
            result = self.run_wrapper(tmp)
            self.assertIn('Missing or empty binder checkpoint', result.stderr)
            for args in (('--num-designs', '0'), ('--noise-scale-ca', 'nan'), ('--sample-id', '../bad')):
                result = self.run_wrapper(tmp, '--stub', *args)
                self.assertNotEqual(result.returncode, 0)


@unittest.skipUnless(os.environ.get('NEXTFLOW_BIN'), 'Set NEXTFLOW_BIN to run Nextflow integration tests.')
class WorkflowTests(unittest.TestCase):
    def test_staging_publication_resume_and_validation(self):
        with tempfile.TemporaryDirectory(prefix='binderflow integration ') as tmp:
            folder = Path(tmp)
            (folder / 'target with spaces.pdb').write_bytes(TARGET.read_bytes())
            sheet = folder / 'samples.csv'
            sheet.write_text('id,target_pdb,contigs,hotspots\n'
                             'one,target with spaces.pdb,C3-23/0 70-100,C3;C10\n'
                             'two,target with spaces.pdb,C3-23/0 70-100,C15;C20\n')
            command = [os.environ['NEXTFLOW_BIN'], 'run', str(ROOT / 'main.nf'),
                       '-profile', 'stub', '-stub-run', '--input', str(sheet),
                       '--stop_after', 'rfdiffusion',
                       '--num_designs', '2', '--design_startnum', '4',
                       '--outdir', str(folder / 'published results'), '-ansi-log', 'false']
            result = subprocess.run(command, cwd=folder, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            for sample in ('one', 'two'):
                output = folder / 'published results' / sample
                self.assertEqual(len(list((output / 'backbones').glob('*.pdb'))), 2)
                self.assertEqual(len(list((output / 'metadata').glob('*.trb'))), 2)
                self.assertTrue(json.loads((output / 'reports/run.json').read_text())['stub'])
                with (output / 'reports/manifest.csv').open() as handle:
                    for row in csv.DictReader(handle):
                        self.assertTrue((output / 'reports' / row['pdb']).is_file())
                        self.assertTrue((output / 'reports' / row['trb']).is_file())
            resumed = subprocess.run(command + ['-resume'], cwd=folder, capture_output=True, text=True)
            self.assertEqual(resumed.returncode, 0, resumed.stdout + resumed.stderr)
            self.assertIn('Cached process', resumed.stdout)
            single_command = command.copy()
            single_command[single_command.index('--num_designs') + 1] = '1'
            single_command[single_command.index('--outdir') + 1] = str(folder / 'single result')
            single = subprocess.run(single_command, cwd=folder, capture_output=True, text=True)
            self.assertEqual(single.returncode, 0, single.stdout + single.stderr)
            for sample in ('one', 'two'):
                self.assertEqual(len(list((folder / 'single result' / sample / 'backbones').glob('*.pdb'))), 1)
            sheet.write_text(sheet.read_text().replace('two,', 'one,'))
            invalid = subprocess.run(command, cwd=folder, capture_output=True, text=True)
            self.assertNotEqual(invalid.returncode, 0)
            self.assertIn('Duplicate sample ID', invalid.stdout + invalid.stderr)


if __name__ == '__main__':
    unittest.main()
