#!/usr/bin/env python3
"""Predict binder:target complexes and retain matched rank-one PDB/score pairs."""
import argparse
import json
from pathlib import Path
import shutil
import sys
from binderflow_utils import chain_data, execute, mock_pdb, sha256, versions, write_json


def validate_prediction(pdb, scores_path, candidate):
    import numpy as np
    _, chains = chain_data(pdb)
    if list(chains) != ['A', 'B']:
        raise ValueError('Expected ColabFold chains A (binder) and B (target), in that order.')
    if chains['A'][0] != candidate['binder_sequence'] or chains['B'][0] != candidate['target_sequence']:
        raise ValueError('Predicted chain sequences do not match the candidate.')
    scores = json.loads(Path(scores_path).read_text())
    n = len(chains['A'][0]) + len(chains['B'][0])
    for key, shape, low, high in [('plddt', (n,), 0, 100), ('pae', (n, n), 0, float('inf')),
                                   ('ptm', (), 0, 1), ('iptm', (), 0, 1)]:
        values = np.asarray(scores.get(key), dtype=float)
        if values.shape != shape or not np.isfinite(values).all() or (values < low).any() or (values > high).any():
            raise ValueError(f'Invalid or missing ColabFold {key}.')
    return scores


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sequences', type=Path, required=True)
    p.add_argument('--executable', default='colabfold_batch')
    p.add_argument('--models', type=Path, required=True)
    p.add_argument('--num-models', type=int, default=1)
    p.add_argument('--recycles', type=int, default=3)
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--stub', action='store_true')
    a = p.parse_args()
    if not 1 <= a.num_models <= 5 or a.recycles < 0 or a.seed < 0:
        raise ValueError('Invalid ColabFold models/recycles/seed settings.')
    candidates = json.loads((a.sequences / 'candidates.json').read_text())
    if not candidates or any(r['stub'] != a.stub for r in candidates):
        raise ValueError('Empty candidate set or mixed real/mock stages.')
    out = Path('predictions')
    out.mkdir()
    shutil.copyfile(a.sequences / 'backbone.pdb', out / 'backbone.pdb')
    command = []
    if not a.stub:
        if not a.models.is_dir():
            raise ValueError('ColabFold model directory missing.')
        command = [a.executable, str(a.sequences / 'complexes.fasta'), str(out / 'raw'),
                   '--model-type', 'alphafold2_multimer_v3', '--msa-mode', 'single_sequence',
                   '--num-models', str(a.num_models), '--num-recycle', str(a.recycles),
                   '--num-seeds', '1', '--random-seed', str(a.seed), '--rank', 'multimer',
                   '--data', str(a.models)]
        execute(command, out / 'inference.log')
    else:
        (out / 'inference.log').write_text('MOCK structures/confidence; ColabFold was not run.\n')
    for record in candidates:
        name = record['candidate_id']
        pdb, scores = out / f'{name}.pdb', out / f'{name}.scores.json'
        if a.stub:
            n = len(record['binder_sequence']) + len(record['target_sequence'])
            mock_pdb(pdb, record['binder_sequence'], record['target_sequence'])
            write_json(scores, dict(plddt=[85.0] * n,
                                    pae=[[0.0 if i == j else 5.0 for j in range(n)] for i in range(n)],
                                    ptm=0.75, iptm=0.7, stub=True))
        else:
            matches = list((out / 'raw').glob(f'{name}_unrelaxed_rank_001_*.pdb'))
            if len(matches) != 1:
                raise ValueError(f'Expected exactly one rank-one prediction for {name}.')
            original = matches[0]
            paired = original.with_name(original.name.replace('_unrelaxed_rank_', '_scores_rank_')).with_suffix('.json')
            if not paired.is_file():
                raise ValueError(f'Missing matching confidence JSON for {name}.')
            shutil.copyfile(original, pdb)
            shutil.copyfile(paired, scores)
        validate_prediction(pdb, scores, record)
        record.update(pdb=pdb.name, scores=scores.name, pdb_sha256=sha256(pdb), scores_sha256=sha256(scores))
    write_json(out / 'candidates.json', candidates)
    write_json(out / 'run.json', dict(stub=a.stub, command=command, msa_mode='single_sequence',
                                    packages=versions(['colabfold', 'jax', 'jaxlib', 'numpy', 'biopython']),
                                    input_sha256=sha256(a.sequences / 'complexes.fasta')))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError) as exc:
        sys.exit(f'BinderFlow ColabFold: {exc}')
