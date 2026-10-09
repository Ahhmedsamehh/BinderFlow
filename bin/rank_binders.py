#!/usr/bin/env python3
"""Extract descriptors and rank candidates within each target using explicit rules."""
import argparse
import csv
import json
import math
from pathlib import Path
import sys
from binderflow_utils import sha256, versions, write_json
from developability import sequence_features, structure_features


def features(directory):
    directory = Path(directory)
    candidates = json.loads((directory / 'candidates.json').read_text())
    rows = []
    for candidate in candidates:
        pdb, scores = directory / candidate['pdb'], directory / candidate['scores']
        if sha256(pdb) != candidate['pdb_sha256'] or sha256(scores) != candidate['scores_sha256']:
            raise ValueError('Prediction hash mismatch.')
        row = dict(candidate)
        row.update(sequence_features(candidate['binder_sequence']))
        row.update(structure_features(pdb, scores, directory / 'backbone.pdb', candidate))
        base = f'{candidate["sample_id"]}/predictions/{candidate["backbone_id"]}'
        row['pdb'] = f'{base}/{pdb.name}'
        row['scores'] = f'{base}/{scores.name}'
        rows.append(row)
    return rows


def rank(rows, config):
    if not rows or len({r['candidate_id'] for r in rows}) != len(rows):
        raise ValueError('No candidates or duplicate candidate IDs.')
    if len({r['stub'] for r in rows}) != 1:
        raise ValueError('Cannot mix mock and real candidates.')
    rules = config['metrics']
    total = sum(rule['weight'] for rule in rules.values())
    if not rules or total <= 0 or any(not math.isfinite(rule['weight']) or rule['weight'] < 0
                                    or rule['direction'] not in ('higher', 'lower') for rule in rules.values()):
        raise ValueError('Invalid ranking weights/directions.')
    ranked = []
    for sample in sorted({r['sample_id'] for r in rows}):
        group = [dict(r) for r in rows if r['sample_id'] == sample]
        if len({r['target_sequence'] for r in group}) != 1:
            raise ValueError('Candidates within a sample must share the target sequence.')
        for row in group:
            failed = []
            for key, limits in config['gates'].items():
                value = row.get(key)
                if value is None or not math.isfinite(value) or ('min' in limits and value < limits['min']) or ('max' in limits and value > limits['max']):
                    failed.append(key)
            row['passes_filters'] = not failed
            row['failed_filters'] = ';'.join(failed)
            row['ranking_score'] = 0.0
        for key, rule in rules.items():
            values = [r.get(key) for r in group]
            if any(v is None or not isinstance(v, (int, float)) or not math.isfinite(v) for v in values):
                raise ValueError(f'Missing/nonfinite ranking metric: {key}')
            for row, value in zip(group, values):
                # Average percentile for ties. A one-candidate set gets 0.5.
                lower = sum(v < value for v in values)
                equal = sum(v == value for v in values)
                percentile = (lower + (equal - 1) / 2) / (len(group) - 1) if len(group) > 1 else 0.5
                if rule['direction'] == 'lower':
                    percentile = 1 - percentile
                row['ranking_score'] += rule['weight'] * percentile / total
        group.sort(key=lambda r: (not r['passes_filters'], -r['ranking_score'], r['candidate_id']))
        for i, row in enumerate(group, 1):
            row['rank_within_target'] = i
        ranked.extend(group)
    return ranked


def write_csv(path, rows):
    with Path(path).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    extract = sub.add_parser('features')
    extract.add_argument('--predictions', type=Path, required=True)
    extract.add_argument('--output', type=Path, required=True)
    aggregate = sub.add_parser('rank')
    aggregate.add_argument('--inputs', type=Path, nargs='+', required=True)
    aggregate.add_argument('--config', type=Path, required=True)
    args = parser.parse_args()
    if args.action == 'features':
        write_json(args.output, features(args.predictions))
    else:
        config = json.loads(args.config.read_text())
        rows = rank([r for path in args.inputs for r in json.loads(path.read_text())], config)
        write_json('ranked_candidates.json', rows)
        write_csv('ranked_candidates.csv', rows)
        Path('ranked_binders.fasta').write_text(''.join(f'>{r["candidate_id"]} target={r["sample_id"]} rank={r["rank_within_target"]} passes_filters={r["passes_filters"]} stub={r["stub"]}\n{r["binder_sequence"]}\n' for r in rows))
        write_json('ranking_method.json', dict(config=config, stub=rows[0]['stub'],
                    packages=versions(['numpy', 'biopython']),
                    note='Within-target heuristic ranking, not a probability of binding or developability. '
                         'Aggregation and solubility are descriptive proxies; no immunogenicity features.',
                    feature_files={p.name: sha256(p) for p in args.inputs}))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError) as exc:
        sys.exit(f'BinderFlow ranking: {exc}')
