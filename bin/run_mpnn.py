#!/usr/bin/env python3
"""Design only the binder chain; preserve the receptor sequence for prediction."""
import argparse
import json
import math
from pathlib import Path
import re
import shutil
import sys
from binderflow_utils import chain_data, execute, mock_pdb, sha256, versions, write_json


def parse_sequences(path, count, length):
    records = []
    for record in Path(path).read_text().split('>')[1:]:
        header, *lines = record.strip().splitlines()
        if not re.search(r'\bsample=\d+', header):
            continue  # ProteinMPNN's first FASTA record is the input sequence.
        sequence = ''.join(lines).strip()
        score = re.search(r'(?:^|,\s*)score=([-+0-9.eE]+)', header)
        if len(sequence) != length or set(sequence) - set('ACDEFGHIKLMNPQRSTVWY'):
            raise ValueError('ProteinMPNN returned a malformed binder sequence.')
        if not score or not math.isfinite(float(score[1])):
            raise ValueError('ProteinMPNN sequence score missing or nonfinite.')
        records.append((sequence, float(score[1])))
    if len(records) != count:
        raise ValueError(f'Expected {count} ProteinMPNN sequences, received {len(records)}.')
    return records


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--backbone', type=Path, required=True)
    p.add_argument('--sample-id', required=True)
    p.add_argument('--backbone-id', required=True)
    p.add_argument('--contigs', required=True)
    p.add_argument('--script', type=Path, required=True)
    p.add_argument('--count', type=int, default=4)
    p.add_argument('--temperature', type=float, default=0.1)
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--stub', action='store_true')
    a = p.parse_args()
    if a.count < 1 or a.seed < 1 or not math.isfinite(a.temperature) or a.temperature <= 0:
        raise ValueError('MPNN count/seed/temperature must be positive and finite.')
    out = Path('sequences')
    out.mkdir()
    if a.stub:
        binder, target = 'AEKLAEKLAEKLAEKLAEKL', 'VTISCTGSSSNIGAGYDVH'
        mock_pdb(out / 'backbone.pdb', binder, target)
        binder_chain, target_chain = 'A', 'B'
        sequences = [(binder[:-1] + 'ALV'[i % 3], 0.5 + i * 0.1) for i in range(a.count)]
        command = []
        (out / 'inference.log').write_text('MOCK sequences; ProteinMPNN was not run.\n')
    else:
        _, chains = chain_data(a.backbone)
        # Pinned RFdiffusion preserves the input target chain ID; the newly
        # generated chain is the other one. Do not blindly assume A/B.
        target_chain = a.contigs[0]
        others = [c for c in chains if c != target_chain]
        if target_chain not in chains or len(chains) != 2 or len(others) != 1:
            raise ValueError('Expected the retained target chain plus exactly one binder chain.')
        binder_chain = others[0]
        binder, target = chains[binder_chain][0], chains[target_chain][0]
        match = re.fullmatch(r'[A-Za-z](\d+)-(\d+)/0 (\d+)(?:-(\d+))?', a.contigs)
        if not match or len(target) != int(match[2]) - int(match[1]) + 1:
            raise ValueError('Backbone target length does not match contigs.')
        if not int(match[3]) <= len(binder) <= int(match[4] or match[3]):
            raise ValueError('Backbone binder length does not match contigs.')
        if not a.script.is_file():
            raise ValueError(f'ProteinMPNN script missing: {a.script}')
        shutil.copyfile(a.backbone, out / 'backbone.pdb')
        command = [sys.executable, str(a.script), '--pdb_path', str(a.backbone),
                   '--pdb_path_chains', binder_chain, '--num_seq_per_target', str(a.count),
                   '--sampling_temp', str(a.temperature), '--seed', str(a.seed),
                   '--model_name', 'v_48_020', '--batch_size', '1', '--out_folder', str(out / 'raw')]
        execute(command, out / 'inference.log')
        sequences = parse_sequences(out / 'raw/seqs' / (a.backbone.stem + '.fa'), a.count, len(binder))
    candidates = []
    for i, (sequence, score) in enumerate(sequences, 1):
        candidates.append(dict(candidate_id=f'{a.backbone_id}__s{i:03d}', sample_id=a.sample_id,
                               backbone_id=a.backbone_id, binder_sequence=sequence, target_sequence=target,
                               reference_binder_chain=binder_chain, reference_target_chain=target_chain,
                               mpnn_score=score, stub=a.stub))
    write_json(out / 'candidates.json', candidates)
    (out / 'complexes.fasta').write_text(''.join(f'>{r["candidate_id"]}\n{r["binder_sequence"]}:{target}\n' for r in candidates))
    (out / 'binders.fasta').write_text(''.join(f'>{r["candidate_id"]}\n{r["binder_sequence"]}\n' for r in candidates))
    write_json(out / 'run.json', dict(stub=a.stub, command=command, seed=a.seed,
                                    packages=versions(['torch', 'numpy', 'biopython']),
                                    input_sha256=sha256(a.backbone),
                                    script_sha256=None if a.stub else sha256(a.script)))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError) as exc:
        sys.exit(f'BinderFlow ProteinMPNN: {exc}')
