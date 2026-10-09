#!/usr/bin/env python3
"""Validate a binder-design task, run RFdiffusion, and audit its outputs."""
import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import re
import subprocess
import sys


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def validate_target(target, contigs, hotspots):
    # This first stage intentionally supports one receptor segment and one binder.
    match = re.fullmatch(r'([A-Za-z])([1-9][0-9]*)-([1-9][0-9]*)/0 ([1-9][0-9]*)(?:-([1-9][0-9]*))?', contigs)
    if not match:
        raise ValueError('Use contigs like C3-23/0 70-100 (one target segment and one binder).')
    chain, start, end, low, high = match.groups()
    start, end, low, high = int(start), int(end), int(low), int(high or low)
    if start > end or low > high:
        raise ValueError('Contig residue and binder-length ranges must be increasing.')
    atoms = {}
    with Path(target).open() as handle:
        for line in handle:
            if line.startswith('ENDMDL'):
                break
            if not line.startswith('ATOM  '):
                continue
            if len(line) < 54:
                raise ValueError('Truncated ATOM record in target PDB.')
            key = (line[21], int(line[22:26]))
            if key[0] == chain and start <= key[1] <= end:
                if line[26].strip():
                    raise ValueError('Insertion-coded target residues are not supported.')
                coordinates = [float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54))]
                if not all(math.isfinite(v) for v in coordinates):
                    raise ValueError('Nonfinite coordinates in target segment.')
                atoms.setdefault(key, set()).add(line[12:16].strip())
    missing = [f'{chain}{i}' for i in range(start, end + 1)
               if not {'N', 'CA', 'C'}.issubset(atoms.get((chain, i), set()))]
    if missing:
        raise ValueError('Target segment lacks N/CA/C atoms: ' + ', '.join(missing[:10]))
    values = hotspots.split(',') if hotspots else []
    if len(values) != len(set(values)):
        raise ValueError('Hotspot residues must be unique.')
    for value in values:
        hit = re.fullmatch(r'([A-Za-z])([1-9][0-9]*)', value)
        if not hit or hit[1] != chain or not start <= int(hit[2]) <= end:
            raise ValueError(f'Hotspot {value!r} is outside the retained target segment.')
    return values


def collect_outputs(directory, sample_id, start, count, stub=False):
    directory = Path(directory)
    expected = {f'{sample_id}_{i}' for i in range(start, start + count)}
    for suffix in ('pdb', 'trb'):
        found = {p.stem for p in directory.glob(f'{sample_id}_*.{suffix}')}
        if found != expected:
            raise ValueError(f'Expected {count} paired designs; {suffix} output names do not match.')
    rows = []
    for i in range(start, start + count):
        name = f'{sample_id}_{i}'
        pdb, trb = directory / f'{name}.pdb', directory / f'{name}.trb'
        if not pdb.stat().st_size or not trb.stat().st_size:
            raise ValueError(f'Empty output for design {name}.')
        if not stub and not any(line.startswith('ATOM  ') for line in pdb.read_text().splitlines()):
            raise ValueError(f'No ATOM records in {pdb.name}.')
        rows.append({'sample_id': sample_id, 'design_id': name,
                     'pdb': f'../backbones/{pdb.name}', 'trb': f'../metadata/{trb.name}',
                     'stub': stub, 'pdb_sha256': sha256(pdb), 'trb_sha256': sha256(trb)})
    with (directory / 'manifest.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample-id', required=True)
    parser.add_argument('--target', type=Path, required=True)
    parser.add_argument('--contigs', required=True)
    parser.add_argument('--hotspots', default='')
    parser.add_argument('--models', type=Path, required=True)
    parser.add_argument('--script', type=Path, required=True)
    parser.add_argument('--num-designs', type=int, required=True)
    parser.add_argument('--design-startnum', type=int, default=0)
    parser.add_argument('--deterministic', choices=('true', 'false'), default='true')
    parser.add_argument('--noise-scale-ca', type=float, default=0)
    parser.add_argument('--noise-scale-frame', type=float, default=0)
    parser.add_argument('--stub', action='store_true')
    args = parser.parse_args(argv)
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', args.sample_id):
        raise ValueError('Invalid sample ID.')
    if args.num_designs < 1 or args.design_startnum < 0:
        raise ValueError('num-designs must be positive; design-startnum must be nonnegative.')
    if any(not math.isfinite(v) or v < 0 for v in (args.noise_scale_ca, args.noise_scale_frame)):
        raise ValueError('Noise scales must be finite and nonnegative.')
    hotspots = validate_target(args.target, args.contigs, args.hotspots)
    checkpoint = args.models / 'Complex_base_ckpt.pt'
    if not args.stub:
        if not args.script.is_file():
            raise ValueError(f'RFdiffusion script not found: {args.script}')
        if not checkpoint.is_file() or checkpoint.stat().st_size == 0:
            raise ValueError(f'Missing or empty binder checkpoint: {checkpoint}')

    out = Path('outputs')
    out.mkdir(exist_ok=True)
    # JSON-quoted Hydra string values and argv execution preserve spaces safely.
    command = [sys.executable, str(args.script.resolve()),
               'inference.input_pdb=' + json.dumps(str(args.target.resolve())),
               'inference.output_prefix=' + json.dumps(str(out / args.sample_id)),
               'inference.model_directory_path=' + json.dumps(str(args.models.resolve())),
               'inference.ckpt_override_path=' + json.dumps(str(checkpoint.resolve())),
               f'inference.num_designs={args.num_designs}',
               f'inference.design_startnum={args.design_startnum}',
               f'inference.deterministic={args.deterministic}',
               'inference.write_trajectory=False',
               'contigmap.contigs=' + json.dumps([args.contigs]),
               'ppi.hotspot_res=' + json.dumps(hotspots),
               f'denoiser.noise_scale_ca={args.noise_scale_ca}',
               f'denoiser.noise_scale_frame={args.noise_scale_frame}']
    versions = {}
    for package in ('torch', 'dgl', 'hydra-core', 'numpy'):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    report = {'stub': args.stub, 'arguments': {k: str(v) if isinstance(v, Path) else v
                                             for k, v in vars(args).items()},
              'command': command, 'python': sys.version, 'packages': versions,
              'target_sha256': sha256(args.target)}
    if not args.stub:
        report.update(script_sha256=sha256(args.script), checkpoint_sha256=sha256(checkpoint))
    (out / 'run.json').write_text(json.dumps(report, indent=2) + '\n')
    if args.stub:
        for i in range(args.design_startnum, args.design_startnum + args.num_designs):
            (out / f'{args.sample_id}_{i}.pdb').write_text('REMARK MOCK OUTPUT - NOT A DESIGNED STRUCTURE\nEND\n')
            (out / f'{args.sample_id}_{i}.trb').write_text('MOCK OUTPUT - NOT RFdiffusion METADATA\n')
        (out / 'inference.log').write_text('STUB: no model was loaded and no inference was performed.\n')
    else:
        with (out / 'inference.log').open('w') as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode:
            print('\n'.join((out / 'inference.log').read_text(errors='replace').splitlines()[-30:]), file=sys.stderr)
            raise ValueError(f'RFdiffusion exited with status {result.returncode}; see outputs/inference.log.')
    collect_outputs(out, args.sample_id, args.design_startnum, args.num_designs, args.stub)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError) as exc:
        sys.exit(f'BinderFlow: {exc}')
