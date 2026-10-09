"""Small shared helpers for BinderFlow's command-line stages."""
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def versions(names):
    result = {}
    for name in names:
        try:
            result[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            result[name] = None
    return result


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def execute(command, logfile):
    with Path(logfile).open('w') as log:
        result = subprocess.run([str(v) for v in command], stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        tail = '\n'.join(Path(logfile).read_text(errors='replace').splitlines()[-20:])
        raise ValueError(f'Command failed ({result.returncode}); see {logfile}\n{tail}')


def chain_data(path):
    import numpy as np
    from Bio.PDB import PDBParser
    from Bio.SeqUtils import seq1
    model = PDBParser(QUIET=True).get_structure('complex', str(path))[0]
    chains = {}
    for chain in model:
        residues = [r for r in chain if r.id[0] == ' ']
        if not residues:
            continue
        if any(not all(a in r for a in ('N', 'CA', 'C')) for r in residues):
            raise ValueError(f'Missing backbone atoms in chain {chain.id}: {path}')
        if any(not np.isfinite(a.coord).all() for r in residues for a in r):
            raise ValueError(f'Nonfinite atom coordinates in chain {chain.id}: {path}')
        sequence = ''.join(seq1(r.resname) for r in residues)
        if set(sequence) - set('ACDEFGHIKLMNPQRSTVWY'):
            raise ValueError(f'Noncanonical residue in chain {chain.id}: {path}')
        chains[chain.id] = (sequence, residues)
    return model, chains


def mock_pdb(path, binder, target):
    """Synthetic geometry for software tests only; never a model prediction."""
    import math
    from Bio.SeqUtils import seq3
    lines = ['REMARK MOCK GEOMETRY - NOT A PREDICTED STRUCTURE']
    serial = 1
    for chain, sequence, offset in [('A', binder, 0), ('B', target, 7)]:
        for i, aa in enumerate(sequence, 1):
            center = (2 * math.cos(i * 1.7), 2 * math.sin(i * 1.7) + offset, i * 1.5)
            for atom, delta, element in [('N', (-1, 0, 0), 'N'), ('CA', (0, 0, 0), 'C'),
                                         ('C', (1, 0, 0), 'C'), ('O', (1, 1, 0), 'O')]:
                x, y, z = [a + b for a, b in zip(center, delta)]
                lines.append(f'ATOM  {serial:5d} {atom:^4s} {seq3(aa).upper():3s} {chain}{i:4d}    '
                             f'{x:8.3f}{y:8.3f}{z:8.3f}{1:6.2f}{85:6.2f}          {element:>2s}')
                serial += 1
        lines.append('TER')
    Path(path).write_text('\n'.join(lines + ['END']) + '\n')
