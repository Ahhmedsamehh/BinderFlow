# Complete workflow: implementation and validation

Local review against BinderFlow main commit
`c35eebaa39cb3d1269138f269f2e2b03787420c2`.

## Implemented

- Parameter-driven container, model directory, resource settings and output location.
- CSV validation, relative target paths and duplicate-ID rejection.
- Target-residue and hotspot checks before inference.
- Explicit binder checkpoint, indexed deterministic seeds and configurable noise.
- Quoted arguments, captured inference logs and error propagation.
- Paired PDB/TRB output checks, manifests and run provenance.
- Pinned RFdiffusion container recipe and official checkpoint download helper.
- A CPU-only mock workflow and automated contract/integration tests.
- ProteinMPNN sequence design that fixes the target chain, with FASTA manifests.
- ColabFold multimer prediction and matched rank-one structures/confidence files.
- Sequence, confidence, SASA, interface, hydrophobic exposure and RMSD features.
- Transparent per-target ranking, editable weights/filters, CSV/JSON/FASTA outputs.
- Full workflow by default, with an optional RFdiffusion-only stop point.
- Source-code fingerprints in task cache keys, including shared Python helpers.

Immunogenicity is excluded. Aggregation and solubility use explicitly labelled
physicochemical proxies; no dedicated aggregation/solubility model is included.

## Local verification

Test environment: macOS arm64, Python 3.9.6, Nextflow 25.10.4 and Temurin Java 17.
Analysis libraries: NumPy 1.26.4 and Biopython 1.85 in an isolated environment.
All 16 automated tests passed. `git diff --check` and download script shell
syntax checks also passed.
The test suite covers malformed parameters, missing weights, hotspot/target
validation, missing or empty outputs, a fake inference engine and its failure
status. Nextflow integration covers two samples, paths with spaces, single and
multiple designs, published manifest links, resume and duplicate-ID rejection.
Additional tests verify target-chain preservation with a fake ProteinMPNN
executable, paired ColabFold output parsing with a fake executable, confidence
dimensions, known hydropathy values, rotation/translation invariant alignment,
loss of interface contacts and buried SASA after chain separation, and ranking
ties/directions/filters. The full integration test runs two targets × two
backbones × two sequences, checks all eight result paths and resumes all 15 tasks.

All mock structures and metadata are labelled as such. No real model weights
were loaded and no designed binders were generated during these checks.

## Remaining GPU validation

The local machine has no Docker/NVIDIA execution environment. Before treating
the workflow as validated for production:

1. Build all three supplied Dockerfiles on the intended Linux/NVIDIA machine.
2. Download the RFdiffusion and ColabFold checkpoints and run one backbone with
   one ProteinMPNN sequence through the Docker profile.
3. Inspect the PDB/TRB, generated sequence, predicted complex, feature table and
   logs; confirm GPU use, target preservation and correct chain assignment.
4. Repeat with `-resume` and confirm that the task is cached.

The official ColabFold 1.5.5/CUDA 11.8 image manifest was checked and its digest
pinned. No container was built or executed locally. Upstream command-line
interfaces were checked against the versions documented in the README.

The example crop and hotspots are technical inputs, not a validated biological
design strategy. The ranking defaults require calibration against experimental
data; this implementation does not establish that any candidate binds or has
acceptable developability.
