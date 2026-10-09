# BinderFlow

A Nextflow workflow that takes a target PDB through **RFdiffusion → ProteinMPNN →
ColabFold/AlphaFold2-Multimer → developability features → ranking**.

ProteinMPNN designs only the generated binder chain. ColabFold predicts each
binder–target complex. The ranker compares its geometry with the designed
backbone and combines explicit confidence, interface and developability criteria.
It exports every candidate and its component features, including failed filters.

![BinderFlow workflow DAG: a sample sheet fans out through RFdiffusion, ProteinMPNN,
ColabFold and developability features once per design, then gathers into a single
ranking task.](docs/workflow-dag.svg)

The workflow fans out once: `rf_diffusion` runs per sample, and its backbones are
split so that every later per-design stage runs once per backbone. The
`--sequences_per_backbone` sequences of a backbone are carried through that single
task together, not as separate tasks. Only `rank_candidates` gathers, comparing
every design in one task.

**Status:** all stages are implemented and tested with mock model engines;
real GPU inference and container builds remain unverified. Candidates are not
experimentally validated binders. Immunogenicity is excluded. Aggregation and
solubility outputs are physicochemical descriptors/proxies, not dedicated
predictor scores. See [feature definitions](docs/FEATURES.md).

## Setup from a clone

This repository tracks only workflow code: no model weights and no third-party
model source. One script fetches everything large, into paths git ignores:

```bash
git clone https://github.com/Ahhmedsamehh/BinderFlow.git
cd BinderFlow
bash scripts/setup.sh
```

That fetches RFdiffusion and ProteinMPNN at the revisions pinned in
`containers/*.Dockerfile`, the RFdiffusion binder checkpoint (about 484 MB) and
the AlphaFold2 multimer-v3 parameters (several GB), then prints the run command
with the right paths filled in. Expect roughly 5 GB of disk. Re-running skips
any step already complete, and `THIRD_PARTY_DIR`, `MODELS_DIR` and
`COLABFOLD_DIR` override where each part lands. The individual scripts
(`setup_third_party.sh`, `download_models.sh`, `download_colabfold_models.sh`)
can also be run on their own; none of them needs Docker.

You still need Nextflow 25.10+ and either Docker or a conda-family tool. See
[Conda/Mamba profile](#condamamba-profile) for GPU hosts without Docker.

## Quick workflow check (no GPU)

Requires Nextflow 25.10 or later, Java 17 or later supported by your Nextflow
version, and Python 3.9–3.12. Install the analysis dependencies in a virtual
environment, then run from this repository:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
nextflow run main.nf -profile stub -stub-run \
  --input sample-data/sample-sheet.csv --num_designs 1 \
  --sequences_per_backbone 2 --outdir stub-results
```

This exercises the complete workflow, including real feature calculations on
synthetic test geometry and ranking. All structures and metadata
are explicitly marked **MOCK** and are unsuitable for scientific analysis.
The `stub` profile requires `-stub-run`.

## Run the complete workflow on a GPU

Use a Linux x86-64 machine with an NVIDIA GPU, a compatible NVIDIA driver, Docker
and the NVIDIA Container Toolkit. The supplied image uses CUDA 11.6.2 and pins
RFdiffusion to commit `86507b6538f51fce57b5a72477165f03999ed7ae`.
ProteinMPNN is pinned to `8907e6671bfbfc92303b5f79c4b5e6ce47cdef57`; ColabFold
uses the official 1.5.5/CUDA 11.8 image pinned by digest. Choose a host driver
compatible with both CUDA runtimes and sufficient GPU memory for your complexes.
Building the images and real inference have not yet been validated in this project;
see [validation status](docs/VALIDATION.md).

```bash
docker build --platform linux/amd64 -f containers/Dockerfile \
  -t binderflow-rfdiffusion:86507b6-cu116 .
docker build --platform linux/amd64 -f containers/ProteinMPNN.Dockerfile \
  -t binderflow-proteinmpnn:8907e66 .
docker build --platform linux/amd64 -f containers/Features.Dockerfile \
  -t binderflow-features:0.2.0 .

bash scripts/download_models.sh "$PWD/models/rfdiffusion"
bash scripts/download_colabfold_models.sh "$PWD/models/colabfold"

nextflow run main.nf -profile docker -params-file params.yaml \
  --models "$PWD/models/rfdiffusion" --colabfold_models "$PWD/models/colabfold"
```

The downloader retrieves the official `Complex_base_ckpt.pt` binder checkpoint
(approximately 484 MB). An existing nonempty file is reused. Each real run records
its SHA-256; this records which weights were used, but does not authenticate an
existing download against an upstream checksum. The separate ColabFold helper
downloads the multimer-v3 model weights, which need several GB of disk space.
ProteinMPNN's weights are included in its repository/image.

For a GPU machine with an existing RFdiffusion installation and its dependencies:

```bash
nextflow run main.nf -profile pod_native -params-file params.yaml \
  --models /absolute/path/to/models \
  --python /absolute/path/to/environment/bin/python \
  --rfdiffusion_script /absolute/path/to/RFdiffusion/scripts/run_inference.py \
  --proteinmpnn_python /absolute/path/to/mpnn/environment/bin/python \
  --proteinmpnn_script /absolute/path/to/ProteinMPNN/protein_mpnn_run.py \
  --colabfold_python /absolute/path/to/colabfold/environment/bin/python \
  --colabfold_executable /absolute/path/to/colabfold/environment/bin/colabfold_batch \
  --colabfold_models /absolute/path/to/colabfold/weights \
  --features_python /absolute/path/to/analysis/environment/bin/python
```

Use the appropriate installed dependencies in each environment. ProteinMPNN
also needs Biopython; the analysis environment uses `requirements.txt`. Download
ColabFold weights with `scripts/download_colabfold_models.sh` into the directory
passed to `--colabfold_models` (it contains a `params/` subdirectory). On either
profile, use `-resume` to reuse successfully cached tasks. Keep the Nextflow work
directory and cache for resuming. `--container_rfdiffusion` can select another
compatible image.

## Conda/Mamba profile

If Docker is unavailable, use Nextflow-managed environments. Run
`bash scripts/setup.sh` first, then:

```bash
nextflow run main.nf -profile micromamba -params-file params.yaml \
  --models "$PWD/models/rfdiffusion" \
  --colabfold_models "$PWD/models/colabfold" \
  --rfdiffusion_script "$PWD/RFdiffusion/scripts/run_inference.py" \
  --proteinmpnn_script "$PWD/ProteinMPNN/protein_mpnn_run.py"
```

Pick the profile matching your tool: `micromamba` sets
`conda.useMicromamba=true`, `mamba` sets `conda.useMamba=true`, and `conda`
uses `conda` directly. Use `micromamba` when `conda` on your `PATH` is a
micromamba shim, because micromamba provides no `bin/activate` for the other
profiles to source. All three disable containers and build the environment
files under `envs/`.

These environments supply dependencies only. RFdiffusion itself and NVIDIA's
SE3Transformer are imported from the checkout that `--rfdiffusion_script` points
into, so that path must be a complete RFdiffusion repository rather than a bare
copy of `run_inference.py`. Environments may take several minutes to solve the
first time; Nextflow caches them for subsequent runs.

ColabFold uses `single_sequence` mode without templates: sequences are not sent
to an external MSA service. This is a screening configuration, not a guarantee
of the best prediction. Default prediction uses one model and three recycles;
increase `--af_num_models` up to 5 for additional model sampling. The rank-one
model under ColabFold's multimer confidence ranking is used for features; raw
outputs are retained. Backbone coordinates are used for comparison, not supplied
as a structure template to ColabFold.

Use `--stop_after rfdiffusion` to run only backbone generation.

## Input

The sample sheet has four required columns:

```csv
id,target_pdb,contigs,hotspots
ctla4,ctla4.pdb,C3-23/0 70-100,C3;C10;C20;C15
```

- `id`: unique letters, numbers, underscores or hyphens; must start with a letter or number.
- `target_pdb`: absolute path, or a path relative to the sample sheet's directory.
- `contigs`: one continuous target segment followed by one binder. `C3-23/0 70-100`
  retains chain C residues 3–23 and generates a binder of 70–100 residues; `/0`
  separates the chains. A fixed binder length, such as `C3-23/0 80`, also works.
- `hotspots`: semicolon-separated residues in the retained target segment, using
  the original PDB numbering. Leave the value empty for no hotspot conditioning.

This version rejects missing backbone N/CA/C atoms, insertion-coded residues in
the selected segment, invalid ranges and hotspots outside the segment. Multiple
target segments and more complex RFdiffusion contigs are not supported yet.

The included CTLA4 crop is a small technical example. Its 21-residue target
segment and chosen hotspots have **not** been established as a biologically
appropriate design target. Choose and inspect a suitable target context before
performing a scientific campaign.

## Settings and outputs

`params.yaml` runs five designs for the example. Useful overrides:

| Parameter | Default | Purpose |
| --- | --- | --- |
| `--num_designs` | 10 (5 in example YAML) | Backbones per sample |
| `--outdir` | `results` | Published output directory |
| `--design_startnum` | 0 | Starting design index |
| `--deterministic` | true | RFdiffusion seeds each design using its index |
| `--noise_scale_ca`, `--noise_scale_frame` | 0 | Denoising noise scales |
| `--cpus`, `--memory` | 2, `8 GB` | CPU resources per task; memory is not GPU VRAM |
| `--max_forks` | 1 | Concurrent tasks; keep at 1 on a single GPU |
| `--sequences_per_backbone` | 4 | ProteinMPNN sequences per backbone |
| `--mpnn_temperature`, `--sequence_seed` | 0.1, 1 | ProteinMPNN sampling |
| `--af_num_models`, `--af_recycles`, `--structure_seed` | 1, 3, 1 | Complex prediction |
| `--af_memory` | `16 GB` | Host RAM for a complex-prediction task |
| `--ranking_config` | `config/ranking.json` | Feature weights and screening gates |

`--models` may also be supplied through the `MODELS_PATH` environment variable.
Zero denoising noise alone does not establish determinism. With deterministic
mode, change `design_startnum` to use different indexed seeds. Results may still
differ across hardware or software versions. Raising `max_forks` does not assign
separate GPUs to tasks. The local executor queue is limited to one task to avoid
overlap between GPU stages. Multi-GPU scheduling is not implemented.

```text
results/
├── ranked_candidates.csv / ranked_candidates.json
├── ranked_binders.fasta
├── ranking_method.json
├── features/<backbone_id>.features.json
└── ctla4/
    ├── backbones/ctla4_0.pdb ...
    ├── metadata/ctla4_0.trb ...
    ├── reports/manifest.csv, run.json, inference.log
    ├── sequences/ctla4_0/   # FASTA files, candidates, original backbone, raw MPNN outputs
    └── predictions/ctla4_0/ # PDB/score pairs, candidates, original backbone, raw AF outputs
```

`manifest.csv` links paired PDB/TRB files using paths relative to the report
directory, with file hashes and a mock-run flag. `run.json` records inputs,
settings, the executed command, package versions and hashes of the target,
inference script and checkpoint (the last two only for real runs). `inference.log`
contains RFdiffusion output. Failed tasks retain logs in their Nextflow work folder.
The workflow rejects missing, empty or incorrectly numbered design pairs.
Sequence and prediction stages have their own `run.json` and `inference.log`.
Ranked-table structure paths are relative to the results directory. Rankings are
**within each target**, with filter-passing designs first; the score is a weighted
average of within-target percentiles, not a probability or an absolute quality
scale. Inspect `failed_filters` and the individual features before selecting a
candidate. All candidates are retained even when none passes the filters.

## Tests

```bash
python3 -m unittest discover -s tests -v

# Also exercise Nextflow staging, publishing, resume and invalid input:
NEXTFLOW_BIN="$(command -v nextflow)" python3 -m unittest discover -s tests -v
```

Tests include fake inference executables, chain preservation, feature geometry,
ranking directions/ties, multi-target execution, output links and resume. They
do not establish model accuracy or GPU compatibility.

## Credits

Backbone generation uses [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion).
The container follows its [upstream Dockerfile](https://github.com/RosettaCommons/RFdiffusion/blob/86507b6538f51fce57b5a72477165f03999ed7ae/docker/Dockerfile).
Sequence design uses [ProteinMPNN](https://github.com/dauparas/ProteinMPNN), complex
prediction uses [ColabFold](https://github.com/sokrypton/ColabFold), and structural
analysis uses [Biopython](https://biopython.org/). Consult these projects for
their licenses, model terms and citations.
