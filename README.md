# field_tomogram_reconstruction

Automated pipeline for **Optical Diffraction Tomography (ODT)** — reconstructing 3D refractive index maps of biological samples from interferometric holographic measurements on the ZPE HPC cluster.

---

## First-time setup

> Do this once per cluster account after cloning the repo.

```bash
# 1. Clone onto the cluster
git clone <repo_url> /beegfs/home/ralajan/matlab/field_tomogram_reconstruction
cd /beegfs/home/ralajan/matlab/field_tomogram_reconstruction

# 2. Install the CLI tool
bash install.sh
source ~/.bashrc
```

<details>
<summary>What install.sh does</summary>

- Creates the scratch directory (`/beegfs/home/ralajan/scratch/tomo_process/`)
- Symlinks `tomo_process` into `~/.local/bin` so it's available system-wide
- Adds `~/.local/bin` to your `PATH` in `~/.bashrc` if not already there
- Verifies that `sbatch` and the network mount are accessible

</details>

<details>
<summary>Configuring cluster settings (config.sh)</summary>

All tunable parameters live in `config.sh`. Edit once if your cluster differs:

| Variable | Default | Meaning |
|---|---|---|
| `GUCK_DIVISION_2_MOUNT` | `/mnt/guck_division2` | Read-only input data mount |
| `ZPE_RESULTS_MOUNT` | `/mnt/ZPE_cluster_results` | Writable results mount |
| `SCRATCH_ROOT` | `/beegfs/home/ralajan/scratch/tomo_process` | Staging area for job data |
| `MAX_PARALLEL_JOBS` | `5` | Max simultaneous SLURM jobs |
| `SLURM_TIME` | `20:00:00` | Walltime per job |
| `MINS_PER_SAMPLE` | `15` | Estimated minutes per sample file (sets chunk size) |
| `SLURM_PARTITION` | `gpu` | SLURM partition for processing jobs |
| `ORCH_PARTITION` | `cpu` | SLURM partition for the orchestrator |

</details>

---

## Processing a dataset

```bash
tomo_process --email "you@institute.de" --path "Members/YourName/experiment_folder"
```

- `--path` is relative to `/mnt/guck_division2/`
- The folder can contain **multiple subdirectories** — all `sample*_Tomog.mat` files found recursively will be processed
- You can **disconnect immediately** after running this command — processing continues on the cluster

**Results appear at:**
```
/mnt/ZPE_cluster_results/Members/YourName/experiment_folder/<subdir>/field_retrieval_zpe_results/
```

<details>
<summary>Additional options</summary>

```bash
tomo_process \
  --email "you@institute.de" \
  --path "Members/YourName/experiment_folder" \
  --max-jobs 3          # parallel SLURM jobs (default: 5)
  --mins-per-sample 10  # override time estimate if your samples are faster/slower
```

**How chunk size is calculated automatically:**

```
max_samples_per_job = floor(SLURM_TIME / MINS_PER_SAMPLE)
                    = floor(20h × 60min / 15min)
                    = 80 samples per job
```

Files from different subdirectories are never mixed in one job — each subdirectory is processed independently, split into multiple jobs if it has more samples than the limit.

</details>

---

## If your SSH disconnects or the login node reboots

The orchestrator process runs on the login node. If it dies mid-run (SSH drop, reboot), **no work is lost** — all state is saved on `/beegfs`. SLURM jobs already submitted will keep running. Just reconnect and resume:

```bash
tomo_process --list                          # find your run_id
tomo_process --resume tomo_20260729_143201   # restart the orchestrator from where it stopped
```

The orchestrator will re-adopt any SLURM jobs still in the queue and continue submitting remaining chunks.

---

## Monitoring

```bash
tomo_process --list                          # list all runs and their status
tomo_process --status tomo_20260729_143201   # detailed status of one run
tomo_process --cancel tomo_20260729_143201   # cancel a run and its active jobs
```

<details>
<summary>Reading the status output</summary>

```
=== Run: tomo_20260729_143201 ===
  Total chunks : 6      ← total number of jobs submitted
  Done         : 4      ← completed and results copied back
  Active       : 2      ← currently running on cluster
  Pending      : 0      ← waiting to be submitted
  Failed       : 0      ← check logs if this is non-zero
```

Full logs are at:
```
/beegfs/home/ralajan/scratch/tomo_process/<run_id>/logs/orchestrator.log
```

Per-job MATLAB output is at:
```
/beegfs/home/ralajan/scratch/tomo_process/<run_id>/logs/tomo_<run_id>_chunk000N.out
```

</details>

---

## Data requirements

Each subdirectory under `--path` must contain:

```
<subdir>/
  bg001_Tomog.mat          ← single shared background hologram stack
  sample001_Tomog.mat      ┐
  sample002_Tomog.mat      ├ one or more sample hologram stacks
  sample003_Tomog.mat      ┘
```

<details>
<summary>What the pipeline produces</summary>

For each subdirectory, results are written to `/mnt/ZPE_cluster_results/<rel_path>/field_retrieval_zpe_results/`:

```
field_retrieval_zpe_results/
  Field_sample001_Tomog.mat    # retPhase, retAmplitude, NA, lambda, res, ZP, f_dx, f_dy
  Field_sample001_Tomog.png    # diagnostic phase overview image
  Tomogram_Field_sample001.mat # Reconimg (3D RI volume), res3, res4, excludeFrame
  Tomogram_Field_sample001.tif # multi-page uint16 TIFF (values × 10000)
  Tomogram_Field_sample001.png # diagnostic orthogonal slice image
  Field_inspection.png         # mean phase and frame-diff curves across all samples
  field_retrieval.log          # Stage 1 log
  tomogram_reconstruction.log  # Stage 2 log
```

</details>

---

<details>
<summary>How it works internally</summary>

```
tomo_process (runs on login node, exits immediately after submitting)
  │
  └─► nohup scripts/orchestrator.sh (stays running on login node as background process)
            │
            │  Groups sample files by subdirectory, auto-sizes chunks from walltime
            │
            ├─ [login node] rsync chunk files: /mnt/guck_division2 → /beegfs/scratch/<run>/<chunk>/
            ├─ sbatch job script (up to 5 jobs in parallel on GPU nodes)
            │      └─ field_Retrieval.m       Stage 1: complex field retrieval
            │      └─ tomogram_Reconstruction.m  Stage 2: 3D RI reconstruction
            │      └─ touch DONE or FAIL marker
            │
            ├─ [login node] rsync results: /beegfs/scratch/<run>/<chunk>/ → /mnt/ZPE_cluster_results
            ├─ rm scratch data for that chunk
            └─ repeat until all chunks done, then send notification email
```

The login node handles all network mount access (`/mnt/...`) because compute nodes do not have those mounts. Only scratch (`/beegfs`) is accessible from compute nodes.

</details>

<details>
<summary>Manual submission (single dataset, advanced)</summary>

If your dataset fits within scratch space and you want to submit directly:

```bash
sbatch scripts/main.sh --data_dir /beegfs/home/ralajan/scratch/my_experiment
```

`scripts/main.sh` runs both stages sequentially on one GPU node. It does not handle copying from/to the network mounts — you must do that manually.

</details>

<details>
<summary>Pipeline stages (technical reference)</summary>

### Stage 1 — Field Retrieval (`src/field_Retrieval.m`)

For each `sample*_Tomog.mat` file:

1. Detect off-axis carrier frequency from background frames via FFT peak
2. Build circular demodulation mask (`mk_ellipse`)
3. Demodulate each frame: shift +1 diffraction order to centre, divide by background
4. Phase pipeline per frame: `angle` → `unwrap2` → `PhiShift` → `phaseCompensation(deg=1)` → cell-exclusion mask → `phaseCompensation(deg=2)`
5. Save `retPhase`, `retAmplitude` and diagnostic PNG

Frames where `unwrap2` fails (e.g. saturated or empty patches) are set to NaN and excluded automatically in Stage 2.

### Stage 2 — Tomogram Reconstruction (`src/tomogram_Reconstruction.m`)

For each `Field_*.mat` output from Stage 1:

1. Detect outlier frames (`excludeFrame`): mean phase > 1.5, frame-to-frame diff > 0.1, NaN frames
2. Build `TomoParam` struct with all optical/geometric parameters
3. `ODTReconstruction` — Ewald sphere mapping into 3D Fourier volume
4. `ODTIteration` (100 iterations, GPU) — enforce `n ≥ n_medium` physical constraint
5. Save `Reconimg` as `.mat`, multi-page `.tif`, and diagnostic `.png`

### Key physical parameters

| Variable | Meaning |
|---|---|
| `lambda` | Illumination wavelength (µm) |
| `NA` | Objective numerical aperture |
| `res` | Camera pixel size (µm) |
| `n_m = 1.337` | Refractive index of culture medium |
| `n_s ≈ 1.377` | Expected sample RI (cells) |
| `ZP` | Hologram FFT size |
| `ZP2 = 512` | Lateral size of 3D Fourier volume |
| `ZP3 = 256` | Axial size of 3D Fourier volume |

</details>
