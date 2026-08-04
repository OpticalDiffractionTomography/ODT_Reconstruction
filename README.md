# ODT_Reconstruction

Automated pipeline for **Optical Diffraction Tomography (ODT)** — reconstructing 3D refractive index maps of biological samples from interferometric holographic measurements on an HPC cluster.

---

## Quick-start checklist

1. [Connect to the cluster](#1-connect-to-the-cluster)
2. [Clone the repository](#2-clone-the-repository)
3. [Configure cluster paths](#3-configure-cluster-paths-configsh)
4. [Install the CLI tool](#4-install-the-cli-tool)
5. [Run your first job](#5-run-your-first-job)

---

## 1. Connect to the cluster

Open a terminal (on Windows: **Command Prompt**, **PowerShell**, or **Windows Terminal**) and SSH into the cluster login node:

```bash
ssh <your_username>@<cluster_login_node>
# Example:
ssh ralajan@zpe.intranet.mpl.mpg.de
```

Enter your password or use your SSH key when prompted. All subsequent steps run inside this SSH session.

<details>
<summary>Tip — stay connected with tmux or screen</summary>  

> ----
> If your SSH session drops, any interactive process you started will die.
> Start a persistent shell before doing the setup steps:
> ```bash
> tmux new -s setup
> # or
> screen -S setup
> ```
> Once the `tomo_process` command is running (step 5), you can disconnect freely — it detaches automatically.
</details>


## 2. Clone the repository

```bash
# Choose a location in your home directory
cd ~
git clone https://github.com/RaghavaAlajangi/ODT_Reconstruction.git
cd ODT_Reconstruction
```

---

## 3. Configure cluster paths (`config.sh`)

Open `config.sh` in a text editor and update the paths to match your cluster:

```bash
nano config.sh    # or: vi config.sh
```

| Variable | Default | What to set it to |
|---|---|---|
| `DATA_MOUNT` | `/mnt/guck_division2` | Path where your raw data is accessible (read-only network mount) |
| `RESULTS_MOUNT` | `/mnt/ZPE_cluster_results` | Path where results should be written (writable network mount) |
| `SCRATCH_ROOT` | `$HOME/scratch/tomo_process` | Fast local scratch for staging; must be reachable from compute nodes |
| `MAX_PARALLEL_JOBS` | `5` | Maximum simultaneous SLURM jobs |
| `SLURM_TIME` | `20:00:00` | Walltime per job (HH:MM:SS) |
| `MINS_PER_SAMPLE` | `15` | Estimated processing time per sample file — used to auto-size job chunks |
| `SLURM_PARTITION` | `gpu` | SLURM partition for GPU processing jobs |

Save and close the file before continuing.

---

## 4. Install the CLI tool

Run the installer **once** per cluster account:

```bash
bash install.sh
source ~/.bashrc
```

Verify it worked:

```bash
tomo_process --help
```

<details>
<summary>What install.sh does</summary>

- Creates the scratch directory (`$SCRATCH_ROOT`)
- Symlinks `tomo_process` into `~/.local/bin` so it is available system-wide
- Adds `~/.local/bin` to your `PATH` in `~/.bashrc` if it is not already there
- Verifies that `sbatch` and the data mount are accessible

</details>

---

## 5. Run your first job

```bash
tomo_process \
  --email "you@institute.de" \
  --path "Members/YourName/experiment_folder"
```

- `--path` is relative to `$DATA_MOUNT/` (as set in `config.sh`)
- The folder can contain **multiple subdirectories** — all `sample*_Tomog.mat` files found recursively will be processed
- **You can close your terminal immediately** after this command — processing continues on the cluster in the background

You will receive an email when processing starts and again when it finishes (or fails).

**Results appear at:**
```
$RESULTS_MOUNT/Members/YourName/experiment_folder/field_retrieval_zpe_results/
```

<details>
<summary>All available options</summary>

```bash
tomo_process \
  --email "you@institute.de" \
  --path "Members/YourName/experiment_folder" \
  --max-jobs 3          # parallel SLURM jobs (default: 5)
  --mins-per-sample 10  # override time estimate if your samples are faster/slower
  --nm 1.340            # refractive index of culture medium (default: 1.337)
```

**How job chunk size is calculated automatically:**

```
max_samples_per_job = floor(SLURM_TIME / MINS_PER_SAMPLE)
                    = floor(20h × 60min / 15min)
                    = 80 samples per job
```

Files from different subdirectories are never mixed in one job — each subdirectory is processed independently, split across multiple jobs if it contains more samples than the per-job limit.

</details>

---

## If your SSH disconnects or the login node reboots

The orchestrator process detaches from your terminal automatically (via `nohup`), so a dropped SSH connection does **not** interrupt it. However, if the login node itself reboots, the orchestrator process will die. In that case, **no work is lost** — all state is saved on scratch and SLURM jobs already submitted keep running. Just reconnect and resume:

```bash
ssh <your_username>@<cluster_login_node>
tomo_process --list                          # find your run_id
tomo_process --resume tomo_20260729_143201   # restart the orchestrator from where it stopped
```

The orchestrator will re-adopt any SLURM jobs still in the queue and continue submitting remaining chunks.

---

## Monitoring a run

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

Full orchestrator log:
```
$SCRATCH_ROOT/<run_id>/logs/orchestrator.log
```

Per-job MATLAB output:
```
$SCRATCH_ROOT/<run_id>/logs/tomo_<run_id>_chunk000N.out
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

For each subdirectory, results are written to `$RESULTS_MOUNT/<rel_path>/field_retrieval_zpe_results/`:

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
            ├─ [login node] rsync chunk files: DATA_MOUNT → SCRATCH_ROOT/<run>/<chunk>/
            ├─ sbatch job script (up to 5 jobs in parallel on GPU nodes)
            │      └─ field_Retrieval.m          Stage 1: complex field retrieval
            │      └─ tomogram_Reconstruction.m  Stage 2: 3D RI reconstruction
            │      └─ touch DONE or FAIL marker
            │
            ├─ [login node] rsync results: SCRATCH_ROOT/<run>/<chunk>/ → RESULTS_MOUNT
            ├─ rm scratch data for that chunk
            └─ repeat until all chunks done, then send notification email
```

The login node handles all network mount access (`DATA_MOUNT` / `RESULTS_MOUNT`) because compute nodes typically do not have those mounts. Only scratch (`SCRATCH_ROOT`) needs to be accessible from compute nodes.

</details>

<details>
<summary>Manual submission (single dataset, advanced)</summary>

If your dataset fits within scratch space and you want to submit directly without the orchestrator:

```bash
sbatch scripts/main.sh --data_dir /path/to/experiment_data
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
| `n_m` | Refractive index of culture medium (default `1.337`; override with `--nm`) |
| `n_s` | Expected sample RI (cells); always `n_m + 0.04` |
| `ZP` | Hologram FFT size |
| `ZP2 = 512` | Lateral size of 3D Fourier volume |
| `ZP3 = 256` | Axial size of 3D Fourier volume |

</details>
