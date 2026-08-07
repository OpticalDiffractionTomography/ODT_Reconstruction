# ODT_Reconstruction

Automated pipeline for **Optical Diffraction Tomography (ODT)** — reconstructing 3D refractive index maps of biological samples from interferometric holographic measurements on an HPC cluster.

> **First time here?** Do the [One-time setup](#one-time-setup-do-once) first (4 steps, ~5 minutes).
> **Already set up?** You only ever need the 2 steps below.

---

## Run jobs (2 steps)

**Step 1 — log in to the cluster:**

```bash
ssh <your_username>@zpe.intranet.mpl.mpg.de
```

**Step 2 — start processing:**

```bash
tomo_process --path "<paste your dataset path here>"
```

Example — copy the folder path from Windows File Explorer and paste it as it is:

```bash
tomo_process --path "U:\Members\YourName\20260715_experiment"
```

Done — you can **close the terminal immediately**; processing continues on the cluster. You will get an email when it starts and when it finishes.

Good to know:

- Windows paths (drive letter, backslashes) are converted automatically — no need to edit them.
- The folder can contain multiple subdirectories; all `sample*_Tomog.mat` files are found and processed.
- **Re-running the same dataset overwrites the old results.**
- Need to change the refractive index of the medium? Add `--nm`, e.g. `--nm 1.340` (default: `1.337`). See [Options](#options).

**Results appear at:**

```
<RESULTS_MOUNT>/Members/YourName/experiment_folder/.../field_retrieval_zpe_results/
```

---

## One-time setup (do once)

These 4 steps are needed **only once** per cluster account. After that, use the 2 steps above.

### 1. Log in to the cluster

Open a terminal (on Windows: **Command Prompt** or **PowerShell**) and connect:

```bash
ssh <your_username>@zpe.intranet.mpl.mpg.de
```

### 2. Get the code

```bash
cd ~
git clone https://github.com/OpticalDiffractionTomography/ODT_Reconstruction.git
cd ODT_Reconstruction
```

### 3. Set your email and check the paths

Open the config file:

```bash
nano config.sh
```

Set your email and check that the two mount paths are correct for your cluster:

```bash
EMAIL="you@institute.de"                    # ← set your email here
DATA_MOUNT="/mnt/guck_division2"            # ← where your raw data is
RESULTS_MOUNT="/mnt/ZPE_cluster_results"    # ← where results are written
```

Save and close: press `Ctrl+O`, then `Enter`, then `Ctrl+X`.

Everything else in the file has working defaults — leave it as is.

### 4. Install the command

```bash
bash install.sh
source ~/.bashrc
```

Check that it works:

```bash
tomo_process --help
```

That's it — go to [Run jobs](#run-jobs-2-steps).

---

## Options

| Option | What it does | Default |
|---|---|---|
| `--nm 1.340` | Refractive index of the culture medium | `1.337` |
| `--email you@lab.de` | Notification email for this run only | `EMAIL` from `config.sh` |
| `--max-jobs 3` | Parallel SLURM jobs | `5` |
| `--mins-per-sample 10` | Time estimate per sample (used to size jobs) | `15` |

Run `tomo_process --help` for the full list.

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

## If your SSH disconnects or the login node reboots

A dropped SSH connection does **not** interrupt processing. If the login node itself reboots, no work is lost — just reconnect and resume:

```bash
ssh <your_username>@zpe.intranet.mpl.mpg.de
tomo_process --list                          # find your run_id
tomo_process --resume tomo_20260729_143201   # continue from where it stopped
```

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

## Advanced

<details>
<summary>How job chunking works</summary>

Jobs are auto-sized from the walltime and the per-sample time estimate:

```
max_samples_per_job = floor(SLURM_TIME / MINS_PER_SAMPLE)
                    = floor(20h × 60min / 15min)
                    = 80 samples per job
```

Files from different subdirectories are never mixed in one job — each subdirectory is processed independently, split across multiple jobs if it contains more samples than the per-job limit.

</details>

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

If a pasted `--path` cannot be found under `DATA_MOUNT`, the script strips the drive letter and drops up to 3 leading folders while searching; if it still fails, it lists the paths it tried and how to fix yours.

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

### Stay connected during setup (tmux)

If your SSH session drops during the one-time setup, interactive commands die with it. Start a persistent shell first:

```bash
tmux new -s setup
```

Once `tomo_process` is running, this is not needed — it detaches automatically.

</details>
