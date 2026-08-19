# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a MATLAB codebase for **Optical Diffraction Tomography (ODT)** — reconstructing 3D refractive index maps of biological samples (e.g., MDCK cells) from interferometric holographic measurements. The pipeline converts raw holographic tomograms into quantitative 3D refractive index volumes.

## Running the Pipeline

This code runs on an HPC cluster (`/beegfs/home/ralajan/matlab/`). There is no build system; scripts are run directly in MATLAB.

**Entry points (in `src/`):**
- [src/field_Retrieval.m](src/field_Retrieval.m) — Stage 1 only
- [src/tomogram_Reconstruction.m](src/tomogram_Reconstruction.m) — Stage 2 only

**SLURM submission (the only way to run):**
```bash
sbatch main.sh --data_dir /beegfs/home/ralajan/matlab/<experiment_date>
```
`main.sh` sets `addpath(genpath('./src'))` and passes `spath` to MATLAB before running each stage script. The stage scripts must **not** contain `clear all`, `addpath`, or a hardcoded `spath` — those are injected by `main.sh`.

`main.sh` requests 1 RTX GPU, 8 CPUs, 64 GB RAM, 24 h walltime and runs both stages sequentially.

**Legacy monolithic script** (all 4 stages, still the authoritative reference for stitching/visualization):
```
field_retrieval_Tomogram_reconstruction.asv
```

GPU (CUDA) and MATLAB R2026a are required. Dependencies:
- `mk_ellipse` — creates circular/elliptical binary masks ([src/mk_ellipse.m](src/mk_ellipse.m))
- `unwrap2` — MATLAB's built-in 2D phase unwrapper (Goldstein branch-cut). It can crash on residue-heavy noisy frames; `field_Retrieval.m` catches this per frame and NaN-fills, and such frames are excluded downstream.
- MATLAB toolboxes required: Image Processing, Parallel Computing (for `gpuArray`/`gather`), Signal Processing (for `xcorr2`)

## Pipeline Architecture

### Stage 1: Field Retrieval ([src/field_Retrieval.m](src/field_Retrieval.m))

Iterates over `batch*/` directories, each containing `bg*_Tomog.mat` (background) and `sample*_Tomog.mat` (sample) files.

- Detects off-axis carrier frequency via FFT peak (`mi`, `mj`); builds circular mask with `mk_ellipse`
- Demodulates each frame: isolates +1 diffraction order, divides by background field
- Per-frame phase pipeline: `angle(Fimg)` → `unwrap2` → `PhiShift` → `phaseCompensation(deg=1)` → `phaseCompensation(deg=2, masked)`
  - Second `phaseCompensation` call uses a cell-exclusion mask from top-hat morphological filtering (`imtophat` + `strel`)
- Saves `retPhase`, `retAmplitude`, `NA`, `lambda`, `res`, `ZP`, `f_dx`, `f_dy` as `Field_*_rev2.mat` in `field_retrieval/`

### Stage 2: Tomogram Reconstruction ([src/tomogram_Reconstruction.m](src/tomogram_Reconstruction.m))

- Loads field retrieval outputs; detects outlier frames into `excludeFrame` (mean phase > 1.5, frame-to-frame diff > 0.1, NaN frames)
- Builds `TomoParam` struct with all optical/geometric parameters (see Key Parameters below)
- `ODTReconstruction` → Ewald sphere mapping → `ORytov` (3D Fourier volume) + initial `Reconimg`
- `ODTIteration` (100 iterations on GPU) → enforces `n >= n_m` physical constraint → final `Reconimg`
- Crops to unpadded region; saves `Reconimg`, `res3`, `res4`, `lambda`, `excludeFrame` to `Tomogram_*.mat`

### Stage 3: Tile Stitching (in `.asv`, lines 237–556)

- v1 (lines 237–348): stitches pre-reconstructed full-volume tomograms
- v2/`_v3` (lines 371–556): reconstructs per-tile first, then stitches
- Bidirectional scan pattern: odd rows left-to-right, even rows right-to-left
- Registration via `xcorr2` on image gradients to find sub-tile offsets

### Stage 4: Visualization (in `.asv`, lines 671–785)

- Orthogonal slices (XZ, YZ, XY) and maximum intensity projections

## Key Functions

| File | Role |
|------|------|
| [src/ODTReconstruction.m](src/ODTReconstruction.m) | Ewald sphere mapping: fills 3D Fourier volume `[ZP2 × ZP2 × ZP3]` from 2D Rytov scattered fields; uses `gpuArray` throughout |
| [src/ODTIteration.m](src/ODTIteration.m) | Iterative projection: enforces `n >= n_m` constraint, replaces measured Fourier voxels from `ORytov` each iteration |
| [src/PhiShift.m](src/PhiShift.m) | Removes linear phase tilt by fitting a plane to 4-pixel-wide border strips |
| [src/phaseCompensation.m](src/phaseCompensation.m) | Least-squares polynomial fit (degree `n`) for phase aberration compensation; optional `mask` arg excludes cell regions from the fit |

## Key Physical Parameters

These appear throughout and must stay consistent. The `TomoParam` struct passed to `ODTReconstruction`/`ODTIteration` bundles all of them:

| Variable | Meaning |
|----------|---------|
| `lambda` | Illumination wavelength (µm) |
| `NA` | Objective numerical aperture |
| `res` | Camera pixel size (µm) |
| `n_m = 1.337` | Refractive index of culture medium |
| `n_s = 1.377` | Expected sample RI (cells) |
| `ZP` | Hologram FFT size (≈ image size; typically `round(1.2 * xx/2) * 2`) |
| `ZP2 = 512` | Lateral size of 3D Fourier volume |
| `ZP3 = 256` | Axial size of 3D Fourier volume |
| `res3`, `res4` | Output lateral/axial voxel sizes |
| `f_dx2`, `f_dy2` | Carrier frequency offset per frame (pixels) |
| `k0_x/y/z` | Illumination wave-vector components per frame |
| `frameList` | Non-excluded frame indices |

## Data Layout

```
/beegfs/home/ralajan/matlab/<experiment_date>/
  batch01_*/
    bg001_Tomog.mat                       # background hologram stack (variable: tomogMap)
    sample001_Tomog.mat                   # sample hologram stack (variable: tomogMap)
    field_retrieval/
      Field_sample001_*_rev2.mat          # retPhase, retAmplitude, NA, lambda, res, ZP, f_dx, f_dy
      Tomogram_Field_*_rev2.mat           # Reconimg, res3, res4, lambda, excludeFrame
      checkpoint_<baseName>_iter_<N>.mat  # progress markers (not auto-resumable)
      stitched_Tomogram.mat               # final stitched volume
```

## Checkpoint System

Field retrieval writes checkpoint markers every 5 frames to `outDir/checkpoint_<baseName>_iter_<N>.mat`. These are informational only — re-running always starts from frame 1 and overwrites previous outputs.

## Refactoring State

All source files live in `src/`. `main.sh` is the sole entry point — it calls `addpath(genpath('./src'))` and injects `spath` before running each stage script. Stage scripts must not set `clear all`, `addpath`, or `spath` themselves.

Stages 3–4 (stitching and visualization) remain only in the `.asv` file and have not yet been split into `src/`. The function `ODTReconstruction_scale` referenced in `test.m` for tiled reconstruction does not exist in this repo.
