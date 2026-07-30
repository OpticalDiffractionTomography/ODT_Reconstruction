# field_tomogram_reconstruction

This is a MATLAB codebase for **Optical Diffraction Tomography (ODT)** — reconstructing 3D refractive index maps of biological samples (e.g., MDCK cells) from interferometric holographic measurements. The pipeline converts raw holographic tomograms into quantitative 3D refractive index volumes.

## Running the Pipeline

This code runs on a HPC cluster (`/beegfs/home/ralajan/matlab/`). There is no build system; scripts are run directly in MATLAB. The main entry point is [field_retrieval_Tomogram_reconstruction.m](field_retrieval_Tomogram_reconstruction.m), which contains the full pipeline as sequential `%%` sections.

To run on the cluster:
```matlab
% Launch MATLAB and run:
run('field_retrieval_Tomogram_reconstruction.m')
```

GPU (CUDA) is required — `gpuArray` is used extensively in `ODTReconstruction` and `ODTIteration`. The `unwrap2` function in [unwrap/](unwrap/) is a compiled MEX binary (`.mexw32`) for 2D phase unwrapping.

## Pipeline Architecture

The pipeline is structured as sequential `%%` sections in [field_retrieval_Tomogram_reconstruction.m](field_retrieval_Tomogram_reconstruction.m):

### Stage 1: Field Retrieval (lines 1–141)
Iterates over `batch*/` directories, each containing `bg*_Tomog.mat` (background) and `sample*_Tomog.mat` (sample) files.

- Finds the off-axis hologram carrier frequency via FFT peak detection
- Applies a circular mask (`mk_ellipse`) to isolate the +1 diffraction order
- Demodulates the complex field by shifting in Fourier space
- Divides sample field by background field to extract relative phase/amplitude
- Calls `PhiShift` → `unwrap2` → `phaseCompensation` to flatten and unwrap phase
- Saves `retPhase` and `retAmplitude` as `Field_*_rev2.mat` in `field_retrieval/` subdirectory

### Stage 2: Tomogram Reconstruction (lines 143–235)
- Loads field retrieval results, filters corrupted/outlier frames into `excludeFrame`
- Calls `ODTReconstruction` to map Rytov fields into 3D Fourier space (the Ewald sphere mapping)
- Calls `ODTIteration` for iterative constraint enforcement (non-negativity of refractive index)
- Saves `Reconimg` (3D refractive index volume) as `Tomogram_*.mat`

### Stage 3: Tile Stitching (lines 237–348 and 371–556)
- Two versions: v1 stitches full-volume tomograms; v2 (`_v3`) reconstructs per-tile then stitches
- Bidirectional scan pattern: odd rows go left-to-right, even rows right-to-left
- Registration uses `xcorr2` on image gradients to find sub-tile offsets

### Stage 4: Visualization (lines 671–785)
- Displays orthogonal slices (XZ, YZ, XY) and maximum intensity projections

## Key Functions

| File | Role |
|------|------|
| [ODTReconstruction.m](ODTReconstruction.m) | Ewald sphere mapping: fills 3D Fourier space from 2D holographic projections using Rytov approximation |
| [ODTIteration.m](ODTIteration.m) | Iterative refinement: enforces `n >= n_m` (refractive index ≥ medium) as a physical constraint |
| [PhiShift.m](PhiShift.m) | Removes linear phase tilt by fitting and subtracting a plane from border pixels |
| [phaseCompensation.m](phaseCompensation.m) | Least-squares polynomial fit (degree `n`) to compensate residual phase aberrations; supports masked fitting to exclude cells |

## Key Physical Parameters

These appear as variables throughout the pipeline and must be consistent:
- `lambda` — illumination wavelength (in µm)
- `NA` — numerical aperture of the objective
- `res` — pixel size in µm
- `n_m = 1.337` — refractive index of culture medium
- `n_s` — expected refractive index of sample (cells ~1.37–1.38)
- `ZP` — zero-padded FFT size (typically `round(1.2 * xx/2) * 2`)

## Data Layout

```
/beegfs/home/ralajan/matlab/<experiment_date>/
  batch01_*/
    bg001_Tomog.mat          # background hologram stack (variable: tomogMap)
    sample001_Tomog.mat      # sample hologram stack (variable: tomogMap)
    field_retrieval/
      Field_sample001_*_rev2.mat   # retPhase, retAmplitude, NA, lambda, res, ZP
      Tomogram_Field_*_rev2.mat    # Reconimg, res3, res4, lambda
      stitched_Tomogram.mat        # final stitched volume
```

## Checkpoint System

Field retrieval saves checkpoint markers every 5 frames to `outDir/checkpoint_<baseName>_iter_<N>.mat`. On restart, these are not automatically used for resuming — re-running overwrites from the beginning.
