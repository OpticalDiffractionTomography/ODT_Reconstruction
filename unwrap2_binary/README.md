# Compiling unwrap2

`unwrap2` is a 2D phase-unwrapping MEX function (Goldstein's branch-cut
algorithm) used by [src/field_Retrieval.m](../src/field_Retrieval.m). The
prebuilt binaries checked in here (`originalUnwrap2/unwrap2/unwrap2.mexw32`,
`originalUnwrap2/unwrap/unwrap2.mexw32`) are old Windows 32-bit builds and
won't load on the HPC cluster (Linux, MATLAB R2023a+). `compile_unwrap2.sh`
rebuilds `unwrap2` from source for the current platform.

## Source layout

```
unwrap2_binary/originalUnwrap2/
  unwrap/unwrap2.c     # MEX entry point (mexFunction), calls into unwrap2/*.c
  unwrap2/*.c, *.h     # supporting sources (gold.c, brcut.c, util.c, ...)
```

## Compiling on the cluster

From the repo root:

```bash
sbatch unwrap2_binary/compile_unwrap2.sh
```

This submits a SLURM job that loads the `matlab/R2023a` module and runs
`mex` against `unwrap/unwrap2.c` plus all supporting sources in `unwrap2/`.
Adjust the `module load matlab/...` line in the script if a different
MATLAB version is required.

Build logs go to `unwrap2_binary/compile_unwrap2_<jobid>.log`. On success
the compiled MEX file is written to `unwrap2_binary/build/unwrap2.mexa64`
(the `mexext` for Linux x86-64).

**Important:** submit the job from the repo root so `$SLURM_SUBMIT_DIR`
resolves to the correct location — the script builds its paths relative to
`$SLURM_SUBMIT_DIR/unwrap2_binary`.

## Compiling interactively (no SLURM)

If MATLAB is already available on `PATH`:

```bash
bash unwrap2_binary/compile_unwrap2.sh
```

## Using the compiled binary

Copy or symlink the resulting `unwrap2.mexa64` into a directory on
MATLAB's path (e.g. into `src/`, or anywhere covered by
`addpath(genpath('./src'))` in `main.sh`) so `unwrap2(...)` calls in
`field_Retrieval.m` resolve to the new build instead of the stale
`.mexw32` binaries.
