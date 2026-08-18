#!/bin/bash
#SBATCH --job-name=compile_unwrap2
#SBATCH --output=compile_unwrap2_%j.log
#SBATCH --error=compile_unwrap2_%j.log
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00

# Compiles the unwrap2 MEX function from source using MATLAB's mex compiler.
#
# Usage:
#   sbatch compile_unwrap2.sh
# or, interactively (no SLURM):
#   bash compile_unwrap2.sh
#
# Source layout expected (relative to this script):
#   originalUnwrap2/unwrap/unwrap2.c   <- MEX entry point (mexFunction)
#   originalUnwrap2/unwrap2/*.c        <- supporting sources (gold.c, brcut.c, etc.)
#   originalUnwrap2/unwrap2/*.h        <- headers included by unwrap2.c

set -euo pipefail

if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    SCRIPT_DIR="${SLURM_SUBMIT_DIR}/unwrap2_binary"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
SRC_MAIN_DIR="${SCRIPT_DIR}/originalUnwrap2/unwrap"
SRC_LIB_DIR="${SCRIPT_DIR}/originalUnwrap2/unwrap2"
OUT_DIR="${SCRIPT_DIR}/build"

mkdir -p "${OUT_DIR}"

module load matlab/R2023a

BUILD_M="${OUT_DIR}/build_unwrap2.m"
cat > "${BUILD_M}" <<EOF
try
    mex('-outdir', '${OUT_DIR}', ...
        '-output', 'unwrap2', ...
        '-I${SRC_LIB_DIR}', ...
        '${SRC_MAIN_DIR}/unwrap2.c', ...
        '${SRC_LIB_DIR}/gold.c', ...
        '${SRC_LIB_DIR}/brcut.c', ...
        '${SRC_LIB_DIR}/dipole.c', ...
        '${SRC_LIB_DIR}/extract.c', ...
        '${SRC_LIB_DIR}/grad.c', ...
        '${SRC_LIB_DIR}/list.c', ...
        '${SRC_LIB_DIR}/maskfat.c', ...
        '${SRC_LIB_DIR}/path.c', ...
        '${SRC_LIB_DIR}/residues.c', ...
        '${SRC_LIB_DIR}/util.c');
    fprintf('unwrap2 MEX build succeeded.\n');
catch err
    fprintf(2, 'unwrap2 MEX build failed: %s\n', err.getReport());
    exit(1);
end
exit(0);
EOF

matlab -batch "run('${BUILD_M}')"

echo "Build complete. Output in: ${OUT_DIR}/"
ls -la "${OUT_DIR}"/unwrap2.*
