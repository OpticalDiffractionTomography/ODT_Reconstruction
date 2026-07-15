#!/bin/bash
#SBATCH -J Matlab_GPU_Job             # Job Name
#SBATCH -o ./out.log                  # Standard output log
#SBATCH -e ./err.log                  # Error log
#SBATCH -D .                          # Working Directory
#SBATCH --partition=gpu               # GPU Partition
#SBATCH --gres=gpu:rtx:1              # Request 1 RTX GPU
#SBATCH --nodes=1                     # Run on a single node
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=8             # CPU cores = parpool workers (field retrieval)
#SBATCH --time=24:00:00               # Max runtime (HH:MM:SS)
#SBATCH --mem=64G                     # RAM: parallel workers each load a tomogMap
#SBATCH --mail-type=all               # Email notifications for start/end/fail
#SBATCH --mail-user=raghava.alajangi@mpzpm.mpg.de

# ---------------------------------------------------------------------------
# Usage:
#   sbatch main.sh --data_dir /beegfs/home/ralajan/matlab/20260203_Cecile_MDCK
#
# Arguments:
#   --data_dir <path>   Experiment folder containing batch*/ subdirectories
# ---------------------------------------------------------------------------

DATA_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --data_dir) DATA_DIR="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if [[ -z "$DATA_DIR" ]]; then
    echo "ERROR: --data_dir is required"; exit 1
fi

module load matlab/R2026a
module load cuda/11.6.0

# Headless rendering: compute nodes have no display; use Qt offscreen plugin
export QT_QPA_PLATFORM=offscreen

REPO_DIR="${SLURM_SUBMIT_DIR}"
SRC_DIR="${REPO_DIR}/src"

run_stage() {
    local stage="$1"
    echo "[$(date '+%H:%M:%S')] Starting: ${stage}"
    # SLURM_CPUS_PER_TASK is already in the environment; MATLAB inherits it.
    matlab -batch "\
        addpath(genpath('${SRC_DIR}')); \
        spath='${DATA_DIR}'; \
        run('${SRC_DIR}/${stage}')" \
        && echo "[$(date '+%H:%M:%S')] Done:     ${stage}" \
        || { echo "[$(date '+%H:%M:%S')] FAILED:   ${stage}"; exit 1; }
}

run_stage field_Retrieval.m
run_stage tomogram_Reconstruction.m
