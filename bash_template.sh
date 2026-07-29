#!/bin/bash
# bash_template.sh — SLURM job template for one processing chunk.
# The orchestrator fills in __PLACEHOLDERS__ and writes a real job script;
# this file is never submitted directly.
#
# Placeholders replaced at submission time:
#   __JOB_NAME__      unique name for this chunk  (e.g. tomo_20260728_chunk003)
#   __LOG_DIR__       directory for out/err logs
#   __EMAIL__         user email for SLURM notifications
#   __DATA_DIR__      chunk staging dir on scratch  (/beegfs/.../scratch/<run>/<chunk>/data)
#   __CHUNK_DONE__    touch-file path; orchestrator watches for this
#   __CHUNK_FAIL__    touch-file path written on failure
#   __PARTITION__     SLURM partition
#   __GRES__          GPU resource string
#   __CPUS__          CPUs per task
#   __MEM__           RAM
#   __TIME__          walltime
#
# NOTE: This job never touches /mnt — all /mnt access happens on the login node
# via the tomo_process orchestrator (rsync in, rsync out).

#SBATCH -J __JOB_NAME__
#SBATCH -o __LOG_DIR__/__JOB_NAME__.out
#SBATCH -e __LOG_DIR__/__JOB_NAME__.err
#SBATCH -D .
#SBATCH --partition=__PARTITION__
#SBATCH --gres=__GRES__
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=__CPUS__
#SBATCH --time=__TIME__
#SBATCH --mem=__MEM__
#SBATCH --mail-type=ALL
#SBATCH --mail-user=__EMAIL__

set -euo pipefail

DATA_DIR="__DATA_DIR__"
CHUNK_DONE="__CHUNK_DONE__"
CHUNK_FAIL="__CHUNK_FAIL__"
REPO_DIR="__REPO_DIR__"
SRC_DIR="${REPO_DIR}/src"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

cleanup_on_fail() {
    log "FAILED — writing fail marker"
    touch "${CHUNK_FAIL}"
    exit 1
}
trap cleanup_on_fail ERR

module load matlab/R2026a
module load cuda/11.6.0

export QT_QPA_PLATFORM=offscreen
unset DISPLAY

run_stage() {
    local stage="$1"
    log "Starting: ${stage}"
    matlab -batch "\
        addpath(genpath('${SRC_DIR}')); \
        spath='${DATA_DIR}'; \
        run('${SRC_DIR}/${stage}')" \
        && log "Done: ${stage}" \
        || { log "FAILED: ${stage}"; return 1; }
}

log "=== Job start: __JOB_NAME__ ==="
log "Data dir : ${DATA_DIR}"

run_stage field_Retrieval.m
run_stage tomogram_Reconstruction.m

# Signal success — the login-node orchestrator will rsync results from
# ${DATA_DIR}/field_retrieval/ back to /mnt/ZPE_cluster_results.
touch "${CHUNK_DONE}"
log "=== Job complete: __JOB_NAME__ ==="
