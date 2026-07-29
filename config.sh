#!/bin/bash
# config.sh — cluster-wide constants for tomo_process
# Edit these values once after cloning; all other scripts source this file.

# ── Cluster paths ────────────────────────────────────────────────────────────

# Root of the read-only network mount on the cluster
ZPE_RESULTS_MOUNT="/mnt/ZPE_results"

# Scratch space for staging data during processing
SCRATCH_ROOT="/beegfs/home/ralajan/scratch/tomo_process"

# Absolute path to this repository on the cluster
REPO_DIR="/beegfs/home/ralajan/matlab/field_tomogram_reconstruction"

# ── Job limits ───────────────────────────────────────────────────────────────

# Maximum number of processing jobs running in parallel
MAX_PARALLEL_JOBS=5

# Target number of sample MAT files per job chunk (auto-sized down if fewer remain)
FILES_PER_CHUNK=10

# How often (seconds) the orchestrator polls for finished jobs
POLL_INTERVAL=60

# ── SLURM resources per processing job ───────────────────────────────────────

SLURM_PARTITION="gpu"
SLURM_GRES="gpu:rtx:1"
SLURM_CPUS=16
SLURM_MEM="128G"
SLURM_TIME="24:00:00"

# ── SLURM resources for the orchestrator job itself ───────────────────────────

ORCH_PARTITION="short"       # low-resource partition; adjust to your cluster's name
ORCH_TIME="72:00:00"         # max orchestrator runtime (must exceed total dataset time)
ORCH_MEM="4G"
