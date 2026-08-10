#!/bin/bash
# config.sh — cluster-wide constants for tomo_process
# Edit these values once after cloning; all other scripts source this file.

# ── User settings ────────────────────────────────────────────────────────────

# Your email for notifications (processing started / finished / failed).
# Set this once; you can still override it per run with --email.
EMAIL=""

# SMTP relay for sending notification emails (host[:port]).
# The login node has no local mail service, so without this, notifications
# fall back to bare SLURM job emails (subject line only, no message body).
# Ask IT for your institute's internal SMTP relay, e.g. "mailhost.mpzpm.mpg.de:25".
SMTP_RELAY=""

# ── Cluster paths ────────────────────────────────────────────────────────────

# Root of the read-only network mount where input data lives.
# Set to the path where your raw tomogram directories are mounted on the cluster.
DATA_MOUNT="/mnt/guck_division2"                    # read-only: input data

# Writable mount where processed results will be written.
# Set to a path accessible from the login node (compute nodes do not need it).
RESULTS_MOUNT="/mnt/ZPE_cluster_results"            # writable: results output

# Scratch space on fast cluster storage for staging data during processing.
# Must be accessible from both the login node and compute nodes.
SCRATCH_ROOT="${HOME}/scratch/tomo_process"

# Absolute path to this repository on the cluster (auto-detected; override if needed)
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || echo "")"

# ── Job limits ───────────────────────────────────────────────────────────────

# Maximum number of processing jobs running in parallel
MAX_PARALLEL_JOBS=5

# How often (seconds) the orchestrator polls for finished jobs
POLL_INTERVAL=60

# ── SLURM resources per processing job ───────────────────────────────────────

SLURM_PARTITION="gpu"
SLURM_GRES="gpu:rtx:1"
SLURM_CPUS=16
SLURM_MEM="128G"

# Maximum walltime per job (HH:MM:SS). Chunk size is derived from this:
#   max_samples_per_job = floor(walltime_hours * 60 / MINS_PER_SAMPLE)
# Leave a safety margin — set to 20h even though the queue allows 24h.
SLURM_TIME="20:00:00"

# Estimated processing time per sample file in minutes (field retrieval + reconstruction).
# Used to auto-size chunks so each job stays within SLURM_TIME.
# Be conservative — better to under-fill than to hit the walltime limit.
MINS_PER_SAMPLE=15

# ── SLURM resources for the orchestrator job itself ───────────────────────────

ORCH_PARTITION="cpu"         # low-resource partition; adjust to your cluster's name
ORCH_TIME="72:00:00"         # max orchestrator runtime (must exceed total dataset time)
ORCH_MEM="4G"
