#!/bin/bash
# orchestrator.sh — runs as a SLURM job; manages the full copy→process→collect loop.
# Submitted by `tomo_process`; never run directly.
#
# Arguments (passed via --export or environment, set by tomo_process):
#   TOMO_INPUT_PATH   path under ZPE_RESULTS_MOUNT to read sample files from
#   TOMO_EMAIL        user email
#   TOMO_RUN_ID       unique run identifier (set by tomo_process)

#SBATCH -J tomo_orchestrator
#SBATCH -o /tmp/tomo_orch_%j.out
#SBATCH -e /tmp/tomo_orch_%j.err
#SBATCH -D .
#SBATCH --partition=ORCH_PARTITION_PLACEHOLDER
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=ORCH_TIME_PLACEHOLDER
#SBATCH --mem=ORCH_MEM_PLACEHOLDER
#SBATCH --mail-type=NONE

set -euo pipefail

# ── Source config ─────────────────────────────────────────────────────────────
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "${REPO_DIR}/config.sh"

# ── Validate environment ──────────────────────────────────────────────────────
: "${TOMO_INPUT_PATH:?TOMO_INPUT_PATH not set}"
: "${TOMO_EMAIL:?TOMO_EMAIL not set}"
: "${TOMO_RUN_ID:?TOMO_RUN_ID not set}"

# ── Directories ───────────────────────────────────────────────────────────────
RUN_SCRATCH="${SCRATCH_ROOT}/${TOMO_RUN_ID}"
STATE_DIR="${RUN_SCRATCH}/state"
LOG_DIR="${RUN_SCRATCH}/logs"
JOBS_DIR="${RUN_SCRATCH}/jobs"       # generated job scripts

SRC_MOUNT="${ZPE_RESULTS_MOUNT}/${TOMO_INPUT_PATH}"
DST_RESULTS="${SRC_MOUNT}/_results"

mkdir -p "${STATE_DIR}" "${LOG_DIR}" "${JOBS_DIR}"

ORCH_LOG="${LOG_DIR}/orchestrator.log"

# ── Logging ───────────────────────────────────────────────────────────────────
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${ORCH_LOG}"; }

send_email() {
    local subject="$1"
    local body="$2"
    echo "${body}" | mail -s "[tomo_process] ${subject}" "${TOMO_EMAIL}" 2>/dev/null || true
}

# ── Discover all sample MAT files ─────────────────────────────────────────────
log "=== tomo_process orchestrator start ==="
log "Run ID    : ${TOMO_RUN_ID}"
log "Input     : ${SRC_MOUNT}"
log "Results   : ${DST_RESULTS}"
log "Scratch   : ${RUN_SCRATCH}"

if [[ ! -d "${SRC_MOUNT}" ]]; then
    log "ERROR: Source path not found: ${SRC_MOUNT}"
    send_email "FAILED — ${TOMO_RUN_ID}" \
        "Source path not found: ${SRC_MOUNT}\nCheck that /mnt/ZPE_results is mounted."
    exit 1
fi

# Collect all sample*_Tomog.mat files, sorted
PENDING_FILE="${STATE_DIR}/pending_files.txt"
CHUNKS_FILE="${STATE_DIR}/chunks.txt"   # format: chunk_id<TAB>file1,file2,...
ACTIVE_FILE="${STATE_DIR}/active_jobs.txt"   # format: chunk_id<TAB>slurm_job_id
DONE_FILE="${STATE_DIR}/done_chunks.txt"
FAIL_FILE="${STATE_DIR}/failed_chunks.txt"

if [[ ! -f "${PENDING_FILE}" ]]; then
    find "${SRC_MOUNT}" -name "sample*_Tomog.mat" | sort > "${PENDING_FILE}"
    TOTAL=$(wc -l < "${PENDING_FILE}")
    log "Found ${TOTAL} sample MAT files"

    if [[ "${TOTAL}" -eq 0 ]]; then
        log "ERROR: No sample*_Tomog.mat files found under ${SRC_MOUNT}"
        send_email "FAILED — ${TOMO_RUN_ID}" \
            "No sample*_Tomog.mat files found under ${SRC_MOUNT}"
        exit 1
    fi

    # Auto-size: aim for FILES_PER_CHUNK files per job; split into chunks
    chunk_idx=0
    chunk_files=()
    while IFS= read -r fpath; do
        chunk_files+=("${fpath}")
        if [[ ${#chunk_files[@]} -ge ${FILES_PER_CHUNK} ]]; then
            chunk_id=$(printf "chunk%04d" ${chunk_idx})
            echo "${chunk_id}"$'\t'"$(IFS=','; echo "${chunk_files[*]}")" >> "${CHUNKS_FILE}"
            chunk_idx=$(( chunk_idx + 1 ))
            chunk_files=()
        fi
    done < "${PENDING_FILE}"
    # Remaining files (last partial chunk)
    if [[ ${#chunk_files[@]} -gt 0 ]]; then
        chunk_id=$(printf "chunk%04d" ${chunk_idx})
        echo "${chunk_id}"$'\t'"$(IFS=','; echo "${chunk_files[*]}")" >> "${CHUNKS_FILE}"
    fi

    NCHUNKS=$(wc -l < "${CHUNKS_FILE}")
    log "Split into ${NCHUNKS} chunks of up to ${FILES_PER_CHUNK} files each"
fi

touch "${ACTIVE_FILE}" "${DONE_FILE}" "${FAIL_FILE}"

send_email "Started — ${TOMO_RUN_ID}" \
    "Processing started.\nInput : ${SRC_MOUNT}\nChunks: $(wc -l < "${CHUNKS_FILE}")\nMax parallel jobs: ${MAX_PARALLEL_JOBS}"

# ── Helper: count active jobs still in SLURM queue ───────────────────────────
count_active_jobs() {
    local count=0
    while IFS=$'\t' read -r cid jid; do
        if squeue -j "${jid}" -h &>/dev/null; then
            count=$(( count + 1 ))
        fi
    done < "${ACTIVE_FILE}"
    echo "${count}"
}

# ── Helper: collect finished jobs ─────────────────────────────────────────────
collect_finished() {
    local new_active=""
    while IFS=$'\t' read -r cid jid; do
        local done_marker="${RUN_SCRATCH}/${cid}/DONE"
        local fail_marker="${RUN_SCRATCH}/${cid}/FAIL"

        if [[ -f "${done_marker}" ]]; then
            log "Chunk ${cid} finished successfully (job ${jid})"
            echo "${cid}" >> "${DONE_FILE}"
            # Remove scratch data for this chunk
            rm -rf "${RUN_SCRATCH:?}/${cid}"
            log "Cleaned scratch for ${cid}"

        elif [[ -f "${fail_marker}" ]]; then
            log "WARNING: Chunk ${cid} FAILED (job ${jid}) — kept in scratch for inspection"
            echo "${cid}" >> "${FAIL_FILE}"
            send_email "Chunk FAILED — ${TOMO_RUN_ID}" \
                "Chunk ${cid} (SLURM job ${jid}) failed.\nLogs: ${LOG_DIR}/${cid}_*.err\nScratch data kept at: ${RUN_SCRATCH}/${cid}"

        else
            # Still running or pending
            new_active+="${cid}"$'\t'"${jid}"$'\n'
        fi
    done < "${ACTIVE_FILE}"
    printf "%s" "${new_active}" > "${ACTIVE_FILE}"
}

# ── Helper: submit one chunk ──────────────────────────────────────────────────
submit_chunk() {
    local chunk_id="$1"
    local files_csv="$2"   # comma-separated list of full paths

    local chunk_scratch="${RUN_SCRATCH}/${chunk_id}"
    local chunk_data="${chunk_scratch}/data"
    local chunk_results="${chunk_scratch}/results"
    mkdir -p "${chunk_data}/field_retrieval" "${chunk_results}"

    log "Copying files for ${chunk_id} ..."
    IFS=',' read -ra flist <<< "${files_csv}"

    # bg01_Tomog.mat is a single common background shared across the whole dataset.
    # It lives in the same directory as the sample files; copy it once per chunk.
    local dataset_dir
    dataset_dir=$(dirname "${flist[0]}")
    local bg_src="${dataset_dir}/bg01_Tomog.mat"
    if [[ -f "${bg_src}" ]]; then
        rsync -a "${bg_src}" "${chunk_data}/" \
            || { log "rsync failed for background ${bg_src}"; return 1; }
        log "Copied background: bg01_Tomog.mat"
    else
        log "ERROR: background file not found: ${bg_src}"
        return 1
    fi

    # Copy the sample files for this chunk
    for fpath in "${flist[@]}"; do
        rsync -a "${fpath}" "${chunk_data}/" \
            || { log "rsync copy failed for ${fpath}"; return 1; }
    done

    # Result destination mirrors the relative path inside ZPE_results
    local rel_path="${TOMO_INPUT_PATH}"
    local result_dest="${DST_RESULTS}/${chunk_id}"

    # Generate job script from template
    local job_script="${JOBS_DIR}/${chunk_id}.sh"
    local job_name="tomo_${TOMO_RUN_ID}_${chunk_id}"

    sed \
        -e "s|__JOB_NAME__|${job_name}|g" \
        -e "s|__LOG_DIR__|${LOG_DIR}|g" \
        -e "s|__EMAIL__|${TOMO_EMAIL}|g" \
        -e "s|__DATA_DIR__|${chunk_data}|g" \
        -e "s|__RESULT_DIR__|${result_dest}|g" \
        -e "s|__CHUNK_DONE__|${chunk_scratch}/DONE|g" \
        -e "s|__CHUNK_FAIL__|${chunk_scratch}/FAIL|g" \
        -e "s|__REPO_DIR__|${REPO_DIR}|g" \
        -e "s|__PARTITION__|${SLURM_PARTITION}|g" \
        -e "s|__GRES__|${SLURM_GRES}|g" \
        -e "s|__CPUS__|${SLURM_CPUS}|g" \
        -e "s|__MEM__|${SLURM_MEM}|g" \
        -e "s|__TIME__|${SLURM_TIME}|g" \
        "${REPO_DIR}/bash_template.sh" > "${job_script}"

    chmod +x "${job_script}"

    local jid
    jid=$(sbatch --parsable "${job_script}")
    log "Submitted ${chunk_id} → SLURM job ${jid}"
    echo "${chunk_id}"$'\t'"${jid}" >> "${ACTIVE_FILE}"
}

# ── Main loop ─────────────────────────────────────────────────────────────────
# Read all chunks into an array; track which ones have been submitted
declare -A submitted_chunks
while IFS=$'\t' read -r cid _; do
    submitted_chunks["${cid}"]=0
done < "${DONE_FILE}"
# Also mark previously active as submitted (resume case)
while IFS=$'\t' read -r cid _; do
    submitted_chunks["${cid}"]=1
done < "${ACTIVE_FILE}"

while true; do
    collect_finished

    # Check if everything is done
    TOTAL_CHUNKS=$(wc -l < "${CHUNKS_FILE}")
    DONE_CHUNKS=$(wc -l < "${DONE_FILE}")
    FAIL_CHUNKS=$(wc -l < "${FAIL_FILE}")
    FINISHED=$(( DONE_CHUNKS + FAIL_CHUNKS ))

    log "Progress: ${DONE_CHUNKS} done, ${FAIL_CHUNKS} failed, ${FINISHED}/${TOTAL_CHUNKS} finished"

    if [[ "${FINISHED}" -ge "${TOTAL_CHUNKS}" ]]; then
        break
    fi

    # Submit new chunks up to MAX_PARALLEL_JOBS
    ACTIVE=$(count_active_jobs)
    SLOTS=$(( MAX_PARALLEL_JOBS - ACTIVE ))

    if [[ "${SLOTS}" -gt 0 ]]; then
        while IFS=$'\t' read -r cid files_csv; do
            [[ "${SLOTS}" -le 0 ]] && break

            # Skip already submitted or failed chunks
            if [[ -n "${submitted_chunks[${cid}]+x}" ]]; then
                continue
            fi
            if grep -qx "${cid}" "${DONE_FILE}" 2>/dev/null; then
                submitted_chunks["${cid}"]=1
                continue
            fi
            if grep -qx "${cid}" "${FAIL_FILE}" 2>/dev/null; then
                submitted_chunks["${cid}"]=1
                continue
            fi

            submit_chunk "${cid}" "${files_csv}" \
                && submitted_chunks["${cid}"]=1 \
                || log "WARNING: submit_chunk failed for ${cid}, will retry next poll"

            SLOTS=$(( SLOTS - 1 ))
        done < "${CHUNKS_FILE}"
    fi

    log "Sleeping ${POLL_INTERVAL}s ..."
    sleep "${POLL_INTERVAL}"
done

# ── Final summary ─────────────────────────────────────────────────────────────
DONE_CHUNKS=$(wc -l < "${DONE_FILE}")
FAIL_CHUNKS=$(wc -l < "${FAIL_FILE}")

log "=== All chunks processed: ${DONE_CHUNKS} succeeded, ${FAIL_CHUNKS} failed ==="

if [[ "${FAIL_CHUNKS}" -gt 0 ]]; then
    FAILED_LIST=$(cat "${FAIL_FILE}")
    send_email "Completed with errors — ${TOMO_RUN_ID}" \
        "Processing complete.\n\nSucceeded: ${DONE_CHUNKS}\nFailed   : ${FAIL_CHUNKS}\n\nFailed chunks:\n${FAILED_LIST}\n\nResults at: ${DST_RESULTS}\nLogs      : ${LOG_DIR}"
else
    send_email "Completed successfully — ${TOMO_RUN_ID}" \
        "All ${DONE_CHUNKS} chunks processed successfully.\n\nResults at: ${DST_RESULTS}\nLogs      : ${LOG_DIR}"
fi

log "Orchestrator done."
