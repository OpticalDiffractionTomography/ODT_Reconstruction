#!/bin/bash
# orchestrator.sh — long-running loop on the login node.
# Launched by tomo_process via nohup; never run directly.
#
# Required environment variables (set by tomo_process via export):
#   TOMO_RUN_ID, TOMO_INPUT_PATH, TOMO_EMAIL,
#   TOMO_CHUNK_SIZE, TOMO_MAX_JOBS, TOMO_REPO_DIR,
#   TOMO_SCRATCH_ROOT, TOMO_DATA_MOUNT, TOMO_RESULTS_MOUNT,
#   TOMO_POLL_INTERVAL, TOMO_SLURM_PARTITION, TOMO_SLURM_GRES,
#   TOMO_SLURM_CPUS, TOMO_SLURM_MEM, TOMO_SLURM_TIME, TOMO_MINS_PER_SAMPLE

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────────
RUN_SCRATCH="${TOMO_SCRATCH_ROOT}/${TOMO_RUN_ID}"
STATE_DIR="${RUN_SCRATCH}/state"
LOG_DIR="${RUN_SCRATCH}/logs"
JOBS_DIR="${RUN_SCRATCH}/jobs"
ORCH_LOG="${LOG_DIR}/orchestrator.log"

SRC_MOUNT="${TOMO_DATA_MOUNT}/${TOMO_INPUT_PATH}"

CHUNKS_FILE="${STATE_DIR}/chunks.txt"
ACTIVE_FILE="${STATE_DIR}/active_jobs.txt"
DONE_FILE="${STATE_DIR}/done_chunks.txt"
FAIL_FILE="${STATE_DIR}/failed_chunks.txt"

# ── Logging ───────────────────────────────────────────────────────────────────
# stdout/stderr are already redirected to ORCH_LOG by nohup in tomo_process;
# write only to stdout here to avoid duplicate lines from tee.
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

send_email() {
    local subject="$1" body="$2"
    local full_subject="[tomo_process] ${subject}"
    local sent=false
    # Redirect both stdout and stderr to suppress dead.letter noise
    if command -v sendmail &>/dev/null; then
        if printf "Subject: %s\nTo: %s\n\n%b\n" \
                "${full_subject}" "${TOMO_EMAIL}" "${body}" \
                | sendmail "${TOMO_EMAIL}" &>/dev/null; then
            sent=true
        fi
    fi
    if [[ "${sent}" == false ]] && command -v mail &>/dev/null; then
        if echo -e "${body}" | mail -s "${full_subject}" "${TOMO_EMAIL}" &>/dev/null; then
            sent=true
        fi
    fi
    # Always write a notification file — readable via --status log
    local note_file="${LOG_DIR}/notification_$(date '+%Y%m%d_%H%M%S').txt"
    printf "Subject: %s\n\n%b\n" "${full_subject}" "${body}" > "${note_file}"
    if [[ "${sent}" == false ]]; then
        log "NOTE: email not sent (no mail relay); notification saved to ${note_file}"
    fi
}

log "=== Orchestrator start (pid $$) ==="
log "Run ID  : ${TOMO_RUN_ID}"
log "Input   : ${SRC_MOUNT}"
log "Results : ${TOMO_RESULTS_MOUNT}/<rel_path>/field_retrieval_zpe_results/ (per source dir)"
log "Scratch : ${RUN_SCRATCH}"

# ── Discover and chunk files (skipped on resume — chunks.txt already exists) ──
# chunks.txt format: chunk_id <TAB> src_dir <TAB> file1,file2,...
# State files on scratch survive login node reboots, so a resumed orchestrator
# simply re-adopts active jobs and continues from where it left off.

if [[ -f "${CHUNKS_FILE}" ]]; then
    NCHUNKS=$(wc -l < "${CHUNKS_FILE}")
    DONE_N=$(wc -l < "${DONE_FILE}")
    FAIL_N=$(wc -l < "${FAIL_FILE}")
    log "RESUMING existing run — ${NCHUNKS} chunks total, ${DONE_N} done, ${FAIL_N} failed"
else
    log "Fresh run — discovering sample files under ${SRC_MOUNT}"
    find "${SRC_MOUNT}" -name "sample*_Tomog.mat" | sort > "${STATE_DIR}/pending_files.txt"
    TOTAL=$(wc -l < "${STATE_DIR}/pending_files.txt")
    log "Found ${TOTAL} sample files"

    if [[ "${TOTAL}" -eq 0 ]]; then
        log "ERROR: no sample*_Tomog.mat files found under ${SRC_MOUNT}"
        send_email "FAILED — ${TOMO_RUN_ID}" "No sample files found under ${SRC_MOUNT}"
        exit 1
    fi

    # Auto-size: max samples per job = floor(walltime_minutes / mins_per_sample)
    walltime_h=$(echo "${TOMO_SLURM_TIME}" | cut -d: -f1)
    walltime_m=$(echo "${TOMO_SLURM_TIME}" | cut -d: -f2)
    walltime_mins=$(( 10#${walltime_h} * 60 + 10#${walltime_m} ))
    max_per_job=$(( walltime_mins / TOMO_MINS_PER_SAMPLE ))
    [[ "${max_per_job}" -lt 1 ]] && max_per_job=1
    log "Walltime ${TOMO_SLURM_TIME} / ${TOMO_MINS_PER_SAMPLE} min per sample → max ${max_per_job} samples per job"

    chunk_idx=0
    current_dir=""
    chunk_files=()

    flush_chunk() {
        [[ ${#chunk_files[@]} -eq 0 ]] && return 0
        local cid
        cid=$(printf "chunk%04d" ${chunk_idx})
        echo "${cid}"$'\t'"${current_dir}"$'\t'"$(IFS=','; echo "${chunk_files[*]}")" >> "${CHUNKS_FILE}"
        log "  ${cid}: ${#chunk_files[@]} files from ${current_dir}"
        chunk_idx=$(( chunk_idx + 1 ))
        chunk_files=()
    }

    while IFS= read -r fpath; do
        local_dir=$(dirname "${fpath}")
        if [[ "${local_dir}" != "${current_dir}" ]]; then
            flush_chunk
            current_dir="${local_dir}"
        fi
        chunk_files+=("${fpath}")
        if [[ ${#chunk_files[@]} -ge ${max_per_job} ]]; then
            flush_chunk
        fi
    done < "${STATE_DIR}/pending_files.txt"
    flush_chunk

    NCHUNKS=$(wc -l < "${CHUNKS_FILE}")
    log "Split into ${NCHUNKS} job(s) — auto-sized to ~${max_per_job} samples each"

    touch "${ACTIVE_FILE}" "${DONE_FILE}" "${FAIL_FILE}"
    send_email "Started — ${TOMO_RUN_ID}" \
        "Processing started.\nInput : ${SRC_MOUNT}\nChunks: ${NCHUNKS}\nMax parallel jobs: ${TOMO_MAX_JOBS}"
fi

# ── Helpers ───────────────────────────────────────────────────────────────────
count_active() {
    local n=0
    [[ -s "${ACTIVE_FILE}" ]] || { echo 0; return 0; }
    while IFS=$'\t' read -r _cid jid; do
        squeue -j "${jid}" -h &>/dev/null 2>&1 && n=$(( n + 1 ))
    done < "${ACTIVE_FILE}"
    echo "${n}"
}

collect_finished() {
    [[ -s "${ACTIVE_FILE}" ]] || return 0
    local still_active=""
    while IFS=$'\t' read -r cid jid; do
        local done_marker="${RUN_SCRATCH}/${cid}/DONE"
        local fail_marker="${RUN_SCRATCH}/${cid}/FAIL"

        if [[ -f "${done_marker}" ]]; then
            # Resolve source dir for this chunk from chunks.txt to build result path
            local src_dir
            src_dir=$(awk -F'\t' -v c="${cid}" '$1==c {print $2}' "${CHUNKS_FILE}")
            local rel_src="${src_dir#${TOMO_DATA_MOUNT}/}"
            local dst_dir="${TOMO_RESULTS_MOUNT}/${rel_src}/field_retrieval_zpe_results"
            log "Chunk ${cid} done — copying results to ${dst_dir}/"
            mkdir -p "${dst_dir}"
            rsync -a "${RUN_SCRATCH}/${cid}/data/field_retrieval/" "${dst_dir}/" \
                && log "  Results copied for ${cid}" \
                || log "WARNING: rsync to /mnt failed for ${cid}"
            rm -rf "${RUN_SCRATCH:?}/${cid}"
            echo "${cid}" >> "${DONE_FILE}"

        elif [[ -f "${fail_marker}" ]]; then
            log "WARNING: Chunk ${cid} FAILED (job ${jid})"
            echo "${cid}" >> "${FAIL_FILE}"
            send_email "Chunk FAILED — ${TOMO_RUN_ID}" \
                "Chunk ${cid} (job ${jid}) failed.\nLogs: ${LOG_DIR}/${cid}.err\nScratch kept at: ${RUN_SCRATCH}/${cid}"
            rm -rf "${RUN_SCRATCH:?}/${cid}"

        else
            still_active+="${cid}"$'\t'"${jid}"$'\n'
        fi
    done < "${ACTIVE_FILE}"
    printf "%s" "${still_active}" > "${ACTIVE_FILE}"
}

submit_chunk() {
    local cid="$1" src_dir="$2" files_csv="$3"
    local chunk_scratch="${RUN_SCRATCH}/${cid}"
    # chunk_data is a flat directory — MATLAB's spath points here directly
    local chunk_data="${chunk_scratch}/data"
    mkdir -p "${chunk_data}/field_retrieval"

    # Copy background file (one per source directory, always present alongside samples)
    local bg_src="${src_dir}/bg001_Tomog.mat"
    if [[ -f "${bg_src}" ]]; then
        rsync -a "${bg_src}" "${chunk_data}/" \
            || { log "ERROR: rsync failed for background ${bg_src}"; return 1; }
        log "  Copied: bg001_Tomog.mat"
    else
        log "ERROR: background not found: ${bg_src}"
        return 1
    fi

    # Copy sample files for this chunk (all from the same source directory)
    IFS=',' read -ra flist <<< "${files_csv}"
    for fpath in "${flist[@]}"; do
        rsync -a "${fpath}" "${chunk_data}/" \
            || { log "ERROR: rsync failed for ${fpath}"; return 1; }
    done
    log "  Copied ${#flist[@]} sample files for ${cid}"

    # Generate job script from template; DATA_DIR = chunk_data (flat, no nesting)
    local job_script="${JOBS_DIR}/${cid}.sh"
    local job_name="tomo_${TOMO_RUN_ID}_${cid}"
    sed \
        -e "s|__JOB_NAME__|${job_name}|g" \
        -e "s|__LOG_DIR__|${LOG_DIR}|g" \
        -e "s|__EMAIL__|${TOMO_EMAIL}|g" \
        -e "s|__DATA_DIR__|${chunk_data}|g" \
        -e "s|__CHUNK_DONE__|${chunk_scratch}/DONE|g" \
        -e "s|__CHUNK_FAIL__|${chunk_scratch}/FAIL|g" \
        -e "s|__REPO_DIR__|${TOMO_REPO_DIR}|g" \
        -e "s|__PARTITION__|${TOMO_SLURM_PARTITION}|g" \
        -e "s|__GRES__|${TOMO_SLURM_GRES}|g" \
        -e "s|__CPUS__|${TOMO_SLURM_CPUS}|g" \
        -e "s|__MEM__|${TOMO_SLURM_MEM}|g" \
        -e "s|__TIME__|${TOMO_SLURM_TIME}|g" \
        "${TOMO_REPO_DIR}/scripts/bash_template.sh" > "${job_script}"
    chmod +x "${job_script}"

    local jid
    jid=$(sbatch --parsable "${job_script}")
    log "Submitted ${cid} → SLURM job ${jid}"
    echo "${cid}"$'\t'"${jid}" >> "${ACTIVE_FILE}"
}

# ── Main loop ─────────────────────────────────────────────────────────────────
log "Entering main loop (bash ${BASH_VERSION})"
declare -A submitted || { log "ERROR: declare -A failed (bash too old?)"; exit 1; }
log "State initialized"

while true; do
    collect_finished

    local_total=$(wc -l < "${CHUNKS_FILE}")
    local_done=$(wc -l < "${DONE_FILE}")
    local_fail=$(wc -l < "${FAIL_FILE}")
    local_finished=$(( local_done + local_fail ))

    log "Progress: ${local_done} done, ${local_fail} failed, ${local_finished}/${local_total} total"

    if [[ "${local_finished}" -ge "${local_total}" ]]; then
        break
    fi

    active=$(count_active)
    slots=$(( TOMO_MAX_JOBS - active ))

    if [[ "${slots}" -gt 0 ]]; then
        while IFS=$'\t' read -r cid src_dir files_csv; do
            [[ "${slots}" -le 0 ]] && break
            [[ -n "${submitted[${cid}]+x}" ]] && continue
            grep -qx "${cid}" "${DONE_FILE}" 2>/dev/null && { submitted["${cid}"]=1; continue; }
            grep -qx "${cid}" "${FAIL_FILE}" 2>/dev/null && { submitted["${cid}"]=1; continue; }

            log "Preparing ${cid} (source: ${src_dir}) ..."
            if submit_chunk "${cid}" "${src_dir}" "${files_csv}"; then
                submitted["${cid}"]=1
                slots=$(( slots - 1 ))
            else
                log "WARNING: submit failed for ${cid}, will retry next poll"
            fi
        done < "${CHUNKS_FILE}"
    fi

    sleep "${TOMO_POLL_INTERVAL}"
done

# ── Final summary ─────────────────────────────────────────────────────────────
local_done=$(wc -l < "${DONE_FILE}")
local_fail=$(wc -l < "${FAIL_FILE}")
log "=== All done: ${local_done} succeeded, ${local_fail} failed ==="

if [[ "${local_fail}" -gt 0 ]]; then
    send_email "Completed with errors — ${TOMO_RUN_ID}" \
        "Processing complete.\nSucceeded: ${local_done}\nFailed: ${local_fail}\n\nFailed chunks:\n$(cat "${FAIL_FILE}")\n\nResults: <src_dir>/field_retrieval_zpe_results/ (per source directory)\nLogs: ${LOG_DIR}"
else
    send_email "Completed successfully — ${TOMO_RUN_ID}" \
        "All ${local_done} chunks processed successfully.\n\nResults: <src_dir>/field_retrieval_zpe_results/ (per source directory)\nLogs: ${LOG_DIR}"
fi
