#!/bin/bash
# orchestrator.sh — long-running loop on the login node.
# Launched by tomo_process via nohup; never run directly.
#
# Required environment variables (set by tomo_process via export):
#   TOMO_RUN_ID, TOMO_INPUT_PATH, TOMO_EMAIL,
#   TOMO_CHUNK_SIZE, TOMO_MAX_JOBS, TOMO_REPO_DIR,
#   TOMO_SCRATCH_ROOT, TOMO_DATA_MOUNT, TOMO_RESULTS_MOUNT,
#   TOMO_POLL_INTERVAL, TOMO_SLURM_PARTITION, TOMO_SLURM_GRES,
#   TOMO_SLURM_CPUS, TOMO_SLURM_MEM, TOMO_SLURM_TIME, TOMO_MINS_PER_SAMPLE,
#   TOMO_NM
# Optional:
#   TOMO_ORCH_PARTITION  cpu partition used for SLURM mail-fallback jobs (default: cpu)
#   TOMO_SMTP_RELAY      SMTP relay host[:port] for direct email via s-nail
#   TOMO_VERSION         pipeline version: default | alice (default: default)

set -euo pipefail

TOMO_VERSION="${TOMO_VERSION:-default}"

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
# Chunks that processed fine but whose results could not be copied to /mnt;
# reported in the final summary email instead of a per-chunk email.
COPYFAIL_FILE="${STATE_DIR}/copyfail_chunks.txt"
# Input files the job skipped as unreadable/corrupt: chunk <TAB> file <TAB> reason
SKIPPED_FILE="${STATE_DIR}/skipped_files.txt"
# Failed chunks whose unfinished samples were re-queued (recovered, not lost)
RETRIED_FILE="${STATE_DIR}/retried_chunks.txt"

# ── Logging ───────────────────────────────────────────────────────────────────
# stdout/stderr are already redirected to ORCH_LOG by nohup in tomo_process;
# write only to stdout here to avoid duplicate lines from tee.
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

fmt_duration() {
    local s=$1
    printf '%dh %02dm' $(( s / 3600 )) $(( (s % 3600) / 60 ))
}

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
    # s-nail can speak SMTP directly to an institute relay — no local MTA
    # needed on the login node. Full subject + body arrives as a normal email.
    if [[ "${sent}" == false && -n "${TOMO_SMTP_RELAY:-}" ]] && command -v mail &>/dev/null; then
        if echo -e "${body}" | mail -s "${full_subject}" \
                -S mta="smtp://${TOMO_SMTP_RELAY}" \
                -S from="${TOMO_EMAIL}" \
                "${TOMO_EMAIL}" &>/dev/null; then
            sent=true
        fi
    fi
    if [[ "${sent}" == false ]] && command -v mail &>/dev/null; then
        if echo -e "${body}" | mail -s "${full_subject}" "${TOMO_EMAIL}" &>/dev/null; then
            sent=true
        fi
    fi
    # Fallback: the login node has no mail relay, but SLURM's own mail
    # (sent by slurmctld) is known to work on this cluster. Submit a
    # zero-work job whose NAME carries the message — the user receives
    # SLURM's END notification with that name in the subject line.
    if [[ "${sent}" == false ]] && command -v sbatch &>/dev/null; then
        local jname="tomo_process ${subject}"
        jname="${jname//—/-}"                       # em-dash → plain dash
        jname="${jname// /_}"                       # spaces → underscores
        jname="$(printf '%s' "${jname}" | tr -cd 'A-Za-z0-9._:-' | cut -c1-120)"
        if sbatch --job-name="${jname}" \
                  --partition="${TOMO_ORCH_PARTITION:-cpu}" \
                  --time=00:02:00 --mem=100M \
                  --output=/dev/null --error=/dev/null \
                  --mail-type=END --mail-user="${TOMO_EMAIL}" \
                  --wrap="true" &>/dev/null; then
            sent=true
            log "Email sent via SLURM notification job: ${jname}"
        fi
    fi
    # Always write a notification file — readable via --status log
    local note_file="${LOG_DIR}/notification_$(date '+%Y%m%d_%H%M%S').txt"
    printf "Subject: %s\n\n%b\n" "${full_subject}" "${body}" > "${note_file}"
    if [[ "${sent}" == false ]]; then
        log "NOTE: email not sent (no mail relay, SLURM mail fallback failed); notification saved to ${note_file}"
    fi
}

log "=== Orchestrator start (pid $$) ==="
log "Run ID  : ${TOMO_RUN_ID}"
log "Input   : ${SRC_MOUNT}"
log "Results : ${TOMO_RESULTS_MOUNT}/<rel_path>/field_retrieval_zpe_results/ (per source dir)"
log "Scratch : ${RUN_SCRATCH}"
log "Version : ${TOMO_VERSION}"

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
    date +%s > "${STATE_DIR}/run_start_epoch"
    # No "started" email from here — the first chunk job carries SLURM
    # --mail-type=BEGIN and notifies the user when processing actually starts.
fi
touch "${COPYFAIL_FILE}" "${SKIPPED_FILE}" "${RETRIED_FILE}"
# Resumed runs started before duration tracking existed: approximate the run
# start with the chunk list's creation time so the summary still has a total.
[[ -f "${STATE_DIR}/run_start_epoch" ]] || stat -c %Y "${CHUNKS_FILE}" > "${STATE_DIR}/run_start_epoch"
DURATIONS_FILE="${STATE_DIR}/chunk_durations.txt"
SUBMIT_TIMES_FILE="${STATE_DIR}/submit_times.txt"
touch "${DURATIONS_FILE}" "${SUBMIT_TIMES_FILE}"

# ── Helpers ───────────────────────────────────────────────────────────────────
# True if the SLURM job is still queued/running. squeue exits 0 with empty
# output for a just-finished job and prints "Invalid job id" once it is purged;
# both mean gone. Any other output (including "unable to contact controller")
# is treated as still active, so a SLURM outage never marks jobs as lost.
job_active() {
    local out
    out=$(squeue -j "$1" -h 2>&1) || true
    [[ -n "${out}" && "${out}" != *"Invalid job id"* ]]
}

count_active() {
    local n=0
    [[ -s "${ACTIVE_FILE}" ]] || { echo 0; return 0; }
    while IFS=$'\t' read -r _cid jid; do
        job_active "${jid}" && n=$(( n + 1 ))
    done < "${ACTIVE_FILE}"
    echo "${n}"
}

# Elapsed time since a chunk was submitted (queue wait + runtime); appended to
# DURATIONS_FILE for the final summary email.
record_chunk_duration() {
    local cid="$1" status="$2"
    local sub_epoch dur_str="unknown"
    sub_epoch=$(awk -F'\t' -v c="${cid}" '$1==c {print $2}' "${SUBMIT_TIMES_FILE}" | tail -1)
    if [[ -n "${sub_epoch}" ]]; then
        dur_str=$(fmt_duration $(( $(date +%s) - sub_epoch )))
    fi
    echo "${cid}: ${dur_str}${status:+ (${status})}" >> "${DURATIONS_FILE}"
    echo "${dur_str}"
}

# Copy a chunk's results to the results mount. Used for finished chunks and, to
# salvage already-processed samples, for failed ones. Returns non-zero (and
# keeps the scratch dir) if the copy fails. Never aborts the caller.
copy_chunk_results() {
    local cid="$1"
    local src_dir
    src_dir=$(awk -F'\t' -v c="${cid}" '$1==c {print $2}' "${CHUNKS_FILE}")
    local rel_src="${src_dir#${TOMO_DATA_MOUNT}/}"
    local dst_dir="${TOMO_RESULTS_MOUNT}/${rel_src}/field_retrieval_zpe_results"
    local res_src="${RUN_SCRATCH}/${cid}/data/field_retrieval"

    if [[ -z "$(ls -A "${res_src}" 2>/dev/null)" ]]; then
        log "  No results to copy for ${cid}"
        return 0
    fi
    log "Copying ${cid} results to ${dst_dir}/"
    mkdir -p "${dst_dir}"
    # Every chunk writes the same log/inspection filenames; suffix them
    # with the chunk id so copies from sibling chunks don't overwrite
    # each other in the shared results directory.
    local f
    for f in field_retrieval.log tomogram_reconstruction.log Field_inspection.png skipped_files.txt; do
        if [[ -f "${res_src}/${f}" ]]; then
            mv "${res_src}/${f}" "${res_src}/${f%.*}_${cid}.${f##*.}"
        fi
    done
    # Input files the job could not read (corrupt MAT etc.) — report them
    local skipped="${res_src}/skipped_files_${cid}.txt"
    if [[ -s "${skipped}" ]]; then
        while IFS=$'\t' read -r sname sstage sreason; do
            log "  WARNING: ${cid} skipped ${sname} (${sstage}): ${sreason}"
            echo "${cid}"$'\t'"${sname}"$'\t'"${sreason}" >> "${SKIPPED_FILE}"
        done < "${skipped}"
    fi
    # --inplace + plain -r: the results mount is SMB/CIFS, which rejects
    # rsync's mkstemp temp files and chown/chmod/utimes for accounts that
    # don't map to the share owner — write directly, data only.
    # --ignore-times: always overwrite existing results from a previous
    # run — the size+mtime quick-check is unreliable on network mounts
    # *.part are in-progress atomic saves from a crashed job — never copy them
    if rsync -r --inplace --ignore-times --exclude='*.part' "${res_src}/" "${dst_dir}/"; then
        log "  Results copied for ${cid}"
        return 0
    fi
    log "WARNING: rsync to /mnt failed for ${cid} — results kept at ${res_src}/"
    return 1
}

# Sample files of a failed chunk that never produced a tomogram and were not
# already reported as unreadable. Prints them comma-separated. Must run after
# copy_chunk_results (which records MATLAB-reported skips) and before the
# chunk's scratch dir is removed.
unfinished_files() {
    local cid="$1" files_csv="$2"
    local res_src="${RUN_SCRATCH}/${cid}/data/field_retrieval"
    local flist fpath base out=()
    IFS=',' read -ra flist <<< "${files_csv}"
    for fpath in "${flist[@]}"; do
        base="$(basename "${fpath}" .mat)"
        [[ -f "${res_src}/Tomogram_Field_${base}.mat" ]] && continue
        compgen -G "${res_src}/Tomogram_Field_${base}_bg_*.mat" > /dev/null && continue
        awk -F'\t' -v c="${cid}" -v b="${base}" '$1==c && index($2, b) {f=1} END {exit !f}' \
            "${SKIPPED_FILE}" && continue
        out+=("${fpath}")
    done
    (IFS=','; echo "${out[*]}")
}

# Re-queue the unfinished samples of a failed chunk as two smaller chunks
# (bisection). A crash on one bad file thus isolates it in log2(N) rounds
# while every good sample still gets processed; a walltime kill shrinks the
# work per job until it fits. A single sample that still fails is reported
# as skipped. Sets REQUEUE_RESULT to: retry | skipped | none.
requeue_unfinished() {
    local cid="$1" src_dir="$2" remaining="$3"
    local rem n half k ncid lo hi
    REQUEUE_RESULT="none"
    [[ -n "${remaining}" ]] || return 0
    IFS=',' read -ra rem <<< "${remaining}"
    n=${#rem[@]}
    if [[ "${n}" -eq 1 ]]; then
        # One extra attempt before giving up, so a one-off node failure or
        # OOM-kill on a lone file is not reported as a bad file.
        local key retry_dir="${STATE_DIR}/single_retry"
        mkdir -p "${retry_dir}"
        key=$(printf '%s' "${rem[0]}" | cksum | cut -d' ' -f1)
        if [[ ! -f "${retry_dir}/${key}" ]]; then
            touch "${retry_dir}/${key}"
            k=$(( $(cat "${STATE_DIR}/retry_counter" 2>/dev/null || echo 0) + 1 ))
            echo "${k}" > "${STATE_DIR}/retry_counter"
            ncid=$(printf "retry%04d" "${k}")
            printf '%s\t%s\t%s\n' "${ncid}" "${src_dir}" "${rem[0]}" >> "${CHUNKS_FILE}"
            echo "${cid} -> ${ncid}" >> "${RETRIED_FILE}"
            log "  ${cid}: $(basename "${rem[0]}") failed alone — one more attempt as ${ncid}"
            REQUEUE_RESULT="retry"
            return 0
        fi
        log "  ${cid}: $(basename "${rem[0]}") failed again on its own — giving up on it"
        printf '%s\t%s\t%s\n' "${cid}" "$(basename "${rem[0]}")" "job crashed, timed out or ran out of memory on this file, twice (see logs)" >> "${SKIPPED_FILE}"
        REQUEUE_RESULT="skipped"
        return 0
    fi
    half=$(( (n + 1) / 2 ))
    for lo in 0 "${half}"; do
        if [[ "${lo}" -eq 0 ]]; then hi=${half}; else hi=$(( n - half )); fi
        k=$(( $(cat "${STATE_DIR}/retry_counter" 2>/dev/null || echo 0) + 1 ))
        echo "${k}" > "${STATE_DIR}/retry_counter"
        ncid=$(printf "retry%04d" "${k}")
        echo "${ncid}"$'\t'"${src_dir}"$'\t'"$(IFS=','; echo "${rem[*]:${lo}:${hi}}")" >> "${CHUNKS_FILE}"
        log "  ${cid}: re-queued ${hi} unfinished sample(s) as ${ncid}"
    done
    echo "${cid} -> bisected into ${n} unfinished sample(s)" >> "${RETRIED_FILE}"
    REQUEUE_RESULT="retry"
}

collect_finished() {
    [[ -s "${ACTIVE_FILE}" ]] || return 0
    local still_active=""
    while IFS=$'\t' read -r cid jid; do
        local done_marker="${RUN_SCRATCH}/${cid}/DONE"
        local fail_marker="${RUN_SCRATCH}/${cid}/FAIL"

        # An orchestrator killed (e.g. by --restart) between removing a chunk's
        # scratch dir and recording its outcome re-adopts it here with nothing
        # left to inspect. Drop it if its outcome was recorded; otherwise count
        # it as failed so the run can still finish.
        if [[ ! -d "${RUN_SCRATCH}/${cid}" ]]; then
            if ! grep -qx "${cid}" "${DONE_FILE}" "${FAIL_FILE}" 2>/dev/null; then
                log "WARNING: scratch dir for ${cid} is gone and no outcome was recorded — counting it as failed"
                echo "${cid}" >> "${FAIL_FILE}"
            fi
            continue
        fi

        if [[ -f "${done_marker}" ]]; then
            local dur
            dur=$(record_chunk_duration "${cid}" "")
            log "Chunk ${cid} finished after ${dur} (queue + run)"
            if copy_chunk_results "${cid}"; then
                rm -rf "${RUN_SCRATCH:?}/${cid}"
            else
                echo "${cid}" >> "${COPYFAIL_FILE}"
            fi
            echo "${cid}" >> "${DONE_FILE}"

        elif [[ -f "${fail_marker}" ]]; then
            record_chunk_duration "${cid}" "failed" > /dev/null
            log "WARNING: Chunk ${cid} FAILED (job ${jid}) — logs: ${LOG_DIR}/${TOMO_RUN_ID}_${cid}.err"
            # Salvage: samples processed before the failure are still valid
            # outputs — don't throw them away along with the failed one.
            local copy_ok=1 files_csv src_dir
            copy_chunk_results "${cid}" || copy_ok=0
            # Then retry whatever never finished, in smaller chunks, so one bad
            # sample (or a walltime/OOM kill) cannot take its neighbours down.
            files_csv=$(awk -F'\t' -v c="${cid}" '$1==c {print $3}' "${CHUNKS_FILE}")
            src_dir=$(awk -F'\t' -v c="${cid}" '$1==c {print $2}' "${CHUNKS_FILE}")
            requeue_unfinished "${cid}" "${src_dir}" "$(unfinished_files "${cid}" "${files_csv}")"
            if [[ "${copy_ok}" -eq 1 ]]; then
                rm -rf "${RUN_SCRATCH:?}/${cid}"
            else
                echo "${cid}" >> "${COPYFAIL_FILE}"
            fi
            if [[ "${REQUEUE_RESULT}" == "none" ]]; then
                echo "${cid}" >> "${FAIL_FILE}"
            else
                # Its remaining work is accounted for (retry chunks or a skipped
                # file in the final summary), so it is not a lost chunk.
                echo "${cid}" >> "${DONE_FILE}"
            fi

        else
            # Neither marker: the job may have been killed outright (SIGKILL,
            # node failure, OOM-kill) before its trap could run. If SLURM no
            # longer knows the job, give the filesystem 2 min to show a late
            # marker, then declare the chunk failed so it is salvaged/retried.
            if ! job_active "${jid}"; then
                local gone="${RUN_SCRATCH}/${cid}/JOB_GONE"
                if [[ ! -f "${gone}" ]]; then
                    touch "${gone}"
                elif [[ -n "$(find "${gone}" -mmin +2 2>/dev/null)" ]]; then
                    log "WARNING: job ${jid} for ${cid} left the queue without DONE/FAIL — treating as failed"
                    touch "${fail_marker}"
                fi
            fi
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

    # Copy background file(s) for this source directory. Background file names
    # vary across experiments (e.g. "bg_Tomog.mat", "bg001_Tomog.mat"), so glob
    # for any "bg*_Tomog.mat" rather than assuming a fixed name.
    local bg_found=0
    for bg_src in "${src_dir}"/bg*_Tomog.mat; do
        [[ -f "${bg_src}" ]] || continue
        rsync -a "${bg_src}" "${chunk_data}/" \
            || { log "ERROR: rsync failed for background ${bg_src}"; return 1; }
        log "  Copied: $(basename "${bg_src}")"
        bg_found=1
    done
    if [[ "${bg_found}" -eq 0 ]]; then
        log "ERROR: no background file (bg*_Tomog.mat) found in: ${src_dir}"
        return 1
    fi

    # Copy sample files for this chunk (all from the same source directory)
    IFS=',' read -ra flist <<< "${files_csv}"
    # A sample that cannot be staged (unreadable/corrupt on the source) is
    # skipped and reported; it must not block the rest of the chunk.
    local staged=0 tiff_idx
    : > "${chunk_data}/tiff_index.txt"
    for fpath in "${flist[@]}"; do
        if rsync -a "${fpath}" "${chunk_data}/"; then
            staged=$(( staged + 1 ))
            # Stable 1-based position of this sample among all samples of its
            # source dir (sorted order), independent of chunking, skips and
            # retries; used for alice TIFF numbering.
            tiff_idx=$(awk -v d="${src_dir}" -v f="${fpath}" \
                '{p=$0; sub(/\/[^\/]*$/, "", p)} p==d {n++; if ($0==f) {print n; exit}}' \
                "${STATE_DIR}/pending_files.txt")
            [[ -n "${tiff_idx}" ]] && \
                printf '%s\t%s\n' "$(basename "${fpath}" .mat)" "${tiff_idx}" >> "${chunk_data}/tiff_index.txt"
        else
            log "WARNING: rsync failed for ${fpath} — skipping this file"
            rm -f "${chunk_data}/$(basename "${fpath}")"
            echo "${cid}"$'\t'"$(basename "${fpath}")"$'\t'"could not be copied from source" >> "${SKIPPED_FILE}"
        fi
    done
    if [[ "${staged}" -eq 0 ]]; then
        log "ERROR: no sample file could be staged for ${cid}"
        return 1
    fi
    log "  Copied ${staged}/${#flist[@]} sample files for ${cid}"

    # Alice-version TIFFs are named <batchName><NNN>_Tomo.tif. The chunk's data
    # dir is a flat staging folder, so pass the source dir's name as the batch
    # name, plus the number of samples from this source dir in earlier chunks
    # so NNN continues across chunks instead of restarting (and overwriting
    # sibling chunks' TIFFs in the shared results dir). tiff_index.txt (written
    # above) gives each sample its exact number; tiff_offset is only the
    # fallback for samples missing from it.
    local batch_name tiff_offset
    batch_name="$(basename "${src_dir}")"
    tiff_offset=$(awk -F'\t' -v c="${cid}" -v d="${src_dir}" \
        '$1==c {exit} $2==d {n+=split($3, a, ",")} END {print n+0}' "${CHUNKS_FILE}")
    # Escape for the bash single-quoted assignment in the template, then for sed
    batch_name="${batch_name//\'/\'\\\'\'}"
    batch_name="$(printf '%s' "${batch_name}" | sed -e 's/[\\&|]/\\&/g')"

    # Generate job script from template; DATA_DIR = chunk_data (flat, no nesting)
    local job_script="${JOBS_DIR}/${cid}.sh"
    local job_name="${TOMO_RUN_ID}_${cid}"
    sed \
        -e "s|__JOB_NAME__|${job_name}|g" \
        -e "s|__LOG_DIR__|${LOG_DIR}|g" \
        -e "s|__DATA_DIR__|${chunk_data}|g" \
        -e "s|__CHUNK_DONE__|${chunk_scratch}/DONE|g" \
        -e "s|__CHUNK_FAIL__|${chunk_scratch}/FAIL|g" \
        -e "s|__REPO_DIR__|${TOMO_REPO_DIR}|g" \
        -e "s|__PARTITION__|${TOMO_SLURM_PARTITION}|g" \
        -e "s|__GRES__|${TOMO_SLURM_GRES}|g" \
        -e "s|__CPUS__|${TOMO_SLURM_CPUS}|g" \
        -e "s|__MEM__|${TOMO_SLURM_MEM}|g" \
        -e "s|__TIME__|${TOMO_SLURM_TIME}|g" \
        -e "s|__NM__|${TOMO_NM}|g" \
        -e "s|__VERSION__|${TOMO_VERSION}|g" \
        -e "s|__BATCH_NAME__|${batch_name}|g" \
        -e "s|__TIFF_OFFSET__|${tiff_offset}|g" \
        "${TOMO_REPO_DIR}/scripts/bash_template.sh" > "${job_script}"
    chmod +x "${job_script}"

    local jid
    jid=$(sbatch --parsable "${job_script}")
    log "Submitted ${cid} → SLURM job ${jid}"
    echo "${cid}"$'\t'"${jid}" >> "${ACTIVE_FILE}"
    # Kept separate from ACTIVE_FILE: tomo_process --cancel parses that file
    # as exactly two tab-separated fields.
    echo "${cid}"$'\t'"$(date +%s)" >> "${SUBMIT_TIMES_FILE}"

    # "Processing started" email: a zero-work job that SLURM holds until the
    # first chunk job actually starts running (--dependency=after). Its name
    # is what appears in the email subject:
    #   Slurm Job_id=... Name=tomo_process_Started_-_<run_id> Began
    # Its job id is saved in state/ so --cancel can kill it (and so a resumed
    # run doesn't email twice).
    if [[ ! -f "${STATE_DIR}/start_notify_jid" ]]; then
        local notify_jid
        if notify_jid=$(sbatch --parsable \
                  --job-name="tomo_process_Started_-_${TOMO_RUN_ID}" \
                  --partition="${TOMO_ORCH_PARTITION:-cpu}" \
                  --time=00:02:00 --mem=100M \
                  --output=/dev/null --error=/dev/null \
                  --dependency="after:${jid}" \
                  --mail-type=BEGIN --mail-user="${TOMO_EMAIL}" \
                  --wrap="true" 2>/dev/null); then
            echo "${notify_jid}" > "${STATE_DIR}/start_notify_jid"
            log "Start-notification job ${notify_jid} submitted (emails when ${cid} begins)"
        else
            log "WARNING: could not submit start-notification job"
        fi
    fi
    return 0
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
copyfail_n=$(wc -l < "${COPYFAIL_FILE}" 2>/dev/null || echo 0)
total_time=$(fmt_duration $(( $(date +%s) - $(cat "${STATE_DIR}/run_start_epoch") )))
log "=== All done in ${total_time}: ${local_done} succeeded, ${local_fail} failed, ${copyfail_n} result-copy failures ==="

timing="\nTotal time: ${total_time}\n\nPer-chunk times (queue + run):\n$(cat "${DURATIONS_FILE}")\n"

retried_n=$(wc -l < "${RETRIED_FILE}")
if [[ "${retried_n}" -gt 0 ]]; then
    timing+="\n${retried_n} chunk failure(s) were retried automatically:\n$(sed 's/^/  /' "${RETRIED_FILE}")\n"
fi
skipped_n=$(wc -l < "${SKIPPED_FILE}")
if [[ "${skipped_n}" -gt 0 ]]; then
    timing+="\nSkipped ${skipped_n} unreadable/corrupt input file(s) (all other samples were processed and copied):\n$(awk -F'\t' '{print "  " $1 ": " $2 " — " $3}' "${SKIPPED_FILE}")\n"
fi

if [[ "${local_fail}" -gt 0 || "${copyfail_n}" -gt 0 || "${skipped_n}" -gt 0 ]]; then
    details=""
    if [[ "${local_fail}" -gt 0 ]]; then
        details+="\nFailed chunks (logs in ${LOG_DIR}):\n$(cat "${FAIL_FILE}")\n"
    fi
    if [[ "${copyfail_n}" -gt 0 ]]; then
        details+="\nChunks processed OK but results could NOT be copied to ${TOMO_RESULTS_MOUNT} (kept on scratch under ${RUN_SCRATCH}):\n$(cat "${COPYFAIL_FILE}")\n"
    fi
    send_email "Completed with errors — ${TOMO_RUN_ID} in ${total_time}" \
        "Processing complete.\nSucceeded: ${local_done}\nFailed: ${local_fail}\nResult-copy failures: ${copyfail_n}\nSkipped input files: ${skipped_n}\n${timing}${details}\nResults: <src_dir>/field_retrieval_zpe_results/ (per source directory)\nLogs: ${LOG_DIR}"
else
    send_email "Completed successfully — ${TOMO_RUN_ID} in ${total_time}" \
        "All ${local_done} chunks processed successfully.\n${timing}\nResults: <src_dir>/field_retrieval_zpe_results/ (per source directory)\nLogs: ${LOG_DIR}"
fi
