#!/bin/bash
# install.sh — one-time setup on the HPC login node.
# Run once after cloning the repository:
#   bash install.sh

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${REPO_DIR}/config.sh"

echo "=== tomo_process installer ==="
echo "Repo : ${REPO_DIR}"
echo ""

# ── 1. Verify mount point ─────────────────────────────────────────────────────
echo "[1/4] Checking data mount (DATA_MOUNT=${DATA_MOUNT}) ..."
if mountpoint -q "${DATA_MOUNT}" 2>/dev/null || [[ -d "${DATA_MOUNT}" ]]; then
    echo "      OK: ${DATA_MOUNT} is accessible"
else
    echo "      WARNING: ${DATA_MOUNT} is not mounted or not accessible."
    echo "      Update DATA_MOUNT in config.sh, then re-run install.sh."
    echo "      (Continuing install anyway.)"
fi

# ── 2. Create scratch root ────────────────────────────────────────────────────
echo "[2/4] Creating scratch root: ${SCRATCH_ROOT} ..."
mkdir -p "${SCRATCH_ROOT}"
echo "      OK"

# ── 3. Add tomo_process to PATH ───────────────────────────────────────────────
echo "[3/4] Installing tomo_process CLI ..."
chmod +x "${REPO_DIR}/tomo_process"

BIN_DIR="${HOME}/.local/bin"
mkdir -p "${BIN_DIR}"

ln -sf "${REPO_DIR}/tomo_process" "${BIN_DIR}/tomo_process"
echo "      Symlink: ${BIN_DIR}/tomo_process -> ${REPO_DIR}/tomo_process"

# Ensure ~/.local/bin is in PATH (add to ~/.bashrc if missing)
if ! grep -q 'PATH.*\.local/bin' "${HOME}/.bashrc" 2>/dev/null; then
    echo '' >> "${HOME}/.bashrc"
    echo '# tomo_process CLI' >> "${HOME}/.bashrc"
    echo 'export PATH="${HOME}/.local/bin:${PATH}"' >> "${HOME}/.bashrc"
    echo "      Added ~/.local/bin to PATH in ~/.bashrc"
else
    echo "      ~/.local/bin already in PATH"
fi

# ── 4. Verify SLURM is available ─────────────────────────────────────────────
echo "[4/4] Checking SLURM ..."
if command -v sbatch &>/dev/null; then
    echo "      OK: sbatch found ($(sbatch --version 2>&1 | head -1))"
else
    echo "      WARNING: sbatch not found. Are you on a SLURM login node?"
fi

echo ""
echo "=== Installation complete ==="
echo ""
echo "Reload your shell or run:"
echo "  source ~/.bashrc"
echo ""
echo "Then process a dataset with:"
echo "  tomo_process --email you@institute.de --path \"experiment_folder/batch_A\""
echo ""
echo "To edit cluster settings (partition, resources, chunk size):"
echo "  ${REPO_DIR}/config.sh"
