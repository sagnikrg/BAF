#!/bin/bash
set -euo pipefail

# -------- Config --------
DATA_DIR="../../data"                           # where the files are
BASE_PREFIX="mbldtc_L8_theta_0.0_"    # filename prefix
EXT="hdf5"                             # extension
START_ITR=1
END_ITR=2000

# Size threshold for corrupted (delete if smaller than this)
# Requirement: "file size is not less than 24 MB" => keep >= 24 MB, delete < 24 MB
THRESHOLD_MB=20
OUTPUT_FILE="small_files.txt"

DRY_RUN=0
if [[ "${1:-}" == "--dry-run" || "${1:-}" == "-n" ]]; then
  DRY_RUN=1
  echo "[DRY-RUN] No files will be modified."
fi

# -------- Helpers --------
bytes_in_mb=$((1024*1024))
THRESHOLD_BYTES=$((THRESHOLD_MB * bytes_in_mb))

filesize_bytes() {
  # cross-platform stat (Linux/macOS)
  if stat --version >/dev/null 2>&1; then
    stat -c%s "$1" 2>/dev/null
  else
    stat -f%z "$1" 2>/dev/null
  fi
}

file_path() {
  local itr="$1"
  printf "%s/%s%d.%s" "$DATA_DIR" "$BASE_PREFIX" "$itr" "$EXT"
}

# -------- 1) Delete small/corrupted files and log their indices --------
: > "$OUTPUT_FILE"
echo "Scanning ${START_ITR}..${END_ITR} and deleting files < ${THRESHOLD_MB} MB..."

for ((itr=START_ITR; itr<=END_ITR; itr++)); do
  f="$(file_path "$itr")"
  [[ -f "$f" ]] || continue

  sz="$(filesize_bytes "$f" || echo 0)"
  if (( sz < THRESHOLD_BYTES )); then
    echo "$itr" >> "$OUTPUT_FILE"
    if (( DRY_RUN )); then
      echo "[DRY-RUN] Would delete: $f (size: $((sz/bytes_in_mb)) MB)"
    else
      echo "Deleting: $f (size: $((sz/bytes_in_mb)) MB)"
      rm -f -- "$f"
    fi
  fi
done

echo "Corrupted indices logged to $OUTPUT_FILE"

# -------- 2) Build reindex plan for remaining files --------
# We will rename good files into a continuous sequence starting at 1.
# To avoid collisions (e.g., renaming 7->1 while 1 still exists), we:
#   (a) First move all files that need renaming to unique temp names
#   (b) Then move temps to their final target names

good_idx=1
declare -a TEMP_SOURCES=()
declare -a FINAL_TARGETS=()

echo "Planning reindexing of remaining good files to a continuous sequence..."

for ((itr=START_ITR; itr<=END_ITR; itr++)); do
  f="$(file_path "$itr")"
  [[ -f "$f" ]] || continue

  final="$(file_path "$good_idx")"
  if [[ "$f" == "$final" ]]; then
    # already in the right place
    (( good_idx++ ))
    continue
  fi

  # Stage 1: move source to a unique temp path (same dir)
  temp="${f}.reseq.$$"
  TEMP_SOURCES+=("$temp")
  FINAL_TARGETS+=("$final")

  if (( DRY_RUN )); then
    echo "[DRY-RUN] Would move to temp: '$f' -> '$temp'"
  else
    mv -f -- "$f" "$temp"
  fi

  (( good_idx++ ))
done

# -------- 3) Apply final renames from temp to target --------
for i in "${!TEMP_SOURCES[@]}"; do
  src="${TEMP_SOURCES[$i]}"
  dst="${FINAL_TARGETS[$i]}"

  if (( DRY_RUN )); then
    echo "[DRY-RUN] Would rename final: '$src' -> '$dst'"
  else
    # Ensure destination directory exists (it does) and overwrite any lingering placeholder
    mv -f -- "$src" "$dst"
    echo "Renamed: $(basename "$dst")"
  fi
done

echo "Reindexing complete."
if (( DRY_RUN )); then
  echo "[DRY-RUN] No changes were made."
fi
