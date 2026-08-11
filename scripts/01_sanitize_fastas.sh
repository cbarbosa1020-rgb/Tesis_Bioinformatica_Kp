#!/usr/bin/env bash
# ==============================================================================
# 01_sanitize_fastas.sh
# Estandarizacion de encabezados FASTA para genomas de K. pneumoniae
# ==============================================================================

set -euo pipefail

RAW_DIR="$HOME/data/fasta_raw"
CLEAN_DIR="$HOME/data/fasta_clean"
THREADS=8

mkdir -p "$CLEAN_DIR"

echo "==> Iniciando sanitización de FASTAs con $THREADS hilos..."

sanitize_single_fasta() {
    local fasta_file="$1"
    local base_name
    base_name=$(basename "$fasta_file")
    local sample_id="${base_name%.*}"
    local output_file="$CLEAN_DIR/${sample_id}.fasta"

    awk '/^>/ {print $1; next} {print}' "$fasta_file" > "$output_file"
}

export CLEAN_DIR
export -f sanitize_single_fasta

find "$RAW_DIR" -type f \( -name "*.fasta" -o -name "*.fna" -o -name "*.fa" \) | \
    parallel -j "$THREADS" --progress sanitize_single_fasta {}

echo "==> Sanitización completada."
