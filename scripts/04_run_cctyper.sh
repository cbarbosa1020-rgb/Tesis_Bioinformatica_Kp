#!/usr/bin/env bash
# ==============================================================================
# 04_run_cctyper.sh
# Minería del Spacerome con cctyper v1.8.0
# ==============================================================================

set -euo pipefail

PASS_DIR="$HOME/data/fasta_pass"
OUTPUT_DIR="$HOME/data/cctyper_results"
THREADS=8

mkdir -p "$OUTPUT_DIR"

echo "==> Iniciando minería del Spacerome con cctyper v1.8.0..."
echo "==> Genomas a procesar: $PASS_DIR"
echo "==> Destino: $OUTPUT_DIR"

run_sample() {
    local fasta_path="$1"
    local sample_id
    sample_id=$(basename "$fasta_path" .fasta)
    local out_dir="$HOME/data/cctyper_results/$sample_id"

    if [ ! -d "$out_dir" ]; then
        source "$HOME/miniforge3/bin/activate" bioinfo
        cctyper "$fasta_path" "$out_dir" --threads 1 --no_plot > /dev/null 2>&1
    fi
}

export -f run_sample

find "$PASS_DIR" -type f -name "*.fasta" | \
    parallel -j "$THREADS" --progress run_sample {}

echo "==> Minería de CRISPRCasTyper completada."
