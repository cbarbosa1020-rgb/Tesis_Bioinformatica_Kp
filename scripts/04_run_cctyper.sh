#!/usr/bin/env bash
# ==============================================================================
# 04_run_cctyper.sh
# Minería del Spacerome con cctyper v1.8.0
# ==============================================================================

set -euo pipefail

PASS_DIR="$HOME/data/fasta_pass"
OUTPUT_DIR="$HOME/data/cctyper_results"
PARALLEL_BIN="$HOME/miniforge3/envs/bioinfo/bin/parallel"
CCTYPER_BIN="$HOME/miniforge3/envs/bioinfo/bin/cctyper"
DB_DIR="$HOME/miniforge3/envs/bioinfo/cct_data"
THREADS=8

mkdir -p "$OUTPUT_DIR"

echo "==> Iniciando minería del Spacerome con cctyper v1.8.0..."

export CCTYPER_BIN DB_DIR OUTPUT_DIR

find "$PASS_DIR" -type f -name "*.fasta" | \
    "$PARALLEL_BIN" -j "$THREADS" --progress '
        sample=$(basename {} .fasta)
        out="'"$OUTPUT_DIR"'/$sample"
        if [ ! -f "$out/cas_operons.tab" ]; then
            '"$CCTYPER_BIN"' {} "$out" --db '"$DB_DIR"' --threads 1 --no_plot > /dev/null 2>&1
        fi
    '

echo "==> Minería de CRISPRCasTyper completada."
