#!/usr/bin/env bash
# ==============================================================================
# 04_run_cctyper.sh
# Minería del Spacerome con cctyper v1.8.0
# ==============================================================================

set -euo pipefail

ENV_DIR="$HOME/miniforge3/envs/bioinfo"
ENV_BIN="$ENV_DIR/bin"

# Exportar el PATH globalmente para que cctyper y sus subprocesos (prodigal/hmmsearch) respondan
export PATH="$ENV_BIN:$PATH"

PASS_DIR="$HOME/data/fasta_pass"
OUTPUT_DIR="$HOME/data/cctyper_results"
PARALLEL_BIN="$ENV_BIN/parallel"
CCTYPER_BIN="$ENV_BIN/cctyper"
DB_DIR="$ENV_DIR/cct_data"
THREADS=8

mkdir -p "$OUTPUT_DIR"

echo "==> Iniciando minería del Spacerome con cctyper v1.8.0..."

find "$PASS_DIR" -type f -name "*.fasta" | \
    "$PARALLEL_BIN" -j "$THREADS" '
        sample=$(basename {} .fasta)
        out="'"$OUTPUT_DIR"'/$sample"
        if [ ! -f "$out/hmmer.log" ]; then
            '"$CCTYPER_BIN"' {} "$out" --db '"$DB_DIR"' --threads 1 --no_plot > /dev/null 2>&1
        fi
    '

echo "==> Minería de CRISPRCasTyper completada."
