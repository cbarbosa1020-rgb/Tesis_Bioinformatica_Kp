#!/usr/bin/env bash
# ==============================================================================
# 04_run_cctyper.sh
# Caracterización masiva del Spacerome con CRISPRCasTyper v1.8.0
# ==============================================================================

set -euo pipefail

PASS_DIR="$HOME/data/fasta_pass"
OUTPUT_DIR="$HOME/data/cctyper_results"
THREADS=8

mkdir -p "$OUTPUT_DIR"

echo "==> Iniciando minería del Spacerome con cctyper v1.8.0..."
echo "==> Genomas a procesar: $PASS_DIR"
echo "==> Destino: $OUTPUT_DIR"

# Ejecución directa en paralelo llamando a cctyper por archivo
find "$PASS_DIR" -type f -name "*.fasta" | \
    parallel -j "$THREADS" --progress '
        sample_id=$(basename {} .fasta)
        out_dir="'"$OUTPUT_DIR"'/$sample_id"
        if [ ! -d "$out_dir" ]; then
            cctyper {} "$out_dir" --threads 1 --no_plot > /dev/null 2>&1
        fi
    '

echo "==> Minería de CRISPRCasTyper completada."
