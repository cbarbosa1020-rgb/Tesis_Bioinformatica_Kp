#!/usr/bin/env bash
# ==============================================================================
# 04_run_cctyper.sh
# Mineria del Spacerome con cctyper v1.8.0 usando ruta explicita del entorno
# ==============================================================================

set -euo pipefail

PASS_DIR="$HOME/data/fasta_pass"
OUTPUT_DIR="$HOME/data/cctyper_results"
CCTYPER_BIN="$HOME/miniforge3/envs/bioinfo/bin/cctyper"
THREADS=8

mkdir -p "$OUTPUT_DIR"

echo "==> Iniciando mineria del Spacerome con cctyper v1.8.0..."
echo "==> Genomas a procesar: $PASS_DIR"
echo "==> Destino: $OUTPUT_DIR"

export OUTPUT_DIR
export CCTYPER_BIN

find "$PASS_DIR" -type f -name "*.fasta" | \
    parallel -j "$THREADS" --progress '
        sample_id=$(basename {} .fasta)
        out_dir="'"$OUTPUT_DIR"'/$sample_id"
        if [ ! -d "$out_dir" ]; then
            '"$CCTYPER_BIN"' {} "$out_dir" --threads 1 --no_plot > /dev/null 2>&1
        fi
    '

echo "==> Mineria de CRISPRCasTyper completada."
