#!/usr/bin/env bash
# ==============================================================================
# 04_run_cctyper.sh
# Caracterizacion masiva del Spacerome y subtipos Cas con CRISPRCasTyper v1.8.0
# ==============================================================================

set -euo pipefail

PASS_DIR="$HOME/data/fasta_pass"
OUTPUT_DIR="$HOME/data/cctyper_results"
THREADS=8

mkdir -p "$OUTPUT_DIR"

echo "==> Iniciando mineración del Spacerome con cctyper v1.8.0..."
echo "==> Genomas a procesar: $PASS_DIR"
echo "==> Destino: $OUTPUT_DIR"

run_cctyper_single() {
    local fasta_file="$1"
    local base_name
    base_name=$(basename "$fasta_file")
    local sample_id="${base_name%.*}"
    local sample_out="$OUTPUT_DIR/${sample_id}"

    # Evitar re-procesar si ya existe la carpeta de la muestra
    if [ ! -d "$sample_out" ]; then
        cctyper "$fasta_file" "$sample_out" --threads 1 --no_plot > /dev/null 2>&1
    fi
}

export OUTPUT_DIR
export -f run_cctyper_single

# Ejecucion en paralelo a 8 hilos
find "$PASS_DIR" -type f -name "*.fasta" | \
    parallel -j "$THREADS" --progress run_cctyper_single {}

echo "==> Minería de CRISPRCasTyper completada."
