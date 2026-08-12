#!/usr/bin/env bash
# ==============================================================================
# 02_run_quast.sh
# Ejecución masiva por lotes de QUAST v5.2.0 para 11,745 genomas
# ==============================================================================

set -euo pipefail

CLEAN_DIR="$HOME/data/fasta_clean"
OUTPUT_DIR="$HOME/data/quast_results"
THREADS=8
BATCH_SIZE=500

mkdir -p "$OUTPUT_DIR"

echo "==> Iniciando análisis con QUAST v5.2.0..."
echo "==> Muestras de origen: $CLEAN_DIR"
echo "==> Carpeta de salida: $OUTPUT_DIR"

# Obtener la lista completa de archivos FASTA
mapfile -t ALL_FASTAS < <(find "$CLEAN_DIR" -type f -name "*.fasta")
TOTAL_FILES=${#ALL_FASTAS[@]}

echo "Total de archivos a procesar: $TOTAL_FILES"

# Procesar en lotes para cuidar el uso de memoria RAM
BATCH_NUM=1
for ((i=0; i<TOTAL_FILES; i+=BATCH_SIZE)); do
    BATCH_FILES=("${ALL_FASTAS[@]:i:BATCH_SIZE}")
    BATCH_OUT="$OUTPUT_DIR/batch_${BATCH_NUM}"
    
    echo "--> Procesando Lote $BATCH_NUM ($((i+1)) a $((i+${#BATCH_FILES[@]})) de $TOTAL_FILES)..."
    
    quast.py "${BATCH_FILES[@]}" \
        -o "$BATCH_OUT" \
        -t "$THREADS" \
        --min-contig 0 \
        --no-plots \
        --no-html \
        --silent
    
    ((BATCH_NUM++))
done

echo "==> Análisis de QUAST completado en todos los lotes."
