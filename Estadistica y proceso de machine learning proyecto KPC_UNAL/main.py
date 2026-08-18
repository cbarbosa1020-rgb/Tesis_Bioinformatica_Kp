"""
Proyecto: KPC_UNAL
Script: main.py
Descripción: Pipeline orquestador maestro para el análisis estadístico y
             modelado de Machine Learning del espacioroma de K. pneumoniae.
"""

import os
import sys
import argparse

# Asegurar importación de módulos en src/
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(CURRENT_DIR, "src"))

from 01_feature_engineering import parse_crispr_spacers, parse_cas_subtypes, merge_features_and_metadata
from 02_preprocessing_splits import run_preprocessing_pipeline
from 03_train_models import run_training_pipeline
from 04_evaluation_metrics import run_evaluation_pipeline
from 05_feature_importance import run_feature_importance_pipeline

PRIORITY_ANTIBIOTICS = [
    "meropenem", "imipenem", "ertapenem",
    "ceftazidime", "ceftriaxone", "cefepime",
    "ceftazidime_avibactam", "piperacillin_tazobactam",
    "colistin", "amikacin", "ciprofloxacin"
]


def main():
    parser = argparse.ArgumentParser(description="Pipeline ML - Espacioroma KPC UNAL")
    parser.add_argument(
        "--antibiotic", 
        type=str, 
        default="all", 
        help="Antibiótico a evaluar (ej: 'meropenem') o 'all' para el panel completo."
    )
    parser.add_argument(
        "--min_spacer_freq", 
        type=int, 
        default=5, 
        help="Frecuencia mínima de aparición de un espaciador para incluirlo en la matriz."
    )
    parser.add_argument(
        "--skip_feature_engineering", 
        action="store_true", 
        help="Omitir el paso 01 si la matriz ya fue construida."
    )

    args = parser.parse_args()

    raw_dir = os.path.join(CURRENT_DIR, "data", "raw")
    processed_dir = os.path.join(CURRENT_DIR, "data", "processed")
    spacers_path = os.path.join(raw_dir, "spacers.tsv")
    cas_path = os.path.join(raw_dir, "cas_subtypes.tsv")
    metadata_path = os.path.join(raw_dir, "mic_metadata.tsv")
    genomic_matrix_path = os.path.join(processed_dir, "genomic_matrix.parquet")

    antibiotics_to_run = PRIORITY_ANTIBIOTICS if args.antibiotic == "all" else [args.antibiotic.lower()]

    print("\n=======================================================")
    print("      PIPELINE MACHINE LEARNING - PROYECTO KPC UNAL    ")
    print("=======================================================\n")

    # Paso 1: Feature Engineering
    if not args.skip_feature_engineering:
        if os.path.exists(spacers_path) and os.path.exists(cas_path) and os.path.exists(metadata_path):
            print("[PASO 1] Construyendo matriz de características genómicas...")
            df_spacers = parse_crispr_spacers(spacers_path, min_occurrence=args.min_spacer_freq)
            df_cas = parse_cas_subtypes(cas_path)
            X, y = merge_features_and_metadata(df_spacers, df_cas, metadata_path, PRIORITY_ANTIBIOTICS)
            X.to_parquet(genomic_matrix_path)
            print(f"[✓] Matriz guardada en: {genomic_matrix_path}")
        else:
            print("[!] Archivos en data/raw no detectados. Asumiendo matrices procesadas existentes.")

    # Pasos 2 al 5 por cada antibiótico
    for ab in antibiotics_to_run:
        print(f"\n>>> INICIANDO PROCESAMIENTO COMPLETO PARA: {ab.upper()} <<<")
        try:
            # 2. Preprocesamiento & Splits
            run_preprocessing_pipeline(genomic_matrix_path, metadata_path, ab, processed_dir)

            # 3. Entrenamiento & Hiperparámetros
            run_training_pipeline(ab, CURRENT_DIR)

            # 4. Evaluación Epidemiológica y Curvas ROC
            run_evaluation_pipeline(ab, CURRENT_DIR)

            # 5. Biomarcadores y SHAP
            run_feature_importance_pipeline(ab, CURRENT_DIR)

        except Exception as e:
            print(f"[ERROR] Fallo en el procesamiento de {ab}: {e}")

    print("\n[✓] Ejecución finalizada. Revisa la carpeta 'results/' para figuras y tablas.")


if __name__ == "__main__":
    main()
