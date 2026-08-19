"""
Proyecto: KPC_UNAL
Script: main.py
Descripción: Pipeline orquestador maestro para el modelado de Machine Learning
             del espacioroma de K. pneumoniae.
"""

import os
import sys
import argparse
import importlib

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(CURRENT_DIR, "src")
if SRC_DIR not in sys.path:
    sys.path.append(SRC_DIR)

mod_01 = importlib.import_module("01_feature_engineering")
mod_02 = importlib.import_module("02_preprocessing_splits")
mod_03 = importlib.import_module("03_train_models")
mod_04 = importlib.import_module("04_evaluation_metrics")
mod_05 = importlib.import_module("05_feature_importance")

parse_crispr_spacers = mod_01.parse_crispr_spacers
parse_cas_subtypes = mod_01.parse_cas_subtypes
merge_features_and_metadata = mod_01.merge_features_and_metadata
run_preprocessing_pipeline = mod_02.run_preprocessing_pipeline
run_training_pipeline = mod_03.run_training_pipeline
run_evaluation_pipeline = mod_04.run_evaluation_pipeline
run_feature_importance_pipeline = mod_05.run_feature_importance_pipeline

PRIORITY_ANTIBIOTICS = [
    "meropenem", "imipenem", "ceftazidime", 
    "ciprofloxacin", "cefotaxime", "gentamicin", 
    "ampicillin", "norfloxacin", "trimethoprim_sulfamethoxazole"
]


def main():
    parser = argparse.ArgumentParser(description="Pipeline ML - Espacioroma KPC UNAL")
    parser.add_argument(
        "--antibiotic", 
        type=str, 
        default="all", 
        help="Antibiótico a evaluar o 'all' para el panel completo."
    )
    parser.add_argument(
        "--min_spacer_freq", 
        type=int, 
        default=5, 
        help="Frecuencia mínima de aparición del espaciador."
    )
    parser.add_argument(
        "--skip_feature_engineering", 
        action="store_true", 
        help="Omitir el paso 01 si las matrices ya existen."
    )

    args = parser.parse_args()

    raw_dir = os.path.join(CURRENT_DIR, "data", "raw")
    processed_dir = os.path.join(CURRENT_DIR, "data", "processed")
    os.makedirs(processed_dir, exist_ok=True)

    spacers_path = os.path.join(raw_dir, "spacers.tsv")
    cas_path = os.path.join(raw_dir, "cas_subtypes.tsv")
    metadata_path = os.path.join(raw_dir, "mic_metadata.tsv")
    
    genomic_matrix_path = os.path.join(processed_dir, "genomic_matrix.parquet")
    phenotypes_matrix_path = os.path.join(processed_dir, "phenotypes_matrix.parquet")

    antibiotics_to_run = PRIORITY_ANTIBIOTICS if args.antibiotic == "all" else [args.antibiotic.lower()]

    print("\n=======================================================")
    print("      PIPELINE MACHINE LEARNING - PROYECTO KPC UNAL    ")
    print("=======================================================\n")

    # 1. Feature Engineering & Matrices
    if not args.skip_feature_engineering:
        if os.path.exists(spacers_path) and os.path.exists(cas_path) and os.path.exists(metadata_path):
            print("[PASO 1] Construyendo matrices genómicas y fenotípicas...")
            df_spacers = parse_crispr_spacers(spacers_path, min_occurrence=args.min_spacer_freq)
            df_cas = parse_cas_subtypes(cas_path)
            X, y = merge_features_and_metadata(df_spacers, df_cas, metadata_path, PRIORITY_ANTIBIOTICS)
            X.to_parquet(genomic_matrix_path)
            y.to_parquet(phenotypes_matrix_path)
            print(f"[✓] Matriz genómica guardada en: {genomic_matrix_path}")
            print(f"[✓] Matriz fenotípica guardada en: {phenotypes_matrix_path}")
        else:
            print("[!] Faltan archivos en data/raw.")

    # 2 al 5 por cada antibiótico
    for ab in antibiotics_to_run:
        print(f"\n>>> INICIANDO MODELADO PARA: {ab.upper()} <<<")
        try:
            if not os.path.exists(genomic_matrix_path) or not os.path.exists(phenotypes_matrix_path):
                print("[!] No existen matrices procesadas.")
                continue

            run_preprocessing_pipeline(genomic_matrix_path, phenotypes_matrix_path, ab, processed_dir)
            run_training_pipeline(ab, CURRENT_DIR)
            run_evaluation_pipeline(ab, CURRENT_DIR)
            run_feature_importance_pipeline(ab, CURRENT_DIR)

        except Exception as e:
            print(f"[ERROR] Fallo en {ab}: {e}")

    print("\n[✓] Pipeline ejecutado.")


if __name__ == "__main__":
    main()
