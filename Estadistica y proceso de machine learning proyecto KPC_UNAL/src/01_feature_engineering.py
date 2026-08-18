"""
Proyecto: KPC_UNAL
Script: 01_feature_engineering.py
Descripción: Procesa las salidas tabulares de CRISPRCasTyper y las combina
             con los metadatos fenotípicos de MIC para generar la matriz de diseño (X)
             y los vectores de etiquetas (y) por antibiótico.
"""

import os
import argparse
import pandas as pd
import numpy as np


def parse_crispr_spacers(spacers_path: str, min_occurrence: int = 5) -> pd.DataFrame:
    """
    Lee la tabla de espaciadores de CRISPRCasTyper y crea una matriz binaria (0/1).
    Filtra espaciadores con frecuencia menor a min_occurrence.
    """
    print(f"[*] Cargando espaciadores desde: {spacers_path}")
    spacers_df = pd.read_csv(spacers_path, sep="\t")

    sample_col = "sample" if "sample" in spacers_df.columns else "genome_id"
    spacer_col = "spacer_id" if "spacer_id" in spacers_df.columns else "Spacer"

    spacer_counts = spacers_df.groupby(spacer_col)[sample_col].nunique()
    valid_spacers = spacer_counts[spacer_counts >= min_occurrence].index
    print(f"[+] Espaciadores totales: {len(spacer_counts)} | Filtrados (>= {min_occurrence}): {len(valid_spacers)}")

    filtered_df = spacers_df[spacers_df[spacer_col].isin(valid_spacers)]

    spacer_matrix = pd.crosstab(
        index=filtered_df[sample_col], columns=filtered_df[spacer_col]
    ).clip(upper=1).astype(np.int8)

    return spacer_matrix.add_prefix("Spacer_")


def parse_cas_subtypes(cas_path: str) -> pd.DataFrame:
    """
    Lee la tabla de subtipos Cas de CRISPRCasTyper y genera la matriz binaria correspondiente.
    """
    print(f"[*] Cargando subtipos Cas desde: {cas_path}")
    cas_df = pd.read_csv(cas_path, sep="\t")

    sample_col = "sample" if "sample" in cas_df.columns else "genome_id"
    subtype_col = "subtype" if "subtype" in cas_df.columns else "cas_subtype"

    cas_matrix = pd.crosstab(
        index=cas_df[sample_col], columns=cas_df[subtype_col]
    ).clip(upper=1).astype(np.int8)

    return cas_matrix.add_prefix("CasType_")


def merge_features_and_metadata(
    spacer_matrix: pd.DataFrame,
    cas_matrix: pd.DataFrame,
    metadata_path: str,
    target_antibiotics: list,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Une las matrices genómicas con las etiquetas 
MIC binarizadas.
    """
    print(f"[*] Integrando con metadatos fenotípicos: {metadata_path}")
    meta_df = pd.read_csv(metadata_path, sep="\t", index_col=0)

    X_genomic = spacer_matrix.join(cas_matrix, how="outer").fillna(0).astype(np.int8)
    common_samples = X_genomic.index.intersection(meta_df.index)
    print(f"[+] Genomas coincidentes: {len(common_samples)}")

    X_final = X_genomic.loc[common_samples]
    y_final = meta_df.loc[common_samples, target_antibiotics]

    return X_final, y_final


if __name__ == "__main__":
    PRIORITY_ANTIBIOTICS = [
        "meropenem", "imipenem", "ertapenem",
        "ceftazidime", "ceftriaxone", "cefepime",
        "ceftazidime_avibactam", "piperacillin_tazobactam",
        "colistin", "amikacin", "ciprofloxacin"
    ]
    print("[✓] Modulo 01 cargado correctamente.")

