"""
Proyecto: KPC_UNAL
Script: 01_feature_engineering.py
Descripción: Construcción de matriz sparse genómica e integración con metadatos CMI.
             Extrae y normaliza identificadores taxonómicos tipo '573.XXXXX'.
"""

import re
import pandas as pd
import numpy as np


def clean_id(val) -> str:
    """Extrae el ID genómico unificado estándar (573.XXXXX)."""
    s = str(val).strip()
    match = re.search(r'(573\.\d+)', s)
    if match:
        return match.group(1)
    
    # Fallback para identificadores por accession o contig
    s = s.replace("accn|", "").replace("accn_", "")
    if '_' in s:
        s = s.split('_')[0]
    return s


def parse_crispr_spacers(spacers_file: str, min_occurrence: int = 5) -> pd.DataFrame:
    """Carga y binariza presencia de clusters de espaciadores CRISPR."""
    print(f"[*] Cargando espaciadores desde: {spacers_file}")
    sep = ',' if spacers_file.endswith('.csv') else '\t'
    df = pd.read_csv(spacers_file, sep=sep, dtype=str)
    
    col_sample = [c for c in df.columns if c.lower() in ['genome_id', 'sample', 'isolate']][0]
    col_spacer = [c for c in df.columns if c.lower() in ['spacer_id', 'spacer', 'cluster', 'cluster_id']][0]
    
    df['genome_id'] = df[col_sample].apply(clean_id)
    df['spacer_id'] = df[col_spacer].astype(str)
    
    spacer_counts = df['spacer_id'].value_counts()
    frequent_spacers = spacer_counts[spacer_counts >= min_occurrence].index
    df_filtered = df[df['spacer_id'].isin(frequent_spacers)]
    
    print(f"[+] Espaciadores totales: {len(spacer_counts)} | Filtrados (>= {min_occurrence}): {len(frequent_spacers)}")
    
    matrix = pd.crosstab(df_filtered['genome_id'], df_filtered['spacer_id'])
    matrix = (matrix > 0).astype(int)
    matrix.columns = [f"Spacer_{c}" for c in matrix.columns]
    return matrix


def parse_cas_subtypes(cas_file: str) -> pd.DataFrame:
    """Carga y binariza subtipos Cas de CCTyper."""
    print(f"[*] Cargando subtipos Cas desde: {cas_file}")
    sep = ',' if cas_file.endswith('.csv') else '\t'
    df = pd.read_csv(cas_file, sep=sep, dtype=str)
    
    col_sample = [c for c in df.columns if c.lower() in ['genome_id', 'sample', 'isolate']][0]
    col_cas = [c for c in df.columns if c.lower() in ['cas_subtype', 'subtype', 'prediction', 'type']][0]
    
    df['genome_id'] = df[col_sample].apply(clean_id)
    df['cas_subtype'] = df[col_cas].astype(str)
    
    matrix = pd.crosstab(df['genome_id'], df['cas_subtype'])
    matrix = (matrix > 0).astype(int)
    matrix.columns = [f"Cas_{c}" for c in matrix.columns]
    return matrix


def merge_features_and_metadata(
    df_spacers: pd.DataFrame, 
    df_cas: pd.DataFrame, 
    metadata_file: str, 
    priority_antibiotics: list
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Integra características genómicas con la matriz fenotípica pivotada."""
    print(f"[*] Integrando con metadatos fenotípicos: {metadata_file}")
    
    meta = pd.read_csv(metadata_file, dtype=str)
    meta.columns = [c.strip().lower() for c in meta.columns]
    
    col_genome = [c for c in meta.columns if c in ['genome_id', 'sample', 'isolate']][0]
    meta['genome_id'] = meta[col_genome].apply(clean_id)
    
    col_ab = [c for c in meta.columns if 'antibiotic' in c][0]
    col_res = [c for c in meta.columns if 'resistant' in c or 'phenotype' in c or 'measurement' in c][0]
    
    meta[col_ab] = meta[col_ab].str.lower().str.replace('-', '_').str.replace(' ', '_').str.replace('/', '_')
    
    def binarize_phenotype(v):
        if pd.isna(v): return np.nan
        v_str = str(v).strip().lower()
        if v_str in ['resistant', 'r', '1', '1.0', 'true']: return 1
        if v_str in ['susceptible', 'sensitive', 's', '0', '0.0', 'false']: return 0
        try:
            val_float = float(v_str.replace('<=', '').replace('>=', '').replace('<', '').replace('>', ''))
            return 1 if val_float >= 4.0 else 0
        except ValueError:
            return np.nan

    meta['target'] = meta[col_res].apply(binarize_phenotype)
    
    # Pivotar formato largo a ancho
    meta_pivot = meta.pivot_table(index='genome_id', columns=col_ab, values='target', aggfunc='max')
    
    # Consolidar matriz genómica
    X_all = df_spacers.join(df_cas, how='outer').fillna(0).astype(int)
    
    common_samples = X_all.index.intersection(meta_pivot.index)
    print(f"[+] Genomas coincidentes (Genómica + Fenotipo): {len(common_samples)}")

    X_final = X_all.loc[common_samples]
    y_final = meta_pivot.loc[common_samples]
    
    return X_final, y_final


if __name__ == "__main__":
    print("[✓] Modulo 01 (Feature Engineering) optimizado.")
