"""
Proyecto: KPC_UNAL
Script: 02_preprocessing_splits.py
Descripción: Control de calidad de la matriz de espaciadores, filtrado por varianza
             y división estratificada (Train 80% / Test 20%) para cada antibiótico.
"""

import os
import joblib
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import VarianceThreshold


def filter_low_variance_features(X: pd.DataFrame, threshold: float = 0.005) -> pd.DataFrame:
    """
    Elimina características cuasi-constantes (varianza < threshold).
    Para variables binarias, Var = p * (1 - p).
    Un threshold de 0.005 remueve variables con p < 0.005 o p > 0.995.
    """
    selector = VarianceThreshold(threshold=threshold)
    selector.fit(X)
    
    retained_cols = X.columns[selector.get_support()]
    print(f"[+] Variables iniciales: {X.shape[1]} | Retenidas post-filtro varianza: {len(retained_cols)}")
    return X[retained_cols]


def create_stratified_partitions(
    X: pd.DataFrame, 
    y: pd.Series, 
    test_size: float = 0.20, 
    random_state: int = 42
) -> tuple:
    """
    Genera particiones balanceadas de Train y Test manteniendo la proporción de resistencia.
    Descarta genomas sin fenotipo etiquetado (NaN).
    """
    # Filtrar muestras sin dato fenotípico
    valid_idx = y.dropna().index
    X_clean = X.loc[valid_idx]
    y_clean = y.loc[valid_idx].astype(int)
    
    print(f"[*] Muestras con fenotipo: {len(y_clean)}")
    print(f"    - Resistentes (1): {y_clean.sum()} ({y_clean.mean()*100:.2f}%)")
    print(f"    - Sensibles   (0): {len(y_clean) - y_clean.sum()} ({(1 - y_clean.mean())*100:.2f}%)")
    
    X_train, X_test, y_train, y_test = train_test_split(
        X_clean, 
        y_clean, 
        test_size=test_size, 
        stratify=y_clean, 
        random_state=random_state
    )
    
    return X_train, X_test, y_train, y_test


def run_preprocessing_pipeline(
    feature_matrix_path: str,
    metadata_path: str,
    target_antibiotic: str,
    output_dir: str
):
    """
    Ejecuta el pipeline de limpieza y partición guardando los arrays serializados (.joblib).
    """
    print(f"\n================ Procesando: {target_antibiotic.upper()} ================")
    if feature_matrix_path.endswith('.parquet'):
        X = pd.read_parquet(feature_matrix_path)
    else:
        X = pd.read_csv(feature_matrix_path, index_col=0)
        
    meta = pd.read_csv(metadata_path, sep='\t', index_col=0)
    
    if target_antibiotic not in meta.columns:
        raise ValueError(f"El antibiótico '{target_antibiotic}' no existe en la tabla de metadatos.")
        
    y = meta[target_antibiotic]
    
    # 1. Filtrado de baja varianza
    X_filtered = filter_low_variance_features(X)
    
    # 2. Split estratificado
    X_train, X_test, y_train, y_test = create_stratified_partitions(X_filtered, y)
    
    # 3. Guardar particiones en disco
    ab_out_dir = os.path.join(output_dir, target_antibiotic)
    os.makedirs(ab_out_dir, exist_ok=True)
    
    joblib.dump(X_train, os.path.join(ab_out_dir, "X_train.joblib"))
    joblib.dump(X_test, os.path.join(ab_out_dir, "X_test.joblib"))
    joblib.dump(y_train, os.path.join(ab_out_dir, "y_train.joblib"))
    joblib.dump(y_test, os.path.join(ab_out_dir, "y_test.joblib"))
    
    print(f"[✓] Particiones guardadas exitosamente en: {ab_out_dir}")


if __name__ == "__main__":
    print("[✓] Modulo 02 (Preprocessing & Stratified Splits) listo para integrarse.")

