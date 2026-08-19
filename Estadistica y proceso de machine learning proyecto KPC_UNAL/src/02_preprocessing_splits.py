"""
Proyecto: KPC_UNAL
Script: 02_preprocessing_splits.py
Descripción: Filtrado de baja varianza y división estratificada (80% Train / 20% Test)
             a partir de las matrices genómica y fenotípica consolidadas.
"""

import os
import joblib
import pandas as pd
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import train_test_split


def run_preprocessing_pipeline(
    genomic_matrix_path: str,
    phenotypes_matrix_path: str,
    antibiotic: str,
    output_dir: str,
    test_size: float = 0.2,
    random_state: int = 42
):
    print(f"\n================ PREPROCESAMIENTO: {antibiotic.upper()} ================")
    
    X = pd.read_parquet(genomic_matrix_path)
    y_df = pd.read_parquet(phenotypes_matrix_path)

    if antibiotic not in y_df.columns:
        raise ValueError(f"El antibiótico '{antibiotic}' no está en las columnas disponibles: {list(y_df.columns)}")

    # Filtrar aislados con dato fenotípico válido para este antibiótico
    y = y_df[antibiotic].dropna()
    common_idx = X.index.intersection(y.index)

    X_ab = X.loc[common_idx]
    y_ab = y.loc[common_idx].astype(int)

    print(f"[+] Aislados evaluados para {antibiotic.upper()}: {len(common_idx)}")
    print(f"    - Resistentes (1): {(y_ab == 1).sum()} ({(y_ab == 1).mean() * 100:.1f}%)")
    print(f"    - Sensibles   (0): {(y_ab == 0).sum()} ({(y_ab == 0).mean() * 100:.1f}%)")

    if len(y_ab.unique()) < 2:
        raise ValueError(f"Se requieren ambas clases (0 y 1) para entrenar. Clases presentes: {y_ab.unique()}")

    # Filtrar variables constantes / sin varianza
    selector = VarianceThreshold(threshold=0.0)
    X_filtered = selector.fit_transform(X_ab)
    retained_features = X_ab.columns[selector.get_support()]
    X_filtered_df = pd.DataFrame(X_filtered, index=X_ab.index, columns=retained_features)
    
    print(f"[+] Predictores con varianza informativa: {X_filtered_df.shape[1]} de {X_ab.shape[1]}")

    # Split Estratificado
    X_train, X_test, y_train, y_test = train_test_split(
        X_filtered_df, y_ab, test_size=test_size, stratify=y_ab, random_state=random_state
    )

    ab_output_dir = os.path.join(output_dir, antibiotic)
    os.makedirs(ab_output_dir, exist_ok=True)

    joblib.dump(X_train, os.path.join(ab_output_dir, "X_train.joblib"))
    joblib.dump(X_test, os.path.join(ab_output_dir, "X_test.joblib"))
    joblib.dump(y_train, os.path.join(ab_output_dir, "y_train.joblib"))
    joblib.dump(y_test, os.path.join(ab_output_dir, "y_test.joblib"))

    print(f"[✓] Splits guardados en {ab_output_dir} (Train: {len(X_train)}, Test: {len(X_test)})")


if __name__ == "__main__":
    print("[✓] Modulo 02 (Preprocessing & Stratified Splits) listo.")
