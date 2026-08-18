"""
Proyecto: KPC_UNAL
Script: 05_feature_importance.py
Descripción: Extracción de biomarcadores genómicos (espaciadores y subtipos Cas)
             mediante Odds Ratios (Elastic Net) y Tree SHAP Values (XGBoost).
"""

import os
import joblib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shap


def explain_elastic_net(model, feature_names: list, top_n: int = 20) -> pd.DataFrame:
    """
    Extrae coeficientes y calcula Odds Ratios (OR = exp(beta)) para Elastic Net.
    OR > 1: Predictor de Resistencia.
    OR < 1: Predictor de Susceptibilidad (Inmunidad CRISPR contra MGEs/plásmidos).
    """
    coefs = model.coef_[0]
    df = pd.DataFrame({
        "Feature": feature_names,
        "Coefficient": coefs,
        "Odds_Ratio": np.exp(coefs)
    })
    
    # Filtrar solo características con coeficientes distintos de cero (seleccionadas por L1)
    df_nonzero = df[df["Coefficient"] != 0].copy()
    df_nonzero["Abs_Effect"] = df_nonzero["Coefficient"].abs()
    
    return df_nonzero.sort_values(by="Abs_Effect", ascending=False).head(top_n).drop(columns=["Abs_Effect"])


def explain_xgboost_shap(model, X_test: pd.DataFrame, figures_dir: str, antibiotic: str, top_n: int = 20) -> pd.DataFrame:
    """
    Calcula valores SHAP exactos usando TreeExplainer para capturar interacciones no lineales.
    """
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X_test)

    # SHAP Summary Plot
    plt.figure(figsize=(10, 8))
    shap.summary_plot(shap_values, X_test, max_display=top_n, show=False)
    plt.title(f"SHAP Values (Impacto en Resistencia) - {antibiotic.upper()}", fontsize=12)
    plt.tight_layout()
    plot_path = os.path.join(figures_dir, f"shap_summary_{antibiotic}.png")
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close()

    # Tabla resumen por importancia media absoluta
    mean_abs_shap = np.abs(shap_values).mean(axis=0)
    shap_df = pd.DataFrame({
        "Feature": X_test.columns,
        "Mean_Abs_SHAP": mean_abs_shap
    }).sort_values(by="Mean_Abs_SHAP", ascending=False).head(top_n)

    return shap_df


def run_feature_importance_pipeline(antibiotic: str, base_dir: str, top_n: int = 25):
    """
    Ejecuta el análisis de importancia de características y exporta tablas y gráficos.
    """
    print(f"\n================ EXTRACCIÓN DE BIOMARCADORES: {antibiotic.upper()} ================")
    data_dir = os.path.join(base_dir, "data", "processed", antibiotic)
    models_dir = os.path.join(base_dir, "results", "models", antibiotic)
    tables_dir = os.path.join(base_dir, "results", "tables")
    figures_dir = os.path.join(base_dir, "results", "figures", antibiotic)

    os.makedirs(tables_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)

    X_train = joblib.load(os.path.join(data_dir, "X_train.joblib"))
    X_test = joblib.load(os.path.join(data_dir, "X_test.joblib"))

    en_model = joblib.load(os.path.join(models_dir, "elastic_net.joblib"))
    xgb_model = joblib.load(os.path.join(models_dir, "xgboost.joblib"))

    feature_names = list(X_train.columns)

    # 1. Odds Ratios (Elastic Net)
    en_biomarkers = explain_elastic_net(en_model, feature_names, top_n=top_n)
    en_csv = os.path.join(tables_dir, f"biomarkers_elasticnet_{antibiotic}.csv")
    en_biomarkers.to_csv(en_csv, index=False)
    print(f"[✓] Top Odds Ratios guardados en: {en_csv}")

    # 2. SHAP Values (XGBoost)
    xgb_biomarkers = explain_xgboost_shap(xgb_model, X_test, figures_dir, antibiotic, top_n=top_n)
    xgb_csv = os.path.join(tables_dir, f"biomarkers_shap_{antibiotic}.csv")
    xgb_biomarkers.to_csv(xgb_csv, index=False)
    print(f"[✓] Top SHAP Features guardados en: {xgb_csv}")


if __name__ == "__main__":
    print("[✓] Modulo 05 (Biological Feature Importance & SHAP) listo para integrarse.")
