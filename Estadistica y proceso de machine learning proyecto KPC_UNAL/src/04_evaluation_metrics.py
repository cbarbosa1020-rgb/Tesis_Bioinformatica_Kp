"""
Proyecto: KPC_UNAL
Script: 04_evaluation_metrics.py
Descripción: Evaluación epidemiológica y bioestadística de modelos en el set de prueba.
             Genera tablas CSV con métricas y figuras de Curvas ROC y Matrices de Confusión.
"""

import os
import joblib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import (
    roc_auc_score,
    roc_curve,
    confusion_matrix,
    f1_score,
    accuracy_score
)


def evaluate_models(
    models_dict: dict, 
    X_test: pd.DataFrame, 
    y_test: pd.Series, 
    antibiotic: str,
    results_dir: str
) -> pd.DataFrame:
    """
    Calcula métricas clave en el test set y exporta gráficos diagnósticos.
    """
    metrics_list = []
    
    figures_dir = os.path.join(results_dir, "figures", antibiotic)
    tables_dir = os.path.join(results_dir, "tables")
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(tables_dir, exist_ok=True)

    plt.figure(figsize=(8, 6))

    for name, model in models_dict.items():
        y_pred = model.predict(X_test)
        
        # Probabilidad de clase resistente (1)
        if hasattr(model, "predict_proba"):
            y_proba = model.predict_proba(X_test)[:, 1]
            auc_val = roc_auc_score(y_test, y_proba)
            fpr, tpr, _ = roc_curve(y_test, y_proba)
            plt.plot(fpr, tpr, lw=2, label=f"{name} (AUC = {auc_val:.3f})")
        else:
            y_proba = None
            auc_val = np.nan

        # Matriz de confusión
        tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
        
        # Métricas epidemiológicas
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0  # Recall clínico
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0  # Capacidad de descartar sensibles
        f1 = f1_score(y_test, y_pred, zero_division=0)
        acc = accuracy_score(y_test, y_pred)

        metrics_list.append({
            "Antibiotic": antibiotic,
            "Model": name,
            "AUC_ROC": round(auc_val, 4),
            "Sensitivity": round(sensitivity, 4),
            "Specificity": round(specificity, 4),
            "F1_Score": round(f1, 4),
            "Accuracy": round(acc, 4),
            "TP": tp, "FP": fp, "TN": tn, "FN": fn
        })

        # Graficar Matriz de Confusión individual
        cm = np.array([[tn, fp], [fn, tp]])
        plt_cm = plt.figure(figsize=(5, 4))
        sns.heatmap(cm, annot=True, fmt="d", cmap="Blues", cbar=False,
                    xticklabels=["Sensible (0)", "Resistente (1)"],
                    yticklabels=["Sensible (0)", "Resistente (1)"])
        plt.title(f"Matriz de Confusión - {name} ({antibiotic.upper()})")
        plt.xlabel("Predicción")
        plt.ylabel("Fenotipo Real")
        plt.tight_layout()
        plt.savefig(os.path.join(figures_dir, f"confusion_matrix_{name}.png"), dpi=300)
        plt.close(plt_cm)

    # Finalizar Curva ROC comparativa
    plt.plot([0, 1], [0, 1], color="grey", linestyle="--", lw=1.5)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("1 - Especificidad (Tasa Falsos Positivos)")
    plt.ylabel("Sensibilidad (Tasa Verdaderos Positivos)")
    plt.title(f"Curvas ROC Comparativas - {antibiotic.upper()}")
    plt.legend(loc="lower right")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "roc_curves_comparison.png"), dpi=300)
    plt.close()

    metrics_df = pd.DataFrame(metrics_list)
    out_csv = os.path.join(tables_dir, f"metrics_{antibiotic}.csv")
    metrics_df.to_csv(out_csv, index=False)
    print(f"[✓] Métricas guardadas en: {out_csv}")
    
    return metrics_df


def run_evaluation_pipeline(antibiotic: str, base_dir: str):
    """
    Carga datos y modelos para ejecutar la evaluación diagnóstica completa.
    """
    print(f"\n================ EVALUANDO MODELOS: {antibiotic.upper()} ================")
    data_dir = os.path.join(base_dir, "data", "processed", antibiotic)
    models_dir = os.path.join(base_dir, "results", "models", antibiotic)
    results_dir = os.path.join(base_dir, "results")

    X_test = joblib.load(os.path.join(data_dir, "X_test.joblib"))
    y_test = joblib.load(os.path.join(data_dir, "y_test.joblib"))

    models = {
        "Elastic_Net": joblib.load(os.path.join(models_dir, "elastic_net.joblib")),
        "XGBoost": joblib.load(os.path.join(models_dir, "xgboost.joblib")),
        "Decision_Tree": joblib.load(os.path.join(models_dir, "decision_tree.joblib"))
    }

    metrics_df = evaluate_models(models, X_test, y_test, antibiotic, results_dir)
    print("\nResumen de Desempeño:")
    print(metrics_df[["Model", "AUC_ROC", "Sensitivity", "Specificity", "F1_Score"]].to_string(index=False))


if __name__ == "__main__":
    print("[✓] Modulo 04 (Model Evaluation & ROC Diagnostics) listo para integrarse.")
