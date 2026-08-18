"""
Proyecto: KPC_UNAL
Script: 03_train_models.py
Descripción: Entrenamiento y optimización de hiperparámetros vía Stratified 5-Fold CV
             para Elastic Net, XGBoost y Decision Tree.
"""

import os
import joblib
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
import xgboost as xgb


def train_elastic_net(X_train: pd.DataFrame, y_train: pd.Series, cv) -> GridSearchCV:
    """
    Optimiza Regresión Logística con regularización Elastic Net (L1 + L2).
    """
    print("[*] Optimizando Elastic Net...")
    model = LogisticRegression(
        penalty="elasticnet",
        solver="saga",
        max_iter=3000,
        class_weight="balanced",
        random_state=42
    )
    param_grid = {
        "C": [0.01, 0.1, 1.0, 10.0],
        "l1_ratio": [0.1, 0.5, 0.7, 0.9]
    }
    grid = GridSearchCV(model, param_grid, cv=cv, scoring="roc_auc", n_jobs=-1, verbose=1)
    grid.fit(X_train, y_train)
    print(f"    - Mejor AUC-ROC (CV): {grid.best_score_:.4f} | Params: {grid.best_params_}")
    return grid.best_estimator_


def train_xgboost(X_train: pd.DataFrame, y_train: pd.Series, cv) -> GridSearchCV:
    """
    Optimiza XGBoost Classifier con ponderación de clases positivas.
    """
    print("[*] Optimizando XGBoost...")
    neg_count = (y_train == 0).sum()
    pos_count = (y_train == 1).sum()
    scale_weight = neg_count / max(pos_count, 1)

    model = xgb.XGBClassifier(
        scale_pos_weight=scale_weight,
        eval_metric="logloss",
        random_state=42,
        tree_method="hist"
    )
    param_grid = {
        "max_depth": [3, 5, 7],
        "learning_rate": [0.01, 0.05, 0.1],
        "n_estimators": [100, 200],
        "subsample": [0.8]
    }
    grid = GridSearchCV(model, param_grid, cv=cv, scoring="roc_auc", n_jobs=-1, verbose=1)
    grid.fit(X_train, y_train)
    print(f"    - Mejor AUC-ROC (CV): {grid.best_score_:.4f} | Params: {grid.best_params_}")
    return grid.best_estimator_


def train_decision_tree(X_train: pd.DataFrame, y_train: pd.Series, cv) -> GridSearchCV:
    """
    Optimiza Árbol de Decisión como baseline interpretable.
    """
    print("[*] Optimizando Decision Tree...")
    model = DecisionTreeClassifier(class_weight="balanced", random_state=42)
    param_grid = {
        "max_depth": [3, 5, 8, 12],
        "min_samples_leaf": [5, 10, 20]
    }
    grid = GridSearchCV(model, param_grid, cv=cv, scoring="roc_auc", n_jobs=-1, verbose=1)
    grid.fit(X_train, y_train)
    print(f"    - Mejor AUC-ROC (CV): {grid.best_score_:.4f} | Params: {grid.best_params_}")
    return grid.best_estimator_


def run_training_pipeline(antibiotic: str, base_dir: str):
    """
    Carga datos particionados, entrena todos los modelos y guarda los estimadores serializados.
    """
    print(f"\n================ ENTRENANDO MODELOS: {antibiotic.upper()} ================")
    data_dir = os.path.join(base_dir, "data", "processed", antibiotic)
    models_dir = os.path.join(base_dir, "results", "models", antibiotic)
    os.makedirs(models_dir, exist_ok=True)

    X_train = joblib.load(os.path.join(data_dir, "X_train.joblib"))
    y_train = joblib.load(os.path.join(data_dir, "y_train.joblib"))

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    trained_models = {
        "elastic_net": train_elastic_net(X_train, y_train, cv),
        "xgboost": train_xgboost(X_train, y_train, cv),
        "decision_tree": train_decision_tree(X_train, y_train, cv)
    }

    for name, model in trained_models.items():
        save_path = os.path.join(models_dir, f"{name}.joblib")
        joblib.dump(model, save_path)
        print(f"[✓] Modelo {name} guardado en: {save_path}")


if __name__ == "__main__":
    print("[✓] Modulo 03 (Model Training & Hyperparameter Tuning) listo para integrarse.")
