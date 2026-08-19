"""
Proyecto: KPC_UNAL
Script: 03_train_models.py
Descripción: Entrenamiento con balanceo de clases (class_weight='balanced' y scale_pos_weight)
             y optimización por 5-Fold Stratified CV para mitigar desbalance clínico.
"""

import os
import warnings
import joblib
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from xgboost import XGBClassifier
from sklearn.model_selection import GridSearchCV, StratifiedKFold

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)


def train_elastic_net(X_train, y_train):
    print("[*] Optimizando Elastic Net (Balanced)...")
    param_grid = {
        'C': [0.01, 0.1, 1.0, 10.0],
        'l1_ratio': [0.1, 0.3, 0.5, 0.7]
    }
    base_model = LogisticRegression(
        solver='saga', 
        penalty='elasticnet', 
        class_weight='balanced', 
        max_iter=2000, 
        random_state=42
    )
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    grid = GridSearchCV(base_model, param_grid, cv=cv, scoring='roc_auc', n_jobs=-1)
    grid.fit(X_train, y_train)
    print(f"    - Mejor AUC-ROC (CV): {grid.best_score_:.4f} | Params: {grid.best_params_}")
    return grid.best_estimator_


def train_xgboost(X_train, y_train):
    print("[*] Optimizando XGBoost (Balanced scale_pos_weight)...")
    n_neg = (y_train == 0).sum()
    n_pos = (y_train == 1).sum()
    scale_weight = float(n_neg) / float(n_pos) if n_pos > 0 else 1.0

    param_grid = {
        'n_estimators': [50, 100],
        'max_depth': [3, 5],
        'learning_rate': [0.01, 0.1],
        'subsample': [0.8, 1.0]
    }
    base_model = XGBClassifier(
        eval_metric='logloss', 
        scale_pos_weight=scale_weight, 
        random_state=42
    )
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    grid = GridSearchCV(base_model, param_grid, cv=cv, scoring='roc_auc', n_jobs=-1)
    grid.fit(X_train, y_train)
    print(f"    - Mejor AUC-ROC (CV): {grid.best_score_:.4f} | Params: {grid.best_params_}")
    return grid.best_estimator_


def train_decision_tree(X_train, y_train):
    print("[*] Optimizando Decision Tree (Balanced)...")
    param_grid = {
        'max_depth': [3, 5, 7],
        'min_samples_leaf': [2, 5, 10]
    }
    base_model = DecisionTreeClassifier(class_weight='balanced', random_state=42)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    grid = GridSearchCV(base_model, param_grid, cv=cv, scoring='roc_auc', n_jobs=-1)
    grid.fit(X_train, y_train)
    print(f"    - Mejor AUC-ROC (CV): {grid.best_score_:.4f} | Params: {grid.best_params_}")
    return grid.best_estimator_


def run_training_pipeline(antibiotic: str, base_dir: str):
    print(f"\n================ ENTRENANDO MODELOS: {antibiotic.upper()} ================")
    data_dir = os.path.join(base_dir, "data", "processed", antibiotic)
    models_dir = os.path.join(base_dir, "results", "models", antibiotic)
    os.makedirs(models_dir, exist_ok=True)

    X_train = joblib.load(os.path.join(data_dir, "X_train.joblib"))
    y_train = joblib.load(os.path.join(data_dir, "y_train.joblib"))

    en_model = train_elastic_net(X_train, y_train)
    xgb_model = train_xgboost(X_train, y_train)
    dt_model = train_decision_tree(X_train, y_train)

    joblib.dump(en_model, os.path.join(models_dir, "elastic_net.joblib"))
    joblib.dump(xgb_model, os.path.join(models_dir, "xgboost.joblib"))
    joblib.dump(dt_model, os.path.join(models_dir, "decision_tree.joblib"))

    print(f"[✓] Modelos balanceados guardados en: {models_dir}")


if __name__ == "__main__":
    print("[✓] Modulo 03 optimizado con balanceo de clases.")
