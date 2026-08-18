# Estadística y Proceso de Machine Learning - Proyecto KPC UNAL

Pipeline integral de Machine Learning aplicado a la genómica de *Klebsiella pneumoniae* para la caracterización del espacioroma (espaciadores CRISPR y subtipos Cas) y su capacidad predictiva sobre fenotipos de resistencia antimicrobiana (AMR / MIC).

## Estructura del Proyecto

- `src/01_feature_engineering.py`: Extracción de espaciadores y subtipos Cas hacia matrices binarias sparse.
- `src/02_preprocessing_splits.py`: Filtrado de baja varianza y división estratificada (80% Train / 20% Test).
- `src/03_train_models.py`: Entrenamiento y optimización (5-Fold CV) para Elastic Net, XGBoost y Decision Trees.
- `src/04_evaluation_metrics.py`: Cálculo de AUC-ROC, Sensibilidad, Especificidad, F1-Score y generación de gráficos diagnósticos.
- `src/05_feature_importance.py`: Extracción de Odds Ratios y Tree SHAP values para identificación de biomarcadores genéticos.
- `main.py`: Orquestador principal de ejecución parametrizable.

## Uso

Ejecución para un antibiótico específico:
```bash
python main.py --antibiotic meropenem
