"""
Proyecto: KPC_UNAL
Script: 06_generate_figures.py
Descripción: Genera figuras integradoras y paneles multipanel de alta resolución
             (300 DPI) para publicación y tesis.
"""

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
sns.set_context("paper", font_scale=1.2)


def plot_overall_performance(tables_dir: str, figures_dir: str):
    """Compila y grafica el desempeño de todos los antibióticos evaluados."""
    metric_files = glob.glob(os.path.join(tables_dir, "metrics_*.csv"))
    if not metric_files:
        print("[!] No se encontraron tablas de métricas.")
        return

    all_metrics = []
    for f in metric_files:
        ab_name = os.path.basename(f).replace("metrics_", "").replace(".csv", "")
        df = pd.read_csv(f)
        df['Antibiotic'] = ab_name.capitalize()
        all_metrics.append(df)

    df_total = pd.concat(all_metrics, ignore_index=True)
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
    
    metrics = [('AUC_ROC', 'AUC-ROC'), ('Sensitivity', 'Sensibilidad'), ('Specificity', 'Especificidad')]
    palette = sns.color_palette("muted")

    for i, (col, title) in enumerate(metrics):
        sns.barplot(
            data=df_total, 
            x='Antibiotic', 
            y=col, 
            hue='Model', 
            ax=axes[i], 
            palette=palette
        )
        axes[i].set_title(title, fontweight='bold', fontsize=13)
        axes[i].set_ylim(0, 1.05)
        axes[i].set_ylabel("Puntuación" if i == 0 else "")
        axes[i].set_xlabel("")
        axes[i].tick_params(axis='x', rotation=45)
        axes[i].axhline(0.5, color='gray', linestyle='--', alpha=0.6)

    plt.tight_layout()
    out_path = os.path.join(figures_dir, "fig1_comparative_performance.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[✓] Figura 1 (Desempeño Comparativo) guardada en: {out_path}")


def plot_forest_odds_ratios(tables_dir: str, figures_dir: str):
    """Forest plot de Odds Ratios para Imipenem y Meropenem (Protectores vs Riesgo)."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=False)
    
    for idx, ab in enumerate(['imipenem', 'meropenem']):
        or_path = os.path.join(tables_dir, f"biomarkers_elasticnet_{ab}.csv")
        if not os.path.exists(or_path):
            continue

        df_or = pd.read_csv(or_path).head(10)
        df_or = df_or.sort_values(by='Odds_Ratio', ascending=True)
        
        colors = ['#2b8cbe' if val < 1.0 else '#de2d26' for val in df_or['Odds_Ratio']]
        
        y_pos = range(len(df_or))
        axes[idx].hlines(y=y_pos, xmin=1.0, xmax=df_or['Odds_Ratio'], color=colors, alpha=0.8, linewidth=2.5)
        axes[idx].scatter(df_or['Odds_Ratio'], y_pos, color=colors, s=70, zorder=3)
        axes[idx].axvline(1.0, color='black', linestyle='--', alpha=0.7)
        axes[idx].set_yticks(y_pos)
        axes[idx].set_yticklabels(df_or['Feature'])
        axes[idx].set_xlabel("Odds Ratio (OR)")
        axes[idx].set_title(f"Biomarcadores Elastic Net: {ab.upper()}", fontweight='bold', fontsize=12)

    plt.tight_layout()
    out_path = os.path.join(figures_dir, "fig2_carbapenems_biomarkers_forest.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[✓] Figura 2 (Forest Plot Odds Ratios) guardada en: {out_path}")


def plot_shap_summary(tables_dir: str, figures_dir: str):
    """Grafica la importancia global de características por SHAP."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=False)
    
    for idx, ab in enumerate(['imipenem', 'meropenem']):
        shap_path = os.path.join(tables_dir, f"biomarkers_shap_{ab}.csv")
        if not os.path.exists(shap_path):
            continue

        df_shap = pd.read_csv(shap_path).head(10)
        df_shap = df_shap.sort_values(by='Mean_Abs_SHAP', ascending=True)

        axes[idx].barh(df_shap['Feature'], df_shap['Mean_Abs_SHAP'], color='#756bb1', edgecolor='black', alpha=0.85)
        axes[idx].set_xlabel("Mean |SHAP Value| (Impacto en la predicción)")
        axes[idx].set_title(f"Top Predictores XGBoost: {ab.upper()}", fontweight='bold', fontsize=12)

    plt.tight_layout()
    out_path = os.path.join(figures_dir, "fig3_shap_importance_panel.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[✓] Figura 3 (Panel SHAP) guardada en: {out_path}")


def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    tables_dir = os.path.join(base_dir, "results", "tables")
    figures_dir = os.path.join(base_dir, "results", "figures")
    os.makedirs(figures_dir, exist_ok=True)

    print("\n================ GENERACIÓN DE FIGURAS PARA TESIS ================")
    plot_overall_performance(tables_dir, figures_dir)
    plot_forest_odds_ratios(tables_dir, figures_dir)
    plot_shap_summary(tables_dir, figures_dir)
    print("[✓] Todas las figuras fueron exportadas exitosamente en 300 DPI.")


if __name__ == "__main__":
    main()
