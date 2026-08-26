import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Configuración estética
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['figure.dpi'] = 300

tables_dir = os.path.join("results", "tables")
fig_dir = os.path.join("results", "figures")
os.makedirs(fig_dir, exist_ok=True)

antibioticos = [
    ("imipenem", "Imipenem"),
    ("meropenem", "Meropenem"),
    ("ceftazidime", "Ceftazidima"),
    ("ciprofloxacin", "Ciprofloxacino"),
    ("gentamicin", "Gentamicina"),
    ("trimethoprim_sulfamethoxazole", "TMP-SMX")
]

fig, axes = plt.subplots(3, 2, figsize=(14, 16), sharey=False)
axes = axes.flatten()

for idx, (ab_key, ab_name) in enumerate(antibioticos):
    ax = axes[idx]
    file_path = os.path.join(tables_dir, f"biomarkers_elasticnet_{ab_key}.csv")
    
    if os.path.exists(file_path):
        df = pd.read_csv(file_path)
        col_feature = 'Feature' if 'Feature' in df.columns else df.columns[0]
        col_or = 'Odds_Ratio' if 'Odds_Ratio' in df.columns else 'Odds Ratio'
        
        # Tomar los 10 biomarcadores con mayor desviación de OR = 1
        df['abs_dev'] = (np.log(df[col_or])).abs()
        df_top = df.sort_values(by='abs_dev', ascending=True).tail(10)
        
        y_pos = np.arange(len(df_top))
        colors = ['#1f77b4' if or_val < 1.0 else '#d62728' for or_val in df_top[col_or]]
        
        ax.hlines(y=y_pos, xmin=1.0, xmax=df_top[col_or], color=colors, alpha=0.85, linewidth=2.2)
        ax.scatter(df_top[col_or], y_pos, color=colors, s=55, zorder=3)
        ax.axvline(1.0, color='black', linestyle='--', linewidth=1.0, alpha=0.7)
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(df_top[col_feature], fontsize=9)
        ax.set_title(f"Elastic Net: {ab_name}", fontsize=11, fontweight='bold')
        ax.set_xlabel("Odds Ratio (OR)", fontsize=9)
        ax.grid(axis='x', linestyle=':', alpha=0.6)
    else:
        ax.text(0.5, 0.5, f"Sin datos: {ab_name}", ha='center', va='center', fontsize=10)

plt.tight_layout()
output_path = os.path.join(fig_dir, "fig2_panel_6_antibioticos_forest.png")
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"[OK] Figura 2 generada exitosamente en: {output_path}")
