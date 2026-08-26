import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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
    file_path = os.path.join(tables_dir, f"biomarkers_shap_{ab_key}.csv")
    
    if os.path.exists(file_path):
        df = pd.read_csv(file_path)
        col_feature = 'Feature' if 'Feature' in df.columns else df.columns[0]
        # Detectar columna de valor SHAP
        shap_cols = [c for c in df.columns if 'shap' in c.lower() or 'importance' in c.lower() or 'mean' in c.lower()]
        col_val = shap_cols[0] if len(shap_cols) > 0 else df.columns[1]
        
        df[col_val] = pd.to_numeric(df[col_val], errors='coerce').fillna(0)
        df_top = df.sort_values(by=col_val, ascending=True).tail(10)
        
        y_pos = np.arange(len(df_top))
        
        ax.barh(y_pos, df_top[col_val], color='#7570b3', alpha=0.85, edgecolor='black', linewidth=0.8)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(df_top[col_feature], fontsize=9)
        ax.set_title(f"XGBoost SHAP: {ab_name}", fontsize=11, fontweight='bold')
        ax.set_xlabel("Mean |SHAP Value| (Impacto en la predicción)", fontsize=9)
        ax.grid(axis='x', linestyle=':', alpha=0.6)
    else:
        ax.text(0.5, 0.5, f"Sin datos: {ab_name}", ha='center', va='center', fontsize=10)

plt.tight_layout()
output_path = os.path.join(fig_dir, "fig3_panel_6_antibioticos_shap.png")
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"[OK] Figura 3 generada exitosamente en: {output_path}")
