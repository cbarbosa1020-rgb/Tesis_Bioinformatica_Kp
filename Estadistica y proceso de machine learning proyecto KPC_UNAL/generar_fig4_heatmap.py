import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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

records = []
for ab_key, ab_name in antibioticos:
    file_path = os.path.join(tables_dir, f"biomarkers_elasticnet_{ab_key}.csv")
    if os.path.exists(file_path):
        df = pd.read_csv(file_path)
        col_feat = 'Feature' if 'Feature' in df.columns else df.columns[0]
        col_or = 'Odds_Ratio' if 'Odds_Ratio' in df.columns else 'Odds Ratio'
        for _, row in df.iterrows():
            records.append({
                'Feature': str(row[col_feat]),
                'Antibiotico': ab_name,
                'Odds_Ratio': float(row[col_or])
            })

df_all = pd.DataFrame(records)

# Mapeo de dianas plasmídicas desde la tabla maestra
master_path = os.path.join(tables_dir, "master_biomarkers_annotated.csv")
plasmid_map = {}
if os.path.exists(master_path):
    df_m = pd.read_csv(master_path)
    feat_cols = [c for c in df_m.columns if 'feature' in c.lower() or 'cluster' in c.lower()]
    target_cols = [c for c in df_m.columns if any(k in c.lower() for k in ['target', 'plasmid', 'diana', 'accn'])]
    if feat_cols and target_cols:
        f_col, t_col = feat_cols[0], target_cols[0]
        for _, r in df_m.iterrows():
            f = str(r[f_col])
            t = str(r[t_col])
            if pd.notna(t) and t != "" and t != "nan":
                short_t = t.split('/')[-1].replace('K. pneumoniae plasmid ', '').replace('plasmid ', '')[:25]
                plasmid_map[f] = short_t

# Matriz pivote
pivot_or = df_all.pivot_table(index='Feature', columns='Antibiotico', values='Odds_Ratio').fillna(1.0)
log2_or = np.log2(pivot_or)

# Seleccionar los 20 biomarcadores con mayor dispersión
max_dev = log2_or.abs().max(axis=1)
top_features = max_dev.sort_values(ascending=False).head(20).index
matrix_to_plot = log2_or.loc[top_features]

# Anotar con diana plasmídica en el eje Y
new_index = [f"{feat} [{plasmid_map[feat]}]" if feat in plasmid_map else feat for feat in matrix_to_plot.index]
matrix_to_plot.index = new_index

# Graficar Heatmap
plt.figure(figsize=(11, 12))
cmap = sns.diverging_palette(240, 10, as_cmap=True)

sns.heatmap(
    matrix_to_plot,
    cmap=cmap,
    center=0.0,
    annot=True,
    fmt=".2f",
    cbar_kws={'label': r'$\log_2(\mathrm{Odds\ Ratio})$  [<0 Protector (Azul) | >0 Riesgo (Rojo)]'},
    linewidths=0.8,
    linecolor='white'
)

plt.title("Mapa de Calor: Efecto Clínico del Espacioroma y Dianas Plasmídicas", fontsize=12, fontweight='bold', pad=15)
plt.xlabel("Antimicrobiano", fontsize=11, labelpad=10)
plt.ylabel("Biomarcador CRISPR / [Diana Plasmídica]", fontsize=11)
plt.xticks(rotation=30, ha='right', fontsize=10)
plt.yticks(fontsize=9)

plt.tight_layout()
output_path = os.path.join(fig_dir, "fig4_heatmap_biomarkers_plasmids.png")
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"[OK] Figura 4 generada exitosamente en: {output_path}")
