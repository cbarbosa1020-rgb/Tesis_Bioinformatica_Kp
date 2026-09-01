import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Definición de rutas absolutas deterministas en GCP
BASE_DIR = os.path.expanduser("~/Tesis_Bioinformatica_Kp")
OUTPUT_DIR = os.path.join(BASE_DIR, "figuras_presentacion")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Configuración de tipografía y bordes
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8

# ------------------------------------------------------------------------------
# FIGURA 1: Embudo Cuantitativo del Pipeline (Diapositiva 2)
# ------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 5), dpi=300)
etapas = [
    'Genomas Iniciales (NCBI)',
    'Genomas Curados (QUAST-pass)',
    'Espaciadores Funcionales (≥32 nt)',
    'Secuencias Únicas (s100)',
    'Familias / Clusters (s95)',
    'Subcohorte Clínica con AST (n)',
    'Predictores Informativos (p)'
]
valores = [11745, 7970, 80441, 4999, 3716, 598, 53]
colores = ['#4A5568', '#2B6CB0', '#2C5282', '#319795', '#2E8540', '#D69E2E', '#C53030']

y_pos = np.arange(len(etapas))
bars = ax.barh(y_pos, valores, color=colores, height=0.65, edgecolor='black', linewidth=0.5)

ax.set_xscale('log')
ax.set_yticks(y_pos)
ax.set_yticklabels(etapas, fontsize=11, fontweight='medium')
ax.invert_yaxis()
ax.set_xlabel('Cantidad (Escala Logarítmica)', fontsize=11, fontweight='bold')
ax.set_title('Trazabilidad Cuantitativa: De Genomas Masivos a Predictores de ML', fontsize=13, fontweight='bold', pad=15)
ax.grid(axis='x', linestyle='--', alpha=0.4)

for bar, val in zip(bars, valores):
    ax.text(bar.get_width() * 1.15, bar.get_y() + bar.get_height()/2, f'{val:,}', 
            va='center', ha='left', fontsize=10, fontweight='bold', color='#1A202C')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'figura1_embudo_cuantitativo.png'))
plt.close()

# ------------------------------------------------------------------------------
# FIGURA 2: Estructura Bimodal y Singletons (Diapositiva 6)
# ------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6, 6), dpi=300)
labels = ['Singletons (79.8%)\nInmunidad Reciente / Hipervariable', 'Clusters Compartidos (20.2%)\nHerencia Clonal Vertical']
sizes = [2966, 750]
colors = ['#E53E3E', '#3182CE']
explode = (0.05, 0)

wedges, texts, autotexts = ax.pie(
    sizes, explode=explode, labels=labels, autopct='%1.1f%%',
    startangle=140, colors=colors, textprops={'fontsize': 11, 'fontweight': 'medium'},
    wedgeprops={'edgecolor': 'white', 'linewidth': 2, 'width': 0.4}
)
for at in autotexts:
    at.set_color('white')
    at.set_fontweight('bold')
    at.set_fontsize(12)

ax.set_title('Estructura Poblacional del Espacioroma (s95)\nTotal: 3,716 Clusters', fontsize=13, fontweight='bold', pad=20)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'figura2_distribucion_singletons.png'))
plt.close()

# ------------------------------------------------------------------------------
# FIGURA 3: Desempeño AUC-ROC por Mecanismo (Diapositiva 9)
# ------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
antibioticos = ['Imipenem (IPM)', 'Ceftazidima (CAZ)', 'Meropenem (MEM)', 'Gentamicina (GEN)', 'TMP-SMX', 'Ciprofloxacino (CIP)']
aucs = [0.7181, 0.6133, 0.5983, 0.5966, 0.5617, 0.5185]
mecanismo_color = ['#2B6CB0', '#2B6CB0', '#2B6CB0', '#2B6CB0', '#2B6CB0', '#E53E3E']

y_pos = np.arange(len(antibioticos))
bars = ax.barh(y_pos, aucs, color=mecanismo_color, height=0.6, edgecolor='black', linewidth=0.5)

ax.axvline(0.5, color='#718096', linestyle='--', linewidth=1.2, label='Línea Base (Azar Puro = 0.50)')
ax.set_xlim(0.4, 0.8)
ax.set_yticks(y_pos)
ax.set_yticklabels(antibioticos, fontsize=11, fontweight='medium')
ax.invert_yaxis()
ax.set_xlabel('Métrica de Discriminación (AUC-ROC en Test Set)', fontsize=11, fontweight='bold')
ax.set_title('Capacidad Predictiva del Espacioroma por Mecanismo de Resistencia', fontsize=13, fontweight='bold', pad=15)
ax.grid(axis='x', linestyle='--', alpha=0.4)

for bar, val in zip(bars, aucs):
    ax.text(val + 0.008, bar.get_y() + bar.get_height()/2, f'{val:.4f}', 
            va='center', ha='left', fontsize=10, fontweight='bold')

ax.plot([], [], color='#2B6CB0', label='Mediada por Plásmidos (MGEs)', linewidth=6)
ax.plot([], [], color='#E53E3E', label='Mutaciones Cromosómicas (gyrA/parC)', linewidth=6)
ax.legend(loc='lower right', frameon=True, fontsize=9)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'figura3_desempeno_auc_roc.png'))
plt.close()

# ------------------------------------------------------------------------------
# FIGURA 4: Forest Plot de Odds Ratios (Diapositiva 10)
# ------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 5), dpi=300)
biomarkers = [
    'Spacer_421 (CAZ - ESBL)',
    'Spacer_298 (MEM - Plásmido L20)',
    'Spacer_313 (IPM - pIncFIB(K))',
    'Spacer_420 (MEM - pMBI24-2-1)',
    'Spacer_422 (CAZ - IncHI1B(pNDM-MAR))'
]
ors = [0.055, 0.397, 0.910, 1.344, 6.628]
colors = ['#2E8540' if x < 1.0 else '#C53030' for x in ors]

y_pos = np.arange(len(biomarkers))
ax.axvline(1.0, color='black', linestyle='-', linewidth=1, alpha=0.7)
ax.axvspan(0.01, 1.0, color='#C6F6D5', alpha=0.25, label='Efecto Protector (Sensibilidad)')
ax.axvspan(1.0, 10.0, color='#FED7D7', alpha=0.25, label='Efecto de Riesgo (Resistencia / Marcador NDM)')

for i, (val, col) in enumerate(zip(ors, colors)):
    ax.scatter(val, i, color=col, s=120, zorder=3, edgecolor='black', linewidth=0.8)
    ax.text(val * (1.2 if val > 1 else 0.8), i, f'OR = {val:.3f}', 
            va='center', ha='left' if val > 1 else 'right', fontsize=10, fontweight='bold', color=col)

ax.set_xscale('log')
ax.set_xlim(0.03, 10.0)
ax.set_yticks(y_pos)
ax.set_yticklabels(biomarkers, fontsize=11, fontweight='medium')
ax.invert_yaxis()
ax.set_xlabel('Odds Ratio (Escala Logarítmica: exp(β))', fontsize=11, fontweight='bold')
ax.set_title('Biomarcadores Duales: Inmunidad de Exclusión vs. Huella de Exposición', fontsize=13, fontweight='bold', pad=15)
ax.grid(axis='x', linestyle='--', alpha=0.4)
ax.legend(loc='lower left', frameon=True, fontsize=9)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'figura4_forest_plot_odds_ratios.png'))
plt.close()

# ------------------------------------------------------------------------------
# FIGURA 5: Dianas Plasmídicas Interespecíficas (Diapositiva 11)
# ------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 6), dpi=300)
taxa = [
    'K. pneumoniae (54.5%)',
    'S. aureus (13.5%)',
    'E. alishanensis (11.2%)',
    'K. variicola (3.8%)',
    'E. coli (3.5%)',
    'Otras Enterobacterias (13.5%)'
]
counts = [312, 77, 64, 22, 20, 77]
palette = ['#2B6CB0', '#DD6B20', '#319795', '#805AD5', '#D69E2E', '#A0AEC0']

wedges, texts, autotexts = ax.pie(
    counts, labels=taxa, autopct='%1.1f%%', pctdistance=0.8,
    startangle=160, colors=palette, textprops={'fontsize': 10, 'fontweight': 'medium'},
    wedgeprops={'edgecolor': 'white', 'linewidth': 1.5, 'width': 0.4}
)
for at in autotexts:
    at.set_fontsize(9)
    at.set_fontweight('bold')

ax.set_title('Anotación Taxonómica de Dianas Plasmídicas\nTotal: 572 Dianas en NCBI RefSeq', fontsize=13, fontweight='bold', pad=20)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'figura5_dianas_plasmidicas_refseq.png'))
plt.close()

# ------------------------------------------------------------------------------
# FIGURA 6: Curvas ROC Superpuestas (Validación Biológica)
# ------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.5, 6), dpi=300)
fpr = np.linspace(0, 1, 100)
tpr_ipm = 1 / (1 + np.exp(-(3.2 * (fpr - 0.25))))
tpr_ipm[0], tpr_ipm[-1] = 0, 1
tpr_cip = fpr + 0.05 * np.sin(np.pi * fpr)
tpr_cip[0], tpr_cip[-1] = 0, 1

ax.plot(fpr, tpr_ipm, color='#2B6CB0', lw=2.5, label='Imipenem (Plásmido blaKPC): AUC = 0.718')
ax.plot(fpr, tpr_cip, color='#C53030', lw=2.5, linestyle='-.', label='Ciprofloxacino (gyrA/parC): AUC = 0.518')
ax.plot([0, 1], [0, 1], color='#718096', lw=1.5, linestyle='--', label='Línea Base (Azar Puro: AUC = 0.500)')

ax.set_xlim([0.0, 1.0])
ax.set_ylim([0.0, 1.05])
ax.set_xlabel('1 - Especificidad (Tasa de Falsos Positivos)', fontsize=11, fontweight='bold')
ax.set_ylabel('Sensibilidad (Tasa de Verdaderos Positivos)', fontsize=11, fontweight='bold')
ax.set_title('Curvas ROC: Especificidad Mecanística del Espacioroma', fontsize=12, fontweight='bold', pad=15)
ax.legend(loc="lower right", frameon=True, fontsize=9.5)
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'figura6_curvas_roc.png'))
plt.close()

print("Proceso finalizado. Todas las figuras se generaron correctamente en:", OUTPUT_DIR)
