#!/usr/bin/env python3
# ==============================================================================
# 03_filter_quast.py
# Consolidador y filtro de QC para reportes de QUAST (K. pneumoniae)
# Criterios: 5.3 MB <= Genome size <= 5.7 MB, GC ~57%, N50 >= 50kb, Contigs <= 150
# ==============================================================================

from pathlib import Path
import pandas as pd

QUAST_DIR = Path.home() / "data" / "quast_results"
OUTPUT_DIR = Path.home() / "data" / "qc_pass_fastas"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 1. Cargar y concatenar todos los archivos transposed_report.tsv
tsv_files = sorted(QUAST_DIR.glob("batch_*/transposed_report.tsv"))
print(f"==> Encontrados {len(tsv_files)} reportes de lotes QUAST.")

df_list = []
for tsv in tsv_files:
    df_batch = pd.read_csv(tsv, sep="\t")
    df_list.append(df_batch)

full_df = pd.concat(df_list, ignore_index=True)
print(f"Total de genomas analizados: {len(full_df)}")

# 2. Renombrar columnas clave de QUAST para facilitar el acceso
full_df = full_df.rename(columns={
    "Assembly": "Sample",
    "Total length": "Total_length",
    "# contigs": "Contigs",
    "GC (%)": "GC_content"
})

# Convertir tipos de datos
full_df["Total_length"] = pd.to_numeric(full_df["Total_length"], errors="coerce")
full_df["N50"] = pd.to_numeric(full_df["N50"], errors="coerce")
full_df["Contigs"] = pd.to_numeric(full_df["Contigs"], errors="coerce")
full_df["GC_content"] = pd.to_numeric(full_df["GC_content"], errors="coerce")

# Guardar la tabla completa consolida previa a filtros
full_df.to_csv(QUAST_DIR / "quast_consolidated_summary.tsv", sep="\t", index=False)

# 3. Aplicar Filtros de Calidad
cond_size = (full_df["Total_length"] >= 5300000) & (full_df["Total_length"] <= 5700000)
cond_gc = (full_df["GC_content"] >= 56.0) & (full_df["GC_content"] <= 58.5)
cond_n50 = full_df["N50"] >= 50000
cond_contigs = full_df["Contigs"] <= 150

pass_df = full_df[cond_size & cond_gc & cond_n50 & cond_contigs].copy()
fail_df = full_df[~(cond_size & cond_gc & cond_n50 & cond_contigs)].copy()

print("\n--- RESUMEN DE CONTROL DE CALIDAD ---")
print(f" Genomas aprobados (PASS): {len(pass_df)} ({len(pass_df)/len(full_df)*100:.2f}%)")
print(f" Genomas descartados (FAIL): {len(fail_df)} ({len(fail_df)/len(full_df)*100:.2f}%)")

# Guardar listas de muestras aprobadas y reprobadas
pass_df.to_csv(QUAST_DIR / "qc_passed_samples.tsv", sep="\t", index=False)
fail_df.to_csv(QUAST_DIR / "qc_failed_samples.tsv", sep="\t", index=False)

# Exportar lista de IDs aprobados a un archivo plano de texto
with open(QUAST_DIR / "passed_samples_list.txt", "w") as f:
    for sample in pass_df["Sample"]:
        f.write(f"{sample}.fasta\n")

print(f"\n Archivos consolidados guardados en: {QUAST_DIR}")
