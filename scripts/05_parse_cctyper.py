#!/usr/bin/env python3
# ==============================================================================
# 05_parse_cctyper.py
# Consolidación y filtrado de operones CRISPR-Cas Tipo I-E / I-E*
# ==============================================================================

import os
import pandas as pd
from pathlib import Path

RESULTS_DIR = Path.home() / "data" / "cctyper_results"
OUTPUT_DIR = Path.home() / "data" / "cctyper_summary"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

operons_list = []
spacers_list = []
crisprs_list = []

print("==> Consolidando resultados de CRISPRCasTyper para los 7,970 genomas...")

for sample_dir in RESULTS_DIR.iterdir():
    if not sample_dir.is_dir() or sample_dir.name == "test_run":
        continue

    sample_id = sample_dir.name

    # 1. Consolidar Operones Cas
    operons_file = sample_dir / "cas_operons.tab"
    if operons_file.exists() and operons_file.stat().st_size > 0:
        try:
            df_op = pd.read_csv(operons_file, sep="\t")
            df_op["Sample"] = sample_id
            operons_list.append(df_op)
        except Exception:
            pass

    # 2. Consolidar Espaciadores
    spacers_file = sample_dir / "spacers.tab"
    if spacers_file.exists() and spacers_file.stat().st_size > 0:
        try:
            df_sp = pd.read_csv(spacers_file, sep="\t")
            df_sp["Sample"] = sample_id
            spacers_list.append(df_sp)
        except Exception:
            pass

    # 3. Consolidar Arreglos CRISPR
    crisprs_file = sample_dir / "crisprs.tab"
    if crisprs_file.exists() and crisprs_file.stat().st_size > 0:
        try:
            df_cr = pd.read_csv(crisprs_file, sep="\t")
            df_cr["Sample"] = sample_id
            crisprs_list.append(df_cr)
        except Exception:
            pass

# Guardar resumen de operones y filtrar Tipo I-E / I-E*
if operons_list:
    df_all_operons = pd.concat(operons_list, ignore_index=True)
    df_all_operons.to_csv(OUTPUT_DIR / "all_cas_operons.tsv", sep="\t", index=False)
    
    # Filtrado estricto Tipo I-E e I-E*
    df_ie = df_all_operons[df_all_operons["Prediction"].isin(["I-E", "I-E*"])]
    df_ie.to_csv(OUTPUT_DIR / "type_IE_operons.tsv", sep="\t", index=False)

    print(f"\n Operones Cas totales detectados: {len(df_all_operons)}")
    print(f" Operones Tipo I-E / I-E* identificados: {len(df_ie)}")
else:
    print("\n No se encontraron operones Cas completos en la cohorte.")

# Guardar resumen de espaciadores
if spacers_list:
    df_all_spacers = pd.concat(spacers_list, ignore_index=True)
    df_all_spacers.to_csv(OUTPUT_DIR / "all_spacers.tsv", sep="\t", index=False)
    print(f" Total de espaciadores extraídos: {len(df_all_spacers)}")

# Guardar resumen de arreglos CRISPR
if crisprs_list:
    df_all_crisprs = pd.concat(crisprs_list, ignore_index=True)
    df_all_crisprs.to_csv(OUTPUT_DIR / "all_crisprs.tsv", sep="\t", index=False)
    print(f" Total de arreglos CRISPR identificados: {len(df_all_crisprs)}")

print(f"\n==> Resumen completado. Archivos guardados en: {OUTPUT_DIR}")
