#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path

SPACEROME_DIR = Path.home() / "data" / "spacerome"
RAW_FASTA = SPACEROME_DIR / "all_spacers_raw.fasta"

REP_100 = SPACEROME_DIR / "spacers_rep_s100.fasta"
UC_100 = SPACEROME_DIR / "spacers_clusters_s100.uc"

REP_95 = SPACEROME_DIR / "spacers_rep_s95.fasta"
UC_95 = SPACEROME_DIR / "spacers_clusters_s95.uc"

THREADS = 32

def run_cmd(cmd, desc):
    print(f"\n==> {desc}...")
    print(f"Comando: {' '.join(cmd)}")
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"[ERROR] Falló el comando:\n{res.stderr}", file=sys.stderr)
        sys.exit(1)
    print(res.stderr) # Vsearch imprime estadisticas en stderr

def count_fasta_records(fasta_path):
    with open(fasta_path, "r") as f:
        return sum(1 for line in f if line.startswith(">"))

def main():
    if not RAW_FASTA.exists():
        print(f"[ERROR] No se encuentra {RAW_FASTA}", file=sys.stderr)
        sys.exit(1)
        
    total_raw = count_fasta_records(RAW_FASTA)
    print(f"==> Iniciando clustering del Spacerome sobre {total_raw} secuencias brutas.")

    # 1. Dereplicacion al 100% de identidad
    cmd_100 = [
        "vsearch",
        "--derep_fulllength", str(RAW_FASTA),
        "--strand", "both",
        "--threads", str(THREADS),
        "--output", str(REP_100),
        "--uc", str(UC_100),
        "--sizeout"
    ]
    run_cmd(cmd_100, "Ejecutando dereplicación exacta al 100% (VSEARCH)")
    total_s100 = count_fasta_records(REP_100)

    # 2. Clustering al 95% de identidad (usando los representantes s100 como base para máxima velocidad)
    cmd_95 = [
        "vsearch",
        "--cluster_fast", str(REP_100),
        "--id", "0.95",
        "--strand", "both",
        "--threads", str(THREADS),
        "--centroids", str(REP_95),
        "--uc", str(UC_95),
        "--sizeout"
    ]
    run_cmd(cmd_95, "Ejecutando clustering al 95% de identidad (VSEARCH)")
    total_s95 = count_fasta_records(REP_95)

    print("\n=======================================================")
    print("           RESUMEN DEL SPACEROME GENERADO              ")
    print("=======================================================")
    print(f"Total espaciadores brutos:          {total_raw:,}")
    print(f"Espaciadores únicos exactos (100%):  {total_s100:,}  (Reducción: {100 - (total_s100/total_raw*100):.1f}%)")
    print(f"Clusters representativos (95%):      {total_s95:,}  (Reducción: {100 - (total_s95/total_raw*100):.1f}%)")
    print("=======================================================")
    print(f"Archivos generados en {SPACEROME_DIR}:")
    print(f" - {REP_100.name} / {UC_100.name}")
    print(f" - {REP_95.name} / {UC_95.name}")

if __name__ == "__main__":
    main()
