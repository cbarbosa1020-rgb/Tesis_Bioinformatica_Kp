#!/usr/bin/env python3
import os
import sys
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

# Configuración de rutas
PASS_DIR = Path.home() / "data" / "fasta_pass"
OUT_DIR = Path.home() / "data" / "cctyper_results"
DB_DIR = Path("/root/miniforge3/envs/bioinfo/cct_data")
CCTYPER_BIN = "/root/miniforge3/envs/bioinfo/bin/cctyper"
NUM_WORKERS = 8

OUT_DIR.mkdir(parents=True, exist_ok=True)

# Asegurar PATH dentro de los procesos de Python
env = os.environ.copy()
env["PATH"] = f"/root/miniforge3/envs/bioinfo/bin:{env.get('PATH', '')}"

def process_genome(fasta_path):
    sample_id = fasta_path.stem
    sample_out = OUT_DIR / sample_id
    
    # Si ya se procesó previamente con éxito, no repetir
    if sample_out.exists() and any(sample_out.iterdir()):
        return f"[SKIP] {sample_id} ya procesado."

    cmd = [
        CCTYPER_BIN,
        str(fasta_path),
        str(sample_out),
        "--db", str(DB_DIR),
        "--threads", "1",
        "--no_plot"
    ]
    
    res = subprocess.run(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    if res.returncode != 0:
        return f"[ERROR] {sample_id}: {res.stderr.strip()}"
    return f"[OK] {sample_id}"

def main():
    fasta_files = sorted(list(PASS_DIR.glob("*.fasta")))
    total = len(fasta_files)
    print(f"==> Iniciando minería de {total} genomas con {NUM_WORKERS} hilos...")

    completed = 0
    errors = 0

    with ProcessPoolExecutor(max_workers=NUM_WORKERS) as executor:
        futures = {executor.submit(process_genome, f): f for f in fasta_files}
        
        for future in as_completed(futures):
            completed += 1
            result = future.result()
            
            if "[ERROR]" in result:
                errors += 1
                print(result, flush=True)
            
            if completed % 100 == 0 or completed == total:
                print(f"Progreso: [{completed}/{total}] ({completed/total*100:.1f}%) | Errores: {errors}", flush=True)

    print(f"\n==> Proceso finalizado. Total procesados: {completed} | Total errores: {errors}")

if __name__ == "__main__":
    main()
