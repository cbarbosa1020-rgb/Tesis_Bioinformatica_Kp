#!/usr/bin/env python3
import os
import sys
from pathlib import Path
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

RESULTS_DIR = Path.home() / "data" / "cctyper_results"
OUT_DIR = Path.home() / "data" / "spacerome"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MASTER_FASTA = OUT_DIR / "all_spacers_raw.fasta"
METADATA_TSV = OUT_DIR / "spacers_metadata.tsv"
NUM_WORKERS = 32

def parse_genome_spacers(genome_dir):
    genome_id = genome_dir.name
    spacers_dir = genome_dir / "spacers"
    records = []
    
    if not spacers_dir.exists() or not spacers_dir.is_dir():
        return records
        
    fa_files = list(spacers_dir.glob("*.fa"))
    for fa_file in fa_files:
        array_name = fa_file.stem
        try:
            with open(fa_file, "r") as f:
                header = None
                seq_lines = []
                for line in f:
                    line = line.strip()
                    if line.startswith(">"):
                        if header and seq_lines:
                            seq = "".join(seq_lines).upper()
                            sp_id = header[1:].split()[0]
                            full_header = f"{genome_id}|{array_name}|{sp_id}"
                            records.append({
                                "header": full_header,
                                "genome_id": genome_id,
                                "array_id": array_name,
                                "spacer_id": sp_id,
                                "length": len(seq),
                                "sequence": seq
                            })
                        header = line
                        seq_lines = []
                    else:
                        seq_lines.append(line)
                        
                if header and seq_lines:
                    seq = "".join(seq_lines).upper()
                    sp_id = header[1:].split()[0]
                    full_header = f"{genome_id}|{array_name}|{sp_id}"
                    records.append({
                        "header": full_header,
                        "genome_id": genome_id,
                        "array_id": array_name,
                        "spacer_id": sp_id,
                        "length": len(seq),
                        "sequence": seq
                    })
        except Exception:
            continue
            
    return records

def main():
    genome_dirs = [d for d in RESULTS_DIR.iterdir() if d.is_dir()]
    total = len(genome_dirs)
    print(f"==> Extrayendo espaciadores de {total} genomas utilizando {NUM_WORKERS} hilos...")
    
    all_spacers = []
    completed = 0
    
    with ProcessPoolExecutor(max_workers=NUM_WORKERS) as executor:
        futures = {executor.submit(parse_genome_spacers, d): d for d in genome_dirs}
        for future in as_completed(futures):
            completed += 1
            res = future.result()
            if res:
                all_spacers.extend(res)
            
            if completed % 1000 == 0 or completed == total:
                print(f"Progreso: [{completed}/{total}] genomas revisados | Espaciadores acumulados: {len(all_spacers)}", flush=True)

    print(f"\n==> Total espaciadores extraídos: {len(all_spacers)}")
    if not all_spacers:
        print("[!] No se encontraron espaciadores.")
        return

    print("==> Escribiendo FASTA maestro...")
    with open(MASTER_FASTA, "w") as f_out:
        for sp in all_spacers:
            f_out.write(f">{sp['header']}\n{sp['sequence']}\n")
            
    print("==> Escribiendo tabla de metadatos...")
    df_meta = pd.DataFrame(all_spacers)
    df_meta = df_meta.drop(columns=["sequence"])
    df_meta.to_csv(METADATA_TSV, sep="\t", index=False)
    
    print(f"==> Generado exitosamente en:")
    print(f"    - FASTA: {MASTER_FASTA}")
    print(f"    - TSV:   {METADATA_TSV}")

if __name__ == "__main__":
    main()
