#!/usr/bin/env python3
import os
import sys
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from concurrent.futures import ProcessPoolExecutor, as_completed

RESULTS_DIR = Path.home() / "data" / "cctyper_results"
OUT_DIR = Path.home() / "data" / "spacerome"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MASTER_FASTA = OUT_DIR / "all_spacers_raw.fasta"
METADATA_TSV = OUT_DIR / "spacers_metadata.tsv"
NUM_WORKERS = 32

def parse_genome_spacers(genome_dir):
    genome_id = genome_dir.name
    spacers_tab = genome_dir / "spacers.tab"
    
    if not spacers_tab.exists() or spacers_tab.stat().st_size == 0:
        return []
    
    records = []
    try:
        df = pd.read_csv(spacers_tab, sep="\t")
        if df.empty or "Spacer" not in df.columns:
            return []
        
        for _, row in df.iterrows():
            array_id = row.get("CRISPR", "Array")
            spacer_num = row.get("Spacer", "0")
            seq = str(row.get("Sequence", "")).strip().upper()
            
            if not seq or set(seq) - set("ACGTN"):
                continue
            
            header = f"{genome_id}|{array_id}|spacer_{spacer_num}"
            records.append({
                "header": header,
                "genome_id": genome_id,
                "array_id": array_id,
                "spacer_num": spacer_num,
                "length": len(seq),
                "sequence": seq
            })
    except Exception as e:
        return []
    
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
                print(f"Progreso: [{completed}/{total}] genomas revisados | Espaciadores extraídos: {len(all_spacers)}", flush=True)

    print(f"\n==> Escribiendo FASTA maestro y tabla de metadatos...")
    
    # Exportar FASTA maestro
    with open(MASTER_FASTA, "w") as f_out:
        for sp in all_spacers:
            f_out.write(f">{sp['header']}\n{sp['sequence']}\n")
            
    # Exportar Metadata TSV
    df_meta = pd.DataFrame(all_spacers)
    df_meta.drop(columns=["sequence"], inplace=True)
    df_meta.to_csv(METADATA_TSV, sep="\t", index=False)
    
    print(f"==> Proceso completado exitosamente.")
    print(f"    - FASTA generado: {MASTER_FASTA} (Total: {len(all_spacers)} secuencias)")
    print(f"    - Tabla metadatos: {METADATA_TSV}")

if __name__ == "__main__":
    main()
