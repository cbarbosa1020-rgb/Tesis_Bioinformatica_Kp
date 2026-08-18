"""
Proyecto: KPC_UNAL
Script: 00_prepare_raw_data.py
Descripción: Parsea la salida de clusterización de VSEARCH (spacers_clusters_s95.uc)
             para generar el archivo tabular estándar 'spacers.tsv'.
"""

import os
import pandas as pd

def parse_vsearch_uc(uc_path: str, output_path: str):
    """
    Parsea el archivo .uc de VSEARCH y mapea cada genoma con su cluster s95.
    """
    print(f"[*] Parseando clusters VSEARCH desde: {uc_path}")
    
    data = []
    with open(uc_path, 'r') as f:
        for line in f:
            if line.startswith(('S', 'H')):
                parts = line.strip().split('\t')
                cluster_id = f"Cluster_{parts[1]}"
                query_label = parts[8]
                
                # Extrae el identificador del genoma
                genome_id = query_label.split('_')[0] if '_' in query_label else query_label
                
                data.append({
                    "genome_id": genome_id,
                    "spacer_id": cluster_id
                })
                
    df = pd.DataFrame(data).drop_duplicates()
    df.to_csv(output_path, sep='\t', index=False)
    print(f"[✓] 'spacers.tsv' generado con éxito ({len(df)} registros) en: {output_path}")

if __name__ == "__main__":
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    RAW_DIR = os.path.join(BASE_DIR, "data", "raw")
    
    uc_file = os.path.join(RAW_DIR, "spacers_clusters_s95.uc")
    out_file = os.path.join(RAW_DIR, "spacers.tsv")
    
    if os.path.exists(uc_file):
        parse_vsearch_uc(uc_file, out_file)
    else:
        print(f"[!] No se encontró {uc_file} en {RAW_DIR}")
