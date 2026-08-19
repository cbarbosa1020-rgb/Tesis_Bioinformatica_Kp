"""
Proyecto: KPC_UNAL
Script: 07_extract_top_spacers_fasta.py
Descripción: Extrae las secuencias FASTA representativas de los clusters de espaciadores 
             biomarcadores (protectores y de riesgo) para análisis BLASTn contra plásmidos.
"""

import os
import glob
import pandas as pd


def get_top_cluster_ids(tables_dir: str, top_n: int = 15) -> set:
    """Recopila los IDs de clusters más relevantes a partir de tablas OR y SHAP."""
    cluster_ids = set()

    for f in glob.glob(os.path.join(tables_dir, "biomarkers_*.csv")):
        df = pd.read_csv(f)
        feat_col = 'Feature' if 'Feature' in df.columns else df.columns[0]
        
        spacers = df[df[feat_col].str.startswith("Spacer_")][feat_col].head(top_n)
        for sp in spacers:
            # Convierte formato 'Spacer_Cluster_281' -> '281' o 'Cluster_281'
            c_num = sp.replace("Spacer_Cluster_", "").replace("Spacer_", "")
            cluster_ids.add(c_num)
            cluster_ids.add(f"Cluster_{c_num}")

    print(f"[+] Total de clusters biomarcadores únicos identificados: {len(cluster_ids)}")
    return cluster_ids


def extract_fasta_sequences(rep_fasta_path: str, target_clusters: set, output_fasta_path: str):
    """Filtra y guarda las secuencias FASTA de los clusters diana."""
    extracted = 0
    writing = False
    current_header = ""

    os.makedirs(os.path.dirname(output_fasta_path), exist_ok=True)

    with open(rep_fasta_path, 'r') as infile, open(output_fasta_path, 'w') as outfile:
        for line in infile:
            if line.startswith(">"):
                header = line.strip()
                # Verificar coincidencia por número de cluster o ID de cabecera
                writing = any(cid in header for cid in target_clusters)
                if writing:
                    outfile.write(f"{header}\n")
                    extracted += 1
            elif writing:
                outfile.write(line)

    print(f"[✓] Se extrajeron {extracted} secuencias FASTA en: {output_fasta_path}")


def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    tables_dir = os.path.join(base_dir, "results", "tables")
    raw_dir = os.path.join(base_dir, "data", "raw")
    output_dir = os.path.join(base_dir, "results", "blast")

    rep_fasta = os.path.join(raw_dir, "spacers_rep_s95.fasta")
    out_fasta = os.path.join(output_dir, "top_spacers_biomarkers.fasta")

    if not os.path.exists(rep_fasta):
        print(f"[!] No se encontró {rep_fasta}")
        return

    print("\n================ EXTRACCIÓN DE FASTAS PARA BLASTN ================")
    top_clusters = get_top_cluster_ids(tables_dir)
    extract_fasta_sequences(rep_fasta, top_clusters, out_fasta)


if __name__ == "__main__":
    main()
