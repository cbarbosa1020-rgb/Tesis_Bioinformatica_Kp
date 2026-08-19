"""
Proyecto: KPC_UNAL
Script: 08_blast_spacers_targets.py
Descripción: Alineamiento blastn-short multihilo de espaciadores biomarcadores
             contra la base de datos de plásmidos curada de NCBI RefSeq.
"""

import os
import subprocess
import pandas as pd


def run_blastn_short(query_fasta: str, db_path: str, output_tsv: str):
    """Ejecuta blastn-short optimizado para espaciadores CRISPR."""
    print(f"[*] Ejecutando blastn-short contra la base de datos: {db_path}")
    
    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"
    
    cmd = [
        "blastn",
        "-task", "blastn-short",
        "-query", query_fasta,
        "-db", db_path,
        "-out", output_tsv,
        "-outfmt", outfmt,
        "-evalue", "1e-3",
        "-perc_identity", "85",
        "-qcov_hsp_perc", "80",
        "-max_target_seqs", "5",
        "-num_threads", "4"
    ]
    subprocess.run(cmd, check=True)
    print(f"[✓] Alineamientos crudos guardados en: {output_tsv}")


def parse_and_summarize_hits(blast_tsv: str, summary_output: str):
    """Procesa los alineamientos y genera el resumen de dianas plasmídicas."""
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
            "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"]
    
    if not os.path.exists(blast_tsv) or os.path.getsize(blast_tsv) == 0:
        print("[!] No se detectaron alineamientos significativos con los umbrales actuales.")
        return

    df = pd.read_csv(blast_tsv, sep='\t', names=cols)
    
    # Extraer identificador simplificado
    df['genome_origin'] = df['qseqid'].apply(lambda x: x.split('|')[0] if '|' in x else x)
    
    # Filtrar mejor hit por espaciador
    top_hits = df.sort_values(by=['bitscore', 'pident'], ascending=[False, False]).groupby('qseqid').first().reset_index()
    
    # Reordenar columnas para legibilidad
    cols_order = ['genome_origin', 'qseqid', 'pident', 'length', 'evalue', 'bitscore', 'sseqid', 'stitle']
    top_hits = top_hits[[c for c in cols_order if c in top_hits.columns]]
    
    top_hits.to_csv(summary_output, sep='\t', index=False)
    print(f"[✓] Resumen de plásmidos diana ({len(top_hits)} hits únicos) guardado en: {summary_output}")


def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    query_fasta = os.path.join(base_dir, "results", "blast", "top_spacers_biomarkers.fasta")
    db_path = os.path.expanduser("~/databases/plasmids/refseq_plasmids_db")
    blast_output = os.path.join(base_dir, "results", "blast", "blast_hits_raw.tsv")
    summary_output = os.path.join(base_dir, "results", "blast", "plasmid_targets_summary.tsv")

    if not os.path.exists(query_fasta):
        print(f"[!] No se encontró el archivo de espaciadores: {query_fasta}")
        return

    print("\n================ ANÁLISIS BLASTN CONTRA PLÁSMIDOS NCBI REFSEQ ================")
    run_blastn_short(query_fasta, db_path, blast_output)
    parse_and_summarize_hits(blast_output, summary_output)


if __name__ == "__main__":
    main()
