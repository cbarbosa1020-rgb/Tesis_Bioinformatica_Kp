"""
Proyecto: KPC_UNAL
Script: 09_map_biomarkers_to_plasmids.py
Descripción: Normaliza cabeceras ignorando el sufijo ';size=X' para vincular
             perfectamente los clusters de Machine Learning con sus dianas plasmídicas.
"""

import os
import glob
import pandas as pd


def parse_vsearch_uc(uc_path: str) -> pd.DataFrame:
    """Mapea cada secuencia (limpia de ';size=') con su respectivo Cluster_X."""
    mapping = []
    with open(uc_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            # 'S' = Semilla/Centroide, 'H' = Miembro del cluster
            if parts[0] in ['S', 'H']:
                cluster_id = f"Cluster_{parts[1]}"
                clean_seq_id = parts[8].split(';')[0]
                mapping.append({
                    'Feature': f"Spacer_{cluster_id}",
                    'clean_seq_id': clean_seq_id
                })
                
    df_uc = pd.DataFrame(mapping).drop_duplicates()
    print(f"[✓] Mapeo VSEARCH: {len(df_uc)} secuencias asociadas a clusters.")
    return df_uc


def load_biomarkers(tables_dir: str) -> pd.DataFrame:
    """Carga y consolida los Odds Ratios (Elastic Net) y SHAP (XGBoost)."""
    records = []

    # Cargar Odds Ratios
    for f in glob.glob(os.path.join(tables_dir, "biomarkers_elasticnet_*.csv")):
        ab = os.path.basename(f).replace("biomarkers_elasticnet_", "").replace(".csv", "")
        df = pd.read_csv(f)
        for _, row in df.iterrows():
            feat = str(row['Feature'])
            if feat.startswith("Spacer_"):
                records.append({
                    'Antibiotic': ab.capitalize(),
                    'Feature': feat,
                    'Odds_Ratio': float(row['Odds_Ratio']) if pd.notna(row.get('Odds_Ratio')) else None
                })

    # Cargar SHAP
    for f in glob.glob(os.path.join(tables_dir, "biomarkers_shap_*.csv")):
        ab = os.path.basename(f).replace("biomarkers_shap_", "").replace(".csv", "")
        df = pd.read_csv(f)
        for _, row in df.iterrows():
            feat = str(row['Feature'])
            if feat.startswith("Spacer_"):
                records.append({
                    'Antibiotic': ab.capitalize(),
                    'Feature': feat,
                    'SHAP_Value': float(row['Mean_Abs_SHAP']) if pd.notna(row.get('Mean_Abs_SHAP')) else None
                })

    df_bio = pd.DataFrame(records)
    if df_bio.empty:
        return pd.DataFrame()

    pivoted = df_bio.groupby(['Antibiotic', 'Feature']).agg({
        'Odds_Ratio': 'first',
        'SHAP_Value': 'first'
    }).reset_index()

    return pivoted


def build_master_table(tables_dir: str, blast_file: str, uc_file: str, output_csv: str):
    """Cruza biomarcadores, secuencias y alineamientos BLASTn."""
    df_bio = load_biomarkers(tables_dir)
    df_uc = parse_vsearch_uc(uc_file)
    
    # Cargar BLAST y limpiar cabecera
    df_blast = pd.read_csv(blast_file, sep='\t')
    df_blast['clean_seq_id'] = df_blast['qseqid'].apply(lambda x: str(x).split(';')[0])

    # 1. Cruzar biomarcadores con secuencias de cada cluster
    merged_uc = pd.merge(df_bio, df_uc, on='Feature', how='inner')

    # 2. Cruzar con resultados BLAST
    merged_all = pd.merge(merged_uc, df_blast, on='clean_seq_id', how='left')

    # 3. Categorización de efecto clínico
    def get_role(or_val):
        if pd.isna(or_val):
            return "Predictor Relevante (SHAP)"
        return "Protector (Sensibilidad)" if float(or_val) < 1.0 else "Riesgo (Resistencia)"

    merged_all['Clinical_Effect'] = merged_all['Odds_Ratio'].apply(get_role)

    # 4. Agrupar y ordenar priorizando hits de BLAST reales
    merged_all['has_hit'] = merged_all['stitle'].notna()
    sorted_df = merged_all.sort_values(
        by=['Antibiotic', 'Feature', 'has_hit', 'bitscore', 'pident'],
        ascending=[True, True, False, False, False]
    )

    master = sorted_df.groupby(['Antibiotic', 'Feature']).first().reset_index()

    cols_order = [
        'Antibiotic', 'Feature', 'Clinical_Effect', 'Odds_Ratio', 'SHAP_Value',
        'pident', 'length', 'evalue', 'bitscore', 'stitle'
    ]
    cols_final = [c for c in cols_order if c in master.columns]
    master = master[cols_final].sort_values(by=['Antibiotic', 'Odds_Ratio'], ascending=[True, True])

    master.to_csv(output_csv, index=False)
    print(f"\n[✓] Tabla Maestra guardada exitosamente en:\n    {output_csv}")

    print("\n=== TOP BIOMARCADORES Y DIANAS PLASMÍDICAS ANOTADAS ===")
    hits_found = master[master['stitle'].notna()]
    if not hits_found.empty:
        print(hits_found[['Antibiotic', 'Feature', 'Clinical_Effect', 'Odds_Ratio', 'pident', 'stitle']].head(15).to_string(index=False))
    else:
        print(master[['Antibiotic', 'Feature', 'Clinical_Effect', 'Odds_Ratio']].head(15).to_string(index=False))


def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    tables_dir = os.path.join(base_dir, "results", "tables")
    raw_dir = os.path.join(base_dir, "data", "raw")
    blast_summary = os.path.join(base_dir, "results", "blast", "plasmid_targets_summary.tsv")
    uc_file = os.path.join(raw_dir, "spacers_clusters_s95.uc")
    output_master = os.path.join(tables_dir, "master_biomarkers_annotated.csv")

    print("\n================ MAPEO DE BIOMARCADORES Y PLÁSMIDOS ================")
    build_master_table(tables_dir, blast_summary, uc_file, output_master)


if __name__ == "__main__":
    main()
