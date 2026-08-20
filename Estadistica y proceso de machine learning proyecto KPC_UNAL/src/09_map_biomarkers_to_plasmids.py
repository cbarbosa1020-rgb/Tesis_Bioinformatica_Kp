"""
Proyecto: KPC_UNAL
Script: 09_map_biomarkers_to_plasmids.py
Descripción: Cruza los clusters de espaciadores biomarcadores (Odds Ratios y SHAP)
             con los alineamientos BLASTn para generar la tabla maestra final.
"""

import os
import glob
import pandas as pd


def load_biomarkers(tables_dir: str) -> pd.DataFrame:
    """Consolida las métricas de importancia (OR y SHAP) de todos los antibióticos."""
    records = []
    
    # 1. Cargar Odds Ratios de Elastic Net
    for f in glob.glob(os.path.join(tables_dir, "biomarkers_elasticnet_*.csv")):
        ab = os.path.basename(f).replace("biomarkers_elasticnet_", "").replace(".csv", "")
        df = pd.read_csv(f)
        for _, row in df.iterrows():
            feat = str(row['Feature'])
            if feat.startswith("Spacer_"):
                records.append({
                    'Antibiotic': ab.capitalize(),
                    'Feature': feat,
                    'Odds_Ratio': row.get('Odds_Ratio', None),
                    'Coefficient': row.get('Coefficient', None),
                    'Metric_Type': 'ElasticNet_OR'
                })

    # 2. Cargar SHAP values de XGBoost
    for f in glob.glob(os.path.join(tables_dir, "biomarkers_shap_*.csv")):
        ab = os.path.basename(f).replace("biomarkers_shap_", "").replace(".csv", "")
        df = pd.read_csv(f)
        for _, row in df.iterrows():
            feat = str(row['Feature'])
            if feat.startswith("Spacer_"):
                records.append({
                    'Antibiotic': ab.capitalize(),
                    'Feature': feat,
                    'SHAP_Value': row.get('Mean_Abs_SHAP', None),
                    'Metric_Type': 'XGBoost_SHAP'
                })

    df_bio = pd.DataFrame(records)
    if df_bio.empty:
        return pd.DataFrame()

    # Consolidar métricas por antibiótico y espaciador
    pivoted = df_bio.groupby(['Antibiotic', 'Feature']).agg({
        'Odds_Ratio': 'first',
        'Coefficient': 'first',
        'SHAP_Value': 'first'
    }).reset_index()
    
    return pivoted


def map_to_blast_hits(df_biomarkers: pd.DataFrame, blast_file: str, spacers_file: str, output_csv: str):
    """Mapea los clusters con las secuencias y hits de plásmidos anotados."""
    if not os.path.exists(blast_file) or not os.path.exists(spacers_file):
        print(f"[!] Archivos faltantes: Verifique existencia de {blast_file} o {spacers_file}")
        return

    df_blast = pd.read_csv(blast_file, sep='\t')
    df_sp = pd.read_csv(spacers_file, sep='\t', dtype=str)

    # Normalizar nombres de columnas
    col_spacer = [c for c in df_sp.columns if 'spacer' in c.lower() or 'cluster' in c.lower()][0]
    col_seq = [c for c in df_sp.columns if 'id' in c.lower() or 'seq' in c.lower() or 'accn' in c.lower()][0]
    
    df_sp['Feature'] = "Spacer_" + df_sp[col_spacer].astype(str)
    
    # Cruce de Cluster con la secuencia correspondiente
    merged_bio_sp = pd.merge(df_biomarkers, df_sp, on='Feature', how='inner')
    
    # Cruce con los hits de BLAST
    final_df = pd.merge(merged_bio_sp, df_blast, left_on=col_seq, right_on='qseqid', how='left')
    
    # Clasificación de rol funcional clínico
    def get_role(or_val):
        if pd.isna(or_val): return "Predictor Relevante (SHAP)"
        try:
            val = float(or_val)
            return "Protector (Sensibilidad)" if val < 1.0 else "Riesgo (Resistencia)"
        except (ValueError, TypeError):
            return "Indeterminado"

    final_df['Clinical_Effect'] = final_df['Odds_Ratio'].apply(get_role)
    
    cols_export = [
        'Antibiotic', 'Feature', 'Clinical_Effect', 'Odds_Ratio', 'SHAP_Value', 
        'pident', 'length', 'evalue', 'sseqid', 'stitle'
    ]
    cols_avail = [c for c in cols_export if c in final_df.columns]
    
    result = final_df[cols_avail].drop_duplicates().sort_values(
        by=['Antibiotic', 'Odds_Ratio'], ascending=[True, True]
    )
    
    result.to_csv(output_csv, index=False)
    print(f"\n[✓] Tabla Maestra guardada exitosamente en:\n    {output_csv}")
    
    print("\n=== MUESTRA DE LA TABLA MAESTRA (Primeros 10 Registros) ===")
    print(result[['Antibiotic', 'Feature', 'Clinical_Effect', 'Odds_Ratio', 'stitle']].head(10).to_string(index=False))


def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    tables_dir = os.path.join(base_dir, "results", "tables")
    raw_dir = os.path.join(base_dir, "data", "raw")
    blast_summary = os.path.join(base_dir, "results", "blast", "plasmid_targets_summary.tsv")
    spacers_file = os.path.join(raw_dir, "spacers.tsv")
    output_master = os.path.join(tables_dir, "master_biomarkers_annotated.csv")

    print("\n================ MAPEO DE BIOMARCADORES Y PLÁSMIDOS ===============")
    df_bio = load_biomarkers(tables_dir)
    map_to_blast_hits(df_bio, blast_summary, spacers_file, output_master)


if __name__ == "__main__":
    main()
