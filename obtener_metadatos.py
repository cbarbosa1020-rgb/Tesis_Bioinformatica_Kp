import requests
import pandas as pd
import os

URL_GENOME = "https://www.bv-brc.org/api/genome/"
URL_AMR = "https://www.bv-brc.org/api/genome_amr/"

HEADERS = {
    "Content-Type": "application/rqlquery+x-www-form-urlencoded",
    "Accept": "application/json"
}

print("-> Paso 1: Extrayendo Genomas de Alta Continuidad (Complete y WGS < 100 contigs)...")

QUERY_GENOME = (
    "eq(taxon_id,573)&"
    "eq(genome_quality,Good)&"
    "in(genome_status,(Complete,WGS))&"
    "select(genome_id,genome_name,genome_status,chromosomes,contigs,n50)&"
    "limit(25000)"
)

try:
    res_genome = requests.post(URL_GENOME, data=QUERY_GENOME, headers=HEADERS)
    if res_genome.status_code != 200:
        print(f"Error en Paso 1: {res_genome.status_code}")
        exit()
        
    df_genome = pd.DataFrame(res_genome.json())
    
    # Filtro estricto local: mantener 'Complete' O 'WGS' con menos de 100 contigs
    df_genome['contigs'] = pd.to_numeric(df_genome['contigs'], errors='coerce').fillna(1)
    df_alta_calidad = df_genome[
        (df_genome['genome_status'] == 'Complete') | 
        ((df_genome['genome_status'] == 'WGS') & (df_genome['contigs'] < 100))
    ].copy()
    
    lista_ids = df_alta_calidad['genome_id'].unique().tolist()
    print(f"-> OK: Se identificaron {len(lista_ids)} genomas de alta calidad estructural.")
    
    print("-> Paso 2: Descargando TODOS los registros AMR (sin filtros de texto restrictivos)...")
    
    chunk_size = 200
    df_amr_consolidado = []
    
    for i in range(0, len(lista_ids), chunk_size):
        chunk = lista_ids[i:i+chunk_size]
        ids_str = ",".join(chunk)
        
        # Eliminamos el filtro de "MIC" y el de los antibióticos para traer TODO el bloque crudo
        QUERY_AMR = (
            f"in(genome_id,({ids_str}))&"
            "select(genome_id,antibiotic,laboratory_typing,measurement_value,measurement_sign,resistant_phenotype)&"
            "limit(25000)"
        )
        
        res_amr = requests.post(URL_AMR, data=QUERY_AMR, headers=HEADERS)
        if res_amr.status_code == 200 and res_amr.json():
            df_amr_consolidado.extend(res_amr.json())
            
    if not df_amr_consolidado:
        print("-> Alerta: Definitivamente no hay datos AMR en BV-BRC para este grupo de genomas.")
        exit()
        
    df_amr_final = pd.DataFrame(df_amr_consolidado)
    print(f"-> OK: Se extrajeron {len(df_amr_final)} registros AMR crudos.")
    
    print("-> Paso 3: Filtrando localmente por tu panel clínico...")
    antibioticos_interes = [
        'ampicillin', 'trimethoprim', 'trimethoprim/sulfamethoxazole', 'sulfamethoxazole/trimethoprim',
        'imipenem', 'meropenem', 'cefotaxime', 'ceftazidime', 'cefalexin',
        'ciprofloxacin', 'norfloxacin', 'gentamicin'
    ]
    
    # Asegurar que no hayan valores nulos que rompan el filtro
    if 'antibiotic' in df_amr_final.columns:
        df_amr_final['antibiotic'] = df_amr_final['antibiotic'].fillna('')
        df_amr_final['antibiotic_lower'] = df_amr_final['antibiotic'].str.lower()
        df_filtrado_local = df_amr_final[df_amr_final['antibiotic_lower'].isin(antibioticos_interes)].copy()
    else:
        df_filtrado_local = pd.DataFrame()
        
    print(f"-> OK: Nos quedamos con {len(df_filtrado_local)} registros que coinciden con tu panel de antibióticos.")
    
    print("-> Paso 4: Vinculando datos estructurales y fenotípicos...")
    if len(df_filtrado_local) > 0:
        matriz_metadatos = df_filtrado_local.merge(df_alta_calidad, on='genome_id', how='inner')
        ruta_salida = os.path.expanduser("~/Documents/Proyecto_Tesis/metadatos/matriz_completa_cmi.csv")
        matriz_metadatos.to_csv(ruta_salida, index=False)
        print(f"-> ¡Éxito absoluto! Tu matriz ha sido guardada en:\n   {ruta_salida}")
    else:
        print("-> Error: Tras el filtro, no quedaron antibióticos de tu panel. Revisa los nombres.")

except Exception as e:
    print(f"-> Ocurrió un error en la ejecución: {str(e)}")
