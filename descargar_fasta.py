import os
import pandas as pd
import requests
import time

# 1. Leer los IDs de tu archivo de metadatos
ruta_csv = os.path.expanduser("~/Documents/Proyecto_Tesis/metadatos/matriz_completa_cmi.csv")
df = pd.read_csv(ruta_csv)

# Sacamos solo los IDs únicos para no descargar el mismo genoma dos veces
genomas = df['genome_id'].unique()

# 2. Carpeta donde se guardarán los archivos .fasta
carpeta_salida = os.path.expanduser("~/Documents/Proyecto_Tesis/secuencias")
print(f"-> Iniciando descarga de {len(genomas)} genomas. Esto tomará tiempo...")

# 3. Descargar uno por uno
for i, genoma_id in enumerate(genomas):
    archivo_salida = os.path.join(carpeta_salida, f"{genoma_id}.fasta")
    
    # Si el archivo ya existe, nos lo saltamos (ideal por si se corta el internet)
    if os.path.exists(archivo_salida):
        continue
        
    # URL especial de BV-BRC para descargar secuencias de ADN
    url = f"https://www.bv-brc.org/api/genome_sequence/?eq(genome_id,{genoma_id})&limit(25000)"
    headers = {"Accept": "application/dna+fasta"}
    
    try:
        respuesta = requests.get(url, headers=headers)
        if respuesta.status_code == 200:
            # Guardamos el texto descargado en un archivo .fasta
            with open(archivo_salida, "w") as archivo:
                archivo.write(respuesta.text)
            print(f"[{i+1}/{len(genomas)}] Descargado: {genoma_id}.fasta")
        else:
            print(f"[{i+1}/{len(genomas)}] Error al descargar {genoma_id}")
    except Exception as e:
        print(f"Error de conexión en {genoma_id}: {e}")
        
    # Pausa pequeña para no saturar al servidor
    time.sleep(0.5) 

print("-> ¡Descarga completada!")
