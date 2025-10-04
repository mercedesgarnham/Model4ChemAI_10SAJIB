# Model4ChemAI

## Descripción
Esta herramienta web permite **cargar, analizar y visualizar datos de moléculas** de manera interactiva.  
Se puede usar con archivos CSV locales que contengan columnas `smiles` y `activity`.

---

## Funcionalidades principales

### 1. Carga de datos
- Subir archivos CSV con moléculas (`smiles`) y actividades (`activity`).

### 2. Análisis
- Conteo de moléculas por clase (`activo` / `inactivo`).
- Ratio de clases (balance de dataset).
- Cálculo de similitud intra-clase e inter-clase usando fingerprints y Tanimoto.
- Histogramas y heatmaps de similitudes moleculares.

### 3. Visualización
- Gráfico de barras de distribución de clases.
- Heatmaps y histogramas de similitud para cada clase.

### 4. Sugerencia de modelo
- Basada en tamaño del dataset, balance de clases y similitud intra-clase.

---

## Cómo usar la herramienta

1. Abrir la aplicación web en la URL: (https://model4chemai.streamlit.app/)[https://model4chemai.streamlit.app/]
.
2. Navegar entre las pestañas:
   - **Carga de datos**: subir tu CSV.
   - **Análisis**: ver conteo de clases, ratio, histogramas y heatmaps.
   - **Visualización**: gráficos interactivos de distribución y similitudes.
   - **Sugerencia de modelo**: obtener recomendaciones de modelo de ML según tus datos.
3. Interactuar con los gráficos y métricas directamente desde la web.

---

## Autoría

- Mercedes Didier Garnham (IIB - EByN - UNSAM / CONICET)  
- Con la colaboración de Nicolás Perez Mauad (IQUIMEFA - UBA / CONICET)  
- Proyecto desarrollado durante la **hackatón 10SAJIB** organizada por **RSG Argentina** (https://rsg-argentina.netlify.app/talk/10sajib/)[https://rsg-argentina.netlify.app/talk/10sajib/]

---

## Requisitos técnicos
- No requiere instalación; funciona desde la interfaz web.
- Se recomienda usar navegador actualizado para una mejor experiencia.
