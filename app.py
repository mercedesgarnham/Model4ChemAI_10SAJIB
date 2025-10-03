import streamlit as st
import pandas as pd
import requests
from rdkit import Chem

st.title("App Molecular con Pestañas y Carga Condicional")

tab1, tab2, tab3 = st.tabs(["Carga de datos", "Análisis", "Visualización"])

with tab1:
    st.header("Carga tus datos")

    modo_carga = st.radio("Selecciona cómo cargar los datos:", ["Subir CSV", "Descargar desde ChEMBL"])

    df = None

    if modo_carga == "Subir CSV":
        uploaded_file = st.file_uploader("Subí tu archivo CSV", type=["csv"])
        st.write('El archivo .csv debe tener una columna "SMILES" y una columna "activity".')
        if uploaded_file is not None:
            try:
                df = pd.read_csv(uploaded_file)

                # Pasar los titulos de las columnas a strings y en minuscula
                df.columns = df.columns.str.lower()

                # Aceptar algunas posibles variaciones en los nombres de las columnas
                rename_headers = {
                    "actividad": "activity",
                    "molecula": "smiles",
                    "moleculas": "smiles",
                    "smile" : "smiles",
                }

                df.rename(columns=rename_headers, inplace=True)

                # Chequear que  las columnas sean correctas
                required_columns = ["smiles", "activity"]
                missing = [col for col in required_columns if col not in df.columns]

                if missing:
                    st.error(f"⚠️ Faltan columnas requeridas: {missing}")
                
                else:
                    st.success("Archivo subido correctamente ✅")

                    # Chequear que los SMILES sean validos
                    invalid_smiles = []
                    for i, smi in enumerate(df["smiles"]):
                        if Chem.MolFromSmiles(str(smi)) is None:
                            invalid_smiles.append(i)
                    if invalid_smiles:
                        st.error(f"⚠️ Se encontraron {len(invalid_smiles)} SMILES invalidos (filas: {invalid_smiles[:]})")
                    else:
                        st.success("Todos los SMILES son válidos ✅")

                        # Aceptar valores posibles en la columna "activity"
                        df["activity"] = df["activity"].astype(str).str.lower()
                        rename_activity = {
                            "a": "activo",
                            "active": "activo",
                            "act": "activo",
                            "1": "activo",
                            "i": "inactivo",
                            "inactive": "inactivo",
                            "inact": "inactivo",
                            "0": "activo"
                        }
                        df["activity"] = df["activity"].replace(rename_activity)

                        # Chequear validez de todos los valores en "activity"
                        valid_values_act = ["activo", "inactivo"]
                        invalid_rows_act = df[~df["activity"].isin(valid_values_act)]

                        if not invalid_rows_act.empty:
                            st.error(f'⚠️ Se encontraron valores inválidos en "activity" en las filas: {invalid_rows_act.index.tolist()}. Los valores deberían ser "activos" o "inactivos".')
                        else:
                            st.success('Todos los valores de "activity" son válidos ✅')
                            st.success('Archivo aceptado correctamente ✅')
                            df.head()

            except Exception as e:
                st.error(f"⚠️ Error al procesar el archivo: {e}")



    elif modo_carga == "Descargar desde ChEMBL":
        st.write("Descargando datos desde la API de ChEMBL...")
        
        target_id = st.text_input("Introduce target_chembl_id (ej: CHEMBL25):", "CHEMBL25")
        st.write("Falta generar el código para obtener los datos de ChEMBL")
        
with tab2:
    st.header("Análisis")
    st.write("Falta generar el código de análisis de los datos cargados.")
    with st.expander("Métricas básicas"):
        st.write("Cantidad total, activos, inactivos...")

with tab3:
    st.header("Visualización")
    st.write("Falta generar el código de gráficos y visualizaciones.")
    with st.expander("Gráficos"):
        st.write("Visualizaciones de la distribución...")

