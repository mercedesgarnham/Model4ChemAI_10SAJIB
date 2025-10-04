import streamlit as st
import pandas as pd
import requests
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

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
                            "0": "inactivo"
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

    with st.expander("Métricas básicas"):
        st.subheader("Distribución de clases:")

        if df is None:
            st.warning("⚠️ Cargá el archivo .csv para comenzar")

        else:
            # Contar cantidad de moléculas en cada clase
            count = df["activity"].value_counts()
            st.write("Cantidad de moléculas por clase")
            st.write(count)

            # gráfico de barras
            fig, ax = plt.subplots(facecolor="none")
            count.plot(kind="bar", ax=ax, color=["#4CAF50", "#F44336"])
            ax.set_facecolor("none")
            ax.yaxis.label.set_color("white")
            ax.tick_params(colors="white")
            ax.xaxis.label.set_color("white")
            ax.title.set_color("white")
            for spine in ax.spines.values():
                spine.set_color("white")
            ax.set_xlabel("Clase")
            ax.set_ylabel("Cantidad")
            ax.set_title("Número de moléculas por clase")
            st.pyplot(fig)

            # grafico de torta
            fig, ax = plt.subplots(facecolor="none")
            ax.set_facecolor("none")
            wedges, texts, autotexts = ax.pie(
                count,
                labels=count.index,
                autopct="%1.1f%%",
                startangle=90,
                colors=["#4CAF50", "#F44336"],
            )
            for t in texts + autotexts:
                t.set_color("white")
            ax.set_ylabel("")  # sacar label "activity"
            ax.set_title("Porcentaje de moléculas por clase", color="white")
            st.pyplot(fig, transparent=True)

    with st.expander("Distribución de similitud"):
        st.header("Heatmap:")

        if df is None:
            st.warning("⚠️ Cargá el archivo .csv para comenzar")

        else:
            # funcion para generar fingerprints
            def fingerprint(smiles):
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
                return None

            # funcion para plottear heatmap de similaridad
            def heatmap_plot(df):
                # separar activos e inactivos
                activos = df[df["activity"] == "activo"]["smiles"].tolist()
                inactivos = df[df["activity"] == "inactivo"]["smiles"].tolist()

                # calcular fingerprints para cada smiles
                fps_act = [fingerprint(smi) for smi in activos]
                fps_inact = [fingerprint(smi) for smi in inactivos]

                # funcion para generar matriz de similitud
                def sim_matrix(group1, group2):
                    sims = np.zeros((len(group1), len(group2)))
                    for i, fp1 in enumerate(group1):
                        for j, fp2 in enumerate(group2):
                            sims[i, j] = DataStructs.TanimotoSimilarity(fp1, fp2)
                    return sims
                
                # calcular matrices (act vs act, inact vs inact, act vs inact)
                mat_aa = sim_matrix(fps_act, fps_act)
                mat_ii = sim_matrix(fps_inact, fps_inact)
                mat_ai = sim_matrix(fps_act, fps_inact)

                # plottear heatmaps
                fig, axs = plt.subplots(1, 3, figsize=(18, 5), facecolor="none")
                plt.subplots_adjust(wspace=0.3)
                for ax, mat, title in zip(
                    axs,
                    [mat_aa, mat_ai, mat_ii],
                    ["Activos vs Activos", "Inactivos vs Inactivos", "Activos vs Inactivos"]
                ):
                    sns.heatmap(mat, ax=ax, cmap="viridis", cbar=False, square=True)
                    ax.set_title(title, color="white")
                    ax.set_facecolor("none")
                    ax.tick_params(colors="white")
                    for spine in ax.spines.values():
                        spine.set_color("white")
                return fig

            fig = heatmap_plot(df)
            st.pyplot(fig, transparent=True)

                        


with tab3:
    st.header("Visualización")
    st.write("Falta generar el código de gráficos y visualizaciones.")
    with st.expander("Gráficos"):
        st.write("Visualizaciones de la distribución...")

