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
        st.write("Falta generar el código para obtener los datos de csv")

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

