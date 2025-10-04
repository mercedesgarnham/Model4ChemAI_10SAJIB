import streamlit as st
from utils.data_loader import load_csv
from utils.chembl_api import download_chembl_bioassay_data
from utils.data_validator import (
    validate_required_columns,
    validate_smiles,
    clean_activity_column,
    validate_activity_values
)
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

st.title("Model4ChemAI")

tab1, tab2, tab3, tab4, tab5 = st.tabs(["Carga de datos", "Comparación entre clases", "Visualización intra clases", "Sugerencia de modelo", "Documentación"])

df = None

with tab1:
    st.header("Carga tus datos")

    modo_carga = st.radio("Selecciona cómo cargar los datos:", ["Subir CSV"])

    if modo_carga == "Subir CSV":
        uploaded_file = st.file_uploader("Subí tu archivo CSV", type=["csv"])
        st.write('El archivo .csv debe tener una columna "SMILES" y una columna "activity".')
        if uploaded_file is not None:
            df, error = load_csv(uploaded_file)

            if error:
                st.error(f"⚠️ Error al cargar el archivo: {error}")
            elif df is not None:
                # Verificar columnas requeridas
                missing = validate_required_columns(df)
                if missing:
                    st.error(f"⚠️ Faltan columnas requeridas: {missing}")
                else:
                    st.success("Archivo subido correctamente ✅")

                    # Validar SMILES
                    invalid_smiles = validate_smiles(df)
                    if invalid_smiles:
                        st.error(f"⚠️ SMILES inválidos en filas: {invalid_smiles}")
                    else:
                        st.success("Todos los SMILES son válidos ✅")

                        # Limpiar y validar activity
                        df = clean_activity_column(df)
                        invalid_activity_rows = validate_activity_values(df)
                        if invalid_activity_rows:
                            st.error(f"⚠️ Valores inválidos en 'activity' en filas: {invalid_activity_rows}")
                        else:
                            st.success("Todos los valores de 'activity' son válidos ✅")
                            st.success("Archivo aceptado correctamente ✅")
                            st.dataframe(df.head())
       
with tab2:
    st.header("Comparación entre clases")

    with st.expander("Métricas básicas"):
        st.subheader("Distribución de clases:")

        if df is None:
            st.warning("⚠️ Cargá el archivo .csv para comenzar")

        else:
            total = len(df)
            st.info(f"📊 Cantidad total de moléculas: **{total}**")

            # Contar cantidad de moléculas en cada clase
            count = df["activity"].value_counts()

            n_activos = count.get("activo", 0)
            n_inactivos = count.get("inactivo", 0)

            if n_inactivos > 0:
                ratio = n_activos / n_inactivos
                st.info(f"⚖️ Ratio Activos/Inactivos: **{ratio:.2f}**")
            else:
                st.info("⚖️ No hay moléculas inactivas para calcular el ratio.")

            # gráfico de barras
            # detectar tema actual
            theme_base = st.get_option("theme.base")  # devuelve "dark" o "light"

            if theme_base == "dark":
                bar_color = "lightgrey"
                text_color = "white"
                edge_color = "white"
            else:
                bar_color = "grey"
                text_color = "black"
                edge_color = "black"

            # gráfico de barras adaptativo
            fig, ax = plt.subplots(facecolor="none")

            count.plot(kind="bar", ax=ax, color=bar_color, edgecolor=edge_color)
            
            y_max = max(count)  # valor máximo de las barras
            ax.set_ylim(0, y_max * 1.15)  # 15% más alto que el máximo

            # valores arriba de las barras
            for p in ax.patches:
                ax.text(
                    p.get_x() + p.get_width() / 2,
                    p.get_height() + 0.1,
                    int(p.get_height()),
                    ha="center", va="bottom", color=text_color, fontsize=11, fontweight="bold"
                )

            # estilo adaptativo
            ax.set_facecolor("none")
            ax.tick_params(colors=text_color, labelsize=10)
            ax.xaxis.label.set_color(text_color)
            ax.yaxis.label.set_color(text_color)
            ax.title.set_color(text_color)
            for spine in ax.spines.values():
                spine.set_color(text_color)

            ax.set_xlabel("Clase", fontsize=12)
            ax.set_ylabel("Cantidad", fontsize=12)
            ax.set_title("Número de moléculas por clase", fontsize=14, fontweight="bold")

            st.pyplot(fig)


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

                # calcular promedios de similitudes
                mean_aa = mat_aa[np.triu_indices(len(fps_act), k=1)].mean() if len(fps_act) > 1 else 0
                mean_ii = mat_ii[np.triu_indices(len(fps_inact), k=1)].mean() if len(fps_inact) > 1 else 0
                mean_ai = mat_ai.mean() if len(fps_act) > 0 and len(fps_inact) > 0 else 0

                st.write(f"🔬 Similitud promedio Activos-Activos: **{mean_aa:.2f}**")
                st.write(f"🔬 Similitud promedio Inactivos-Inactivos: **{mean_ii:.2f}**")
                st.write(f"🔬 Similitud promedio Activos-Inactivos: **{mean_ai:.2f}**")

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
    st.header("Visualización intra clases")

    if df is None:
        st.warning("⚠️ Primero cargá los datos")
    else:
        # Función para obtener fingerprints
        def fingerprint(smiles):
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            return None

        # Función para generar matriz de similitud
        def sim_matrix(fps):
            n = len(fps)
            sims = np.zeros((n, n))
            for i in range(n):
                for j in range(n):
                    sims[i, j] = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            return sims

        # ---------- EXPANDER ACTIVOS ----------
        with st.expander("Métricas básicas – Activos"):
            activos = df[df["activity"] == "activo"]["smiles"].tolist()
            fps_act = [fingerprint(smi) for smi in activos if fingerprint(smi) is not None]

            if len(fps_act) > 1:
                mat_aa = sim_matrix(fps_act)

                # Histograma
                sims_flat = mat_aa[np.triu_indices(len(fps_act), k=1)]
                fig, ax = plt.subplots()
                ax.hist(sims_flat, bins=20, color="#4CAF50", alpha=0.7)
                ax.set_xlabel("Similitud Tanimoto")
                ax.set_ylabel("Frecuencia")
                ax.set_title("Distribución de similitudes (Activos)")
                st.pyplot(fig)

                # Heatmap
                fig, ax = plt.subplots(figsize=(6, 5))
                sns.heatmap(mat_aa, cmap="viridis", square=True, cbar=True, ax=ax)
                ax.set_title("Heatmap de similitudes (Activos)")
                st.pyplot(fig)
            else:
                st.info("No hay suficientes moléculas activas para calcular similitudes.")

        # ---------- EXPANDER INACTIVOS ----------
        with st.expander("Métricas básicas – Inactivos"):
            inactivos = df[df["activity"] == "inactivo"]["smiles"].tolist()
            fps_inact = [fingerprint(smi) for smi in inactivos if fingerprint(smi) is not None]

            if len(fps_inact) > 1:
                mat_ii = sim_matrix(fps_inact)

                # Histograma
                sims_flat = mat_ii[np.triu_indices(len(fps_inact), k=1)]
                fig, ax = plt.subplots()
                ax.hist(sims_flat, bins=20, color="#F44336", alpha=0.7)
                ax.set_xlabel("Similitud Tanimoto")
                ax.set_ylabel("Frecuencia")
                ax.set_title("Distribución de similitudes (Inactivos)")
                st.pyplot(fig)

                # Heatmap
                fig, ax = plt.subplots(figsize=(6, 5))
                sns.heatmap(mat_ii, cmap="viridis", square=True, cbar=True, ax=ax)
                ax.set_title("Heatmap de similitudes (Inactivos)")
                st.pyplot(fig)
            else:
                st.info("No hay suficientes moléculas inactivas para calcular similitudes.")

with tab4:
    st.header("Sugerencia de modelo")

    if df is None:
        st.warning("⚠️ Primero cargá los datos")
    else:
        # ---- Métricas de dataset ----
        total = len(df)
        count = df["activity"].value_counts()
        n_activos = count.get("activo", 0)
        n_inactivos = count.get("inactivo", 0)

        st.write(f"**Total de moléculas:** {total}")
        st.write(f"**Activos:** {n_activos}")
        st.write(f"**Inactivos:** {n_inactivos}")

        # ---- Función para fingerprints y matriz ----
        def fingerprint(smiles):
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            return None

        def avg_intra_similarity(subset):
            fps = [fingerprint(smi) for smi in subset if fingerprint(smi) is not None]
            if len(fps) > 1:
                sims = []
                for i in range(len(fps)):
                    for j in range(i + 1, len(fps)):
                        sims.append(DataStructs.TanimotoSimilarity(fps[i], fps[j]))
                return np.mean(sims)
            return 0

        avg_activos = avg_intra_similarity(df[df["activity"] == "activo"]["smiles"].tolist())
        avg_inactivos = avg_intra_similarity(df[df["activity"] == "inactivo"]["smiles"].tolist())

        st.write(f"**Similitud promedio (activos):** {avg_activos:.2f}")
        st.write(f"**Similitud promedio (inactivos):** {avg_inactivos:.2f}")

        # ---- Reglas simples (umbrales random por ahora) ----
        sugerencia = "Baseline"

        if total < 200:
            sugerencia = "Baseline"
        elif abs(n_activos - n_inactivos) > total * 0.4:
            sugerencia = "Random Forest"
        elif (avg_activos + avg_inactivos) / 2 > 0.6:
            sugerencia = "XGBoost"
        else:
            sugerencia = "GNN"

        st.success(f"📌 Sugerencia de modelo: **{sugerencia}**")
with tab5:
    st.header("Documentación de la herramienta")

    st.markdown("""
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

    1. Abrir la aplicación web en la URL proporcionada.
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
    """)