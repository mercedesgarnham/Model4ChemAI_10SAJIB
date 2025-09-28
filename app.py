import streamlit as st
import pandas as pd
import requests

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
    st.write("Falta generar el código de análisis de los datos cargados.")
    with st.expander("Métricas básicas"):
        st.write("Cantidad total, activos, inactivos...")

with tab3:
    st.header("Visualización")
    st.write("Falta generar el código de gráficos y visualizaciones.")
    with st.expander("Gráficos"):
        st.write("Visualizaciones de la distribución...")

