import streamlit as st
import pandas as pd
import plotly.express as px
import re
import os

st.set_page_config(page_title="VCF Dashboard", layout="wide")
st.title("🧬 Variant Analysis Dashboard")

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
VCF_PATH = os.path.join(BASE_DIR, "calls", "all.vcf")
# ============ LEER VCF ============
@st.cache_data
def load_vcf(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 8:
                continue
            
            chrom, pos, _, ref, alt, qual, filt, info = cols[:8]
            
            # extraer DP del campo INFO
            dp_match = re.search(r"DP=(\d+)", info)
            dp = int(dp_match.group(1)) if dp_match else 0
            
            # tipo de variant
            if len(ref) != len(alt):
                vtype = "INDEL"
            else:
                vtype = "SNP"
            
            rows.append({
                "CHROM": chrom,
                "POS": int(pos),
                "REF": ref,
                "ALT": alt,
                "QUAL": float(qual) if qual != "." else 0,
                "FILTER": filt,
                "DP": dp,
                "TYPE": vtype
            })
    
    df = pd.DataFrame(rows)
    
    # simular grupos responders / non-responders
    import random
    random.seed(42)
    df["RESPONSE"] = random.choices(
        ["Responder", "Non-responder"], 
        k=len(df)
    )
    return df

df = load_vcf(VCF_PATH)

# ============ SIDEBAR FILTROS ============
st.sidebar.header("Filtros")

qual_min = st.sidebar.slider(
    "QUAL mínimo", 
    min_value=0, 
    max_value=100, 
    value=0
)

dp_min = st.sidebar.slider(
    "DP mínimo (cobertura)", 
    min_value=0, 
    max_value=10, 
    value=0
)

tipo = st.sidebar.multiselect(
    "Tipo de variante",
    options=["SNP", "INDEL"],
    default=["SNP", "INDEL"]
)

respuesta = st.sidebar.multiselect(
    "Grupo de respuesta",
    options=["Responder", "Non-responder"],
    default=["Responder", "Non-responder"]
)

# aplicar filtros
df_f = df[
    (df["QUAL"] >= qual_min) &
    (df["DP"] >= dp_min) &
    (df["TYPE"].isin(tipo)) &
    (df["RESPONSE"].isin(respuesta))
]

# ============ METRICAS ============
col1, col2, col3, col4 = st.columns(4)
col1.metric("Total variantes", len(df_f))
col2.metric("SNPs", len(df_f[df_f["TYPE"] == "SNP"]))
col3.metric("INDELs", len(df_f[df_f["TYPE"] == "INDEL"]))
col4.metric("QUAL promedio", f"{df_f['QUAL'].mean():.1f}")

st.divider()

# ============ GRAFICOS ============
col_a, col_b = st.columns(2)

with col_a:
    st.subheader("Distribución de calidad (QUAL)")
    fig1 = px.histogram(
        df_f, x="QUAL",
        color="RESPONSE",
        color_discrete_map={
            "Responder": "#1D9E75",
            "Non-responder": "#D85A30"
        },
        barmode="overlay",
        opacity=0.7
    )
    st.plotly_chart(fig1, use_container_width=True)

with col_b:
    st.subheader("Variantes por posición")
    fig2 = px.scatter(
        df_f, x="POS", y="QUAL",
        color="RESPONSE",
        symbol="TYPE",
        color_discrete_map={
            "Responder": "#1D9E75",
            "Non-responder": "#D85A30"
        },
        hover_data=["REF", "ALT", "DP"]
    )
    st.plotly_chart(fig2, use_container_width=True)

col_c, col_d = st.columns(2)

with col_c:
    st.subheader("Cobertura (DP) por variante")
    fig3 = px.box(
        df_f, x="RESPONSE", y="DP",
        color="RESPONSE",
        color_discrete_map={
            "Responder": "#1D9E75",
            "Non-responder": "#D85A30"
        }
    )
    st.plotly_chart(fig3, use_container_width=True)

with col_d:
    st.subheader("SNPs vs INDELs por grupo")
    counts = df_f.groupby(
        ["RESPONSE", "TYPE"]
    ).size().reset_index(name="count")
    fig4 = px.bar(
        counts, x="RESPONSE", y="count",
        color="TYPE", barmode="group"
    )
    st.plotly_chart(fig4, use_container_width=True)

st.divider()

# ============ TABLA ============
st.subheader("Tabla de variantes")
st.dataframe(
    df_f[["CHROM","POS","REF","ALT","QUAL","DP","TYPE","RESPONSE"]],
    use_container_width=True
)