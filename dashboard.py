"""
dashboard.py - Interactive Streamlit dashboard for Teiko Bio clinical trial analysis.
Run with: streamlit run dashboard.py

"""

import sqlite3
import os
import pandas as pd
import numpy as np
from scipy import stats
import streamlit as st
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import warnings
warnings.filterwarnings("ignore")

# Config
DB_PATH = "cell_counts.db"
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
POP_LABELS = {
    "b_cell": "B Cell", "cd8_t_cell": "CD8 T Cell",
    "cd4_t_cell": "CD4 T Cell", "nk_cell": "NK Cell", "monocyte": "Monocyte",
}

RESPONDER_COLOR    = "#2E86AB"
NONRESPONDER_COLOR = "#E84855"
ACCENT            = "#9D2A88"
BG                = "#0F1923"

st.set_page_config(
    page_title="Teiko Bio · Immune Cell Dashboard",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
  [data-testid="stAppViewContainer"] { background: #000000; color: #F4F4F4; }
  [data-testid="stSidebar"] { background: #898989; border-right: 1px solid #E0E0E0; }
  [data-testid="stSidebar"] * { color: #000000 !important; }
  h1, h2, h3 { color: #F4F4F4 !important; }
  .metric-card {
    background: #162230; border-radius: 10px; padding: 1rem 1.4rem;
    border: 1px solid #1E2F3F; margin-bottom: 0.5rem;
  }
  .metric-val { font-size: 2rem; font-weight: 700; color: #F4F4F4; }
  .metric-lbl { font-size: 0.85rem; color: #8ba0b4; margin-top: -4px; }
  .sig-yes { color: #4ECDC4; font-weight: 700; }
  .sig-no  { color: #8ba0b4; }
  [data-testid="stDataFrame"] { background: #162230; }
  div[data-testid="stSelectbox"] label { color: #000000 !important; }
  div[data-testid="stMultiSelect"] label { color: #000000 !important; }
  .stTabs [data-baseweb="tab"] { color: #8ba0b4; }
  .stTabs [aria-selected="true"] { color: #F4F4F4 !important; border-bottom-color: #4ECDC4 !important; }
  span[data-baseweb="tag"] { background-color: #000000 !important; }
  span[data-baseweb="tag"] span { color: #F4F4F4 !important; }
  [data-testid="stSidebar"] [data-baseweb="tag"] { background-color: #000000 !important; }
  [data-testid="stSidebar"] [data-baseweb="tag"] span { color: #F4F4F4 !important; }

</style>
""", unsafe_allow_html=True)




@st.cache_resource
def get_conn():
    return sqlite3.connect(DB_PATH, check_same_thread=False)


@st.cache_data
def load_frequency_table():
    conn = get_conn()
    query = """
        SELECT s.sample_id AS sample, sub.project_id,
               sub.condition, sub.treatment, sub.response,
               sub.sex, sub.age, s.sample_type,
               s.time_from_treatment_start,
               cc.population, cc.count
        FROM cell_counts cc
        JOIN samples s  ON cc.sample_id  = s.sample_id
        JOIN subjects sub ON s.subject_id = sub.subject_id
    """
    df = pd.read_sql_query(query, conn)
    totals = df.groupby("sample")["count"].sum().rename("total_count")
    df = df.join(totals, on="sample")
    df["percentage"] = (df["count"] / df["total_count"] * 100).round(3)
    return df


@st.cache_data
def load_part4_data():
    conn = get_conn()
    q = """
        SELECT s.sample_id, sub.subject_id, sub.project_id,
               sub.condition, sub.treatment, sub.response, sub.sex, sub.age,
               s.sample_type, s.time_from_treatment_start
        FROM samples s
        JOIN subjects sub ON s.subject_id = sub.subject_id
        WHERE sub.condition='melanoma' AND sub.treatment='miraclib'
          AND s.sample_type='PBMC' AND s.time_from_treatment_start=0
    """
    return pd.read_sql_query(q, conn)


def make_plotly_layout(title="", height=420):
    return dict(
        title=dict(text=title, font=dict(color="#1A1A2E", size=14)),
        paper_bgcolor="#FFFFFF", plot_bgcolor="#F8FAFB",
        font=dict(color="#1A1A2E"),
        height=height,
        xaxis=dict(gridcolor="#E0E0E0", zeroline=False, tickfont=dict(color="#1A1A2E"), title_font=dict(color="#1A1A2E")),
        yaxis=dict(gridcolor="#E0E0E0", zeroline=False, tickfont=dict(color="#1A1A2E"), title_font=dict(color="#1A1A2E")),
        legend=dict(bgcolor="rgba(0,0,0,0)", bordercolor="#E0E0E0", font=dict(color="#1A1A2E")),
        margin=dict(l=50, r=20, t=50, b=50),
    )

# Slidebar
st.sidebar.markdown("# Teiko Bio Dashboard")
st.sidebar.markdown("---")

df_all = load_frequency_table()

all_projects   = sorted(df_all["project_id"].unique())
all_conditions = sorted(df_all["condition"].dropna().unique())
all_treatments = sorted(df_all["treatment"].dropna().unique())
all_sampletypes = sorted(df_all["sample_type"].dropna().unique())

sel_projects   = st.sidebar.multiselect("Projects",   all_projects,   default=all_projects)
sel_conditions = st.sidebar.multiselect("Conditions", all_conditions, default=all_conditions)
sel_treatments = st.sidebar.multiselect("Treatments", all_treatments, default=all_treatments)
sel_sampletypes = st.sidebar.multiselect("Sample Type", all_sampletypes, default=all_sampletypes)

df_filtered = df_all[
    df_all["project_id"].isin(sel_projects) &
    df_all["condition"].isin(sel_conditions) &
    df_all["treatment"].isin(sel_treatments) &
    df_all["sample_type"].isin(sel_sampletypes)
]

n_samples  = df_filtered["sample"].nunique()
n_subjects = df_all[df_all["sample"].isin(df_filtered["sample"].unique())]["sample"].nunique()

st.sidebar.markdown("---")
st.sidebar.markdown(f"**{n_samples:,}** samples selected")
st.sidebar.markdown(f"**{df_filtered['project_id'].nunique()}** projects")

# Header
st.title("Immune Cell Population Dashboard")
st.markdown("*Clinical Trial Analysis - Miraclib Drug Candidate*")

col1, col2, col3, col4 = st.columns(4)
metrics = [
    (col1, str(df_filtered["sample"].nunique()),     "Samples"),
    (col2, str(df_filtered["project_id"].nunique()), "Projects"),
    (col3, str(len(sel_conditions)),                 "Conditions"),
    (col4, str(len(sel_treatments)),                 "Treatments"),
]
for col, val, lbl in metrics:
    with col:
        st.markdown(f'<div class="metric-card"><div class="metric-val">{val}</div><div class="metric-lbl">{lbl}</div></div>', unsafe_allow_html=True)

st.markdown("---")

# Tabs
tab1, tab2, tab3, tab4 = st.tabs([" Overview", " Frequencies", " Statistics", " Subset Analysis"])

# Tab 1: Overview

with tab1:
    st.subheader("Data Overview")
    c1, c2 = st.columns(2)

    with c1:
        cond_counts = df_filtered.drop_duplicates("sample").groupby("condition").size().reset_index(name="n")
        fig = go.Figure(go.Bar(
            x=cond_counts["condition"], y=cond_counts["n"],
            marker_color=ACCENT, opacity=0.85,
        ))
        fig.update_layout(**make_plotly_layout("Samples by Condition"), showlegend=False)
        st.plotly_chart(fig, use_container_width=True)

    with c2:
        treat_counts = df_filtered.drop_duplicates("sample").groupby("treatment").size().reset_index(name="n")
        fig = go.Figure(go.Bar(
            x=treat_counts["treatment"], y=treat_counts["n"],
            marker_color=RESPONDER_COLOR, opacity=0.85,
        ))
        fig.update_layout(**make_plotly_layout("Samples by Treatment"), showlegend=False)
        st.plotly_chart(fig, use_container_width=True)

    st.subheader("Mean Cell Population Composition")
    mean_pct = df_filtered.groupby("population")["percentage"].mean().reset_index()
    mean_pct["population_label"] = mean_pct["population"].map(POP_LABELS)
    fig = go.Figure(go.Pie(
        labels=mean_pct["population_label"],
        values=mean_pct["percentage"].round(2),
        hole=0.42,
        marker=dict(colors=["#2E86AB","#F34343","#4ECD70","#F7B731","#9D2A88"]),
        textfont=dict(color="#000000"),
    ))
    fig.update_layout(**make_plotly_layout("Average Relative Frequency (all selected samples)", height=380))
    st.plotly_chart(fig, use_container_width=True)

# Tab 2: Frequency Table
with tab2:
    st.subheader("Part 2 · Relative Frequencies per Sample")

    display_df = df_filtered[["sample","total_count","population","count","percentage"]].copy()
    display_df["population"] = display_df["population"].map(POP_LABELS)
    display_df = display_df.sort_values(["sample","population"]).reset_index(drop=True)

    st.dataframe(display_df, use_container_width=True, height=400)
    st.download_button("Download CSV", display_df.to_csv(index=False),
                       "frequency_table.csv", "text/csv")

    st.subheader("Distribution of Relative Frequencies")
    pop_sel = st.selectbox("Population", POPULATIONS, format_func=lambda x: POP_LABELS[x])
    pop_data = df_filtered[df_filtered["population"] == pop_sel]

    fig = go.Figure()
    groups = [("yes", "Responders", RESPONDER_COLOR), ("no", "Non-Responders", NONRESPONDER_COLOR)]
    for resp_val, label, color in groups:
        vals = pop_data[pop_data["response"] == resp_val]["percentage"]
        if len(vals):
            fig.add_trace(go.Histogram(x=vals, name=label, marker_color=color, opacity=0.7, nbinsx=30))
    fig.update_layout(**make_plotly_layout(f"{POP_LABELS[pop_sel]} — Frequency Distribution"))
    st.plotly_chart(fig, use_container_width=True)

# Tab 3: Statistics
with tab3:
    st.subheader("Part 3 · Responders vs Non-Responders · Melanoma · Miraclib · PBMC")

    miracle_df = df_all[
        (df_all["condition"]   == "melanoma") &
        (df_all["treatment"]   == "miraclib") &
        (df_all["sample_type"] == "PBMC") &
        (df_all["response"].isin(["yes","no"]))
    ]

    # Boxplot
    fig = make_subplots(rows=1, cols=5, subplot_titles=[POP_LABELS[p] for p in POPULATIONS])
    for i, pop in enumerate(POPULATIONS, 1):
        pop_df = miracle_df[miracle_df["population"] == pop]
        for resp_val, label, color in [("yes","Responders",RESPONDER_COLOR),
                                        ("no","Non-Responders",NONRESPONDER_COLOR)]:
            vals = pop_df[pop_df["response"] == resp_val]["percentage"]
            fig.add_trace(go.Box(
                y=vals, name=label, marker_color=color,
                boxmean=True, showlegend=(i == 1),
                legendgroup=label,
            ), row=1, col=i)

    fig.update_layout(
        **make_plotly_layout("Cell Frequencies: Responders vs Non-Responders", height=480),
        boxmode="group",
    )
    for i in range(1, 6):
        fig.update_yaxes(gridcolor="#000000", row=1, col=i)
        fig.update_xaxes(showticklabels=False, row=1, col=i)
    st.plotly_chart(fig, use_container_width=True)

    # Stats table
    results = []
    for pop in POPULATIONS:
        pop_df = miracle_df[miracle_df["population"] == pop]
        resp    = pop_df[pop_df["response"] == "yes"]["percentage"]
        nonresp = pop_df[pop_df["response"] == "no"]["percentage"]
        if len(resp) < 2 or len(nonresp) < 2:
            continue
        stat, pval = stats.mannwhitneyu(resp, nonresp, alternative="two-sided")
        results.append({
            "Population":          POP_LABELS[pop],
            "Resp Median %":       round(resp.median(), 2),
            "Non-Resp Median %":   round(nonresp.median(), 2),
            "n Responders":        len(resp),
            "n Non-Responders":    len(nonresp),
            "Mann-Whitney U":      round(stat, 1),
            "p-value":             round(pval, 6),
            "Significant p<0.05":  "YES" if pval < 0.05 else " no",
        })

    stats_df = pd.DataFrame(results)
    st.dataframe(stats_df, use_container_width=True)

    sig_pops = [r["Population"] for r in results if "YES" in r["Significant p<0.05"]]
    if sig_pops:
        st.markdown(f'<div style="background:#898989; color:#F4F4F4; padding:0.75rem 1rem; border-radius:6px; font-weight:600;">Significant differences found in: {", ".join(sig_pops)}</div>', unsafe_allow_html=True)
    else:
        st.info("No populations showed significant differences at p<0.05.")

    st.caption(
        "Mann-Whitney U test (two-sided). Compares relative frequencies of each immune cell "
        "population between responders and non-responders in melanoma patients treated with miraclib (PBMC only)."
    )

# Tab 4: Subset Analysis
with tab4:
    st.subheader("Part 4 -Melanoma - Miraclib - PBMC - Baseline (time=0)")

    part4_df = load_part4_data()
    st.markdown(f"**Total baseline samples:** {len(part4_df)}")

    c1, c2, c3 = st.columns(3)

    with c1:
        st.markdown("#### Samples per Project")
        by_proj = part4_df.groupby("project_id").size().reset_index(name="n_samples")
        st.dataframe(by_proj, use_container_width=True, hide_index=True)

    with c2:
        st.markdown("#### Subjects by Response")
        subj_uniq = part4_df.drop_duplicates("subject_id")
        by_resp = subj_uniq["response"].value_counts().reset_index()
        by_resp.columns = ["response", "n_subjects"]
        by_resp["response"] = by_resp["response"].map({"yes":"Responder","no":"Non-Responder"})
        st.dataframe(by_resp, use_container_width=True, hide_index=True)

    with c3:
        st.markdown("#### Subjects by Sex")
        by_sex = subj_uniq["sex"].value_counts().reset_index()
        by_sex.columns = ["sex", "n_subjects"]
        by_sex["sex"] = by_sex["sex"].map({"M":"Male","F":"Female"})
        st.dataframe(by_sex, use_container_width=True, hide_index=True)

    # Average B cells
    conn = get_conn()
    b_query = """
        SELECT cc.count AS b_cell_count, s.sample_id
        FROM cell_counts cc
        JOIN samples s ON cc.sample_id = s.sample_id
        WHERE cc.population = 'b_cell'
    """
    b_df = pd.read_sql_query(b_query, conn)
    male_resp = part4_df[(part4_df["sex"] == "M") & (part4_df["response"] == "yes")]
    merged_b = male_resp.merge(b_df, on="sample_id")
    avg_b = merged_b["b_cell_count"].mean()

    st.markdown("---")
    st.markdown("#### Average B Cells -Melanoma Males -Responders -Baseline")
    st.markdown(
        f'<div class="metric-card" style="max-width:300px">'
        f'<div class="metric-val">{avg_b:.2f}</div>'
        f'<div class="metric-lbl">cells -n={len(merged_b)} samples</div>'
        f'</div>',
        unsafe_allow_html=True,
    )

    # Visualize B cells distribution
    fig = go.Figure()
    male_resp_vals = merged_b["b_cell_count"]
    male_nonresp = part4_df[(part4_df["sex"] == "M") & (part4_df["response"] == "no")]
    merged_b_nr = male_nonresp.merge(b_df, on="sample_id")

    fig.add_trace(go.Box(y=male_resp_vals, name="Responders", marker_color=RESPONDER_COLOR, boxmean=True))
    fig.add_trace(go.Box(y=merged_b_nr["b_cell_count"], name="Non-Responders", marker_color=NONRESPONDER_COLOR, boxmean=True))
    fig.update_layout(**make_plotly_layout("B Cell Counts -Melanoma Males -Baseline -Miraclib", height=380))
    st.plotly_chart(fig, use_container_width=True)
