import sqlite3
import pandas as pd
import plotly.express as px
import streamlit as st
from scipy import stats

DB_FILE = "cell_count.db"
CELL_TYPES = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

import os
if not os.path.exists(DB_FILE):
    import load_data
    load_data.main()

st.set_page_config(page_title="Loblaw Bio Dashboard", layout="wide")
st.title("Cell Count Dashboard")


@st.cache_data
def load_data():
    conn = sqlite3.connect(DB_FILE)
    df = pd.read_sql("SELECT * FROM samples", conn)
    conn.close()

    df["total_count"] = df[CELL_TYPES].sum(axis=1)
    for cell in CELL_TYPES:
        df[f"{cell}_pct"] = df[cell] / df["total_count"] * 100
    return df


df = load_data()

tab1, tab2, tab3 = st.tabs([
    "Part 2: Frequencies",
    "Part 3: Responders vs Non-responders",
    "Part 4: Baseline Subset",
])

with tab1:
    st.header("Cell type frequencies per sample")

    long_df = df.melt(
        id_vars=["sample", "total_count"],
        value_vars=CELL_TYPES,
        var_name="population",
        value_name="count",
    )
    long_df["percentage"] = (long_df["count"] / long_df["total_count"] * 100).round(4)
    long_df = long_df[["sample", "total_count", "population", "count", "percentage"]]

    all_samples = sorted(long_df["sample"].unique())
    chosen = st.multiselect("Filter samples", all_samples, default=all_samples[:5])

    if chosen:
        st.dataframe(long_df[long_df["sample"].isin(chosen)], use_container_width=True)
    else:
        st.dataframe(long_df, use_container_width=True, height=500)

with tab2:
    st.header("Responders vs non-responders")
    st.caption("Melanoma + miraclib + PBMC")

    sub = df[
        (df["condition"] == "melanoma")
        & (df["treatment"] == "miraclib")
        & (df["sample_type"] == "PBMC")
        & (df["response"].isin(["yes", "no"]))
    ].copy()

    plot_df = sub.melt(
        id_vars=["sample", "response"],
        value_vars=[f"{c}_pct" for c in CELL_TYPES],
        var_name="population",
        value_name="percentage",
    )
    plot_df["population"] = plot_df["population"].str.replace("_pct", "")
    plot_df["Response"] = plot_df["response"].map({"yes": "Responder", "no": "Non-responder"})

    fig = px.box(
        plot_df,
        x="population", y="percentage", color="Response",
        category_orders={"population": CELL_TYPES},
        labels={"percentage": "Percentage (%)", "population": "Cell type"},
    )
    st.plotly_chart(fig, use_container_width=True)

    st.subheader("T-test results")
    results = []
    for cell_type in CELL_TYPES:
        r = sub[sub["response"] == "yes"][f"{cell_type}_pct"]
        nr = sub[sub["response"] == "no"][f"{cell_type}_pct"]
        _, p = stats.ttest_ind(r, nr)
        results.append({
            "population": cell_type,
            "responder_mean_%": round(r.mean(), 2),
            "non_responder_mean_%": round(nr.mean(), 2),
            "p_value": round(p, 4),
            "significant": "yes" if p < 0.05 else "no",
        })
    st.dataframe(pd.DataFrame(results), use_container_width=True, hide_index=True)

with tab3:
    st.header("Baseline subset")
    st.caption("Melanoma + miraclib + PBMC at time=0")

    baseline = df[
        (df["condition"] == "melanoma")
        & (df["treatment"] == "miraclib")
        & (df["sample_type"] == "PBMC")
        & (df["time_from_treatment_start"] == 0)
    ]

    col1, col2, col3 = st.columns(3)

    with col1:
        st.subheader("Samples per project")
        st.dataframe(
            baseline.groupby("project").size().reset_index(name="n_samples"),
            hide_index=True,
        )

    subjects = baseline[["subject", "response", "sex"]].drop_duplicates()

    with col2:
        st.subheader("By response")
        by_resp = subjects["response"].value_counts().reset_index()
        by_resp.columns = ["response", "n_subjects"]
        st.dataframe(by_resp, hide_index=True)

    with col3:
        st.subheader("By sex")
        by_sex = subjects["sex"].value_counts().reset_index()
        by_sex.columns = ["sex", "n_subjects"]
        st.dataframe(by_sex, hide_index=True)

    male_resp = baseline[(baseline["sex"] == "M") & (baseline["response"] == "yes")]
    avg = male_resp["b_cell"].mean()
    st.metric(
        "Avg B cells (males, responders, baseline)",
        f"{avg:.2f}",
        help=f"n = {len(male_resp)}",
    )
