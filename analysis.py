import sqlite3
import os
import pandas as pd
from scipy import stats

DB_FILE = "cell_count.db"
OUTPUT_DIR = "outputs"
CELL_TYPES = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def get_frequency_table(conn):
    """Part 2 — percentages per sample."""
    df = pd.read_sql("SELECT sample, b_cell, cd8_t_cell, cd4_t_cell, nk_cell, monocyte FROM samples", conn)
    df["total_count"] = df[CELL_TYPES].sum(axis=1)

    long_df = df.melt(
        id_vars=["sample", "total_count"],
        value_vars=CELL_TYPES,
        var_name="population",
        value_name="count",
    )
    long_df["percentage"] = (long_df["count"] / long_df["total_count"] * 100).round(4)
    return long_df[["sample", "total_count", "population", "count", "percentage"]]


def compare_responders(conn, freq_df):
    """Part 3 — t-test per cell type, melanoma + miraclib + PBMC only."""
    meta = pd.read_sql("SELECT sample, condition, treatment, response, sample_type FROM samples", conn)
    df = freq_df.merge(meta, on="sample")

    df = df[
        (df["condition"] == "melanoma")
        & (df["treatment"] == "miraclib")
        & (df["sample_type"] == "PBMC")
        & (df["response"].isin(["yes", "no"]))
    ]

    results = []
    for cell_type in CELL_TYPES:
        responders = df[(df["population"] == cell_type) & (df["response"] == "yes")]["percentage"]
        non_responders = df[(df["population"] == cell_type) & (df["response"] == "no")]["percentage"]
        _, p_value = stats.ttest_ind(responders, non_responders)
        results.append({
            "population": cell_type,
            "responder_mean_%": round(responders.mean(), 2),
            "non_responder_mean_%": round(non_responders.mean(), 2),
            "p_value": round(p_value, 4),
            "significant": "yes" if p_value < 0.05 else "no",
        })
    return pd.DataFrame(results)


def baseline_subset(conn):
    """Part 4 — baseline melanoma + miraclib + PBMC."""
    df = pd.read_sql("""
        SELECT * FROM samples
        WHERE condition = 'melanoma'
          AND treatment = 'miraclib'
          AND sample_type = 'PBMC'
          AND time_from_treatment_start = 0
    """, conn)

    by_project = df.groupby("project").size().reset_index(name="n_samples")

    subjects = df[["subject", "response", "sex"]].drop_duplicates()
    by_response = subjects["response"].value_counts().reset_index()
    by_response.columns = ["response", "n_subjects"]
    by_sex = subjects["sex"].value_counts().reset_index()
    by_sex.columns = ["sex", "n_subjects"]

    male_responders = df[(df["sex"] == "M") & (df["response"] == "yes")]
    avg_b_cells = round(male_responders["b_cell"].mean(), 2)

    return by_project, by_response, by_sex, avg_b_cells


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    conn = sqlite3.connect(DB_FILE)

    freq = get_frequency_table(conn)
    freq.to_csv(f"{OUTPUT_DIR}/frequency_table.csv", index=False)
    print(f"Part 2: frequency_table.csv ({len(freq)} rows)")

    stats_table = compare_responders(conn, freq)
    stats_table.to_csv(f"{OUTPUT_DIR}/responder_stats.csv", index=False)
    print("\nPart 3:")
    print(stats_table.to_string(index=False))

    by_project, by_response, by_sex, avg_b = baseline_subset(conn)
    by_project.to_csv(f"{OUTPUT_DIR}/subset_by_project.csv", index=False)
    by_response.to_csv(f"{OUTPUT_DIR}/subset_by_response.csv", index=False)
    by_sex.to_csv(f"{OUTPUT_DIR}/subset_by_sex.csv", index=False)

    print("\nPart 4:")
    print(by_project.to_string(index=False))
    print()
    print(by_response.to_string(index=False))
    print()
    print(by_sex.to_string(index=False))
    print(f"\nAvg B cells (melanoma males, responders, baseline): {avg_b}")

    conn.close()


if __name__ == "__main__":
    main()
