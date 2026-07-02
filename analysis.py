"""
analysis.py - Parts 2, 3, 4 analysis pipeline.
Reads from the SQLite database, produces tables and plots.
"""

import sqlite3
import os
import pandas as pd 
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings("ignore")

DB_PATH = "cell_counts.db"
OUTPUT_DIR = "outputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
POP_LABELS = {
    "b_cell": "B Cell",
    "cd8_t_cell": "CD8 T Cell",
    "cd4_t_cell": "CD4 T Cell",
    "nk_cell": "NK Cell",
    "monocyte": "Monocyte",
}

# Colors
RESPONDER_COLOR    = "#2E86AB"
NONRESPONDER_COLOR = "#E84855"
BG_COLOR           = "#FFFFFF"
PANEL_COLOR        = "#F8FAFB"
TEXT_COLOR         = "#1A1A2E"
ACCENT_COLOR       = "#E84855"
GRID_COLOR         = "#E0E0E0"



def get_conn() -> sqlite3.Connection:
    return sqlite3.connect(DB_PATH)

# Part 2: Frequency summary table

def part2_frequency_table(conn: sqlite3.Connection) -> pd.DataFrame:
    """Return relative frequency of each cell population per sample."""
    query = """
        SELECT
            s.sample_id  AS sample,
            cc.population,
            cc.count
        FROM cell_counts cc
        JOIN samples s ON cc.sample_id = s.sample_id
    """
    df = pd.read_sql_query(query, conn)

    totals = df.groupby("sample")["count"].sum().rename("total_count")
    df = df.join(totals, on="sample")
    df["percentage"] = (df["count"] / df["total_count"] * 100).round(4)
    df = df[["sample", "total_count", "population", "count", "percentage"]]
    df = df.sort_values(["sample", "population"]).reset_index(drop=True)

    out_path = os.path.join(OUTPUT_DIR, "part2_frequency_table.csv")
    df.to_csv(out_path, index=False)
    print(f"[Part 2] Frequency table saved → {out_path}  ({len(df):,} rows)")
    return df

# Part 3: Statistical analysis - melanoma + miraclib + PBMC

def part3_stats(conn: sqlite3.Connection, freq_df: pd.DataFrame) -> pd.DataFrame:
    """Compare cell frequencies in responders vs non-responders."""

    # Pull relevant metadata
    meta_query = """
        SELECT s.sample_id, s.sample_type, s.time_from_treatment_start,
               sub.condition, sub.treatment, sub.response
        FROM samples s
        JOIN subjects sub ON s.subject_id = sub.subject_id
        WHERE sub.condition = 'melanoma'
          AND sub.treatment = 'miraclib'
          AND s.sample_type = 'PBMC'
          AND sub.response IN ('yes', 'no')
    """
    meta = pd.read_sql_query(meta_query, conn)
    merged = freq_df.merge(meta, left_on="sample", right_on="sample_id", how="inner")

    results = []
    for pop in POPULATIONS:
        pop_df = merged[merged["population"] == pop]
        resp    = pop_df[pop_df["response"] == "yes"]["percentage"]
        nonresp = pop_df[pop_df["response"] == "no"]["percentage"]
        stat, pval = stats.mannwhitneyu(resp, nonresp, alternative="two-sided")
        results.append({
            "population": pop,
            "responder_median_%":    round(resp.median(),    2),
            "nonresponder_median_%": round(nonresp.median(), 2),
            "n_responders":          len(resp),
            "n_nonresponders":       len(nonresp),
            "mannwhitney_U":         round(stat, 1),
            "p_value":               round(pval, 6),
            "significant (p<0.05)":  "YES" if pval < 0.05 else "no",
        })

    stats_df = pd.DataFrame(results)
    out_path = os.path.join(OUTPUT_DIR, "part3_stats.csv")
    stats_df.to_csv(out_path, index=False)
    print(f"[Part 3] Statistics saved → {out_path}")
    print(stats_df.to_string(index=False))

# Boxplot

    fig, axes = plt.subplots(1, 5, figsize=(18, 7))
    fig.patch.set_facecolor(BG_COLOR)
    fig.suptitle(
        "Melanoma · Miraclib · PBMC\nCell Population Frequencies: Responders vs Non-Responders",
        color=TEXT_COLOR, fontsize=14, fontweight="bold", y=1.01
    )

    for ax, pop in zip(axes, POPULATIONS):
        pop_df = merged[merged["population"] == pop]
        resp_vals    = pop_df[pop_df["response"] == "yes"]["percentage"].values
        nonresp_vals = pop_df[pop_df["response"] == "no"]["percentage"].values

        ax.set_facecolor(PANEL_COLOR)
        for spine in ax.spines.values():
            spine.set_color(GRID_COLOR)

        bp = ax.boxplot(
            [resp_vals, nonresp_vals],
            patch_artist=True,
            widths=0.55,
            medianprops=dict(color="white", linewidth=2.5),
            whiskerprops=dict(color=TEXT_COLOR, linewidth=1.2),
            capprops=dict(color=TEXT_COLOR, linewidth=1.5),
            flierprops=dict(marker="o", markerfacecolor=TEXT_COLOR,
                            markersize=3, alpha=0.5, linestyle="none"),
        )
        bp["boxes"][0].set_facecolor(RESPONDER_COLOR)
        bp["boxes"][0].set_alpha(0.85)
        bp["boxes"][1].set_facecolor(NONRESPONDER_COLOR)
        bp["boxes"][1].set_alpha(0.85)

        # Jitter overlay
        for i, (vals, color) in enumerate([(resp_vals, RESPONDER_COLOR),
                                            (nonresp_vals, NONRESPONDER_COLOR)], 1):
            jitter = np.random.uniform(-0.12, 0.12, size=len(vals))
            ax.scatter(np.full(len(vals), i) + jitter, vals,
                       color=color, s=12, alpha=0.4, zorder=3)

        # p-value annotation
        row = stats_df[stats_df["population"] == pop].iloc[0]
        pval = row["p_value"]
        sig  = row["significant (p<0.05)"]
        pstr = f"p={pval:.4f}" if pval >= 0.0001 else f"p<0.0001"
        color = ACCENT_COLOR if sig == "YES" else TEXT_COLOR
        ax.set_title(f"{POP_LABELS[pop]}\n{pstr}", color=color,
                     fontsize=10, fontweight="bold")

        ax.set_xticks([1, 2])
        ax.set_xticklabels(["Resp", "Non-R"], color=TEXT_COLOR, fontsize=9)
        ax.tick_params(colors=TEXT_COLOR, labelsize=8)
        ax.set_ylabel("Relative Frequency (%)", color=TEXT_COLOR, fontsize=8)
        ax.yaxis.label.set_color(TEXT_COLOR)
        ax.grid(axis="y", color=GRID_COLOR, linewidth=0.6, alpha=0.7)

    legend_handles = [
        mpatches.Patch(facecolor=RESPONDER_COLOR,    label="Responders"),
        mpatches.Patch(facecolor=NONRESPONDER_COLOR, label="Non-Responders"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=2,
               frameon=False, fontsize=10, labelcolor=TEXT_COLOR,
               bbox_to_anchor=(0.5, -0.04))

    plt.tight_layout()
    plot_path = os.path.join(OUTPUT_DIR, "part3_boxplot.png")
    fig.savefig(plot_path, dpi=150, bbox_inches="tight", facecolor=BG_COLOR)
    plt.close()
    print(f"[Part 3] Boxplot saved → {plot_path}")

    return stats_df

# Part 4: Subset analysis - melanoma, miraclib, PBMC, baseline

def part4_subset(conn: sqlite3.Connection) -> None:
    base_query = """
        SELECT
            s.sample_id, s.sample_type, s.time_from_treatment_start,
            sub.subject_id, sub.project_id, sub.condition,
            sub.treatment, sub.response, sub.sex, sub.age
        FROM samples s
        JOIN subjects sub ON s.subject_id = sub.subject_id
        WHERE sub.condition  = 'melanoma'
          AND sub.treatment  = 'miraclib'
          AND s.sample_type  = 'PBMC'
          AND s.time_from_treatment_start = 0
    """
    df = pd.read_sql_query(base_query, conn)

    print(f"\n[Part 4] Melanoma · miraclib · PBMC · baseline samples: {len(df)}")

    # Samples per project
    by_project = df.groupby("project_id").size().reset_index(name="n_samples")
    print("\nSamples per project:")
    print(by_project.to_string(index=False))
    by_project.to_csv(os.path.join(OUTPUT_DIR, "part4_by_project.csv"), index=False)

    # Responders vs non-responders (unique subjects)
    subj = df.drop_duplicates("subject_id")
    by_response = subj["response"].value_counts().reset_index()
    by_response.columns = ["response", "n_subjects"]
    print("\nSubjects by response:")
    print(by_response.to_string(index=False))
    by_response.to_csv(os.path.join(OUTPUT_DIR, "part4_by_response.csv"), index=False)

    # Males vs females
    by_sex = subj["sex"].value_counts().reset_index()
    by_sex.columns = ["sex", "n_subjects"]
    print("\nSubjects by sex:")
    print(by_sex.to_string(index=False))
    by_sex.to_csv(os.path.join(OUTPUT_DIR, "part4_by_sex.csv"), index=False)

    # Average B cells: melanoma, male, responder, time=0
    male_resp = df[(df["sex"] == "M") & (df["response"] == "yes")]
    b_cells_query = """
        SELECT cc.count AS b_cell_count, s.sample_id
        FROM cell_counts cc
        JOIN samples s ON cc.sample_id = s.sample_id
        WHERE cc.population = 'b_cell'
    """
    b_cells = pd.read_sql_query(b_cells_query, conn)
    merged = male_resp.merge(b_cells, on="sample_id")
    avg_b = merged["b_cell_count"].mean()
    print(f"\nAverage B cells (melanoma males, responders, time=0): {avg_b:.2f}")
    with open(os.path.join(OUTPUT_DIR, "part4_avg_bcell.txt"), "w") as f:
        f.write(f"Average B cells (melanoma males, miraclib, responders, baseline): {avg_b:.2f}\n")
        f.write(f"N samples: {len(merged)}\n")

if __name__=="__main__":
    conn = get_conn()
    freq_df = part2_frequency_table(conn)
    part3_stats(conn, freq_df)
    part4_subset(conn)
    conn.close()
    print("\nAll outputs saved to ./outputs/")