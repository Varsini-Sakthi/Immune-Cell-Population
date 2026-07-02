# Teiko Bio - Immune Cell Population Analysis

A clinical trial analysis pipeline and interactive dashboard built for Bob Loblaw at Loblaw Bio, examining how the drug candidate **miraclib** affects immune cell populations in melanoma and carcinoma patients.

---

## Quick Start (GitHub Codespaces)

```bash
make setup      # Install all dependencies
make pipeline   # Run full data pipeline (Parts 1–4)
make dashboard  # Launch interactive dashboard
```

The dashboard will open at `http://localhost:8501`.

---

## Project Structure

```
teiko/
├── cell-count.csv          # Source data
├── load_data.py            # Part 1: Database initialization and CSV loading
├── analysis.py             # Parts 2–4: Frequency table, statistics, subset analysis
├── dashboard.py            # Interactive Streamlit dashboard
├── requirements.txt        # Python dependencies
├── Makefile                # setup / pipeline / dashboard targets
├── outputs/                # Generated tables and plots (after make pipeline)
│   ├── part2_frequency_table.csv
│   ├── part3_stats.csv
│   ├── part3_boxplot.png
│   ├── part4_by_project.csv
│   ├── part4_by_response.csv
│   ├── part4_by_sex.csv
│   └── part4_avg_bcell.txt
└── README.md
```

---

## Database Schema

```
projects (project_id PK)
    │
    └── subjects (subject_id PK, project_id FK,
    │            condition, age, sex, treatment, response)
    │
    └── samples (sample_id PK, subject_id FK,
    │           sample_type, time_from_treatment_start)
    │
    └── cell_counts (id PK, sample_id FK,
                    population, count)
```

### Design Rationale

**Normalization over flat storage.** The source CSV stores subject metadata redundantly across every sample row. The schema separates `projects`, `subjects`, `samples`, and `cell_counts` so each fact is stored exactly once. This means updating a subject's metadata requires editing a single row rather than thousands.

**Long-format cell counts.** Rather than five separate columns (`b_cell`, `cd8_t_cell`, etc.), cell populations are stored as rows in `cell_counts(sample_id, population, count)`. This makes it trivial to add new cell types without schema changes and enables clean aggregate queries using `GROUP BY population`.

**Scaling to hundreds of projects and thousands of samples.** The schema handles this by design. Projects and subjects each live in their own table. Indexes on `cell_counts(sample_id)`, `samples(subject_id)`, and `subjects(project_id)` keep lookups fast. For very large datasets, partitioning `cell_counts` by `project_id` in a production database like PostgreSQL would allow per-project scans to skip entire partitions.

**For diverse analytics.** The long-format cell_counts table naturally supports:
- Per-population aggregations (median, mean, distribution)
- Joining against any subject-level metadata (treatment, response, sex, condition)
- Time-series queries (`GROUP BY time_from_treatment_start`)
- Cross-project comparisons without schema changes

---

## Code Structure

### `load_data.py`
Initializes the SQLite database with the normalized schema and loads all rows from `cell-count.csv`. Run directly with `python load_data.py`. Creates `cell_counts.db` in the root directory.

### `analysis.py`
Runs the full analytical pipeline across Parts 2–4:
- **Part 2:** Computes relative frequency of each cell population per sample and saves to CSV.
- **Part 3:** Filters to melanoma patients on miraclib (PBMC only), runs Mann-Whitney U tests comparing responders vs non-responders for each cell population, and generates a boxplot.
- **Part 4:** Queries the database for melanoma PBMC baseline samples on miraclib, summarizes by project, response, and sex, and calculates the average B cell count for male responders.

### `dashboard.py`
Interactive Streamlit dashboard with four tabs and real-time sidebar filters. Uses Plotly for all charts. Reads directly from the SQLite database.

---

## Dashboard Description

This interactive dashboard provides a comprehensive view of immune cell population dynamics in a clinical trial evaluating the drug candidate miraclib. Built with Streamlit and Plotly, it allows researchers to explore and filter data across projects, conditions, treatments, and sample types in real time.

The dashboard is organized into four tabs:

**Overview** - High-level summary of the dataset including sample counts by condition and treatment, and the average immune cell composition across all selected samples shown as a donut chart.

**Frequencies** - The full relative frequency table from Part 2, searchable and downloadable as CSV. Includes a histogram showing the distribution of any selected cell population, split by responder status.

**Statistics** - Side-by-side boxplots comparing cell population frequencies between responders and non-responders in melanoma patients treated with miraclib (PBMC samples only). Includes Mann-Whitney U test results for each population, highlighting statistically significant differences (p < 0.05).

**Subset Analysis** - Focused view of melanoma patients on miraclib at baseline (time = 0), broken down by project, response status, and sex. Displays the average B cell count for male responders and a comparative boxplot.

Use the sidebar filters to narrow the analysis by project, condition, treatment, or sample type. All charts update dynamically based on your selection.

---

## Key Findings

### Part 3 - Statistical Analysis
Mann-Whitney U tests (two-sided) were run on all five cell populations comparing responders vs non-responders in melanoma patients treated with miraclib (PBMC samples only):

| Population  | Responder Median % | Non-Responder Median % | p-value | Significant |
|-------------|-------------------|------------------------|---------|-------------|
| B Cell      | 9.43              | 9.79                   | 0.0557  | No          |
| CD8 T Cell  | 24.73             | 24.60                  | 0.6391  | No          |
| CD4 T Cell  | 30.22             | 29.66                  | 0.0133  | **YES**     |
| NK Cell     | 14.51             | 14.80                  | 0.1211  | No          |
| Monocyte    | 19.61             | 19.94                  | 0.1631  | No          |

**CD4 T Cell** showed a statistically significant difference in relative frequency between responders and non-responders (p = 0.013), suggesting it may be a potential predictor of treatment response to miraclib.

### Part 4 - Subset Analysis
- **656** melanoma PBMC baseline samples treated with miraclib
- **prj1:** 384 samples - **prj3:** 272 samples
- **331** responders - **325** non-responders
- **344** males - **312** females
- **Average B cells (melanoma males, responders, time=0): 10401.28**

---

## Dashboard Link

**Live Dashboard:** https://immune-cell-population.streamlit.app/

---

## Author

Varsini Sakthivadivel Ramasamy

MS Bioinformatics

Johns Hopkins University

