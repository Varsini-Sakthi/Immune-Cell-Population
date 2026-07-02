import sqlite3
import csv
import os

DB_PATH = "cell_counts.db"
CSV_PATH = "cell-count.csv"

def init_db(conn: sqlite3.Connection) -> None:
    conn.executescript("""
        PRAGMA foreign_keys = ON;

        CREATE TABLE IF NOT EXISTS projects (
            project_id TEXT PRIMARY KEY
        );

        CREATE TABLE IF NOT EXISTS subjects (
            subject_id   TEXT PRIMARY KEY,
            project_id   TEXT NOT NULL REFERENCES projects(project_id),
            condition    TEXT,
            age          INTEGER,
            sex          TEXT,
            treatment    TEXT,
            response     TEXT  -- 'yes', 'no', or NULL for healthy/no treatment
        );
        CREATE TABLE IF NOT EXISTS samples (
            sample_id                  TEXT PRIMARY KEY,
            subject_id                 TEXT NOT NULL REFERENCES subjects(subject_id),
            sample_type                TEXT,
            time_from_treatment_start  INTEGER
        );
        CREATE TABLE IF NOT EXISTS cell_counts (
            id          INTEGER PRIMARY KEY AUTOINCREMENT,
            sample_id   TEXT NOT NULL REFERENCES samples(sample_id),
            population  TEXT NOT NULL,
            count       INTEGER NOT NULL
        );

        CREATE INDEX IF NOT EXISTS idx_cell_counts_sample ON cell_counts(sample_id);
        CREATE INDEX IF NOT EXISTS idx_samples_subject    ON samples(subject_id);
        CREATE INDEX IF NOT EXISTS idx_subjects_project   ON subjects(project_id);
    """)
    conn.commit()

def load_csv(conn: sqlite3.Connection, csv_path: str) -> None:
    POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

    projects_seen = set()
    subjects_seen = set()
    samples_seen  = set()

    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        cell_rows = []

        for row in reader:
            project_id = row["project"]
            subject_id = row["subject"]
            sample_id  = row["sample"]

            if project_id not in projects_seen:
                conn.execute(
                    "INSERT OR IGNORE INTO projects VALUES (?)", (project_id,)
                )
                projects_seen.add(project_id)

            if subject_id not in subjects_seen:
                conn.execute(
                                       """INSERT OR IGNORE INTO subjects
                       (subject_id, project_id, condition, age, sex, treatment, response)
                       VALUES (?, ?, ?, ?, ?, ?, ?)""",
                    (
                        subject_id,
                        project_id,
                        row["condition"],
                        int(row["age"]) if row["age"] else None,
                        row["sex"],
                        row["treatment"] if row["treatment"] else None,
                        row["response"]  if row["response"]  else None,
                    ),
                )
                subjects_seen.add(subject_id)

            if sample_id not in samples_seen:
                conn.execute(
                    """INSERT OR IGNORE INTO samples
                       (sample_id, subject_id, sample_type, time_from_treatment_start)
                       VALUES (?, ?, ?, ?)""",
                    (
                        sample_id,
                        subject_id,
                        row["sample_type"],
                        int(row["time_from_treatment_start"])
                        if row["time_from_treatment_start"] != ""
                        else None,
                    ),
                )
                samples_seen.add(sample_id)

            for pop in POPULATIONS:
                cell_rows.append((sample_id, pop, int(row[pop])))

        conn.executemany(
            "INSERT INTO cell_counts (sample_id, population, count) VALUES (?, ?, ?)",
            cell_rows,
        )

    conn.commit()
    print(f"Loaded {len(projects_seen)} projects, {len(subjects_seen)} subjects, "
          f"{len(samples_seen)} samples.")
    
if __name__ =="__main__":
    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)
    conn = sqlite3.connect(DB_PATH)
    init_db(conn)
    load_csv(conn, CSV_PATH)
    conn.close()
    print(f"Database created: {DB_PATH}") 
            