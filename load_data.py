import sqlite3
import pandas as pd

CSV_FILE = "cell-count.csv"
DB_FILE = "cell_count.db"


def main():
    df = pd.read_csv(CSV_FILE)
    conn = sqlite3.connect(DB_FILE)
    df.to_sql("samples", conn, if_exists="replace", index=False)
    conn.close()
    print(f"Loaded {len(df)} rows into {DB_FILE}")


if __name__ == "__main__":
    main()
