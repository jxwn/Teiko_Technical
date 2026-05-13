# Teiko Technical

## Installation

```bash
make setup      # install packages
make pipeline   # load data + run analyses
make dashboard  # open the dashboard
```
Secondary Installation

```bash
pip install -r requirements.txt
python load_data.py
python analysis.py
python -m streamlit run dashboard.py
```

## Database and Scaling

Created a SQLite table that mirrors the CSV columns. Kept it simple by keeping it on one table since the dataset wasn't too large.

If this had to scale to hundreds of projects and thousands of samples I'd split it up into separate tables and add indexes on the columns used for filtering.

## Overview

I used 3 files to make each piece do one job, which in my opinion makes the code easier to read.
I chose to use pandas since it keeps the code shorter and the analysis in the assignment was mostly filtering and grouping.
Lastly I used streamlit for the dashboard because in my classes we used it and it was the simplest way I know of to make a dashboard.

## Dashboard Link
https://john-lam-teiko.streamlit.app/