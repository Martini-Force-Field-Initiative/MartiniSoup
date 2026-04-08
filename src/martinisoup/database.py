import io

import pandas as pd
import requests

DATABASE_URL = (
    "https://raw.githubusercontent.com/Martini-Force-Field-Initiative/"
    "M3-Metabolome/refs/heads/main/misc/database.csv"
)


def load_metabolite_classes(source: str = DATABASE_URL) -> dict:
    """Return a {resname: class} mapping from a URL or a local CSV path.

    Defaults to fetching the M3-Metabolome database from GitHub.
    """
    if source.startswith("http"):
        response = requests.get(source)
        response.raise_for_status()
        df = pd.read_csv(io.StringIO(response.content.decode('utf-8')),
                         usecols=['Metabolite name', 'resname', 'class'],
                         index_col='resname')
    else:
        df = pd.read_csv(source,
                         usecols=['Metabolite name', 'resname', 'class'],
                         index_col='resname')
    return df.to_dict()['class']
