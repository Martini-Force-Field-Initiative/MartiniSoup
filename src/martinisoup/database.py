import io

import pandas as pd
import requests

DATABASE_URL = (
    "https://raw.githubusercontent.com/Martini-Force-Field-Initiative/"
    "M3-Metabolome/refs/heads/main/misc/database.csv"
)


def load_metabolite_classes(url: str | None = None, local: str | None = None) -> dict:
    """Return a {resname: class} mapping from either a local CSV or a URL.

    If neither is provided, the M3-Metabolome database is fetched from GitHub.
    """
    if local:
        df = pd.read_csv(local,
                         usecols=['Metabolite name', 'resname', 'class'],
                         index_col='resname')
    else:
        source = url or DATABASE_URL
        response = requests.get(source)
        response.raise_for_status()
        df = pd.read_csv(io.StringIO(response.content.decode('utf-8')),
                         usecols=['Metabolite name', 'resname', 'class'],
                         index_col='resname')
    return df.to_dict()['class']
