# utils/data_loader.py
import pandas as pd

def load_csv(file):
    try:
        df = pd.read_csv(file)
        df.columns = df.columns.str.lower()

        rename_headers = {
            "actividad": "activity",
            "molecula": "smiles",
            "moleculas": "smiles",
            "smile": "smiles",
        }

        df.rename(columns=rename_headers, inplace=True)

        return df, None

    except Exception as e:
        return None, str(e)
