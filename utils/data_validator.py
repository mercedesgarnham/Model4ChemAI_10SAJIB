# utils/data_validator.py
from rdkit import Chem

def validate_required_columns(df, required=["smiles", "activity"]):
    missing = [col for col in required if col not in df.columns]
    return missing

def validate_smiles(df):
    invalid_indices = []
    for i, smi in enumerate(df["smiles"]):
        if Chem.MolFromSmiles(str(smi)) is None:
            invalid_indices.append(i)
    return invalid_indices

def clean_activity_column(df):
    df["activity"] = df["activity"].astype(str).str.lower()
    rename_activity = {
        "a": "activo",
        "active": "activo",
        "act": "activo",
        "1": "activo",
        "i": "inactivo",
        "inactive": "inactivo",
        "inact": "inactivo",
        "0": "inactivo"
    }
    df["activity"] = df["activity"].replace(rename_activity)
    return df

def validate_activity_values(df, valid=["activo", "inactivo"]):
    invalid_rows = df[~df["activity"].isin(valid)]
    return invalid_rows.index.tolist()
