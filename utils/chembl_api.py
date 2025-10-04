# utils/chembl_api.py
from chembl_webresource_client.new_client import new_client
import pandas as pd

def download_chembl_bioassay_data(target_chembl_id, max_records=1000):
    """
    Descarga moléculas y actividades bioquímicas para un target usando chembl_webresource_client.

    Args:
        target_chembl_id (str): ID del target, ej: 'CHEMBL25'.
        max_records (int): máximo número de moléculas a traer.

    Returns:
        pd.DataFrame con columnas ['smiles', 'activity'] (activity en float), o (None, error_str).
    """
    try:
        activity = new_client.activity
        # Filtramos por target_chembl_id
        res = activity.filter(target_chembl_id=target_chembl_id).only(
            ['canonical_smiles', 'standard_value', 'standard_type']
        )

        # res es un iterable (generator). Lo convertimos a lista y limitamos max_records
        records = list(res)[:max_records]

        if not records:
            return None, "No se encontraron actividades para ese target."

        data = []
        for r in records:
            smiles = r.get('canonical_smiles')
            value = r.get('standard_value')

            # Convertimos standard_value a float si es posible
            try:
                value = float(value)
            except (ValueError, TypeError):
                value = None

            if smiles and value is not None:
                data.append({
                    'smiles': smiles,
                    'activity': value,
                    'standard_type': r.get('standard_type')
                })

        if not data:
            return None, "No se encontraron datos válidos."

        df = pd.DataFrame(data)
        return df, None

    except Exception as e:
        return None, f"Error descargando datos de ChEMBL: {str(e)}"
