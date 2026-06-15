import json
import numpy as np
from collections import defaultdict
import logging

logger = logging.getLogger("energy_distribution")

def extract_energy_data_from_master_json(master_json_path: str):
    """Извлекает энергии по типам из master_data.json."""
    try:
        with open(master_json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except Exception as e:
        logger.error(f"Ошибка загрузки {master_json_path}: {e}")
        return {}
    
    energy_types_data = defaultdict(lambda: {'ref': [], 'alt': [], 'diff': []})
    
    for entry in data:
        energy_type = entry.get('energy_type')
        if not energy_type:
            continue
            
        try:
            ref_val = float(entry.get('ref_energy', 0))
            alt_val = float(entry.get('alt_energy', 0))
            energy_types_data[energy_type]['ref'].append(ref_val)
            energy_types_data[energy_type]['alt'].append(alt_val)
            energy_types_data[energy_type]['diff'].append(ref_val - alt_val)
        except (ValueError, TypeError):
            continue
    
    # Конвертируем в numpy arrays
    result = {}
    for energy_type, energies in energy_types_data.items():
        if energies['ref']:
            result[energy_type] = {
                'ref': np.array(energies['ref']),
                'alt': np.array(energies['alt']),
                'diff': np.array(energies['diff'])
            }
    
    return result