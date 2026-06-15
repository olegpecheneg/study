import logging
from pathlib import Path

def setup_logger(name: str = "energy_distribution", log_file: Path = None):
    """Настраивает логгер."""
    logger = logging.getLogger(name)
    
    if not logger.handlers:
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(console_handler)
        
        if log_file:
            file_handler = logging.FileHandler(log_file, encoding='utf-8')
            file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            logger.addHandler(file_handler)
        
        logger.setLevel(logging.INFO)
    
    return logger

def get_mutation_id(master_json_path: str):
    """Извлекает ID мутации из имени файла."""
    from pathlib import Path
    json_path = Path(master_json_path)
    stem = json_path.stem
    
    if stem.endswith('_master_data'):
        return stem.replace('_master_data', '')
    elif '_' in stem:
        return stem.split('_')[0]
    else:
        return stem or "unknown_mutation"