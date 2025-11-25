"""Инициализация логирования и matplotlib."""

import logging
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt


def setup_matplotlib() -> None:
    """Настройка бэкенда matplotlib и шрифтов по умолчанию."""
    matplotlib.use('Agg')
    plt.rcParams['font.family'] = 'DejaVu Sans'


def setup_logger(name: str = "visualization") -> logging.Logger:
    """Настройка и возврат логгера для модуля.

    Args:
        name: Название логгера.

    Returns:
        Сконфигурированный logger.
    """
    logger = logging.getLogger(name)
    if not logger.hasHandlers():
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(
            logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        )
        
        # Основной лог файл
        base_dir = Path(__file__).resolve().parent.parent
        log_path = base_dir / 'visualization.log'
        file_handler = logging.FileHandler(log_path, encoding='utf-8')
        file_handler.setFormatter(
            logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        )
        
        logger.addHandler(stream_handler)
        logger.addHandler(file_handler)
        logger.setLevel(logging.INFO)
        
        # Отдельный файл только для ошибок
        error_log_path = base_dir / 'errors.log'
        error_handler = logging.FileHandler(error_log_path, encoding='utf-8')
        error_handler.setLevel(logging.ERROR)
        error_handler.setFormatter(
            logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        )
        logger.addHandler(error_handler)
    
    return logger
