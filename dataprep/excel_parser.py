"""Загрузка и парсинг данных из Excel файлов."""

from typing import Dict, List, Optional, Tuple, Any
from pathlib import Path
import pandas as pd
from logging import getLogger

logger = getLogger(__name__)


def load_excel_file(excel_path: Path, sheet_name: str = 0) -> Optional[pd.DataFrame]:
    """Загрузить данные из Excel файла.

    Args:
        excel_path: Путь к Excel файлу.
        sheet_name: Название или индекс листа.

    Returns:
        DataFrame с данными или None в случае ошибки.
    """
    try:
        excel_path = Path(excel_path)
        if not excel_path.exists():
            logger.error(f"Excel файл не найден: {excel_path}")
            return None
        
        df = pd.read_excel(excel_path, sheet_name=sheet_name)
        logger.info(f"Загружено {len(df)} строк из {excel_path}")
        return df
    except Exception as e:
        logger.error(f"Ошибка загрузки Excel файла {excel_path}: {e}")
        return None


def get_sheet_names(excel_path: Path) -> List[str]:
    """Получить список названий листов в Excel файле.

    Args:
        excel_path: Путь к Excel файлу.

    Returns:
        Список названий листов.
    """
    try:
        excel_path = Path(excel_path)
        with pd.ExcelFile(excel_path) as xls:
            return xls.sheet_names
    except Exception as e:
        logger.error(f"Ошибка получения списка листов: {e}")
        return []


def load_all_sheets(excel_path: Path) -> Dict[str, pd.DataFrame]:
    """Загрузить все листы из Excel файла.

    Args:
        excel_path: Путь к Excel файлу.

    Returns:
        Словарь {название_листа: DataFrame}.
    """
    try:
        excel_path = Path(excel_path)
        sheets = pd.read_excel(excel_path, sheet_name=None)
        logger.info(f"Загружено {len(sheets)} листов из {excel_path}")
        return sheets
    except Exception as e:
        logger.error(f"Ошибка загрузки всех листов: {e}")
        return {}


def extract_column_values(
    df: pd.DataFrame,
    column_name: str,
    skip_nan: bool = True
) -> List[Any]:
    """Извлечить значения из столбца DataFrame.

    Args:
        df: DataFrame.
        column_name: Название столбца.
        skip_nan: Пропускать ли NaN значения.

    Returns:
        Список значений.
    """
    try:
        if column_name not in df.columns:
            logger.warning(f"Столбец '{column_name}' не найден. Доступные: {list(df.columns)}")
            return []
        
        values = df[column_name]
        if skip_nan:
            values = values.dropna()
        
        return values.tolist()
    except Exception as e:
        logger.error(f"Ошибка извлечения значений из столбца '{column_name}': {e}")
        return []


def filter_dataframe_by_column_value(
    df: pd.DataFrame,
    column_name: str,
    value: Any,
    exact_match: bool = True
) -> pd.DataFrame:
    """Отфильтровать DataFrame по значению в столбце.

    Args:
        df: DataFrame.
        column_name: Название столбца для фильтрации.
        value: Значение для поиска.
        exact_match: Точное ли совпадение или содержание.

    Returns:
        Отфильтрованный DataFrame.
    """
    try:
        if column_name not in df.columns:
            logger.warning(f"Столбец '{column_name}' не найден")
            return df.copy()
        
        if exact_match:
            return df[df[column_name] == value].copy()
        else:
            # Для строковых значений используем contains
            if df[column_name].dtype == 'object':
                return df[df[column_name].astype(str).str.contains(str(value), case=False, na=False)].copy()
            else:
                return df[df[column_name] == value].copy()
    except Exception as e:
        logger.error(f"Ошибка фильтрации DataFrame: {e}")
        return df.copy()


def merge_dataframes(
    dfs: List[pd.DataFrame],
    on: Optional[str] = None,
    how: str = 'outer'
) -> Optional[pd.DataFrame]:
    """Объединить несколько DataFrame.

    Args:
        dfs: Список DataFrames.
        on: Столбец для слияния (если None, используется индекс).
        how: Тип слияния ('inner', 'outer', 'left', 'right').

    Returns:
        Объединённый DataFrame или None.
    """
    try:
        if not dfs:
            logger.warning("Список DataFrame пуст")
            return None
        
        if len(dfs) == 1:
            return dfs[0].copy()
        
        result = dfs[0]
        for df in dfs[1:]:
            if on:
                result = pd.merge(result, df, on=on, how=how)
            else:
                result = pd.concat([result, df], axis=0, ignore_index=True)
        
        logger.info(f"Объединено {len(dfs)} DataFrame")
        return result
    except Exception as e:
        logger.error(f"Ошибка объединения DataFrame: {e}")
        return None


def validate_required_columns(
    df: pd.DataFrame,
    required_columns: List[str]
) -> Tuple[bool, List[str]]:
    """Проверить наличие необходимых столбцов в DataFrame.

    Args:
        df: DataFrame.
        required_columns: Список необходимых столбцов.

    Returns:
        Кортеж (все ли столбцы присутствуют, список отсутствующих столбцов).
    """
    missing = [col for col in required_columns if col not in df.columns]
    
    if missing:
        logger.warning(f"Отсутствуют столбцы: {missing}. Доступные: {list(df.columns)}")
        return False, missing
    
    return True, []


def save_dataframe_to_excel(
    df: pd.DataFrame,
    output_path: Path,
    sheet_name: str = 'Sheet1',
    include_index: bool = False
) -> bool:
    """Сохранить DataFrame в Excel файл.

    Args:
        df: DataFrame для сохранения.
        output_path: Путь к выходному файлу.
        sheet_name: Название листа в Excel.
        include_index: Включать ли индекс.

    Returns:
        True если успешно, False иначе.
    """
    try:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        df.to_excel(output_path, sheet_name=sheet_name, index=include_index, engine='openpyxl')
        logger.info(f"DataFrame сохранён в Excel: {output_path}")
        return True
    except Exception as e:
        logger.error(f"Ошибка сохранения DataFrame в Excel: {e}")
        return False
