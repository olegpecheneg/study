import pandas as pd
import numpy as np
from typing import List, Tuple
import os

def find_top_fdr_rows(
    file_path: str, 
    search_criteria: List[Tuple[str, str, Tuple[int, int]]]
) -> None:
    """
    Ищет топ-3 УНИКАЛЬНЫХ ПО ПОЗИЦИИ строк по максимальному FDR для каждого набора критериев.
    Сохраняет результаты в CSV в той же папке.
    Сравнение аллелей НЕЗАВИСИМО от регистра.
    
    Args:
        file_path: путь к Excel файлу
        search_criteria: список критериев в формате [(ref, alt, (start, end)), ...]
    """
    
    # 1. Читаем Excel файл
    df = pd.read_excel(file_path)
    
    # 2. Конвертируем нужные колонки в правильные типы
    # Position -> int
    df['Position'] = pd.to_numeric(df['Position'], errors='coerce').astype('Int64')
    
    # FDR -> float (заменяем запятые на точки, если есть)
    if df['FDR'].dtype == object:
        df['FDR'] = (
            df['FDR']
            .astype(str)
            .str.replace(',', '.', regex=False)
            .str.replace(' ', '', regex=False)
        )
    df['FDR'] = pd.to_numeric(df['FDR'], errors='coerce')
    
    # 3. Создаем копии колонок с аллелями в верхнем регистре для сравнения
    # НЕ меняем исходные колонки, создаем временные
    df['Allele1_upper'] = df['Allele1'].astype(str).str.upper()
    df['Allele2_upper'] = df['Allele2'].astype(str).str.upper()
    df['Minor allele_upper'] = df['Minor allele'].astype(str).str.upper()
    
    # 4. Подготовка списка для результатов
    all_results = []
    
    # 5. Поиск для каждого набора критериев
    for ref, alt, (start_pos, end_pos) in search_criteria:
        # Приводим критерии поиска к верхнему регистру
        ref_upper = ref.upper()
        alt_upper = alt.upper()
        
        # Фильтрация по позиции
        mask_pos = (df['Position'] >= start_pos) & (df['Position'] <= end_pos)
        
        # Фильтрация по аллелям (НЕЗАВИСИМО ОТ РЕГИСТРА)
        mask_alleles = (
            (df['Minor allele_upper'] == alt_upper) & 
            ((df['Allele1_upper'] == ref_upper) | (df['Allele2_upper'] == ref_upper))
        )
        
        # Применяем все фильтры
        filtered = df[mask_pos & mask_alleles].copy()
        
        if not filtered.empty:
            # УДАЛЯЕМ ДУБЛИКАТЫ ПО ПОЗИЦИИ, оставляя строку с МАКСИМАЛЬНЫМ FDR для каждой позиции
            # Сначала сортируем по FDR по убыванию, чтобы при удалении дубликатов оставались строки с max FDR
            filtered_sorted = filtered.sort_values('FDR', ascending=False)
            # Удаляем дубликаты по позиции, оставляя первую (с максимальным FDR)
            filtered_unique = filtered_sorted.drop_duplicates(subset=['Position'], keep='first')
            
            # Берем топ-3 УНИКАЛЬНЫХ ПО ПОЗИЦИИ строк
            top_rows = filtered_unique.head(3)
            
            # Добавляем в результаты
            for _, row in top_rows.iterrows():
                all_results.append({
                    'pos': row['Position'],
                    'ref': ref_upper,
                    'alt': alt_upper
                })
    
    # 6. Сохраняем результаты в CSV
    if all_results:
        results_df = pd.DataFrame(all_results)
        
        # Формируем путь для сохранения
        dir_path = os.path.dirname(file_path)
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        output_path = os.path.join(dir_path, f"{base_name}_results.csv")
        
        results_df.to_csv(output_path, index=False)
        print(f"Результаты сохранены в: {output_path}")
        print(f"Найдено уникальных строк: {len(results_df)}")
        
        # Выводим краткую статистику
        print("\nСтатистика по критериям:")
        for i, (ref, alt, (start, end)) in enumerate(search_criteria):
            count = len([r for r in all_results if r['ref'] == ref.upper() and r['alt'] == alt.upper()])
            print(f"Критерий {i+1}: ref={ref}, alt={alt}, pos=({start}-{end}) -> найдено: {count} уникальных позиций")
        
        return output_path
    else:
        print("Не найдено ни одной строки по заданным критериям.")
        return None

# Пример использования:
if __name__ == "__main__":
    # Пример вызова функции
    search_criteria_example = [
        ("G", "A", (7970, 8342)),
        ("T", "C", (8376, 8567)),
        ("C", "T", (12606, 12801)),
        # Добавьте свои критерии здесь
    ]
    
    # Вызов функции с вашим файлом
    result = find_top_fdr_rows("D:/pythonProject/MitoFragility/DataPreparing/raw_data/MitoPhewas_associations.xlsx", search_criteria_example)