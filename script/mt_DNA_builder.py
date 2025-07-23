from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
import pandas as pd
from pathlib import Path
import random
import os
import sys
import logging

# Настройка логирования
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Добавляем путь к модулю для импорта функций
# sys.path.append("D:/pythonProject/MitoFragility/DataPreparing/plots")
from scatter_plus_n_std import parse_construct_id, calculate_arm_ranges

def csv_constructor(excel_path: Path, output_path: Path) -> pd.DataFrame:
    """Создает CSV с SNV из XLSX"""
    # Загрузка данных
    snv_df = pd.read_excel(excel_path)
    
    # Фильтрация по FDR
    filtered_snv_df = snv_df[snv_df['FDR'] < 0.056].copy()
    
    # Определение референсного и альтернативного аллелей
    ref_alleles = []
    alt_alleles = []
    
    for _, row in filtered_snv_df.iterrows():
        # Минорный аллель всегда alt_allele
        alt_allele = row['Minor allele']
        alt_alleles.append(alt_allele)
        
        # Определение референсного аллеля:
        if row['Allele1'] == alt_allele:
            ref_alleles.append(row['Allele2'])
        elif row['Allele2'] == alt_allele:
            ref_alleles.append(row['Allele1'])
        else:
            raise ValueError(
                f"Ошибка в позиции {row['Position']}: "
                f"Минорный аллель '{alt_allele}' не совпадает "
                f"ни с Allele1 ({row['Allele1']}), ни с Allele2 ({row['Allele2']})"
            )
    
    # Добавляем новые столбцы
    filtered_snv_df['ref_allele'] = ref_alleles
    filtered_snv_df['alt_allele'] = alt_alleles
    
    # Сохранение нужных столбцов
    final_df = filtered_snv_df[['Position', 'ref_allele', 'alt_allele']].rename(
        columns={'Position': 'position'}
    )
    
    final_df.to_csv(output_path, index=False)
    logger.info(f"Сохранено {len(final_df)} уникальных SNV в {output_path}")
    return final_df

def get_covered_positions(ref_constructs_dir: str) -> set:
    """Возвращает все позиции, покрытые конструктами референса"""
    # Проверка существования директории
    if not os.path.exists(ref_constructs_dir):
        logger.error(f"Директория с конструктами не существует: {ref_constructs_dir}")
        return set()
    
    # Проверка что это директория
    if not os.path.isdir(ref_constructs_dir):
        logger.error(f"Путь не является директорией: {ref_constructs_dir}")
        return set()
    
    logger.info(f"Сканирую директорию с конструктами: {ref_constructs_dir}")
    
    # 1. Загружаем информацию о конструктах референса
    ref_constructs = []
    file_count = 0
    
    # Получаем список файлов
    files = os.listdir(ref_constructs_dir)
    logger.info(f"Найдено {len(files)} файлов в директории")
    
    for filename in files:
        if filename.endswith("-EF.csv"):
            file_count += 1
            filepath = os.path.join(ref_constructs_dir, filename)
            try:
                df = pd.read_csv(filepath)
                construct_count = len(df)
                ref_constructs.extend(df['ConstructID'].tolist())
                logger.debug(f"Файл {filename}: загружено {construct_count} конструктов")
            except Exception as e:
                logger.error(f"Ошибка чтения файла {filename}: {str(e)}")
    
    logger.info(f"Проверено {file_count} файлов с конструктами")
    logger.info(f"Всего загружено {len(ref_constructs)} конструктов референса")
    
    # 2. Собираем все позиции, которые покрыты конструктами референса
    covered_positions = set()
    processed_constructs = 0
    
    if not ref_constructs:
        logger.warning("Не найдено ни одного конструкта референса!")
        return covered_positions
    
    for construct_id in ref_constructs:
        arm_size, center, arm3_start, arm4_start = parse_construct_id(construct_id)
        if None in (arm_size, center, arm3_start, arm4_start):
            continue
        
        try:
            arm_ranges = calculate_arm_ranges(arm_size, center, arm3_start, arm4_start)
            for start, end in arm_ranges:
                # Добавляем +1 чтобы включить конечную позицию
                covered_positions.update(range(start, end + 1))
            processed_constructs += 1
        except Exception as e:
            logger.error(f"Ошибка обработки конструкта {construct_id}: {str(e)}")
    
    logger.info(f"Обработано {processed_constructs} конструктов")
    
    if covered_positions:
        logger.info(f"Найдено {len(covered_positions)} позиций, покрытых конструктами")
        logger.info(f"Диапазон покрытия: от {min(covered_positions)} до {max(covered_positions)}")
    else:
        logger.warning("Не найдено ни одной покрытой позиции!")
    
    return covered_positions

def apply_snvs(ref_record, snv_df, log_path: Path, covered_positions: set, num: int) -> SeqRecord:
    """Применяет 2 случайные SNV к референсной последовательности, гарантируя их присутствие в конструктах референса"""
    original_seq = str(ref_record.seq).upper()
    mutable_seq = MutableSeq(original_seq)
    
    # 1. Фильтруем SNV: только те, что покрыты конструктами
    filtered_snvs = []
    for _, row in snv_df.iterrows():
        position = int(row['position'])
        if position in covered_positions:
            filtered_snvs.append({
                'position': position,
                'ref_allele': row['ref_allele'].upper(),
                'alt_allele': row['alt_allele'].upper()
            })
    
    logger.info(f"Найдено {len(filtered_snvs)} SNV, покрытых конструктами")
    
    # 2. Случайный выбор 2 уникальных позиций из покрытых конструктами
    unique_positions = list({snv['position'] for snv in filtered_snvs})
    if len(unique_positions) >= 2:
        selected_positions = random.sample(unique_positions, 2)
    else:
        selected_positions = unique_positions
    
    logger.info(f"Выбрано позиций для мутации: {selected_positions}")
    
    # 3. Применение выбранных SNV
    applied_count = 0
    mismatch_log = []
    selected_snvs = []

    for position in selected_positions:
        position_snvs = [snv for snv in filtered_snvs if snv['position'] == position]
        idx = position - 1
        
        if idx >= len(original_seq):
            mismatch_log.append({
                'position': position,
                'original_base': None,
                'ref_allele': '|'.join([s['ref_allele'] for s in position_snvs]),
                'alt_allele': '|'.join([s['alt_allele'] for s in position_snvs]),
                'status': 'SKIPPED',
                'notes': 'Position out of sequence bounds'
            })
            continue
        
        current_base = original_seq[idx]
        applied = False
        
        for snv in position_snvs:
            ref_allele = snv['ref_allele']
            alt_allele = snv['alt_allele']
            
            if ref_allele == alt_allele:
                continue
            
            if current_base == ref_allele:
                mutable_seq[idx] = alt_allele
                applied_count += 1
                applied = True
                selected_snvs.append(snv)
                mismatch_log.append({
                    'position': position,
                    'original_base': current_base,
                    'ref_allele': ref_allele,
                    'alt_allele': alt_allele,
                    'status': 'APPLIED',
                    'notes': ''
                })
                break
        
        if not applied:
            notes = []
            for snv in position_snvs:
                if current_base == snv['alt_allele']:
                    notes.append('ALT allele already present')
                else:
                    notes.append(f"Expected ref: {snv['ref_allele']}, found: {current_base}")
            
            mismatch_log.append({
                'position': position,
                'original_base': current_base,
                'ref_allele': '|'.join([s['ref_allele'] for s in position_snvs]),
                'alt_allele': '|'.join([s['alt_allele'] for s in position_snvs]),
                'status': 'SKIPPED',
                'notes': '; '.join(notes)
            })
    
    # 4. Сохранение лога
    if mismatch_log:
        mismatch_df = pd.DataFrame(mismatch_log)
        mismatch_df.to_csv(log_path, index=False)
        logger.info(f"Лог мутаций сохранён в {log_path}")
    
    # Статистика
    logger.info(f"\nСтатистика применения SNV для последовательности #{num}:")
    logger.info(f"Всего SNV покрытых конструктами: {len(filtered_snvs)}")
    logger.info(f"Уникальных позиций: {len(unique_positions)}")
    logger.info(f"Выбрано позиций: {len(selected_positions)}")
    logger.info(f"Успешно применено: {applied_count}")
    
    return SeqRecord(
        seq=mutable_seq,
        id=f"custom_mtDNA_{num}",
        description=f"Modified from {ref_record.id} | Applied {applied_count} of {len(selected_positions)} selected SNVs"
    )

def main(num: int):
    # Конфигурация
    SNV_CSV_PATH = Path("D:/pythonProject/MitoFragility/DataPreparing/snv_csv/snvs.csv")
    INPUT_FASTA = Path("D:/pythonProject/MitoFragility/DataPreparing/sequences/ref_seq/Homo_sapiens_assembly38.chrM.fasta")
    XLSX_PATH = Path("D:/pythonProject/MitoFragility/DataPreparing/raw_data/MitoPhewas_associations.xlsx")
    LOG_PATH = Path(f"D:/pythonProject/MitoFragility/DataPreparing/snv_log/snv_log_{num}.csv")
    OUTPUT_FASTA = Path(f"D:/pythonProject/MitoFragility/DataPreparing/sequences/relative_seq/test_individual_{num+4}.fasta")
    REF_CONSTRUCTS_DIR = "D:/pythonProject/MitoFragility/MitoFragilityScore/Energies/SEQ-g38_Mt-Short_Test"
    
    # 1. Создаем CSV с SNV (только при первом запуске)
    if num == 0 and not SNV_CSV_PATH.exists():
        snv_df = csv_constructor(XLSX_PATH, SNV_CSV_PATH)
    else:
        snv_df = pd.read_csv(SNV_CSV_PATH)
    
    # 2. Загружаем референсную последовательность
    ref_record = SeqIO.read(INPUT_FASTA, "fasta")
    logger.info(f"Загружена референсная последовательность: {ref_record.id}")
    logger.info(f"Длина: {len(ref_record.seq)} bp")
    
    # 3. Получаем покрытые позиции
    covered_positions = get_covered_positions(REF_CONSTRUCTS_DIR)
    
    # 4. Применяем SNV
    custom_record = apply_snvs(ref_record, snv_df, LOG_PATH, covered_positions, num)
    
    # 5. Сохраняем результат
    SeqIO.write(custom_record, OUTPUT_FASTA, "fasta")
    logger.info(f"Результат сохранен в {OUTPUT_FASTA}")
    logger.info(f"ID: {custom_record.id}")
    logger.info(f"Описание: {custom_record.description}")

if __name__ == "__main__":
    # Создаем 5 вариантов последовательностей
    for i in range(5):
        logger.info(f"\n{'='*50}")
        logger.info(f"Создание последовательности #{i}")
        logger.info(f"{'='*50}")
        main(i)