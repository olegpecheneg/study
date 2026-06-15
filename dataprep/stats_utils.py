import numpy as np
import pandas as pd
from scipy.stats import genpareto, wilcoxon, rankdata, anderson, norm
from statsmodels.stats.multitest import multipletests
import logging

logger = logging.getLogger("energy_distribution")




def apply_fdr_fast(p_values, alpha=0.05):
    """
    Быстрая FDR коррекция Бенджамини-Хохберга
    
    Parameters:
    -----------
    p_values : array-like
        Массив p-values
    alpha : float
        Уровень значимости FDR (по умолчанию 0.05)
    
    Returns:
    --------
    reject : np.ndarray[bool]
        Булев массив: True для значимых гипотез
    q_values : np.ndarray
        Скорректированные p-values (q-values)
    threshold : float
        Порог для исходных p-values
    """
    # FDR коррекция
    reject, q_values, _, _ = multipletests(p_values, alpha=alpha, method='fdr_bh')
    
    # Находим порог для исходных p-values
    if np.any(reject):
        # Максимальное исходное p-value среди значимых
        threshold = np.max(p_values[reject])
    else:
        threshold = 0
    log_treshold = -np.log10(threshold) if threshold > 0 else 0.0
    return reject, q_values, threshold, log_treshold

def log_transform(values: np.ndarray, base: float = 10.0) -> np.ndarray:
    """Выполняет логарифмическое преобразование значений с заданным основанием."""
    with np.errstate(divide='ignore'):
        log_values = -np.log(values) / np.log(base)
    log_values[np.isinf(log_values)] = 0.0  # Заменяем -inf на 0
    return log_values

def calculate_z_score_p_value(diffs: np.ndarray):
    """
    Устойчивый Z-score на основе медианы и MAD
    Не требует предположения о нормальности
    
    Parameters:
    -----------
    energies_diff : np.ndarray
        Массив РАЗНОСТЕЙ энергий (alt - ref) для каждого конструкта
    
    Returns:
    --------
    p_values : np.ndarray
        p-values в том же порядке, что и входные разности
    """
    # Медиана вместо среднего
    median = np.median(diffs)
    
    # Median Absolute Deviation (MAD) вместо стандартного отклонения
    mad = np.median(np.abs(diffs - median))
    
    # Преобразование MAD в оценку стандартного отклонения
    # Для нормального распределения: σ = MAD / 0.6745
    # Это делает Z-scores сопоставимыми с классическими
    if mad == 0:
        mad = 1e-10  # защита от деления на 0
    
    robust_std = mad / 0.6745
    
    # Устойчивый Z-score
    z_robust = (diffs - median) / robust_std
    
    # P-values (все еще из нормального распределения, но данные устойчивее)
    p_values = 2 * (1 - norm.cdf(np.abs(z_robust)))
    
    return z_robust, p_values

def calculate_anderson(energies: np.ndarray):
    # Добавляем тест нормальности anderson-darling
    normality_text = ""
    if len(energies) <= 5000:  # anderson работает до 5000 точек
        result = anderson(energies)
        anderson_statistic = result.statistic
        anderson_critical_values = result.critical_values
        anderson_significance_level = result.significance_level
        
        if anderson_statistic > anderson_critical_values[2]:  # 5% уровень значимости
            normality_text = f"Данные НЕ нормальны (stat={anderson_statistic:.3f} > crit={anderson_critical_values[2]:.3f} при α=5%)"
        else:
            normality_text = f"Данные нормальны (stat={anderson_statistic:.3f} ≤ crit={anderson_critical_values[2]:.3f} при α=5%)"
        anderson_p = anderson_p_value(anderson_statistic, anderson_critical_values, anderson_significance_level)
    return anderson_p, normality_text

def anderson_p_value(statistic, critical_values, significance_levels, ):
    """
    Вычисляет приближенное p-value для теста Андерсона-Дарлинга
    Возвращает: (statistic, p_value_approx, message)
    """
    from scipy.interpolate import interp1d
    
    significance_levels = significance_levels / 100  # преобразуем в доли
    
    # Если статистика меньше минимального критического значения
    if statistic < critical_values[0]:
        p_value = significance_levels[0]  # p > 0.15
        return p_value
    
    # Если статистика больше максимального критического значения
    if statistic > critical_values[-1]:
        p_value = significance_levels[-1]  # p < 0.01
        return p_value
    
    # Интерполяция для нахождения приближенного p-value
    # Переворачиваем массивы для интерполяции (x=critical_values, y=significance_levels)
    f = interp1d(critical_values, significance_levels, kind='linear', fill_value='extrapolate')
    p_value = float(f(statistic))
    
    return p_value

def calculate_percentiles(energies: np.ndarray):
    """Вычисляет 95-й, 99-й, 99.9-й перцентили."""
    if len(energies) < 10:
        return {}
    
    percentiles_dict = {}
    for p in [95, 99, 99.9]:
        value = np.percentile(energies, p)
        percentiles_dict[str(p)] = float(value)
    
    # Проверяем монотонность
    if '95' in percentiles_dict and '99' in percentiles_dict and '99.9' in percentiles_dict:
        if percentiles_dict['99.9'] < percentiles_dict['99']:
            logger.warning(f"99.9 перцентиль ({percentiles_dict['99.9']:.3f}) меньше 99 перцентиля ({percentiles_dict['99']:.3f})!")
    
    return percentiles_dict

def calculate_full_distribution_stats(energies: np.ndarray):
    """Рассчитывает полную статистику распределения."""
    if len(energies) < 10:
        return {}
    
    stats = {
        'count': len(energies),
        'mean': float(np.mean(energies)),
        'std': float(np.std(energies)),
        'min': float(np.min(energies)),
        'max': float(np.max(energies)),
        'median': float(np.median(energies)),
        'skewness': float(pd.Series(energies).skew()),
        'kurtosis': float(pd.Series(energies).kurtosis()),
        'percentiles': calculate_percentiles(energies)
    }
    
    return stats

def fit_gpd_to_left_tail(energies: np.ndarray, threshold_percentile: float = 5):
    """Подгоняет GPD к левому хвосту распределения (низкие значения)."""
    if len(energies) < 50:
        return {'error': 'Недостаточно данных для GPD левого хвоста'}
    
    try:
        # Для левого хвоста берем нижний перцентиль
        threshold = np.percentile(energies, threshold_percentile)
        
        # Превышения порога: порог - значение (делаем положительными)
        exceedances = threshold - energies[energies < threshold]
        
        if len(exceedances) < 10:
            return {
                'error': f'Недостаточно значений в левом хвосте: {len(exceedances)} < 10',
                'threshold': float(threshold),
                'exceedances_count': len(exceedances)
            }
        
        params = genpareto.fit(exceedances, floc=0)
        shape, loc, scale = params
        
        gpd_result = {
            'threshold': float(threshold),
            'threshold_percentile': threshold_percentile,
            'exceedances_count': len(exceedances),
            'shape': float(shape),
            'scale': float(scale),
            'mean_exceedance': float(np.mean(exceedances)),
            'tail_direction': 'left'
        }
        
        # Определяем тип хвоста
        if shape > 0:
            gpd_result['distribution_type'] = 'heavy_tailed'
            gpd_result['tail_behavior'] = 'Тяжелый хвост'
        elif shape == 0:
            gpd_result['distribution_type'] = 'exponential'
            gpd_result['tail_behavior'] = 'Экспоненциальный хвост'
        else:
            gpd_result['distribution_type'] = 'short_tailed'
            gpd_result['tail_behavior'] = 'Короткий хвост'
        
        return gpd_result
        
    except Exception as e:
        logger.error(f"Ошибка подгонки GPD к левому хвосту: {e}")
        return {'error': str(e)}

def fit_gpd_to_tail(energies: np.ndarray, threshold_percentile: float = 95):
    """Подгоняет GPD к хвосту распределения."""
    if len(energies) < 50:
        return {'error': 'Недостаточно данных для GPD'}
    
    try:
        threshold = np.percentile(energies, threshold_percentile)
        exceedances = energies[energies > threshold] - threshold
        
        if len(exceedances) < 10:
            return {
                'error': f'Недостаточно превышений порога: {len(exceedances)} < 10',
                'threshold': float(threshold),
                'exceedances_count': len(exceedances)
            }
        
        params = genpareto.fit(exceedances, floc=0)
        shape, loc, scale = params
        
        gpd_result = {
            'threshold': float(threshold),
            'threshold_percentile': threshold_percentile,
            'exceedances_count': len(exceedances),
            'shape': float(shape),
            'scale': float(scale),
            'mean_exceedance': float(np.mean(exceedances)),
            'exceedances': exceedances.tolist()
        }
        
        if shape > 0:
            gpd_result['distribution_type'] = 'heavy_tailed'
            gpd_result['tail_behavior'] = 'Тяжелый хвост'
        elif shape == 0:
            gpd_result['distribution_type'] = 'exponential'
            gpd_result['tail_behavior'] = 'Экспоненциальный хвост'
        else:
            gpd_result['distribution_type'] = 'short_tailed'
            gpd_result['tail_behavior'] = 'Короткий хвост'
        
        return gpd_result
        
    except Exception as e:
        logger.error(f"Ошибка подгонки GPD: {e}")
        return {'error': str(e)}

def calculate_survival_and_cdf(energies: np.ndarray):
    """Вычисляет функцию выживания и CDF."""
    if len(energies) < 10:
        return {}
    
    sorted_energies = np.sort(energies)
    n = len(sorted_energies)
    
    # Эмпирическая CDF
    cdf = np.arange(1, n + 1) / n
    
    # Функция выживания
    survival = 1 - cdf
    
    return {
        'energies': sorted_energies.tolist(),
        'survival': survival.tolist(),
        'cdf': cdf.tolist()
    }

def calculate_wilcoxon_test(ref_energies: np.ndarray, alt_energies: np.ndarray):
    """
    Рассчитывает Wilcoxon signed-rank test для парных данных.
    
    Внимание! Это тест для СРАВНЕНИЯ ref vs alt, а не для каждой выборки отдельно!
    Возвращает ОДИН p-value для сравнения двух выборок.
    """
    if len(ref_energies) != len(alt_energies):
        return {'error': 'Количество ref и alt точек не совпадает'}
    
    if len(ref_energies) < 10:
        return {'error': 'Недостаточно парных данных (минимум 10 пар)'}
    
    try:
        # 1. Wilcoxon signed-rank test (основной тест)
        stat, p_value = wilcoxon(ref_energies, alt_energies)
        
        # 2. Разности между парами
        differences = ref_energies - alt_energies
        
        # 3. Абсолютные разности для расчета рангов
        abs_differences = np.abs(differences)
        
        # 4. Ранжируем абсолютные разности
        # Используем метод 'average' для одинаковых значений
        ranks = rankdata(abs_differences, method='average')
        
        # 5. C-value на основе рангов
        # Концепция: большая абсолютная разность = большой ранг = маленький c-value
        n = len(ranks)
        max_rank = n
        
        # Преобразуем ранги в c-values:
        # c-value = 1 - (ранг/макс_ранг)
        # Это означает: ранг 1 (самая маленькая разность) → c-value ≈ 1
        # ранг n (самая большая разность) → c-value ≈ 0
        c_values = 1 - (ranks / max_rank)
        
        # Избегаем нулевых c-values (для логарифма)
        c_values = np.where(c_values < 1e-10, 1e-10, c_values)
        
        # 6. Направление эффекта
        effect_sign = np.sign(differences)
        
        # 7. Основные статистики
        ref_mean = float(np.mean(ref_energies))
        alt_mean = float(np.mean(alt_energies))
        mean_diff = float(np.mean(differences))
        
        return {
            # Результаты теста
            'wilcoxon_statistic': float(stat),
            'wilcoxon_p_value': float(p_value),
            'n_pairs': n,
            
            # Данные для визуализации
            'differences': differences.tolist(),
            'abs_differences': abs_differences.tolist(),
            'ranks': ranks.tolist(),
            'c_values': c_values.tolist(),
            'log_c_values': (-np.log10(c_values)).tolist(),
            'effect_sign': effect_sign.tolist(),
            
            # Статистики
            'mean_difference': mean_diff,
            'median_difference': float(np.median(differences)),
            'std_difference': float(np.std(differences)),
            'ref_mean': ref_mean,
            'alt_mean': alt_mean,
            'ref_median': float(np.median(ref_energies)),
            'alt_median': float(np.median(alt_energies)),
            
            # Интерпретация
            'ref_vs_alt': 'ref > alt в среднем' if mean_diff > 0 else 'ref < alt в среднем',
            'is_significant': bool(p_value < 0.05)  # Явно преобразуем в bool
        }
        
    except Exception as e:
        logger.error(f"Ошибка расчета Wilcoxon: {e}")
        return {'error': str(e)}

def calculate_empirical_cvalue(energies: np.ndarray):
    """
    Рассчитывает эмпирический c-value для одиночной выборки.
    Используется для ref или alt, когда нет парных данных.
    
    c-value = 1 - CDF(значение)
    Т.е. доля точек, которые БОЛЬШЕ или РАВНЫ данной.
    """
    if len(energies) < 10:
        return {'error': 'Недостаточно данных'}
    
    try:
        # Сортируем для расчета CDF
        sorted_energies = np.sort(energies)
        n = len(sorted_energies)
        
        # Для каждой исходной точки находим ее c-value
        c_values = np.zeros_like(energies, dtype=float)
        
        for i, energy in enumerate(energies):
            # Количество точек, которые МЕНЬШЕ или РАВНЫ данной
            count_less_equal = np.sum(sorted_energies <= energy)
            
            # Эмпирическая CDF = P(X ≤ energy)
            empirical_cdf = count_less_equal / n
            
            # c-value = 1 - CDF = P(X > energy)
            c_value = 1 - empirical_cdf
            if c_value >= 0.5:
                c_value = empirical_cdf
            
            # Избегаем нулевых c-values
            c_values[i] = max(c_value, 1e-10)
        # P-values (все еще из нормального распределения, но данные устойчивее)
        log_c_values = -np.log10(c_values)

        return {
            'c_values': c_values.tolist(),
            'log_c_values': log_c_values.tolist(),
            'max_log_c': float(np.max(log_c_values)),
        }
        
    except Exception as e:
        logger.error(f"Ошибка расчета эмпирического c-value: {e}")
        return {'error': str(e)}