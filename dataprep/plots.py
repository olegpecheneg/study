import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import probplot, genpareto
from .stats_utils import calculate_full_distribution_stats, fit_gpd_to_tail, calculate_survival_and_cdf, calculate_empirical_cvalue, calculate_anderson, calculate_z_score_p_value, apply_fdr_fast, log_transform

def create_interactive_analysis_plot(energy_type: str, data_type: str, 
                                    energies: np.ndarray, individual_id: str,
                                    paired_data=None, wilcoxon_results=None):
    """Создает полный интерактивный анализ распределения."""
    
    # Рассчитываем статистики
    distribution_stats = calculate_full_distribution_stats(energies)
    percentiles = distribution_stats.get('percentiles', {})
    gpd_result = fit_gpd_to_tail(energies)
    survival_data = calculate_survival_and_cdf(energies)
    
    # Создаем фигуру с subplots
    fig = make_subplots(
        rows=2, cols=3,
        subplot_titles=(
            '1. Полное распределение и перцентили',
            '2. Функция выживания и CDF',
            '3. Хвост распределения и GPD',
            '4. QQ-plot',
            '5. Манхэттен-плот (c-value)',
            '6. Сравнение с GPD'
        ),
        vertical_spacing=0.12,
        horizontal_spacing=0.1
    )
    
    # 1. Полное распределение и перцентили
    fig.add_trace(
        go.Histogram(
            x=energies, nbinsx=30, name='Распределение',
            marker_color='lightblue', opacity=0.7,
            hovertemplate='Энергия: %{x:.2f}<br>Частота: %{y}<extra></extra>'
        ),
        row=1, col=1
    )
    
    # Добавляем перцентили с повернутыми подписями
    colors = {'95': 'orange', '99': 'red', '99.9': 'darkred'}
    y_positions = [0.9, 0.7, 0.5]
    
    for idx, (p_key, p_val) in enumerate(percentiles.items()):
        if idx < len(y_positions):
            fig.add_vline(
                x=p_val, line_dash="dash", line_color=colors.get(p_key, 'black'),
                line_width=2, row=1, col=1
            )
            
            fig.add_annotation(
                x=p_val, y=y_positions[idx] * max(np.histogram(energies, bins=30)[0]),
                text=f'{p_key}%: {p_val:.2f}',
                showarrow=False,
                textangle=-90,
                font=dict(size=10, color=colors.get(p_key, 'black')),
                bgcolor='rgba(255,255,255,0.8)',
                bordercolor=colors.get(p_key, 'black'),
                borderwidth=1,
                row=1, col=1
            )
    
    # 2. Функция выживания и CDF на одном графике
    if survival_data:
        fig.add_trace(
            go.Scatter(
                x=survival_data['energies'], y=survival_data['survival'],
                mode='lines', name='S(x) - Функция выживания',
                line=dict(color='green', width=2),
                hovertemplate='Энергия: %{x:.2f}<br>S(x): %{y:.4f}<extra></extra>'
            ),
            row=1, col=2
        )
        
        fig.add_trace(
            go.Scatter(
                x=survival_data['energies'], y=survival_data['cdf'],
                mode='lines', name='F(x) - CDF',
                line=dict(color='blue', width=2, dash='dash'),
                hovertemplate='Энергия: %{x:.2f}<br>F(x): %{y:.4f}<extra></extra>'
            ),
            row=1, col=2
        )
        
        fig.update_yaxes(range=[0, 1.05], row=1, col=2)
    
    # 3. Хвост распределения и GPD
    if 'threshold' in gpd_result:
        threshold = gpd_result['threshold']
        tail_data = energies[energies > threshold]
        
        if len(tail_data) > 0:
            fig.add_trace(
                go.Histogram(
                    x=tail_data, nbinsx=20, name='Хвост распределения',
                    marker_color='orange', opacity=0.7,
                    hovertemplate='Энергия: %{x:.2f}<br>Частота: %{y}<extra></extra>'
                ),
                row=1, col=3
            )
            
            if 'shape' in gpd_result and 'scale' in gpd_result:
                shape = gpd_result['shape']
                scale = gpd_result['scale']
                
                x_min, x_max = threshold, max(tail_data)
                if x_max > x_min:
                    x_gpd = np.linspace(x_min, x_max, 100)
                    y_gpd = genpareto.pdf(x_gpd - threshold, shape, scale=scale)
                    
                    bin_width = (x_max - x_min) / 20
                    y_gpd_scaled = y_gpd * len(tail_data) * bin_width
                    
                    fig.add_trace(
                        go.Scatter(
                            x=x_gpd, y=y_gpd_scaled, mode='lines',
                            name='GPD плотность', line=dict(color='red', width=2),
                            hovertemplate='Энергия: %{x:.2f}<br>GPD плотность: %{y:.4f}<extra></extra>'
                        ),
                        row=1, col=3
                    )
    
    # 4. QQ-plot
    try:
        (os, qs), (slope, intercept, r) = probplot(energies, dist="norm", fit=True)
        
        fig.add_trace(
            go.Scatter(
                x=os, y=qs, mode='markers', name='QQ-plot',
                marker=dict(color='blue', size=6, opacity=0.7),
                hovertemplate='Теор. квантиль: %{x:.3f}<br>Эмп. квантиль: %{y:.3f}<extra></extra>'
            ),
            row=2, col=1
        )
        
        x_line = np.array([os.min(), os.max()])
        y_line = slope * x_line + intercept
        
        fig.add_trace(
            go.Scatter(
                x=x_line, y=y_line, mode='lines',
                line=dict(color='red', width=2, dash='dash'),
                name='Теоретическая нормаль',
                hovertemplate='Теор. квантиль: %{x:.3f}<br>Ожидаемый квантиль: %{y:.3f}<extra></extra>'
            ),
            row=2, col=1
        )

        anderson_p, normality_text = calculate_anderson(energies)
        
        
        fig.add_annotation(
            x=0.05, y=0.95, xref="x domain", yref="y domain",
            text=f"R² = {r**2:.3f}<br>{normality_text}",
            showarrow=False,
            bgcolor="white",
            bordercolor="red" if anderson_p < 0.05 else "green",
            borderwidth=1,
            font=dict(size=10),
            row=2, col=1
        )
        
    except Exception as e:
        fig.add_annotation(
            x=0.5, y=0.5, xref="x domain", yref="y domain",
            text=f"Ошибка QQ-plot: {e}",
            showarrow=False, row=2, col=1
        )
    
    # 5. Манхэттен-плот
    if data_type == 'diff' and wilcoxon_results and 'error' not in wilcoxon_results:
        # Используем c-values из Wilcoxon теста
        n_pairs = wilcoxon_results['n_pairs']
        log_c_values = wilcoxon_results['log_c_values']
        effect_sign = wilcoxon_results['effect_sign']
        differences = wilcoxon_results['differences']
        
        _, z_score_p_values = calculate_z_score_p_value(energies)
        
        q_reject, q_values, q_threshold, q_log_treshold = apply_fdr_fast(z_score_p_values, alpha=0.05)
        z_score_p_values_log = log_transform(z_score_p_values)

        # Пороги статистической значимости
        normal_threshold = -np.log10(0.05)  # Обычный порог 0.05
        bonferroni_threshold = -np.log10(0.05 / len(energies)) if len(energies) > 0 else 0
        # Цвета по знаку эффекта
        colors = ['red' if sign > 0 else 'blue' for sign in effect_sign]

        # Создаем массив alpha (прозрачности) для каждой точки
        alpha_values = np.zeros_like(z_score_p_values_log, dtype=float)

        # Зоны прозрачности:
        # 1. Ниже normal_threshold: alpha = 0.2
        mask_lowest = z_score_p_values_log < normal_threshold
        alpha_values[mask_lowest] = 0.2

        # 2. Между normal_threshold и q_log_treshold: alpha = 0.5
        mask_mid_low = (z_score_p_values_log >= normal_threshold) & (z_score_p_values_log < q_log_treshold)
        alpha_values[mask_mid_low] = 0.5

        # 3. Между q_log_treshold и bonferroni_threshold: alpha = 0.7
        mask_mid_high = (z_score_p_values_log >= q_log_treshold) & (z_score_p_values_log < bonferroni_threshold)
        alpha_values[mask_mid_high] = 0.7

        # 4. Выше bonferroni_threshold: alpha = 1.0
        mask_highest = z_score_p_values_log >= bonferroni_threshold
        alpha_values[mask_highest] = 1.0


        # Создаем цвета с прозрачностью в формате rgba
        colors_rgba = []
        for sign, alpha in zip(effect_sign, alpha_values):
            if sign > 0:
                # Красный с заданной прозрачностью
                colors_rgba.append(f'rgba(255, 0, 0, {alpha})')
            else:
                # Синий с заданной прозрачностью
                colors_rgba.append(f'rgba(0, 0, 255, {alpha})')

        # Создаем цвет контура с той же прозрачностью
        line_colors = [f'rgba(0, 0, 0, {alpha})' for alpha in alpha_values]

        fig.add_trace(
            go.Scatter(
                x=list(range(n_pairs)),
                y=z_score_p_values_log,
                mode='markers',
                name='p-value (Wilcoxon)',
                marker=dict(
                    size=5,
                    color=colors_rgba,
                    showscale=False,
                    line=dict(width=0.1, color='black')
                ),
                hovertemplate=(
                    'Пара: %{x}<br>'
                    'Разность (ref-alt): %{customdata:.3f}<br>'
                    '-log₁₀(p-value): %{y:.3f}<br>'
                    'p-value: 10^(-%{y:.3f})<extra></extra>'
                ),
                customdata=differences
            ),
            row=2, col=2
        )

            # Линия статистической значимости (обычный порог 0.05)
        fig.add_hline(
                y=normal_threshold,
                line_dash="dash", line_color="red",
                annotation_text=f'Порог 0.05 = {normal_threshold:.3f}',
                annotation_position="top right",
                annotation_font_size=10,
                row=2, col=2
        )
        
            # Линия статистической значимости (поправка Бонферрони)
        fig.add_hline(
                y=bonferroni_threshold,
                line_dash="dash", line_color="red",
                annotation_text=f'Порог Бонферрони = {bonferroni_threshold:.3f}',
                annotation_position="top right",
                annotation_font_size=10,
                row=2, col=2
        )
            # Линия статистической значимости (поправка FDR)
        fig.add_hline(
                y=q_log_treshold,
                line_dash="dash", line_color="red",
                annotation_text=f'Порог FDR = {q_log_treshold:.3f}',
                annotation_position="top right",
                annotation_font_size=10,
                row=2, col=2
        )
        # Добавляем p-value от Wilcoxon теста
        wilcoxon_p = wilcoxon_results['wilcoxon_p_value']
        mean_diff = wilcoxon_results['mean_difference']
        median_diff = wilcoxon_results['median_difference']
        direction = "ref > alt" if mean_diff > 0 else "ref < alt"
        
        # Цвет текста в зависимости от значимости
        text_color = 'red' if wilcoxon_p < 0.05 else 'black'
        border_color = 'red' if wilcoxon_p < 0.05 else 'green'
        
        fig.add_annotation(
            x=0.05, y=0.95, xref="x domain", yref="y domain",
            text=(
                f"Wilcoxon p = {wilcoxon_p:.4f}<br>"
                f"Медианная разность: {median_diff:.3f}<br>"
                f"Направление: {direction}"
            ),
            showarrow=False,
            bgcolor="white",
            bordercolor=border_color,
            borderwidth=2,
            font=dict(size=10, color=text_color),
            row=2, col=2
        )
        
    elif data_type in ['ref', 'alt']:
        # Для ref или alt используем эмпирический c-value
        empirical_result = calculate_empirical_cvalue(energies)
        
        if 'error' not in empirical_result:
            log_c_values = empirical_result['log_c_values']
             # Пороги статистической значимости
            normal_threshold = -np.log10(0.05)  # Обычный порог 0.05
            bonferroni_threshold = -np.log10(0.05 / len(energies)) if len(energies) > 0 else 0

            
            fig.add_trace(
                go.Scatter(
                    x=list(range(len(log_c_values))),
                    y=log_c_values,
                    mode='markers',
                    name=f'c-value ({data_type})',
                    marker=dict(
                        size=5,
                        color=log_c_values,
                        colorscale='Viridis',
                        showscale=False,
                        cmin=0,
                        cmax=empirical_result.get('max_log_c', 5)
                    ),
                    hovertemplate='Индекс: %{x}<br>-log₁₀(c-value): %{y:.3f}<br>Энергия: %{text:.2f}<extra></extra>',
                    text=energies
                ),
                row=2, col=2
            )
            
                    # Линия статистической значимости (обычный порог 0.05)
            fig.add_hline(
                y=normal_threshold,
                line_dash="dash", line_color="red",
                annotation_text=f'Порог 0.05 = {normal_threshold:.3f}',
                annotation_position="top right",
                annotation_font_size=10,
                row=2, col=2
            )
        
            # Линия статистической значимости (поправка Бонферрони)
            fig.add_hline(
                y=bonferroni_threshold,
                line_dash="dash", line_color="red",
                annotation_text=f'Порог Бонферрони = {bonferroni_threshold:.3f}',
                annotation_position="top right",
                annotation_font_size=10,
                row=2, col=2
            )
            
            # Если есть результаты Wilcoxon, показываем их
            if wilcoxon_results and 'error' not in wilcoxon_results:
                wilcoxon_p = wilcoxon_results['wilcoxon_p_value']
                text_color = 'red' if wilcoxon_p < 0.05 else 'black'
                
                fig.add_annotation(
                    x=0.05, y=0.85, xref="x domain", yref="y domain",
                    text=f"Wilcoxon p = {wilcoxon_p:.3f}",
                    showarrow=False,
                    bgcolor="white",
                    bordercolor="black",
                    borderwidth=1,
                    font=dict(size=10, color=text_color),
                    row=2, col=2
                )
    
    # 6. Сравнение с GPD
    if survival_data and 'shape' in gpd_result and 'scale' in gpd_result:
        threshold = gpd_result['threshold']
        shape = gpd_result['shape']
        scale = gpd_result['scale']
        
        tail_mask = np.array(survival_data['energies']) > threshold
        tail_energies = np.array(survival_data['energies'])[tail_mask]
        
        if len(tail_energies) > 0:
            gpd_cdf_tail = genpareto.cdf(tail_energies - threshold, shape, scale=scale)
            fraction_below = np.sum(energies <= threshold) / len(energies)
            gpd_cdf_scaled = fraction_below + (1 - fraction_below) * gpd_cdf_tail
            
            fig.add_trace(
                go.Scatter(
                    x=tail_energies, y=gpd_cdf_scaled,
                    mode='lines', name='GPD CDF',
                    line=dict(color='red', width=2, dash='dash'),
                    hovertemplate='Энергия: %{x:.2f}<br>GPD CDF: %{y:.4f}<extra></extra>'
                ),
                row=2, col=3
            )
            
            empirical_cdf_tail = np.array(survival_data['cdf'])[tail_mask]
            
            fig.add_trace(
                go.Scatter(
                    x=tail_energies, y=empirical_cdf_tail,
                    mode='lines', name='Эмп. CDF',
                    line=dict(color='blue', width=2),
                    hovertemplate='Энергия: %{x:.2f}<br>Эмп. CDF: %{y:.4f}<extra></extra>'
                ),
                row=2, col=3
            )
    
    # Обновление layout
    fig.update_layout(
        title=dict(
            text=f"{individual_id} - {energy_type} ({data_type})<br>"
                 f"Анализ распределения (N={len(energies)})",
            x=0.5, font=dict(size=18)
        ),
        height=1000,
        showlegend=True,
        hovermode='closest',
        plot_bgcolor='white',
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=1.02,
            bgcolor='rgba(255, 255, 255, 0.8)'
        )
    )
    
    # Настройки осей
    fig.update_xaxes(title_text="Энергия (ккал/моль)", row=1, col=1)
    fig.update_yaxes(title_text="Частота", row=1, col=1)
    
    fig.update_xaxes(title_text="Энергия (ккал/моль)", row=1, col=2)
    fig.update_yaxes(title_text="Вероятность", row=1, col=2)
    
    fig.update_xaxes(title_text="Энергия (ккал/моль)", row=1, col=3)
    fig.update_yaxes(title_text="Частота", row=1, col=3)
    
    fig.update_xaxes(title_text="Теоретические квантили", row=2, col=1)
    fig.update_yaxes(title_text="Эмпирические квантили", row=2, col=1)
    
    fig.update_xaxes(title_text="Индекс пары/точки", row=2, col=2)
    fig.update_yaxes(title_text="-log₁₀(c-value)", row=2, col=2)
    
    fig.update_xaxes(title_text="Энергия (ккал/моль)", row=2, col=3)
    fig.update_yaxes(title_text="CDF", row=2, col=3)
    
    return fig

def create_statistics_table(energy_type: str, data_type: str,
                           distribution_stats: dict, 
                           gpd_result_right: dict,  # Правый хвост
                           gpd_result_left: dict,   # Левый хвост
                           survival_data: dict, 
                           wilcoxon_results=None):
    """Создает таблицу со статистикой, включая оба хвоста."""
    
    stats = {
        'Количество точек': distribution_stats.get('count', 0),
        'Среднее': f"{distribution_stats.get('mean', 0):.3f}",
        'Стандартное отклонение': f"{distribution_stats.get('std', 0):.3f}",
        'Минимум': f"{distribution_stats.get('min', 0):.3f}",
        'Максимум': f"{distribution_stats.get('max', 0):.3f}",
        'Медиана': f"{distribution_stats.get('median', 0):.3f}",
        'Асимметрия': f"{distribution_stats.get('skewness', 0):.3f}",
        'Эксцесс': f"{distribution_stats.get('kurtosis', 0):.3f}"
    }
    
    # Перцентили
    percentiles = distribution_stats.get('percentiles', {})
    for p_key in ['95', '99', '99.9']:
        if p_key in percentiles:
            stats[f'{p_key}%-й перцентиль'] = f"{percentiles[p_key]:.3f}"
    
    # GPD статистика для ПРАВОГО хвоста
    if 'shape' in gpd_result_right:
        stats['Правый хвост - порог (95%-й)'] = f"{gpd_result_right['threshold']:.3f}"
        stats['Правый хвост - форма (ξ)'] = f"{gpd_result_right['shape']:.4f}"
        stats['Правый хвост - масштаб (σ)'] = f"{gpd_result_right['scale']:.4f}"
        stats['Правый хвост - точек'] = gpd_result_right['exceedances_count']
        stats['Правый хвост - тип'] = gpd_result_right.get('tail_behavior', 'N/A')
    
    # GPD статистика для ЛЕВОГО хвоста
    if 'shape' in gpd_result_left:
        stats['Левый хвост - порог (5%-й)'] = f"{gpd_result_left['threshold']:.3f}"
        stats['Левый хвост - форма (ξ)'] = f"{gpd_result_left['shape']:.4f}"
        stats['Левый хвост - масштаб (σ)'] = f"{gpd_result_left['scale']:.4f}"
        stats['Левый хвост - точек'] = gpd_result_left['exceedances_count']
        stats['Левый хвост - тип'] = gpd_result_left.get('tail_behavior', 'N/A')
    elif 'error' in gpd_result_left:
        stats['Левый хвост - статус'] = gpd_result_left['error']
    
    # Wilcoxon статистика
    if wilcoxon_results and 'error' not in wilcoxon_results:
        stats['Wilcoxon p-value'] = f"{wilcoxon_results['wilcoxon_p_value']:.4f}"
        stats['Wilcoxon статистика'] = f"{wilcoxon_results['wilcoxon_statistic']:.2f}"
        stats['Количество пар'] = wilcoxon_results['n_pairs']
        stats['Средняя разность (ref-alt)'] = f"{wilcoxon_results['mean_difference']:.3f}"
        stats['Медианная разность'] = f"{wilcoxon_results['median_difference']:.3f}"
        stats['Станд. отклонение разности'] = f"{wilcoxon_results['std_difference']:.3f}"
        stats['Среднее ref'] = f"{wilcoxon_results['ref_mean']:.3f}"
        stats['Среднее alt'] = f"{wilcoxon_results['alt_mean']:.3f}"
        stats['Направление эффекта'] = wilcoxon_results.get('ref_vs_alt', 'N/A')
        stats['Статистическая значимость'] = 'ДА' if wilcoxon_results.get('is_significant', False) else 'НЕТ'
    
    # Создание таблицы
    fig = go.Figure(data=[go.Table(
        header=dict(
            values=['<b>Параметр</b>', '<b>Значение</b>'],
            fill_color='lightgrey', align='center',
            font=dict(size=14, color='black')
        ),
        cells=dict(
            values=[list(stats.keys()), list(stats.values())],
            fill_color='white', align=['left', 'right'],
            font=dict(size=12)
        )
    )])
    
    fig.update_layout(
        title=dict(
            text=f"Статистика {energy_type} ({data_type})",
            x=0.5, font=dict(size=16)
        ),
        height=400 + len(stats) * 25,
        margin=dict(l=20, r=20, t=60, b=20)
    )
    
    return fig