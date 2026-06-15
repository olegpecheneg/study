import pandas as pd
import plotly.io as pio
from pathlib import Path
import json

def generate_html_report(mutation_id, energy_type, data_type,
                        plot_fig, table_fig, energies,
                        distribution_stats, gpd_result_right, gpd_result_left,
                        survival_data, wilcoxon_results, output_dir):  # параметр output_dir
    """Генерирует HTML отчет."""
    
    # Исправляем: используем output_dir вместо output_path
    html_path = output_dir / f"{mutation_id}_{energy_type}_{data_type}_analysis.html"
    
    plot_html = pio.to_html(plot_fig, full_html=False, include_plotlyjs='cdn')
    table_html = pio.to_html(table_fig, full_html=False, include_plotlyjs=False)
    
    # Добавляем информацию о Wilcoxon test в отчет
    wilcoxon_info = ""
    if wilcoxon_results and 'error' not in wilcoxon_results:
        wilcoxon_p = wilcoxon_results['wilcoxon_p_value']
        significance = "значимо" if wilcoxon_p < 0.05 else "не значимо"
        mean_diff = wilcoxon_results['mean_difference']
        direction = "ref > alt" if mean_diff > 0 else "ref < alt"
        
        wilcoxon_info = f"""
        <div class="wilcoxon-info">
            <h4>📊 Wilcoxon Signed-Rank Test (ref vs alt)</h4>
            <p><strong>p-value:</strong> {wilcoxon_p:.4f} ({significance} при α=0.05)</p>
            <p><strong>Статистика:</strong> {wilcoxon_results['wilcoxon_statistic']:.2f}</p>
            <p><strong>Количество пар:</strong> {wilcoxon_results['n_pairs']}</p>
            <p><strong>Средняя разность (ref-alt):</strong> {mean_diff:.3f} ({direction})</p>
            <p><strong>Медианная разность:</strong> {wilcoxon_results['median_difference']:.3f}</p>
            <p><strong>Станд. отклонение разности:</strong> {wilcoxon_results['std_difference']:.3f}</p>
        </div>
        """
    
    # Добавляем информацию о хвостах
    tails_info = ""
    if 'shape' in gpd_result_right or 'shape' in gpd_result_left:
        tails_info = """
        <div class="tails-info" style="background: #f0f7ff; padding: 15px; border-radius: 5px; margin: 20px 0; border-left: 4px solid #4CAF50;">
            <h4>📈 Анализ хвостов распределения</h4>
        """
        
        if 'shape' in gpd_result_right:
            tails_info += f"""
            <p><strong>Правый хвост (значения > {gpd_result_right.get('threshold', 0):.2f}):</strong> 
            {gpd_result_right.get('tail_behavior', 'N/A')} (ξ={gpd_result_right.get('shape', 0):.4f}, 
            {gpd_result_right.get('exceedances_count', 0)} точек)</p>
            """
        
        if 'shape' in gpd_result_left:
            tails_info += f"""
            <p><strong>Левый хвост (значения < {gpd_result_left.get('threshold', 0):.2f}):</strong> 
            {gpd_result_left.get('tail_behavior', 'N/A')} (ξ={gpd_result_left.get('shape', 0):.4f}, 
            {gpd_result_left.get('exceedances_count', 0)} точек)</p>
            """
        
        tails_info += "</div>"
    
    full_html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Анализ мутации {mutation_id} - {energy_type} ({data_type})</title>
        <meta charset="utf-8">
        <style>
            body {{ font-family: Arial; margin: 20px; background: #f5f5f5; }}
            .container {{ max-width: 1400px; margin: auto; background: white; 
                        padding: 20px; border-radius: 10px; box-shadow: 0 0 10px rgba(0,0,0,0.1); }}
            .header {{ text-align: center; margin-bottom: 20px; padding-bottom: 10px; 
                     border-bottom: 2px solid #4CAF50; }}
            .plot-container {{ margin: 20px 0; }}
            .stats-container {{ margin: 20px 0; }}
            .wilcoxon-info {{ 
                background: #e8f4f8; 
                padding: 15px; 
                border-radius: 5px; 
                margin: 20px 0;
                border-left: 4px solid #2196F3;
            }}
            .note {{ 
                background: #fff3cd; 
                padding: 10px; 
                border-radius: 5px; 
                margin: 10px 0;
                border-left: 4px solid #ffc107;
                font-size: 14px;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="header">
                <h1>🧬 Анализ распределения энергий для мутации</h1>
                <h2>Мутация: {mutation_id} | Тип энергии: {energy_type} | Данные: {data_type}</h2>
                <p>Всего точек: {len(energies)} | Дата анализа: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}</p>
            </div>
            
            <div class="note">
                <p><strong>Примечание:</strong> Wilcoxon signed-rank test сравнивает парные данные ref и alt. 
                Один и тот же p-value отображается для всех трех типов данных (ref, alt, diff), 
                так как тест выполняется для сравнения ref vs alt.</p>
            </div>
            
            {wilcoxon_info}
            
            {tails_info}
            
            <div class="plot-container">
                {plot_html}
            </div>
            
            <div class="stats-container">
                {table_html}
            </div>
        </div>
    </body>
    </html>
    """
    
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(full_html)
    
    return html_path

def save_json_results(mutation_id: str, energy_type: str, data_type: str,
                     distribution_stats: dict, gpd_result_right: dict, gpd_result_left: dict,
                     survival_data: dict, wilcoxon_results: dict, output_path: Path):
    """Сохраняет результаты в JSON (с двумя GPD результатами)."""
    
    json_results = {
        'mutation_id': mutation_id,
        'energy_type': energy_type,
        'data_type': data_type,
        'distribution_stats': distribution_stats,
        'gpd_analysis_right': gpd_result_right,  # Правый хвост
        'gpd_analysis_left': gpd_result_left,    # Левый хвост
        'survival_analysis': survival_data,
        'wilcoxon_analysis': wilcoxon_results,
        'generated_at': pd.Timestamp.now().isoformat()
    }
    
    json_path = output_path / f"{mutation_id}_{energy_type}_{data_type}_full_results.json"
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(json_results, f, indent=2, ensure_ascii=False)
    
    return json_path