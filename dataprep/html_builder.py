import json
import pandas as pd
import numpy as np
from pathlib import Path
import plotly.graph_objects as go
import plotly.io as pio
from typing import Dict, List, Optional, Tuple
import logging
from collections import defaultdict
import re

try:
    from .snvs_helper import (
        infer_snv_pair_for_master_json,
        resolve_point_position,
        load_snv_nucleotide_data,
    )
except ImportError:
    from snvs_helper import (
        infer_snv_pair_for_master_json,
        resolve_point_position,
        load_snv_nucleotide_data,
    )

logger = logging.getLogger(__name__)

MT_GENE_RANGES = {
    'ND1': (3307, 4262),
    'ND2': (4470, 5511),
    'CO1': (5904, 7445),
    'CO2': (7586, 8269),
    'ATP8': (8366, 8572),
    'ATP6': (8527, 9207),
    'CO3': (9207, 9990),
    'ND3': (10059, 10404),
    'ND4L': (10470, 10766),
    'ND4': (10760, 12141),
    'ND5': (12337, 14148),
    'ND6': (14149, 14673),
    'CYB': (14747, 15887),
}

def infer_mt_gene(position: Optional[int]) -> str:
    if position is None:
        return 'unknown'
    for gene_name, (start, end) in MT_GENE_RANGES.items():
        if start <= position <= end:
            return gene_name
    return 'non-coding'


def format_nucleotide_transition(position: Optional[int], snv_data: Dict[int, Tuple[str, str]]) -> str:
    if position is None or position == 0:
        return 'N/A'
    allele_pair = snv_data.get(position)
    if allele_pair and allele_pair[0] and allele_pair[1]:
        return f"{allele_pair[0]}→{allele_pair[1]}"
    return 'N/A'


def rgb_to_hex(rgb_tuple: Tuple[float, float, float]) -> str:
    """Конвертация RGB кортежа в HEX строку."""
    return '#%02x%02x%02x' % tuple(int(255 * max(0, min(1, v))) for v in rgb_tuple)

def generate_distinct_colors(n: int) -> List[Tuple[float, float, float]]:
    """Генерация списка визуально различимых RGB цветов."""
    import colorsys
    colors = []
    for i in range(n):
        hue = (i * 0.618033988749895) % 1.0
        saturation = 0.9
        value = 0.9
        rgb = colorsys.hsv_to_rgb(hue, saturation, value)
        colors.append(rgb)
    return colors

def generate_plotly_html_from_master_json(master_json_path: str, output_dir: Optional[str] = None) -> Dict[str, str]:
    """
    Генерация интерактивных графиков Plotly из master_data.json.
    
    Args:
        master_json_path: Путь к файлу master_data.json
        output_dir: Директория для сохранения HTML файлов (если None, используется директория JSON файла)
    
    Returns:
        Словарь с путями к созданным HTML файлам
    """
    
    # Настройка пути вывода
    master_path = Path(master_json_path)
    if output_dir is None:
        output_dir = master_path.parent
    else:
        output_dir = Path(output_dir)
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Загрузка данных
    try:
        with open(master_json_path, 'r', encoding='utf-8') as f:
            all_points = json.load(f)
        
        if not all_points:
            logger.error(f"Файл {master_json_path} пустой")
            return {}
            
        logger.info(f"Загружено {len(all_points)} точек из {master_json_path}")
        
    except Exception as e:
        logger.error(f"Ошибка загрузки {master_json_path}: {e}")
        return {}
    
    # Группировка по типу энергии
    points_by_energy = defaultdict(list)
    for point in all_points:
        energy_type = point.get('energy_type')
        if energy_type:
            points_by_energy[energy_type].append(point)
    
    # Определяем individual_id из первого элемента
    individual_id = all_points[0].get('construct_id', '').split('_')[-1] if all_points else 'unknown'
    
    generated_files = {}
    snp_info = infer_snv_pair_for_master_json(master_json_path)
    master_position = None
    master_ref = None
    master_alt = None
    master_mutation = 'N/A'
    if snp_info is not None:
        position, ref_allele, alt_allele = snp_info
        master_position = position
        master_ref = ref_allele
        master_alt = alt_allele
        master_mutation = f"{ref_allele}→{alt_allele}"

    # SNV map for per-point allele lookup (CSV + FASTA overrides)
    try:
        snv_map = load_snv_nucleotide_data(None, None, None)
    except Exception:
        snv_map = {}
    
    # Создаем график для каждого типа энергии
    for energy_type, points in points_by_energy.items():
        logger.info(f"Создание графика для {energy_type} ({len(points)} точек)")
        
        # Подготовка данных
        ref_energies = []
        alt_energies = []
        snp_values = []
        colors = []
        hover_texts = []
        point_indices = []
        
        snp_colors = {master_position: '#4ECDC4'}
        
        # Подготовка данных для каждой точки
        for idx, point in enumerate(points):
            ref_energies.append(float(point['ref_energy']))
            alt_energies.append(float(point['alt_energy']))
            
            mutation_position = resolve_point_position(point, master_json_path)
            if mutation_position is None:
                mutation_position = master_position
            point['snp'] = mutation_position if mutation_position is not None else None
            if mutation_position is not None:
                snp_values.append(str(mutation_position))
                # try to get allele pair for this mutation position
                allele_pair = None
                if mutation_position == master_position and master_ref and master_alt:
                    allele_pair = (master_ref, master_alt)
                else:
                    try:
                        allele_pair = snv_map.get(mutation_position)
                    except Exception:
                        allele_pair = None

                if allele_pair:
                    a_ref, a_alt = allele_pair
                    tpair = f"{a_ref}{a_alt}"
                    transition = 'transition' if tpair in {'AG','GA','CT','TC'} else 'transversion'
                    color = '#4ECDC4' if transition == 'transition' else '#FF6B6B'
                    transition_text = f"{a_ref}→{a_alt}"
                else:
                    color = snp_colors.get(mutation_position, '#808080')
                    transition_text = master_mutation
            else:
                snp_values.append(None)
                color = '#808080'
            colors.append(color)
            
            # Создание текста для hover
            outlier_status = point.get('outlier', 'Норма')
            gene_label = infer_mt_gene(mutation_position)
            transition_text = master_mutation
            hover_text = (
                f"Точка {idx}<br>"
                f"Конструкт: {point.get('construct_id', 'N/A')}<br>"
                f"Ref: {point['ref_energy']:.2f}, Alt: {point['alt_energy']:.2f}<br>"
                f"SNP: {mutation_position if mutation_position is not None else 'Нет'}<br>"
                f"Позиция: {mutation_position if mutation_position is not None else 'N/A'}<br>"
                f"Переход: {transition_text}<br>"
                f"Ген: {gene_label}<br>"
                f"Статус: {outlier_status}<br>"
                f"Центр: {point.get('center_position', 'N/A')}<br>"
                f"Плечо1: {point.get('arm1_start', 'N/A')}, Плечо2: {point.get('arm2_start', 'N/A')}<br>"
                f"Плечо3: {point.get('arm3_start', 'N/A')}, Плечо4: {point.get('arm4_start', 'N/A')}<br>"
            )
            
            # Добавляем информацию о последовательностях, если есть
            for w in range(1, 7):
                window_key = f'sequence_window{w}'
                if window_key in point and point[window_key] != 'N/A':
                    hover_text += f"Window{w}: {point[window_key]}<br>"
            
            hover_texts.append(hover_text)
            point_indices.append(idx)

        # Создание графика Plotly
        fig = go.Figure()
        fig.update_layout(template='plotly_white')
        fig.update_traces(marker=dict(line=dict(width=0)), marker_pattern_shape='')
        # UNIFIED approach: одна трасса со ВСЕМИ точками (все 4000 доступны для hover)
        # Маркеры различаются по размеру (выбросы больше) и символу (выбросы diamond)
        marker_sizes = [16 if p.get('outlier') in ['Верхний', 'Нижний'] else 8 for p in points]
        marker_symbols = ['diamond' if p.get('outlier') in ['Верхний', 'Нижний'] else 'circle' for p in points]
        
        fig.add_trace(go.Scatter(
            x=ref_energies,
            y=alt_energies,
            mode='markers',
            name='Точки',
            marker=dict(
                color=colors,
                size=marker_sizes,
                symbol=marker_symbols,
                line=dict(width=1, color='Black')
            ),
            hoverinfo='text',
            text=hover_texts,
            customdata=list(range(len(points))),
            showlegend=False
        ))
        
        gene_title = infer_mt_gene(master_position)
        # Диагональная линия
        min_val = min(min(ref_energies), min(alt_energies))
        max_val = max(max(ref_energies), max(alt_energies))
        fig.add_trace(go.Scatter(
            x=[min_val, max_val],
            y=[min_val, max_val],
            mode='lines',
            line=dict(dash='dash', color='black'),
            name='Диагональ',
            hoverinfo='skip'
        ))
        
        # Настройка layout
        n_points = len(points)
        mean_ref = float(np.mean(ref_energies)) if ref_energies else float('nan')
        mean_alt = float(np.mean(alt_energies)) if alt_energies else float('nan')
        mean_delta = mean_alt - mean_ref
        n_outliers_upper = sum(1 for p in points if p.get('outlier') == 'Верхний')
        n_outliers_lower = sum(1 for p in points if p.get('outlier') == 'Нижний')
        n_snp_points = sum(1 for s in snp_values if s)
        snp_list_str = ', '.join(sorted(set(s for s in snp_values if s))) if any(snp_values) else 'None'

        fig.update_layout(
            title=f"Сравнение {energy_type} для индивидуума {individual_id} (ген: {gene_title})<br><sub>Точек: {n_points} | Δ: {mean_delta:.2f} | SNP: {n_snp_points}</sub>",
            xaxis_title='Ref energy (kcal/mol)',
            yaxis_title='Alt energy (kcal/mol)',
            width=1200,
            height=800,
            showlegend=False,
            clickmode='event+select'
        )
        
        # JavaScript + HTML для интерактивности (аналогичный исходному коду)
        modal_js = f"""
        <style>
        /* Top search panel */
        #topSearchPanel {{
            position: absolute;
            z-index: 1002;
            right: 12px;
            top: 30px;
            background: rgba(255,255,255,0.97);
            padding: 10px;
            border-radius: 6px;
            box-shadow: 0 2px 12px rgba(0,0,0,0.12);
            font-family: Arial, sans-serif;
            font-size: 13px;
            width: 360px;
        }}
        
        /* Bottom filters panel */
        #filterPanel {{
            position: absolute;
            z-index: 1002;
            right: 12px;
            bottom: 30px;
            background: rgba(255,255,255,0.97);
            padding: 12px;
            border-radius: 6px;
            box-shadow: 0 2px 12px rgba(0,0,0,0.12);
            font-family: Arial, sans-serif;
            font-size: 13px;
            width: 360px;
            max-height: 48%;
            overflow: auto;
        }}
        
        .panel-row {{ display:flex; align-items:center; gap:8px; margin-top:6px }}
        .panel-row label {{ width:120px; font-size:12px; color:#333 }}
        .panel-row input, .panel-row select, .panel-row textarea {{ padding:6px; border:1px solid #ccc; border-radius:4px; flex:1 }}
        .panel-actions {{ display:flex; gap:8px; margin-top:8px }}
        .panel-actions button {{ padding:6px 8px; border-radius:4px; border:1px solid #999; background:#f7f7f7; cursor:pointer }}
        .panel-actions button:hover {{ background:#eee }}
        
        #resultsSummary {{ margin-top:8px; font-size:12px }}
        #resultsTable {{ margin-top:6px; border-top:1px solid #eee; padding-top:6px }}
        #resultsTable table {{ width:100%; border-collapse:collapse; font-size:12px }}
        #resultsTable th, #resultsTable td {{ padding:6px; border:1px solid #eee; text-align:left }}
        
        .modal {{ display: none; position: fixed; z-index: 1000; left: 0; top: 0; width: 100%; height: 100%; background-color: rgba(0,0,0,0.4); }}
        .modal-content {{ background-color: #fff; margin: 10% auto; padding: 18px; border-radius: 6px; width: 60%; max-width: 700px; max-height: 70vh; overflow-y: auto; box-shadow: 0 6px 18px rgba(0,0,0,0.16) }}
        .close {{ color: #777; float: right; font-size: 22px; font-weight: bold; cursor: pointer; position: sticky; top: 0; background: white; padding: 0 4px; }}
        .close:hover {{ color: #111 }}
        a.result-link {{ color:#0366d6; cursor:pointer; text-decoration:underline }}
        </style>
        
        <div id="topSearchPanel">
            <div style="font-weight:bold; margin-bottom:6px">Быстрый поиск (по всем полям)</div>
            <div class="panel-row">
                <label>Поиск (всё):</label>
                <input id="searchAny" placeholder="введите текст, ищет по всем полям" />
            </div>
            <div class="panel-actions">
                <button id="btnTopFind">Найти</button>
                <button id="btnTopClear">Сброс</button>
            </div>
            <div id="topInfo" style="margin-top:8px;font-size:12px;color:#444">Результатов: <span id="resultsCount">0</span></div>
        </div>
        
        <div id="filterPanel">
            <div style="font-weight:bold">Фильтрация (подробно)</div>
            
            <div class="panel-row">
                <label>ConstructID</label>
                <input id="filterConstruct" placeholder="точное или подстрока" />
            </div>
            <div style="display:flex;gap:8px;margin-top:6px;align-items:center">
                <label style="width:120px">Режим для ConstructID</label>
                <label><input type="radio" name="constructMode" value="exact" checked/> Точное</label>
                <label><input type="radio" name="constructMode" value="substring"/> Подстрока</label>
            </div>
            
            <div class="panel-row">
                <label>SNP (multi)</label>
                <select id="filterSnpMulti" multiple size="4"></select>
            </div>
            <div style="display:flex;gap:8px;margin-top:6px;align-items:center">
                <label style="width:120px">Режим SNP</label>
                <label><input type="radio" name="snpMode" value="exact" checked/> Точное</label>
                <label><input type="radio" name="snpMode" value="substring"/> Подстрока</label>
                <button id="btnSelectAllSnps" style="margin-left:auto">Выбрать все</button>
            </div>
            
            <div class="panel-row">
                <label>Energy min / max</label>
                <input id="filterMinEnergy" placeholder="min" />
                <input id="filterMaxEnergy" placeholder="max" />
            </div>
            
            <div style="margin-top:8px;font-weight:600">Поля плеч (arms)</div>
            <div class="panel-row">
                <label>Arm1</label>
                <input id="filterArm1" placeholder="значение или часть" />
            </div>
            <div class="panel-row">
                <label>Arm2</label>
                <input id="filterArm2" placeholder="значение или часть" />
            </div>
            <div class="panel-row">
                <label>Arm3</label>
                <input id="filterArm3" placeholder="значение или часть" />
            </div>
            <div class="panel-row">
                <label>Arm4</label>
                <input id="filterArm4" placeholder="значение или часть" />
            </div>
            
            <div style="margin-top:8px;font-weight:600">Поля Window (последовательности)</div>
            <div class="panel-row"><label>Window1</label><input id="filterW1" placeholder="часть последовательности" /></div>
            <div class="panel-row"><label>Window2</label><input id="filterW2" placeholder="часть последовательности" /></div>
            <div class="panel-row"><label>Window3</label><input id="filterW3" placeholder="часть последовательности" /></div>
            <div class="panel-row"><label>Window4</label><input id="filterW4" placeholder="часть последовательности" /></div>
            <div class="panel-row"><label>Window5</label><input id="filterW5" placeholder="часть последовательности" /></div>
            <div class="panel-row"><label>Window6</label><input id="filterW6" placeholder="часть последовательности" /></div>
            
            <div class="panel-actions">
                <button id="btnApply">Применить</button>
                <button id="btnClear">Очистить</button>
                <button id="btnDownloadCSV">Сохранить CSV</button>
                <button id="btnDownloadJSON">Сохранить JSON</button>
            </div>
            
            <div id="resultsTable"></div>
        </div>
        
        <div id="dataModal" class="modal"><div class="modal-content"><span class="close" id="modalClose">&times;</span><div id="dataTable"></div></div></div>
        
        <script>
        // Встраиваем данные точек
        window.POINT_DATA = {json.dumps(points, ensure_ascii=False)};
        
        function getPlotDiv() {{ 
            var p = document.getElementById('plotly'); 
            if(!p) p = document.getElementsByClassName('plotly-graph-div')[0]; 
            return p; 
        }}
        
        var ORIGINAL_MARKER_SIZES = [];
        var ORIGINAL_MARKER_COLORS = [];
        var LAST_FILTERED = [];
        
        function storeOriginalMarkerState(){{
            var plot = getPlotDiv(); 
            if(!plot) return; 
            var data = plot.data || [];
            for(var t=0; t<data.length; t++){{ 
                var m = data[t].marker||{{}}; 
                var size = m.size||10; 
                var color = m.color||'#808080';
                if(!Array.isArray(size)){{ 
                    var arr = new Array((data[t].x && data[t].x.length)||0).fill(size); 
                    size = arr; 
                }}
                if(!Array.isArray(color)){{ 
                    var carr = new Array((data[t].x && data[t].x.length)||0).fill(color); 
                    color = carr; 
                }}
                ORIGINAL_MARKER_SIZES[t] = size.slice(); 
                ORIGINAL_MARKER_COLORS[t] = color.slice(); 
            }}
        }}
        
        function resetAllMarkers(){{
            var plot = getPlotDiv(); 
            if(!plot) return; 
            for(var t=0; t<plot.data.length; t++){{ 
                var s = ORIGINAL_MARKER_SIZES[t]||[]; 
                var c = ORIGINAL_MARKER_COLORS[t]||[]; 
                try{{ 
                    Plotly.restyle(plot, {{'marker.size':[s], 'marker.color':[c]}}, [t]); 
                }}catch(e){{}}
            }} 
            document.getElementById('resultsCount').innerText = '0'; 
            document.getElementById('resultsTable').innerHTML = ''; 
            LAST_FILTERED = []; 
        }}
        
        function selectAllSnps(){{
            var sel = document.getElementById('filterSnpMulti'); 
            if(!sel) return; 
            for(var i=0; i<sel.options.length; i++){{
                sel.options[i].selected = true; 
            }}
        }}
        
        function applyFilters(){{
            var construct = document.getElementById('filterConstruct').value.trim();
            var constructMode = document.querySelector('input[name="constructMode"]:checked') ? 
                document.querySelector('input[name="constructMode"]:checked').value : 'exact';
            var snpSelect = document.getElementById('filterSnpMulti');
            var selectedSnps = [];
            if(snpSelect){{
                for(var i=0; i<snpSelect.options.length; i++){{
                    if(snpSelect.options[i].selected) selectedSnps.push(snpSelect.options[i].value);
                }}
            }}
            var snpMode = document.querySelector('input[name="snpMode"]:checked') ? 
                document.querySelector('input[name="snpMode"]:checked').value : 'exact';
            var minE = parseFloat(document.getElementById('filterMinEnergy').value);
            var maxE = parseFloat(document.getElementById('filterMaxEnergy').value);
            var arm1 = document.getElementById('filterArm1').value.trim();
            var arm2 = document.getElementById('filterArm2').value.trim();
            var arm3 = document.getElementById('filterArm3').value.trim();
            var arm4 = document.getElementById('filterArm4').value.trim();
            var w1 = document.getElementById('filterW1').value.trim();
            var w2 = document.getElementById('filterW2').value.trim();
            var w3 = document.getElementById('filterW3').value.trim();
            var w4 = document.getElementById('filterW4').value.trim();
            var w5 = document.getElementById('filterW5').value.trim();
            var w6 = document.getElementById('filterW6').value.trim();
            var any = document.getElementById('searchAny').value.trim().toLowerCase();
            
            var matches = [];
            for(var i=0; i<POINT_DATA.length; i++){{
                var r = POINT_DATA[i];
                var ok = true;
                
                // ConstructID filter
                if(construct){{
                    if(constructMode === 'exact'){{
                        if(String(r.construct_id) !== construct) ok = false;
                    }} else {{
                        if(!(String(r.construct_id||'').toLowerCase().indexOf(construct.toLowerCase()) !== -1)) ok = false;
                    }}
                }}
                
                // SNP multi filter
                if(selectedSnps.length > 0){{
                    var rSnp = r.snp === null ? '' : String(r.snp);
                    var matchedSnp = false;
                    for(var si=0; si<selectedSnps.length; si++){{
                        var sval = selectedSnps[si];
                        if(snpMode === 'exact'){{
                            if(rSnp === sval) matchedSnp = true;
                        }} else {{
                            if(rSnp.toLowerCase().indexOf(sval.toLowerCase()) !== -1) matchedSnp = true;
                        }}
                    }}
                    if(!matchedSnp) ok = false;
                }}
                
                // Energy range
                if(!isNaN(minE) || !isNaN(maxE)){{
                    var lo = isNaN(minE) ? -1e99 : minE;
                    var hi = isNaN(maxE) ? 1e99 : maxE;
                    if(!((r.ref_energy >= lo && r.ref_energy <= hi) || (r.alt_energy >= lo && r.alt_energy <= hi))) ok = false;
                }}
                
                // Arms filters
                if(arm1 && String(r.arm1_start||'').toLowerCase().indexOf(arm1.toLowerCase()) === -1) ok = false;
                if(arm2 && String(r.arm2_start||'').toLowerCase().indexOf(arm2.toLowerCase()) === -1) ok = false;
                if(arm3 && String(r.arm3_start||'').toLowerCase().indexOf(arm3.toLowerCase()) === -1) ok = false;
                if(arm4 && String(r.arm4_start||'').toLowerCase().indexOf(arm4.toLowerCase()) === -1) ok = false;
                
                // Windows filters
                if(w1 && String(r.sequence_window1||'').toLowerCase().indexOf(w1.toLowerCase()) === -1) ok = false;
                if(w2 && String(r.sequence_window2||'').toLowerCase().indexOf(w2.toLowerCase()) === -1) ok = false;
                if(w3 && String(r.sequence_window3||'').toLowerCase().indexOf(w3.toLowerCase()) === -1) ok = false;
                if(w4 && String(r.sequence_window4||'').toLowerCase().indexOf(w4.toLowerCase()) === -1) ok = false;
                if(w5 && String(r.sequence_window5||'').toLowerCase().indexOf(w5.toLowerCase()) === -1) ok = false;
                if(w6 && String(r.sequence_window6||'').toLowerCase().indexOf(w6.toLowerCase()) === -1) ok = false;
                
                // Any-field search
                if(any){{
                    var j = JSON.stringify(r).toLowerCase();
                    if(j.indexOf(any) === -1) ok = false;
                }}
                
                if(ok) matches.push(r);
            }}
            
            LAST_FILTERED = matches.slice();
            highlightFiltered(matches.map(function(x){{ return x.index; }}));
            populateResults(matches);
        }}
        
        function highlightFiltered(indices){{
            var plot = getPlotDiv(); 
            if(!plot) return; 
            for(var t=0; t<plot.data.length; t++){{
                var data = plot.data[t];
                var sizes = ORIGINAL_MARKER_SIZES[t] ? ORIGINAL_MARKER_SIZES[t].slice() : 
                    new Array((data.x && data.x.length)||0).fill(10);
                var colors = ORIGINAL_MARKER_COLORS[t] ? ORIGINAL_MARKER_COLORS[t].slice() : 
                    new Array((data.x && data.x.length)||0).fill('#808080');
                var custom = data.customdata||[];
                for(var j=0; j<custom.length; j++){{
                    if(indices.indexOf(custom[j]) !== -1){{
                        sizes[j] = Math.max(sizes[j]||10, 18);
                        colors[j] = '#FFD54F';
                    }}
                }}
                try{{
                    Plotly.restyle(plot, {{'marker.size':[sizes], 'marker.color':[colors]}}, [t]);
                }}catch(e){{}}
            }}
            if(indices.length > 0){{
                var first = indices[0];
                try{{
                    Plotly.relayout(plot, {{
                        'xaxis.range': [POINT_DATA[first].ref_energy-1, POINT_DATA[first].ref_energy+1],
                        'yaxis.range': [POINT_DATA[first].alt_energy-1, POINT_DATA[first].alt_energy+1]
                    }});
                }}catch(e){{}}
            }}
            document.getElementById('resultsCount').innerText = indices.length;
        }}
        
        function populateResults(matches){{
            var html = '<table><thead><tr><th>#</th><th>Construct</th><th>SNP</th><th>Ref</th><th>Alt</th><th>Energy</th></tr></thead><tbody>';
            for(var i=0; i<matches.length; i++){{
                var r = matches[i];
                html += '<tr><td>'+r.index+'</td><td><a class="result-link" data-idx="'+r.index+'">'+r.construct_id+'</a></td><td>'+
                    (r.snp===null?'':r.snp)+'</td><td>'+r.ref_energy.toFixed(3)+'</td><td>'+r.alt_energy.toFixed(3)+
                    '</td><td>'+r.energy_type+'</td></tr>';
            }}
            html += '</tbody></table>';
            document.getElementById('resultsTable').innerHTML = html;
            
            // Добавляем обработчики кликов на ссылки
            var links = document.querySelectorAll('.result-link');
            for(var j=0; j<links.length; j++){{
                links[j].addEventListener('click', function(e){{
                    var idx = e.target.getAttribute('data-idx');
                    openModalByIndex(parseInt(idx));
                }});
            }}
        }}
        
        function openModalByIndex(idx){{
            var rec = (window.POINT_DATA && window.POINT_DATA[idx]) ? window.POINT_DATA[idx] : null;
            if(rec) openModal(rec);
        }}
        
        function openModal(data){{
            var modal = document.getElementById('dataModal');
            var table = document.getElementById('dataTable');
            var content = '<table style="width:100%; border-collapse:collapse">';
            content += '<tr style="background:#f6f6f6"><th style="padding:8px;border:1px solid #eee">Поле</th>'+
                       '<th style="padding:8px;border:1px solid #eee">Значение</th></tr>';
            for(var k in data){{
                var v = data[k];
                if(Array.isArray(v)) v = v.join(', ');
                content += '<tr><td style="padding:8px;border:1px solid #eee;font-weight:bold">'+k+
                           '</td><td style="padding:8px;border:1px solid #eee">'+String(v)+'</td></tr>';
            }}
            content += '</table>';
            table.innerHTML = content;
            modal.style.display = 'block';
        }}
        
        function closeModal(){{
            document.getElementById('dataModal').style.display = 'none';
        }}
        
        function downloadCSV(){{
            var rows = LAST_FILTERED;
            if(!rows || !rows.length){{
                alert('Нет данных для сохранения');
                return;
            }}
            var keys = Object.keys(rows[0]);
            var csv = keys.join(',') + '\\n';
            for(var i=0; i<rows.length; i++){{
                var vals = keys.map(function(k){{
                    return '"'+String(rows[i][k]).replace(/"/g,'""')+'"';
                }});
                csv += vals.join(',') + '\\n';
            }}
            var blob = new Blob([csv], {{type:'text/csv;charset=utf-8;'}});
            var url = URL.createObjectURL(blob);
            var a = document.createElement('a');
            a.href = url;
            a.download = 'filtered_points_energy.csv';
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            URL.revokeObjectURL(url);
        }}
        
        function downloadJSON(){{
            var rows = LAST_FILTERED;
            if(!rows || !rows.length){{
                alert('Нет данных для сохранения');
                return;
            }}
            var s = JSON.stringify(rows, null, 2);
            var blob = new Blob([s], {{type:'application/json;charset=utf-8;'}});
            var url = URL.createObjectURL(blob);
            var a = document.createElement('a');
            a.href = url;
            a.download = 'filtered_points_energy.json';
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            URL.revokeObjectURL(url);
        }}
        
        function clearFilters(){{
            document.getElementById('filterConstruct').value = '';
            var sel = document.getElementById('filterSnpMulti');
            if(sel){{
                for(var i=0; i<sel.options.length; i++){{
                    sel.options[i].selected = false;
                }}
            }}
            document.getElementById('filterMinEnergy').value = '';
            document.getElementById('filterMaxEnergy').value = '';
            document.getElementById('filterArm1').value = '';
            document.getElementById('filterArm2').value = '';
            document.getElementById('filterArm3').value = '';
            document.getElementById('filterArm4').value = '';
            document.getElementById('filterW1').value = '';
            document.getElementById('filterW2').value = '';
            document.getElementById('filterW3').value = '';
            document.getElementById('filterW4').value = '';
            document.getElementById('filterW5').value = '';
            document.getElementById('filterW6').value = '';
            document.getElementById('searchAny').value = '';
            resetAllMarkers();
        }}
        
        function populateSnpOptions(){{
            try{{
                var sel = document.getElementById('filterSnpMulti');
                if(!sel) return;
                var unique = {{}};
                for(var i=0; i<POINT_DATA.length; i++){{
                    var s = POINT_DATA[i].snp;
                    if(s !== null && s !== undefined && String(s) !== '') unique[String(s)] = true;
                }}
                var keys = Object.keys(unique).sort(function(a,b){{return Number(a)-Number(b);}});
                for(var k=0; k<keys.length; k++){{
                    var opt = document.createElement('option');
                    opt.value = keys[k];
                    opt.text = keys[k];
                    sel.appendChild(opt);
                }}
            }}catch(e){{}}
        }}
        
        // Инициализация после загрузки
        window.addEventListener('load', function(){{
            setTimeout(function(){{
                try{{
                    storeOriginalMarkerState();
                    populateSnpOptions();
                    
                    // Привязка обработчиков событий
                    document.getElementById('btnTopFind').addEventListener('click', applyFilters);
                    document.getElementById('btnTopClear').addEventListener('click', clearFilters);
                    document.getElementById('btnApply').addEventListener('click', applyFilters);
                    document.getElementById('btnClear').addEventListener('click', clearFilters);
                    document.getElementById('btnDownloadCSV').addEventListener('click', downloadCSV);
                    document.getElementById('btnDownloadJSON').addEventListener('click', downloadJSON);
                    document.getElementById('btnSelectAllSnps').addEventListener('click', selectAllSnps);
                    document.getElementById('modalClose').addEventListener('click', closeModal);
                    
                    // Клик по графику
                    var plot = getPlotDiv();
                    if(plot){{
                        plot.on('plotly_click', function(data){{
                            if(data.points && data.points[0]){{
                                var idx = data.points[0].customdata;
                                openModalByIndex(idx);
                            }}
                        }});
                    }}
                    
                    // Закрытие модального окна при клике вне его
                    window.onclick = function(event){{
                        var modal = document.getElementById('dataModal');
                        if(event.target == modal) modal.style.display = 'none';
                    }};
                    
                }}catch(e){{
                    console.error('Ошибка инициализации:', e);
                }}
            }}, 300);
        }});
        </script>
        """
        
        # Добавляем обработчик кликов
        click_handler = """
        <script>
        document.addEventListener('DOMContentLoaded', function() {
            var plot = document.getElementsByClassName('plotly-graph-div')[0];
            if(plot) {
                plot.on('plotly_click', function(data) {
                    if(data.points && data.points[0]) {
                        var idx = data.points[0].customdata;
                        if(window.openModalByIndex) {
                            window.openModalByIndex(idx);
                        }
                    }
                });
            }
        });
        </script>
        """
        
        # Сохраняем HTML файл
        html_path = output_dir / f"{energy_type}_comparison_interactive.html"
        config = {
            'displayModeBar': True,
            'scrollZoom': True,
            'modeBarButtonsToAdd': ['select2d', 'lasso2d'],
        }
        
        try:
            # Генерируем HTML с Plotly
            html_content = pio.to_html(
                fig,
                config=config,
                include_plotlyjs='cdn',
                full_html=True,
                validate=True
            )
            
            # Вставляем JavaScript и данные перед закрывающим тегом body
            html_content = html_content.replace('</body>', f'{modal_js}{click_handler}</body>')
            
            with open(html_path, 'w', encoding='utf-8') as f:
                f.write(html_content)
            
            generated_files[energy_type] = str(html_path)
            logger.info(f"Интерактивный график сохранён: {html_path}")
            
        except Exception as e:
            logger.error(f"Ошибка сохранения HTML для {energy_type}: {e}")
    
    return generated_files


# Функция для запуска восстановления графиков из всех найденных master_data.json
def restore_all_interactive_plots(base_output_dir: str, max_files: Optional[int] = None):
    """
    Восстановление интерактивных графиков из всех найденных master_data.json.

    Args:
        base_output_dir: Базовая директория с результатами (содержит поддиректории test_individual_*)
        max_files: максимальное число файлов для обработки за один запуск.
    """
    base_path = Path(base_output_dir)
    # Find all master_data json files recursively
    master_files = sorted(base_path.rglob('*_master_data.json'))
    if max_files is not None:
        master_files = master_files[:max_files]

    logger.info(f"Найдено {len(master_files)} master_data.json файлов для обработки")

    all_generated = {}
    for idx, master_file in enumerate(master_files, start=1):
        logger.info(f"[{idx}/{len(master_files)}] Обработка: {master_file}")
        try:
            generated = generate_plotly_html_from_master_json(str(master_file))
            all_generated[str(master_file)] = generated
        except Exception as e:
            logger.error(f"Ошибка обработки {master_file}: {e}")

    return all_generated


# Пример использования:
if __name__ == "__main__":
    import sys
    
    # Настройка логирования
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    if len(sys.argv) > 1:
        if sys.argv[1] == '--restore-all':
            base_dir = sys.argv[2] if len(sys.argv) > 2 else None
            max_files = int(sys.argv[3]) if len(sys.argv) > 3 else None
            if base_dir is None:
                raise ValueError('Укажите базовую директорию для --restore-all')
            print(f"Восстановление интерактивных графиков в: {base_dir} (до {max_files if max_files else 'всех'} файлов)")
            result = restore_all_interactive_plots(base_dir, max_files=max_files)
            print(f"Готово: обработано {len(result)} файлов")
        else:
            # Если передали путь к master_data.json как аргумент
            master_json_path = sys.argv[1]
            output_dir = sys.argv[2] if len(sys.argv) > 2 else None
            print(f"Генерация графиков из: {master_json_path}")
            result = generate_plotly_html_from_master_json(master_json_path, output_dir)
            print(f"Создано {len(result)} графиков:")
            for energy_type, path in result.items():
                print(f"  {energy_type}: {path}")
    else:
        # Восстановление всех графиков из указанной директории
        base_dir = input("Введите путь к директории с результатами (например, D:/pythonProject/MitoFragility/DataPreparing/output): ")
        result = restore_all_interactive_plots(base_dir)
        print(f"Готово: обработано {len(result)} файлов")
        print(f"Поиск master_data.json в {base_dir}...")
        
        restore_all_interactive_plots(base_dir)