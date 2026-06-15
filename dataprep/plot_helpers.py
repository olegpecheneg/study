from typing import Any, List, Dict, Tuple, Optional
import numpy as np
from matplotlib.lines import Line2D
from logging import getLogger

logger = getLogger(__name__)


def plot_scatter_points(ax: Any, ref_data: List[float], alt_data: List[float], snp_values: List[Optional[int]], snp_colors: Dict[int, Tuple[float, float, float]], upper_outliers: np.ndarray, lower_outliers: np.ndarray, normal_points: np.ndarray) -> None:
    gray_rgb = (0.5, 0.5, 0.5)
    colors = [snp_colors.get(snp, gray_rgb) for snp in snp_values]
    ax.scatter(
        np.array(ref_data)[normal_points], 
        np.array(alt_data)[normal_points], 
        c=np.array(colors)[normal_points], 
        alpha=0.7, s=80, edgecolor='black', linewidth=0.5
    )
    ax.scatter(
        np.array(ref_data)[upper_outliers], 
        np.array(alt_data)[upper_outliers], 
        c=np.array(colors)[upper_outliers], 
        alpha=0.9, s=150, edgecolor='green', linewidth=3
    )
    ax.scatter(
        np.array(ref_data)[lower_outliers], 
        np.array(alt_data)[lower_outliers], 
        c=np.array(colors)[lower_outliers], 
        alpha=0.9, s=150, edgecolor='red', linewidth=3
    )


def add_diagonal_line(ax: Any, min_e: float, max_e: float) -> None:
    ax.plot([min_e, max_e], [min_e, max_e], 'k--', linewidth=2, alpha=0.7)


def add_outlier_zones(ax: Any, x: np.ndarray, mean_diff: float, std_diff: float, min_e: float, max_e: float) -> None:
    line_upper = x - (mean_diff + 2 * std_diff)
    line_lower = x - (mean_diff - 2 * std_diff)
    ax.fill_between(x, min_e, line_upper, color='red', alpha=0.1, label='Верхние выбросы (+2std)')
    ax.fill_between(x, line_lower, max_e, color='green', alpha=0.1, label='Нижние выбросы (-2std)')


def create_legend_elements(snp_colors: Dict[int, Tuple[float, float, float]]) -> List[Line2D]:
    legend_elements = []
    gray_rgb = (0.5, 0.5, 0.5)
    legend_snps = list(snp_colors.keys())
    if len(legend_snps) > 20:
        legend_elements.append(
            Line2D([0], [0], marker='o', color='w', markerfacecolor=gray_rgb, 
                   markersize=12, label=f'Другие SNP ({len(snp_colors)-20})', markeredgecolor='black')
        )
        legend_snps = legend_snps[:20]
    for snp in legend_snps:
        legend_elements.append(
            Line2D([0], [0], marker='o', color='w', markerfacecolor=snp_colors[snp], 
                   markersize=12, label=f'SNP {snp}', markeredgecolor='black')
        )
    legend_elements.append(
        Line2D([0], [0], marker='o', color='w', markerfacecolor=gray_rgb, 
               markersize=12, label='Без SNP/Ошибка', markeredgecolor='black')
    )
    legend_elements.append(
        Line2D([0], [0], marker='o', color='w', markerfacecolor=gray_rgb, 
               markeredgecolor='green', markersize=12, label='Верхние выбросы (+2std)', linewidth=3)
    )
    legend_elements.append(
        Line2D([0], [0], marker='o', color='w', markerfacecolor=gray_rgb, 
               markeredgecolor='red', markersize=12, label='Нижние выбросы (-2std)', linewidth=3)
    )
    legend_elements.append(
        Line2D([0], [0], color='k', linestyle='--', linewidth=2, label='Диагональ (x=y)')
    )
    return legend_elements
