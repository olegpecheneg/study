from typing import List, Tuple
import numpy as np


def calculate_outlier_stats(ref_data: List[float], alt_data: List[float], std_mult: float = 2.0) -> Tuple[float, float, np.ndarray, np.ndarray, np.ndarray]:
    diff = np.array(ref_data) - np.array(alt_data)
    mean_diff = np.mean(diff)
    std_diff = np.std(diff)
    upper_outliers = diff > mean_diff + std_mult * std_diff
    lower_outliers = diff < mean_diff - std_mult * std_diff
    normal_points = ~(upper_outliers | lower_outliers)
    return mean_diff, std_diff, upper_outliers, lower_outliers, normal_points
