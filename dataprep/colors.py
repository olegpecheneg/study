import colorsys
from typing import List, Tuple


def generate_distinct_colors(n: int) -> List[Tuple[float, float, float]]:
    colors = []
    for i in range(n):
        hue = (i * 0.618033988749895) % 1.0
        saturation = 0.9
        value = 0.9
        rgb = colorsys.hsv_to_rgb(hue, saturation, value)
        colors.append(rgb)
    return colors
