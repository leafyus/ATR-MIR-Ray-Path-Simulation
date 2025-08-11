from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple, Dict

import numpy as np


def normalize(vector: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(vector)
    if norm == 0:
        return vector
    return vector / norm


@dataclass
class Ray:
    origin: np.ndarray
    direction: np.ndarray

    def __post_init__(self) -> None:
        self.origin = np.asarray(self.origin, dtype=float)
        self.direction = normalize(np.asarray(self.direction, dtype=float))


class Polygon:
    def __init__(self, vertices: List[Tuple[float, float]], material_n: float, default_boundary_n: float) -> None:
        self.vertices: List[np.ndarray] = [np.array(v, dtype=float) for v in vertices]
        self.material_n: float = float(material_n)
        # One boundary per edge (same count as vertices)
        self.boundary_media_n: List[float] = [float(default_boundary_n)] * len(self.vertices)

    def set_boundary_medium(self, boundary_index: int, medium_n: float) -> None:
        if 0 <= boundary_index < len(self.boundary_media_n):
            self.boundary_media_n[boundary_index] = float(medium_n)

    def get_segments_and_normals(self) -> List[Dict[str, np.ndarray]]:
        # Assumes vertices are ordered counter-clockwise
        segments: List[Dict[str, np.ndarray]] = []
        num_vertices = len(self.vertices)
        for i in range(num_vertices):
            q1 = self.vertices[i]
            q2 = self.vertices[(i + 1) % num_vertices]
            s = q2 - q1
            # For CCW polygons, outward normal is right-hand (clockwise) rotation by 90 deg
            outward_normal = normalize(np.array([s[1], -s[0]], dtype=float))
            segments.append({
                "q1": q1,
                "q2": q2,
                "normal": outward_normal,
            })
        return segments


