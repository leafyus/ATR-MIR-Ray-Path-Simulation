from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np

from .geometry import Ray, Polygon


def normalize(vector: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(vector)
    if norm == 0:
        return vector
    return vector / norm


def dot(v1: np.ndarray, v2: np.ndarray) -> float:
    return float(np.dot(v1, v2))


def cross2d(v1: np.ndarray, v2: np.ndarray) -> float:
    # 2D scalar cross product (z-component of 3D cross)
    return float(v1[0] * v2[1] - v1[1] * v2[0])


def find_closest_intersection(ray: Ray, polygon: Polygon, epsilon: float = 1e-9) -> Optional[Dict]:
    """Find the closest valid intersection point of a ray with a polygon's edges.

    Returns a dict with keys: intersection_point, segment_index, normal, distance_t
    or None if no valid intersections are found.
    """
    origin = ray.origin
    direction = ray.direction

    closest_t: float = float("inf")
    closest: Optional[Dict] = None

    segments = polygon.get_segments_and_normals()
    for idx, seg in enumerate(segments):
        q1 = seg["q1"]
        q2 = seg["q2"]
        s = q2 - q1

        denom = cross2d(direction, s)
        # Parallel or collinear
        if abs(denom) < epsilon:
            continue

        q1_minus_p = q1 - origin
        t = cross2d(q1_minus_p, s) / denom
        u = cross2d(q1_minus_p, direction) / denom

        if t > epsilon and 0.0 <= u <= 1.0:
            if t < closest_t:
                intersection_point = origin + t * direction
                closest_t = t
                closest = {
                    "intersection_point": intersection_point,
                    "segment_index": idx,
                    "normal": seg["normal"],
                    "distance_t": t,
                }

    return closest


def reflect(direction_in: np.ndarray, normal_outward: np.ndarray) -> np.ndarray:
    # Specular reflection: R = I - 2*(I·N)*N
    return normalize(direction_in - 2.0 * dot(direction_in, normal_outward) * normal_outward)


def _incident_context(direction_in: np.ndarray, normal_outward: np.ndarray, n_inside: float, n_outside: float) -> Tuple[float, float, np.ndarray, float, bool]:
    """Compute incident medium, transmitted medium, normal into incident medium, cos(theta_i), and inside->outside flag."""
    I = normalize(direction_in)
    N = normalize(normal_outward)
    d = dot(I, N)
    if d > 0.0:
        # inside -> outside; flip N to point into incident medium
        n_incident = n_inside
        n_transmitted = n_outside
        n_inc_normal = -N
        inside_to_outside = True
    else:
        # outside -> inside; keep N pointing into incident medium (outside)
        n_incident = n_outside
        n_transmitted = n_inside
        n_inc_normal = N
        inside_to_outside = False
    cos_theta_i = -dot(n_inc_normal, I)
    # clamp for safety
    cos_theta_i = max(-1.0, min(1.0, cos_theta_i))
    return n_incident, n_transmitted, n_inc_normal, cos_theta_i, inside_to_outside


def compute_incident_debug(direction_in: np.ndarray, normal_outward: np.ndarray, n_inside: float, n_outside: float) -> Dict[str, float]:
    I = normalize(direction_in)
    n_inc, n_trn, n_inc_normal, cos_theta_i, inside_to_outside = _incident_context(I, normal_outward, n_inside, n_outside)
    sin_theta_i_sq = max(0.0, 1.0 - cos_theta_i * cos_theta_i)
    sin_theta_i = np.sqrt(sin_theta_i_sq)
    angle_i_deg = float(np.degrees(np.arccos(max(-1.0, min(1.0, cos_theta_i)))))
    critical_deg = None
    if n_inc > n_trn:
        ratio = n_trn / n_inc
        critical_deg = float(np.degrees(np.arcsin(ratio)))
    return {
        "n_incident": float(n_inc),
        "n_transmitted": float(n_trn),
        "cos_theta_i": float(cos_theta_i),
        "sin_theta_i": float(sin_theta_i),
        "angle_i_deg": angle_i_deg,
        "critical_deg": -1.0 if critical_deg is None else float(critical_deg),
        "inside_to_outside": 1.0 if inside_to_outside else 0.0,
    }


def refract(direction_in: np.ndarray, normal_outward: np.ndarray, n_inside: float, n_outside: float) -> Tuple[Optional[np.ndarray], bool, bool]:
    """Return (refracted_direction, is_tir, inside_to_outside).

    Uses a robust TIR check based on incident angle.

    normal_outward must point from inside to outside of the polygon.
    """
    I = normalize(direction_in)
    n_inc, n_trn, n_inc_normal, cos_theta_i, inside_to_outside = _incident_context(I, normal_outward, n_inside, n_outside)

    sin_theta_i_sq = max(0.0, 1.0 - cos_theta_i * cos_theta_i)
    # TIR if going from higher to lower index and exceeding critical angle
    if n_inc > n_trn:
        critical_ratio = n_trn / n_inc
        if sin_theta_i_sq > critical_ratio * critical_ratio:
            return None, True, inside_to_outside

    eta = n_inc / n_trn
    k = 1.0 - eta * eta * (1.0 - cos_theta_i * cos_theta_i)
    if k < 0.0:
        return None, True, inside_to_outside
    T = eta * I + (eta * cos_theta_i - np.sqrt(k)) * n_inc_normal
    return normalize(T), False, inside_to_outside


def calculate_interaction(direction_in: np.ndarray, normal_outward: np.ndarray, n_inside: float, n_outside: float) -> Tuple[np.ndarray, str, bool]:
    """Calculate outgoing direction.

    Returns (new_direction, mode, inside_to_outside) where mode is 'reflect' or 'refract'.
    inside_to_outside is True only for refraction when exiting the crystal.
    """
    refracted, is_tir, inside_to_outside = refract(direction_in, normal_outward, n_inside, n_outside)
    if is_tir:
        return reflect(direction_in, normal_outward), "reflect", False
    assert refracted is not None
    return refracted, "refract", inside_to_outside


