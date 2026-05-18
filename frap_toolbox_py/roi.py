from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Tuple

import numpy as np
from skimage.draw import disk, polygon


@dataclass
class CircularROI:
    """Simple representation of a circular ROI in pixel coordinates."""

    center_x: float
    center_y: float
    radius: float

    def to_mask(self, shape: Tuple[int, int]) -> np.ndarray:
        # MATLAB specifies ROI centres in 1-based image coordinates where pixel
        # centres fall on integer locations. Convert to NumPy's 0-based index
        # space by subtracting one before rasterising the disk so the binary
        # mask matches the original toolbox geometry.
        row = self.center_y - 1.0
        col = self.center_x - 1.0
        rr, cc = disk((row, col), self.radius, shape=shape)
        mask = np.zeros(shape, dtype=bool)
        mask[rr, cc] = True
        return mask


@dataclass
class PolygonROI:
    """Representation of a user-defined polygon ROI."""

    vertices: Iterable[Tuple[float, float]]

    def to_mask(self, shape: Tuple[int, int]) -> np.ndarray:
        vertices = np.asarray(self.vertices, dtype=float)
        rr, cc = polygon(vertices[:, 1], vertices[:, 0], shape=shape)
        mask = np.zeros(shape, dtype=bool)
        mask[rr, cc] = True
        return mask


def adjacent_circle(base_circle: CircularROI, offset_factor: float = 2.5) -> CircularROI:
    """Return an adjacent circular ROI offset along the x-axis."""

    return CircularROI(
        center_x=base_circle.center_x + base_circle.radius * offset_factor,
        center_y=base_circle.center_y,
        radius=base_circle.radius,
    )
