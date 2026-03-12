"""
Parameter data classes for the rainfall-infiltration-runoff model.
"""

from dataclasses import dataclass


@dataclass
class SoilParameters:
    """Physical properties of the soil column."""
    D: float         # Effective soil depth [mm]
    theta_s: float   # Saturated moisture content [m^3/m^3]
    theta_0: float   # Initial moisture content [m^3/m^3]
    K_s: float       # Saturated hydraulic conductivity [mm/h]


@dataclass
class HortonParams:
    """
    Horton's empirical infiltration parameters.
    Equation: I(t) = f_c + (f_0 - f_c) * exp(-k * t)
    """
    f_0: float  # Initial infiltration capacity [mm/h]
    f_c: float  # Steady-state infiltration capacity [mm/h]
    k: float    # Decay constant [1/h]


@dataclass
class GreenAmptParams:
    """
    Green-Ampt physically-based infiltration parameters.
    Equation: I(t) = K_s * (1 + (psi * delta_theta) / F_cum)
    """
    K_s: float  # Saturated hydraulic conductivity [mm/h]
    psi: float  # Wetting front suction head [mm]
