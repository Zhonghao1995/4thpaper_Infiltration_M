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


@dataclass
class RichardsParams:
    """
    Van Genuchten-Mualem parameters for Richards equation solver.

    Constitutive relations:
        Se      = [1 + |alpha*h|^n]^(-m),  m = 1 - 1/n
        theta   = theta_r + (theta_s - theta_r) * Se
        K       = K_s * Se^l * [1 - (1 - Se^(1/m))^m]^2
    """
    theta_r: float          # Residual water content [m^3/m^3]
    alpha: float            # VG inverse air-entry value [1/mm]
    n_vg: float             # VG pore-size distribution index [-]
    l: float = 0.5          # Mualem pore connectivity parameter [-]
    dz: float = 5.0         # Spatial step size [mm]
    max_iter: int = 50      # Max Picard iterations per time step
    tol: float = 1e-6       # Picard convergence tolerance [mm]
