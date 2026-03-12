"""
Horton's empirical infiltration model.

I(t) = f_c + (f_0 - f_c) * exp(-k * t)

Reference: Horton, R.E. (1941). Soil Science Society of America Journal.
"""

import numpy as np
from infiltration_model.core.parameters import HortonParams


def horton_infiltration_rate(t_since_rain: float, params: HortonParams) -> float:
    """
    Compute infiltration rate using Horton's equation.

    Parameters
    ----------
    t_since_rain : float  — Time since rainfall onset [h]
    params : HortonParams — Model parameters (f_0, f_c, k)

    Returns
    -------
    float — Instantaneous infiltration rate [mm/h]
    """
    return params.f_c + (params.f_0 - params.f_c) * np.exp(-params.k * t_since_rain)
