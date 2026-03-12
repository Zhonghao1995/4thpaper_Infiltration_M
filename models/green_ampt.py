"""
Green-Ampt physically-based infiltration model.

I(t) = K_s * (1 + (psi * delta_theta) / F_cum)

Reference: Green, W.H. and Ampt, G.A. (1911). Journal of Agricultural Science.
"""

from infiltration_model.core.parameters import GreenAmptParams


def green_ampt_infiltration_rate(
    F_cum: float,
    theta: float,
    params: GreenAmptParams,
    theta_s: float
) -> float:
    """
    Compute infiltration rate using the Green-Ampt equation.

    Parameters
    ----------
    F_cum   : float — Cumulative infiltration since rainfall onset [mm]
    theta   : float — Current soil moisture [m^3/m^3]
    params  : GreenAmptParams — Model parameters (K_s, psi)
    theta_s : float — Saturated soil moisture [m^3/m^3]

    Returns
    -------
    float — Instantaneous infiltration rate [mm/h]
    """
    delta_theta = theta_s - theta
    F_cum = max(F_cum, 0.1)  # Prevent division by zero at t=0
    return params.K_s * (1.0 + (params.psi * delta_theta) / F_cum)
