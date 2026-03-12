"""
Core simulation engine and verification utilities.
"""

import numpy as np
from typing import Union
from infiltration_model.core.parameters import SoilParameters, HortonParams, GreenAmptParams, RichardsParams
from infiltration_model.models.horton import horton_infiltration_rate
from infiltration_model.models.green_ampt import green_ampt_infiltration_rate
from infiltration_model.models.richards import solve_richards


def simulate(
    precip: np.ndarray,
    dt: float,
    soil: SoilParameters,
    infiltration_model: str,
    infil_params: Union[HortonParams, GreenAmptParams, RichardsParams]
) -> dict:
    """
    Run forward rainfall-infiltration-runoff simulation.

    At each time step:
        1. Compute remaining storage capacity: SC = D * (theta_s - theta)
        2. Compute infiltration rate I from sub-model (Horton or Green-Ampt)
        3. Compute infiltration depth: F = I * Dt (forward Euler)
        4. Apply constraints: F = min(F, P, SC)
        5. Update soil moisture: theta += F / D
        6. Compute runoff: Q = (P - F) / Dt

    If infiltration_model == 'richards', delegates directly to solve_richards.

    Parameters
    ----------
    precip             : Rainfall depth per time step [mm], array of length n
    dt                 : Time step duration [h]
    soil               : Soil physical properties
    infiltration_model : 'horton', 'green_ampt', or 'richards'
    infil_params       : Parameters for the chosen model

    Returns
    -------
    dict with arrays: time, precip_intensity, theta, infil_rate,
                      infil_depth, runoff_intensity, storage_capacity, F_cumulative
    """
    if infiltration_model == 'richards':
        if not isinstance(infil_params, RichardsParams):
            raise TypeError("Expected RichardsParams for 'richards' model")
        return solve_richards(precip, dt, soil, infil_params)

    if infiltration_model not in ('horton', 'green_ampt'):
        raise ValueError(f"Unknown model: '{infiltration_model}'")

    n = len(precip)

    # Output arrays
    time_arr = np.zeros(n)
    theta_arr = np.zeros(n)
    infil_rate_arr = np.zeros(n)
    infil_depth_arr = np.zeros(n)
    runoff_arr = np.zeros(n)
    sc_arr = np.zeros(n)
    F_cum_arr = np.zeros(n)
    precip_int_arr = np.zeros(n)

    # Initial conditions
    theta_current = soil.theta_0
    F_cumulative = 0.0
    t_elapsed = 0.0

    for i in range(n):
        t_elapsed += dt
        time_arr[i] = t_elapsed
        P_i = precip[i]
        precip_int_arr[i] = P_i / dt

        # Step 1: Remaining storage capacity [mm]
        SC = max(soil.D * (soil.theta_s - theta_current), 0.0)
        sc_arr[i] = SC

        # Step 2: Infiltration rate from sub-model (Trapezoidal: Average of start and end timestep rates)
        if infiltration_model == 'horton':
            I_start = horton_infiltration_rate(t_elapsed - dt, infil_params)
            I_end = horton_infiltration_rate(t_elapsed, infil_params)
            I_t = (I_start + I_end) / 2.0  # Average rate for the timestep
        else:
            # For Green-Ampt, we use the rate based on current moisture deficit
            # A true trapezoidal rule for GA requires an implicit solver.
            # We approximate by using the capacity at the start of the timestep.
            I_t = green_ampt_infiltration_rate(
                F_cumulative, theta_current, infil_params, soil.theta_s
            )

        # Step 3: Potential infiltration depth (Trapezoidal integration for Horton)
        F_potential = I_t * dt

        # Step 4: Apply constraints (supply and capacity limits)
        F_actual = max(min(F_potential, P_i, SC), 0.0)
        infil_depth_arr[i] = F_actual
        infil_rate_arr[i] = F_actual / dt  # Effective rate after constraints

        # Update cumulative infiltration
        F_cumulative += F_actual
        F_cum_arr[i] = F_cumulative

        # Step 5: Update soil moisture
        theta_current = np.clip(theta_current + F_actual / soil.D, 0.0, soil.theta_s)
        theta_arr[i] = theta_current

        # Step 6: Surface runoff
        runoff_arr[i] = max((P_i - F_actual) / dt, 0.0)

    return {
        'time': time_arr, 'precip_intensity': precip_int_arr,
        'theta': theta_arr, 'infil_rate': infil_rate_arr,
        'infil_depth': infil_depth_arr, 'runoff_intensity': runoff_arr,
        'storage_capacity': sc_arr, 'F_cumulative': F_cum_arr
    }


def verify_water_balance(precip: np.ndarray, results: dict, dt: float) -> bool:
    """Check mass conservation: Total P = Total F + Total Q."""
    total_precip = np.sum(precip)
    total_infil = np.sum(results['infil_depth'])
    total_runoff = np.sum(results['runoff_intensity'] * dt)
    balance_error = abs(total_precip - total_infil - total_runoff)
    relative_error = balance_error / max(total_precip, 1e-10)

    print("=" * 60)
    print("WATER BALANCE VERIFICATION")
    print("=" * 60)
    print(f"  Total rainfall      : {total_precip:.4f} mm")
    print(f"  Total infiltration  : {total_infil:.4f} mm")
    print(f"  Total runoff        : {total_runoff:.4f} mm")
    print(f"  Balance error       : {balance_error:.6f} mm")
    print(f"  Status              : {'PASS' if relative_error < 1e-6 else 'FAIL'}")
    print("=" * 60)
    return relative_error < 1e-6


def verify_physical_constraints(results: dict, soil: SoilParameters) -> bool:
    """Check theta in [0, theta_s], infiltration >= 0, runoff >= 0."""
    all_pass = True
    print("\n" + "=" * 60)
    print("PHYSICAL CONSTRAINTS VERIFICATION")
    print("=" * 60)

    theta_min, theta_max = np.min(results['theta']), np.max(results['theta'])
    c1 = (theta_min >= -1e-10) and (theta_max <= soil.theta_s + 1e-10)
    print(f"  theta in [0, {soil.theta_s}] : [{theta_min:.6f}, {theta_max:.6f}] -> {'PASS' if c1 else 'FAIL'}")
    all_pass &= c1

    c2 = np.all(results['infil_depth'] >= -1e-10)
    print(f"  Infiltration >= 0   : {'PASS' if c2 else 'FAIL'}")
    all_pass &= c2

    c3 = np.all(results['runoff_intensity'] >= -1e-10)
    print(f"  Runoff >= 0         : {'PASS' if c3 else 'FAIL'}")
    all_pass &= c3

    c4 = np.all(results['storage_capacity'] >= -1e-10)
    print(f"  Storage cap >= 0    : {'PASS' if c4 else 'FAIL'}")
    all_pass &= c4

    print(f"  Overall             : {'ALL PASS' if all_pass else 'SOME FAILED'}")
    print("=" * 60)
    return all_pass
