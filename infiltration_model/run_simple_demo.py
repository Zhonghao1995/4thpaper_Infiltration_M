"""
============================================================================
Rainfall-Infiltration-Runoff Model — Simple 3-Hour Demo Case
============================================================================
A simplified demonstration based on specific user parameters:
    - Time step Dt = 1 hour
    - Soil depth D = 100 mm
    - Initial moisture theta_0 = 0.3
    - Saturated moisture theta_s = 0.45
    - Rainfall duration = 3 hours
    - Rainfall intensity P = [2, 12, 4] mm/h

Usage:
    python -m infiltration_model.run_simple_demo
============================================================================
"""

import numpy as np
from infiltration_model.core import (
    SoilParameters, HortonParams, GreenAmptParams,
    simulate, verify_water_balance, verify_physical_constraints
)
from infiltration_model.visualization import plot_results


def run_simple_demo():
    print("=" * 60)
    print("  SIMPLE 3-HOUR DEMONSTRATION CASE")
    print("=" * 60)

    # 1. Time & Rainfall
    dt = 1.0  # Time step = 1 hour
    precip_intensity = np.array([2.0, 12.0, 4.0])  # mm/h
    # Since dt=1, precip depth [mm] = precip intensity [mm/h]
    precip = precip_intensity * dt

    print(f"\nTime step (Dt) : {dt} hour")
    print(f"Rainfall [mm/h]: {precip_intensity}")
    print(f"Total rainfall : {np.sum(precip):.2f} mm")

    # 2. Soil Parameters
    # K_s is usually equal to Horton's f_c (5.0 mm/h)
    soil = SoilParameters(
        D=100.0,       # Effective depth = 100 mm
        theta_s=0.45,  # Saturated moisture
        theta_0=0.30,  # Initial moisture
        K_s=5.0        # Saturated hydraulic conductivity
    )

    print(f"\nSoil params    : D={soil.D}mm, theta_0={soil.theta_0}, theta_s={soil.theta_s}, K_s={soil.K_s}mm/h")

    # 3. Infiltration Parameters
    # Horton parameters provided by user
    horton_params = HortonParams(f_0=20.0, f_c=5.0, k=1.0)
    
    # Green-Ampt equivalent parameters
    # We use the same K_s=5.0. We need to choose a suction head (psi) that 
    # gives roughly similar behavior. 100-200mm is typical for medium soils.
    ga_params = GreenAmptParams(K_s=5.0, psi=150.0)

    print(f"Horton params  : f_0={horton_params.f_0}, f_c={horton_params.f_c}, k={horton_params.k}")
    print(f"GreenAmpt px   : K_s={ga_params.K_s}, psi={ga_params.psi}mm")

    # 4. Run Horton Simulation
    print("\n" + "-" * 40)
    print("Running HORTON Model...")
    results_horton = simulate(
        precip=precip, dt=dt, soil=soil,
        infiltration_model='horton', infil_params=horton_params
    )
    verify_water_balance(precip, results_horton, dt)

    # 5. Run Green-Ampt Simulation
    print("\n" + "-" * 40)
    print("Running GREEN-AMPT Model...")
    results_ga = simulate(
        precip=precip, dt=dt, soil=soil,
        infiltration_model='green_ampt', infil_params=ga_params
    )
    verify_water_balance(precip, results_ga, dt)

    # 6. Detailed Step-by-Step Output for Horton
    print("\n" + "=" * 60)
    print("STEP-BY-STEP RESULTS (HORTON)")
    print("=" * 60)
    print(f"{'Hour':>4} | {'Rain[mm]':>8} | {'I_cap[mm/h]':>11} | {'F_act[mm]':>9} | {'Q[mm/h]':>8} | {'theta':>6}")
    print("-" * 60)
    for i in range(len(precip)):
        print(f"{i+1:>4} | {precip[i]:>8.2f} | {results_horton['infil_rate'][i]:>11.2f} | "
              f"{results_horton['infil_depth'][i]:>9.2f} | {results_horton['runoff_intensity'][i]:>8.2f} | "
              f"{results_horton['theta'][i]:>6.4f}")

    # 7. Detailed Step-by-Step Output for Green-Ampt
    print("\n" + "=" * 60)
    print("STEP-BY-STEP RESULTS (GREEN-AMPT)")
    print("=" * 60)
    print(f"{'Hour':>4} | {'Rain[mm]':>8} | {'I_cap[mm/h]':>11} | {'F_act[mm]':>9} | {'Q[mm/h]':>8} | {'theta':>6}")
    print("-" * 60)
    for i in range(len(precip)):
        print(f"{i+1:>4} | {precip[i]:>8.2f} | {results_ga['infil_rate'][i]:>11.2f} | "
              f"{results_ga['infil_depth'][i]:>9.2f} | {results_ga['runoff_intensity'][i]:>8.2f} | "
              f"{results_ga['theta'][i]:>6.4f}")

    # 8. Generate Plot
    save_path = "simple_demo_output.png"
    print("\nGenerating comparison plot...")
    plot_results(results_horton, results_ga, dt, soil, save_path=save_path)
    print(f"Saved to: {save_path}")


if __name__ == "__main__":
    run_simple_demo()
