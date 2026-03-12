"""
Demonstration case: 6-hour synthetic storm on silt loam soil.

Usage:
    python -m infiltration_model.run_demo
"""

import numpy as np
from infiltration_model.core import (
    SoilParameters, HortonParams, GreenAmptParams, RichardsParams,
    simulate, verify_water_balance, verify_physical_constraints
)
from infiltration_model.visualization import plot_results


def run_demo():
    """Run a 6-hour storm simulation with both Horton and Green-Ampt models."""

    print("=" * 60)
    print("  RAINFALL-INFILTRATION-RUNOFF MODEL DEMONSTRATION")
    print("=" * 60)

    # Time discretization
    dt = 1.0 / 12.0   # 5-minute steps [h]
    total_hours = 6.0
    n_steps = int(total_hours / dt)

    print(f"\nDt = {dt*60:.1f} min | Duration = {total_hours:.0f} h | Steps = {n_steps}")

    # Rainfall schedule: (duration [h], intensity [mm/h])
    rainfall_schedule = [
        (1.0,  5.0),   # Light onset
        (1.0, 15.0),   # Moderate
        (1.0, 40.0),   # Peak
        (1.0, 20.0),   # Declining
        (1.0,  8.0),   # Trailing
        (1.0,  2.0),   # End
    ]

    # Convert to depth per step [mm]
    precip = np.zeros(n_steps)
    idx = 0
    for dur, intensity in rainfall_schedule:
        for _ in range(int(dur / dt)):
            if idx < n_steps:
                precip[idx] = intensity * dt
                idx += 1

    print(f"Total rainfall = {np.sum(precip):.2f} mm")

    # Soil parameters (silt loam)
    soil = SoilParameters(D=500.0, theta_s=0.45, theta_0=0.20, K_s=6.5)

    # Horton parameters
    horton_params = HortonParams(f_0=25.0, f_c=5.0, k=1.5)

    # Green-Ampt parameters
    ga_params = GreenAmptParams(K_s=6.5, psi=150.0)

    # Richards parameters (Van Genuchten-Mualem for silt loam)
    richards_params = RichardsParams(theta_r=0.067, alpha=0.002, n_vg=1.41, dz=5.0)

    print(f"\nSoil: D={soil.D}mm, theta_s={soil.theta_s}, theta_0={soil.theta_0}, K_s={soil.K_s}mm/h")
    print(f"Horton: f_0={horton_params.f_0}, f_c={horton_params.f_c}, k={horton_params.k}")
    print(f"Green-Ampt: K_s={ga_params.K_s}, psi={ga_params.psi}mm")
    print(f"Richards: theta_r={richards_params.theta_r}, alpha={richards_params.alpha}, n_vg={richards_params.n_vg}, dz={richards_params.dz}mm")

    # Run Horton
    print("\n" + "=" * 60)
    print("Running Horton model...")
    print("=" * 60)
    results_horton = simulate(precip, dt, soil, 'horton', horton_params)
    wb_h = verify_water_balance(precip, results_horton, dt)
    pc_h = verify_physical_constraints(results_horton, soil)

    # Run Green-Ampt
    print("\n" + "=" * 60)
    print("Running Green-Ampt model...")
    print("=" * 60)
    results_ga = simulate(precip, dt, soil, 'green_ampt', ga_params)
    wb_g = verify_water_balance(precip, results_ga, dt)
    pc_g = verify_physical_constraints(results_ga, soil)

    # Run Richards
    print("\n" + "=" * 60)
    print("Running Richards equation (PDE solver)...")
    print("=" * 60)
    results_richards = simulate(precip, dt, soil, 'richards', richards_params)
    wb_r = verify_water_balance(precip, results_richards, dt)
    pc_r = verify_physical_constraints(results_richards, soil)

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY COMPARISON")
    print("=" * 60)
    print(f"{'Metric':<25} {'Horton':>12} {'Green-Ampt':>12} {'Richards':>12}")
    print("-" * 64)
    print(f"{'Total infiltration [mm]':<25} {np.sum(results_horton['infil_depth']):>12.2f} {np.sum(results_ga['infil_depth']):>12.2f} {np.sum(results_richards['infil_depth']):>12.2f}")
    print(f"{'Total runoff [mm]':<25} {np.sum(results_horton['runoff_intensity']*dt):>12.2f} {np.sum(results_ga['runoff_intensity']*dt):>12.2f} {np.sum(results_richards['runoff_intensity']*dt):>12.2f}")
    print(f"{'Peak runoff [mm/h]':<25} {np.max(results_horton['runoff_intensity']):>12.2f} {np.max(results_ga['runoff_intensity']):>12.2f} {np.max(results_richards['runoff_intensity']):>12.2f}")
    print(f"{'Final theta [m^3/m^3]':<25} {results_horton['theta'][-1]:>12.4f} {results_ga['theta'][-1]:>12.4f} {results_richards['theta'][-1]:>12.4f}")
    print(f"{'Water balance':<25} {'PASS' if wb_h else 'FAIL':>12} {'PASS' if wb_g else 'FAIL':>12} {'PASS' if wb_r else 'FAIL':>12}")
    print(f"{'Physical constraints':<25} {'PASS' if pc_h else 'FAIL':>12} {'PASS' if pc_g else 'FAIL':>12} {'PASS' if pc_r else 'FAIL':>12}")

    # Plot
    print("\nGenerating comparison plot...")
    plot_results(results_horton, results_ga, dt, soil, save_path="demo_output.png", results_richards=results_richards)

    print("\n" + "=" * 60)
    print("  DEMONSTRATION COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    run_demo()
