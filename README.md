# Point-Scale Infiltration Model: Horton vs Green-Ampt Comparison

## Overview
A point-scale rainfall-infiltration-runoff simulation comparing Horton's empirical model and the Green-Ampt physically-based model under identical boundary conditions, for quantifying model form uncertainty.

## Model Assumptions
- **Point-scale (1D column)**: No spatial routing, catchment area, or overland flow
- **Single-layer soil**: Uniform soil column tracked by bulk moisture content θ
- **Forward simulation**: Explicit time stepping with trapezoidal integration (Horton)
- **Constraint-based**: F_actual = min(F_potential, P_supply, Storage_capacity)

## Governing Equations

### Horton (Empirical)
```
I(t) = f_c + (f_0 - f_c) * exp(-k * t)
```
Integration: Trapezoidal rule — I_avg = (I(t_start) + I(t_end)) / 2

### Green-Ampt (Physically-based)
```
I(t) = K_s * (1 + ψ * Δθ / F_cum)
```
where Δθ = θ_s - θ (moisture deficit)

## File Structure
```
submission_to_professor/
├── README.md                  ← This file
├── run_demo.py                ← Main entry point (6-hour storm simulation)
├── core/
│   ├── parameters.py          ← Soil & model parameter definitions
│   └── simulation.py          ← Core simulation engine + verification
├── models/
│   ├── horton.py              ← Horton infiltration rate function
│   └── green_ampt.py          ← Green-Ampt infiltration rate function
├── visualization/
│   └── plotting.py            ← 4-panel comparison plot generator
└── results/
    ├── demo_output.png        ← 4-panel comparison plot (6h storm)
    └── Step_by_Step_Results_Demo.xlsx  ← Detailed timestep results
```

## Demo Parameters

| Parameter | Value | Unit |
|-----------|-------|------|
| Δt | 5 min | - |
| Duration | 6 hours | - |
| Soil type | Silt loam | - |
| D (soil depth) | 500 | mm |
| θ_s (saturated) | 0.45 | m³/m³ |
| θ_0 (initial) | 0.20 | m³/m³ |
| K_s | 6.5 | mm/h |
| Horton f₀ | 25.0 | mm/h |
| Horton f_c | 5.0 | mm/h |
| Horton k | 1.5 | 1/h |
| Green-Ampt ψ | 150.0 | mm |

### Rainfall Schedule (6 hours)
| Hour | Intensity (mm/h) |
|------|------------------|
| 1 | 5.0 (light) |
| 2 | 15.0 (moderate) |
| 3 | 40.0 (peak) |
| 4 | 20.0 (declining) |
| 5 | 8.0 (trailing) |
| 6 | 2.0 (end) |

## Key Results

| Metric | Horton | Green-Ampt |
|--------|--------|------------|
| Total infiltration [mm] | 29.97 | 54.72 |
| Total runoff [mm] | 60.03 | 35.28 |
| Peak runoff [mm/h] | 34.76 | 28.04 |
| Final θ [m³/m³] | 0.2599 | 0.3094 |
| Water balance | PASS | PASS |
| Physical constraints | PASS | PASS |

## How to Run
```bash
# From the parent directory containing infiltration_model/
python -m infiltration_model.run_demo
```

> **Note**: This code package is extracted for review purposes. The runnable version requires the full `infiltration_model` package structure with `__init__.py` files.
