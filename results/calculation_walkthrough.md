# Detailed Hand-Verification: First 4 Timesteps (Horton vs Green-Ampt)

## 1. Initial Conditions & Parameters

### Time Discretization

| Parameter | Value | Description |
|-----------|-------|-------------|
| Δt | 1/12 h = 5 min = 0.08333 h | Timestep duration |
| Rainfall intensity (Hour 1) | 5.0 mm/h | Constant during first hour |
| P_i (depth per step) | 5.0 × 1/12 = **0.41667 mm** | P_depth = intensity × Δt |

### Soil Parameters (Silt Loam)

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Effective soil depth | D | 500.0 | mm |
| Saturated moisture content | θ_s | 0.45 | m³/m³ |
| Initial moisture content | θ_0 | 0.20 | m³/m³ |
| Saturated hydraulic conductivity | K_s | 6.5 | mm/h |

### Horton Model Parameters

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Initial infiltration capacity | f₀ | 25.0 | mm/h |
| Steady-state infiltration capacity | f_c | 5.0 | mm/h |
| Decay constant | k | 1.5 | 1/h |

### Green-Ampt Model Parameters

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Saturated hydraulic conductivity | K_s | 6.5 | mm/h |
| Wetting front suction head | ψ | 150.0 | mm |

### Initial State (shared by both models)

| Variable | Initial Value |
|----------|---------------|
| θ (current) | 0.2000 |
| F_cum (cumulative infiltration) | 0.0 mm |
| t_elapsed | 0.0 h |

---

## 2. Simulation Algorithm (per timestep)

```
For each step i:
  1. SC = D × (θ_s − θ)                    ← Remaining storage capacity
  2. Compute infiltration rate I(t)          ← From sub-model
  3. F_potential = I(t) × Δt                 ← Potential infiltration depth
  4. F_actual = max(min(F_pot, P_i, SC), 0)  ← Apply constraints
  5. θ_new = θ + F_actual / D               ← Update soil moisture
  6. Q = (P_i − F_actual) / Δt              ← Surface runoff
```

**Key constraint logic**: `F_actual = min(F_potential, P_i, SC)`
- If **P_i** is smallest → rainfall-limited (all rain infiltrates)
- If **F_potential** is smallest → capacity-limited (excess becomes runoff)
- If **SC** is smallest → storage-limited (soil nearly saturated)

---

## 3. Horton Model — Step-by-Step

### Governing Equation

$$I(t) = f_c + (f_0 - f_c) \cdot e^{-k \cdot t}$$

**Trapezoidal integration**: Average of instantaneous rates at timestep start and end:

$$I_{avg} = \frac{I(t_{start}) + I(t_{end})}{2}$$

---

### Step 1: t = 0.0000 → 0.0833 h

**① Storage capacity**
```
SC = D × (θ_s − θ) = 500 × (0.45 − 0.2000) = 125.0000 mm
```

**② Infiltration rate (Trapezoidal Rule)**
```
I(t_start = 0.0) = 5.0 + (25.0 − 5.0) × exp(−1.5 × 0.0)
                  = 5.0 + 20.0 × 1.000000
                  = 25.0000 mm/h

I(t_end = 0.0833) = 5.0 + 20.0 × exp(−1.5 × 0.0833)
                   = 5.0 + 20.0 × 0.882497
                   = 22.6499 mm/h

I_avg = (25.0000 + 22.6499) / 2 = 23.8250 mm/h
```

**③ Potential infiltration depth**
```
F_potential = 23.8250 × 0.08333 = 1.9854 mm
```

**④ Constraint check** ⭐
```
min(F_pot=1.9854, P_i=0.4167, SC=125.0000) = 0.4167
→ F_actual = 0.4167 mm  ← [RAINFALL-LIMITED: P_i is smallest]
→ I_effective = 0.4167 / 0.08333 = 5.0000 mm/h (equals rainfall intensity)
```

**⑤⑥ State update**
```
F_cum = 0 + 0.4167 = 0.4167 mm
θ = 0.2000 + 0.4167/500 = 0.200833
Q = (0.4167 − 0.4167) / 0.08333 = 0.0000 mm/h  ← No runoff
```

---

### Step 2: t = 0.0833 → 0.1667 h

**① SC**
```
SC = 500 × (0.45 − 0.200833) = 124.5833 mm
```

**② Infiltration rate**
```
I(0.0833) = 5.0 + 20.0 × exp(−1.5 × 0.0833) = 5.0 + 20.0 × 0.882497 = 22.6499 mm/h
I(0.1667) = 5.0 + 20.0 × exp(−1.5 × 0.1667) = 5.0 + 20.0 × 0.778801 = 20.5760 mm/h
I_avg = (22.6499 + 20.5760) / 2 = 21.6130 mm/h
```

**③ F_potential**
```
F_potential = 21.6130 × 0.08333 = 1.8011 mm
```

**④ Constraint check** ⭐
```
min(F_pot=1.8011, P_i=0.4167, SC=124.5833) = 0.4167
→ F_actual = 0.4167 mm ← [RAINFALL-LIMITED]
→ I_effective = 5.0000 mm/h
```

**⑤⑥ State update**
```
F_cum = 0.4167 + 0.4167 = 0.8333 mm
θ = 0.200833 + 0.4167/500 = 0.201667
Q = 0.0000 mm/h
```

---

### Step 3: t = 0.1667 → 0.2500 h

**① SC**
```
SC = 500 × (0.45 − 0.201667) = 124.1667 mm
```

**② Infiltration rate**
```
I(0.1667) = 5.0 + 20.0 × 0.778801 = 20.5760 mm/h
I(0.2500) = 5.0 + 20.0 × 0.687289 = 18.7458 mm/h
I_avg = (20.5760 + 18.7458) / 2 = 19.6609 mm/h
```

**③④ Constraint check**
```
F_potential = 19.6609 × 0.08333 = 1.6384 mm
min(1.6384, 0.4167, 124.1667) = 0.4167 ← [RAINFALL-LIMITED]
I_effective = 5.0000 mm/h
```

**⑤⑥ State update**
```
F_cum = 0.8333 + 0.4167 = 1.2500 mm
θ = 0.201667 + 0.4167/500 = 0.202500
Q = 0.0000 mm/h
```

---

### Step 4: t = 0.2500 → 0.3333 h

**① SC**
```
SC = 500 × (0.45 − 0.202500) = 123.7500 mm
```

**② Infiltration rate**
```
I(0.2500) = 5.0 + 20.0 × 0.687289 = 18.7458 mm/h
I(0.3333) = 5.0 + 20.0 × 0.606531 = 17.1306 mm/h
I_avg = (18.7458 + 17.1306) / 2 = 17.9382 mm/h
```

**③④ Constraint check**
```
F_potential = 17.9382 × 0.08333 = 1.4949 mm
min(1.4949, 0.4167, 123.7500) = 0.4167 ← [RAINFALL-LIMITED]
I_effective = 5.0000 mm/h
```

**⑤⑥ State update**
```
F_cum = 1.2500 + 0.4167 = 1.6667 mm
θ = 0.202500 + 0.4167/500 = 0.203333
Q = 0.0000 mm/h
```

---

## 4. Green-Ampt Model — Step-by-Step

### Governing Equation

$$I(t) = K_s \cdot \left(1 + \frac{\psi \cdot \Delta\theta}{F_{cum}}\right)$$

where $\Delta\theta = \theta_s - \theta$

> [!NOTE]
> When F_cum = 0 (first step), the code uses `F_used = max(F_cum, 0.1)` to prevent division by zero.

---

### Step 1: t = 0.0000 → 0.0833 h

**① SC**
```
SC = 500 × (0.45 − 0.2000) = 125.0000 mm
```

**② Infiltration rate**
```
Δθ = 0.45 − 0.2000 = 0.2500
F_cum = 0.0 → F_used = max(0.0, 0.1) = 0.1  ← [Division-by-zero guard]

I(t) = 6.5 × (1 + 150 × 0.25 / 0.1)
     = 6.5 × (1 + 37.5 / 0.1)
     = 6.5 × (1 + 375.0)
     = 6.5 × 376.0
     = 2444.0000 mm/h  ← Very large (because F_cum ≈ 0)
```

**③ F_potential**
```
F_potential = 2444.0000 × 0.08333 = 203.6667 mm
```

**④ Constraint check** ⭐
```
min(F_pot=203.6667, P_i=0.4167, SC=125.0000) = 0.4167
→ F_actual = 0.4167 mm ← [RAINFALL-LIMITED]
→ I_effective = 5.0000 mm/h
```

**⑤⑥ State update**
```
F_cum = 0 + 0.4167 = 0.4167 mm
θ = 0.2000 + 0.4167/500 = 0.200833
Q = 0.0000 mm/h
```

---

### Step 2: t = 0.0833 → 0.1667 h

**① SC**
```
SC = 500 × (0.45 − 0.200833) = 124.5833 mm
```

**② Infiltration rate**
```
Δθ = 0.45 − 0.200833 = 0.249167
F_cum = 0.4167 → F_used = 0.4167  ← F_cum > 0.1, use actual value

I(t) = 6.5 × (1 + 150 × 0.249167 / 0.4167)
     = 6.5 × (1 + 37.375 / 0.4167)
     = 6.5 × (1 + 89.700)
     = 6.5 × 90.700
     = 589.5500 mm/h
```

**③④ Constraint check**
```
F_potential = 589.5500 × 0.08333 = 49.1292 mm
min(49.1292, 0.4167, 124.5833) = 0.4167 ← [RAINFALL-LIMITED]
I_effective = 5.0000 mm/h
```

**⑤⑥ State update**
```
F_cum = 0.4167 + 0.4167 = 0.8333 mm
θ = 0.200833 + 0.4167/500 = 0.201667
Q = 0.0000 mm/h
```

---

### Step 3: t = 0.1667 → 0.2500 h

**① SC**
```
SC = 500 × (0.45 − 0.201667) = 124.1667 mm
```

**② Infiltration rate**
```
Δθ = 0.45 − 0.201667 = 0.248333
F_used = 0.8333

I(t) = 6.5 × (1 + 150 × 0.248333 / 0.8333)
     = 6.5 × (1 + 37.250 / 0.8333)
     = 6.5 × 45.700
     = 297.0500 mm/h
```

**③④ Constraint check**
```
F_potential = 297.05 × 0.08333 = 24.7542 mm
min(24.7542, 0.4167, 124.1667) = 0.4167 ← [RAINFALL-LIMITED]
I_effective = 5.0000 mm/h
```

**⑤⑥ State update**
```
F_cum = 0.8333 + 0.4167 = 1.2500 mm
θ = 0.201667 + 0.4167/500 = 0.202500
Q = 0.0000 mm/h
```

---

### Step 4: t = 0.2500 → 0.3333 h

**① SC**
```
SC = 500 × (0.45 − 0.202500) = 123.7500 mm
```

**② Infiltration rate**
```
Δθ = 0.45 − 0.202500 = 0.247500
F_used = 1.2500

I(t) = 6.5 × (1 + 150 × 0.2475 / 1.25)
     = 6.5 × (1 + 37.125 / 1.25)
     = 6.5 × 30.700
     = 199.5500 mm/h
```

**③④ Constraint check**
```
F_potential = 199.55 × 0.08333 = 16.6292 mm
min(16.6292, 0.4167, 123.7500) = 0.4167 ← [RAINFALL-LIMITED]
I_effective = 5.0000 mm/h
```

**⑤⑥ State update**
```
F_cum = 1.2500 + 0.4167 = 1.6667 mm
θ = 0.202500 + 0.4167/500 = 0.203333
Q = 0.0000 mm/h
```

---

## 5. Summary Comparison Tables

### Horton Model

| Step | Time [h] | I_start | I_end | I_avg [mm/h] | F_pot [mm] | P_i [mm] | SC [mm] | F_actual | Limiting Factor | θ | Q [mm/h] |
|------|----------|---------|-------|-------------|------------|----------|---------|----------|----------------|-------|----------|
| 1 | 0→0.0833 | 25.0000 | 22.6499 | 23.8250 | 1.9854 | 0.4167 | 125.0000 | 0.4167 | P_i | 0.200833 | 0 |
| 2 | 0.0833→0.1667 | 22.6499 | 20.5760 | 21.6130 | 1.8011 | 0.4167 | 124.5833 | 0.4167 | P_i | 0.201667 | 0 |
| 3 | 0.1667→0.2500 | 20.5760 | 18.7458 | 19.6609 | 1.6384 | 0.4167 | 124.1667 | 0.4167 | P_i | 0.202500 | 0 |
| 4 | 0.2500→0.3333 | 18.7458 | 17.1306 | 17.9382 | 1.4949 | 0.4167 | 123.7500 | 0.4167 | P_i | 0.203333 | 0 |

### Green-Ampt Model

| Step | Time [h] | Δθ | F_used [mm] | I(t) [mm/h] | F_pot [mm] | P_i [mm] | SC [mm] | F_actual | Limiting Factor | θ | Q [mm/h] |
|------|----------|----|------------|-------------|------------|----------|---------|----------|----------------|-------|----------|
| 1 | 0→0.0833 | 0.2500 | 0.1* | 2444.00 | 203.67 | 0.4167 | 125.00 | 0.4167 | P_i | 0.200833 | 0 |
| 2 | 0.0833→0.1667 | 0.2492 | 0.4167 | 589.55 | 49.13 | 0.4167 | 124.58 | 0.4167 | P_i | 0.201667 | 0 |
| 3 | 0.1667→0.2500 | 0.2483 | 0.8333 | 297.05 | 24.75 | 0.4167 | 124.17 | 0.4167 | P_i | 0.202500 | 0 |
| 4 | 0.2500→0.3333 | 0.2475 | 1.2500 | 199.55 | 16.63 | 0.4167 | 123.75 | 0.4167 | P_i | 0.203333 | 0 |

*\* F_used = max(F_cum=0, 0.1) = 0.1, division-by-zero guard*

> [!IMPORTANT]
> **Key observation for all 4 steps**: Both models have infiltration capacity (F_pot) far exceeding the rainfall supply (P_i = 0.4167 mm), so all rainfall infiltrates with zero surface runoff. Both models produce identical results in these first 4 steps. The models diverge only when rainfall intensity increases beyond infiltration capacity (Hour 2: 15 mm/h, Hour 3: 40 mm/h).

---

## 6. Source Code Reference

| Algorithm Step | Code Location |
|----------------|---------------|
| Horton equation | [horton.py](file:///e:/antigravity/2.26%204thpaper/infiltration_model/models/horton.py) `horton_infiltration_rate()` |
| Green-Ampt equation | [green_ampt.py](file:///e:/antigravity/2.26%204thpaper/infiltration_model/models/green_ampt.py) `green_ampt_infiltration_rate()` |
| Simulation loop (6-step algorithm) | [simulation.py](file:///e:/antigravity/2.26%204thpaper/infiltration_model/core/simulation.py) `simulate()` L73–L108 |
