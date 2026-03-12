"""
Richards equation 1D solver for variably saturated vertical flow.

PDE (z positive downward):
    dtheta/dt = d/dz [K(h)(dh/dz - 1)]

Constitutive relations: Van Genuchten-Mualem model.

Boundary conditions:
    Top  (z=0): Switching flux BC (rainfall) / head BC (ponding, h=0)
    Bottom (z=D): Free drainage (dh/dz = 0, gravity-only flux)

Numerical method: Implicit (backward Euler) finite differences with
    Picard iteration for nonlinearity.

References:
    Richards, L.A. (1931). Journal of Applied Physics
    van Genuchten, M.Th. (1980). Soil Science Society of America Journal
"""

import numpy as np
from infiltration_model.core.parameters import SoilParameters, RichardsParams


# ── Van Genuchten constitutive relations (vectorised) ──────────────

def _vg_Se(h, alpha, n_vg, m):
    """Effective saturation Se(h)."""
    return np.where(h >= 0.0, 1.0,
                    (1.0 + (alpha * np.abs(h)) ** n_vg) ** (-m))


def _vg_theta(h, alpha, n_vg, m, theta_r, theta_s):
    """Volumetric water content theta(h)."""
    return theta_r + (theta_s - theta_r) * _vg_Se(h, alpha, n_vg, m)


def _vg_K(h, alpha, n_vg, m, l_param, K_s):
    """Unsaturated hydraulic conductivity K(h)."""
    Se = _vg_Se(h, alpha, n_vg, m)
    Se_safe = np.clip(Se, 1e-15, 1.0 - 1e-15)
    term = (1.0 - (1.0 - Se_safe ** (1.0 / m)) ** m) ** 2
    return np.where(Se >= 1.0, K_s, K_s * Se ** l_param * term)


def _vg_C(h, alpha, n_vg, m, theta_r, theta_s):
    """Specific moisture capacity C(h) = dtheta/dh."""
    abs_h = np.maximum(np.abs(h), 1e-10)
    ah_n = (alpha * abs_h) ** n_vg
    numer = alpha * m * n_vg * (theta_s - theta_r) * (alpha * abs_h) ** (n_vg - 1.0)
    denom = (1.0 + ah_n) ** (m + 1.0)
    return np.where(h >= 0.0, 1e-10, numer / denom)


def _h_from_theta(theta, alpha, n_vg, m, theta_r, theta_s):
    """Invert VG curve to obtain h from theta (scalar)."""
    if theta >= theta_s:
        return 0.0
    Se = np.clip((theta - theta_r) / (theta_s - theta_r), 1e-10, 1.0 - 1e-10)
    return -(1.0 / alpha) * (Se ** (-1.0 / m) - 1.0) ** (1.0 / n_vg)


# ── Thomas algorithm (tridiagonal solver) ──────────────────────────

def _thomas(a, d, c, b):
    """Solve a*x_{i-1} + d*x_i + c*x_{i+1} = b  via Thomas algorithm."""
    n = len(d)
    cp = np.zeros(n)
    dp = np.zeros(n)

    cp[0] = c[0] / d[0]
    dp[0] = b[0] / d[0]
    for i in range(1, n):
        w = d[i] - a[i] * cp[i - 1]
        cp[i] = c[i] / w if i < n - 1 else 0.0
        dp[i] = (b[i] - a[i] * dp[i - 1]) / w

    x = np.zeros(n)
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


# ── Single time-step Picard solver ─────────────────────────────────

def _picard_step(h_old, P_rate, dt, nz, dz,
                 alpha, n_vg, m, theta_r, theta_s, l_param, K_s,
                 max_iter, tol, ponding):
    """Solve one implicit time step with Picard linearisation."""
    h = h_old.copy()

    for _ in range(max_iter):
        h_prev = h.copy()

        # Constitutive coefficients at current iterate
        K_n = _vg_K(h, alpha, n_vg, m, l_param, K_s)
        C_n = np.maximum(_vg_C(h, alpha, n_vg, m, theta_r, theta_s), 1e-10)

        # Interface K (harmonic mean between adjacent nodes)
        K_sum = K_n[:-1] + K_n[1:]
        K_half = np.where(K_sum > 0, 2.0 * K_n[:-1] * K_n[1:] / K_sum, 0.0)

        # Tridiagonal coefficients
        a_vec = np.zeros(nz)
        d_vec = np.zeros(nz)
        c_vec = np.zeros(nz)
        rhs   = np.zeros(nz)

        # Interior nodes j = 1 .. nz-2
        if nz > 2:
            kl = K_half[0:nz - 2]          # interface left of j
            kr = K_half[1:nz - 1]          # interface right of j
            a_vec[1:nz - 1] = -kl / dz
            c_vec[1:nz - 1] = -kr / dz
            d_vec[1:nz - 1] = C_n[1:nz - 1] * dz / dt + (kl + kr) / dz
            rhs[1:nz - 1]   = C_n[1:nz - 1] * dz / dt * h_old[1:nz - 1] + (kl - kr)

        # Top boundary (j=0)
        if ponding:
            # Dirichlet: h[0] = 0 (surface ponded at atmospheric pressure)
            d_vec[0] = 1.0
            c_vec[0] = 0.0
            rhs[0]   = 0.0
        else:
            # Neumann: prescribed flux = rainfall rate
            kr0 = K_half[0] if nz > 1 else K_n[0]
            d_vec[0] = C_n[0] * dz / dt + kr0 / dz
            c_vec[0] = -kr0 / dz
            rhs[0]   = C_n[0] * dz / dt * h_old[0] + P_rate - kr0

        # Bottom boundary (j=nz-1): free drainage  (dh/dz = 0 → q = K)
        kl_bot = K_half[-1] if nz > 1 else K_n[0]
        K_bot  = K_n[-1]
        a_vec[-1] = -kl_bot / dz
        d_vec[-1] = C_n[-1] * dz / dt + kl_bot / dz
        rhs[-1]   = C_n[-1] * dz / dt * h_old[-1] + kl_bot - K_bot

        # Solve
        h = _thomas(a_vec, d_vec, c_vec, rhs)

        # Convergence check
        if np.max(np.abs(h - h_prev)) < tol:
            break

    return h


# ── Public API ─────────────────────────────────────────────────────

def solve_richards(precip, dt, soil, params):
    """
    Solve 1D Richards equation for a rainfall event.

    Returns a dict with the same keys as Horton / Green-Ampt results
    so that all three models can be compared directly.

    Parameters
    ----------
    precip : np.ndarray  Rainfall depth per time step [mm]
    dt     : float       Time step [h]
    soil   : SoilParameters
    params : RichardsParams

    Returns
    -------
    dict  Keys: time, precip_intensity, theta, infil_rate,
          infil_depth, runoff_intensity, storage_capacity, F_cumulative
    """
    n_steps = len(precip)
    nz = int(round(soil.D / params.dz))
    dz = params.dz
    K_s, theta_s = soil.K_s, soil.theta_s
    theta_r = params.theta_r
    alpha, n_vg = params.alpha, params.n_vg
    m = 1.0 - 1.0 / n_vg
    l_param = params.l

    # Uniform initial pressure-head profile
    h_init = _h_from_theta(soil.theta_0, alpha, n_vg, m, theta_r, theta_s)
    h = np.full(nz, h_init, dtype=np.float64)

    # Pre-allocate output
    time_arr      = np.zeros(n_steps)
    precip_int    = np.zeros(n_steps)
    theta_arr     = np.zeros(n_steps)
    infil_rate    = np.zeros(n_steps)
    infil_depth   = np.zeros(n_steps)
    runoff_arr    = np.zeros(n_steps)
    sc_arr        = np.zeros(n_steps)
    F_cum_arr     = np.zeros(n_steps)

    F_cum = 0.0
    t = 0.0

    for step in range(n_steps):
        t += dt
        P_rate = precip[step] / dt          # rainfall intensity [mm/h]
        h_old = h.copy()
        theta_old = _vg_theta(h_old, alpha, n_vg, m, theta_r, theta_s)

        # --- Pass 1: try flux BC (rainfall enters soil) ---
        h_new = _picard_step(h_old, P_rate, dt, nz, dz,
                             alpha, n_vg, m, theta_r, theta_s, l_param, K_s,
                             params.max_iter, params.tol, ponding=False)

        ponding = False
        if h_new[0] > 0.0:
            # Surface pressure > 0 → ponding → re-solve with h[0]=0
            ponding = True
            h_new = _picard_step(h_old, P_rate, dt, nz, dz,
                                 alpha, n_vg, m, theta_r, theta_s, l_param, K_s,
                                 params.max_iter, params.tol, ponding=True)

        h = h_new

        h = h_new

        h = h_new

        # --- Compute Darcy Flux at the surface (Infiltration Rate) ---
        # Evaluate unsaturated hydraulic conductivity at the top interface
        K_top = _vg_K(h[0:2], alpha, n_vg, m, l_param, K_s)
        kr0 = 2.0 * K_top[0] * K_top[1] / (K_top[0] + K_top[1]) if (K_top[0] + K_top[1]) > 0 else 0.0

        if ponding:
            # Under ponding (h[0] = 0), infiltration is limited by soil capacity
            # Darcy flux: q = -K * (dh/dz - 1)
            # Inward infiltration rate I = -q
            I_rate = kr0 * ((h[0] - h[1]) / dz + 1.0)
            
            # Since we assume zero pond depth tracking, if rainfall > infiltration capacity, 
            # the excess is instantly runoff. However, if rainfall < capacity but we remained 
            # in ponding state (numerical overshoot), we limit I_rate to max capacity or P_rate.
            # To be strictly physical without tracking pond head, I_rate cannot exceed what the
            # soil can physically suck in, which is exactly the Darcy flux computed above.
            I_rate = max(I_rate, 0.0)

            # Important: if the computed Darcy flux > P_rate, and there is no ponded water 
            # available from previous steps (since we don't track it), the actual infiltration 
            # cannot exceed the rainfall rate. This prevents the "spikes" when rainfall drops.
            I_rate = min(I_rate, P_rate)

        else:
            # If not ponded, all rainfall infiltrates exactly (Neumann boundary was met)
            I_rate = P_rate

        # --- Compute Runoff ---
        # Runoff is exactly the rainfall that couldn't infiltrate
        Q_rate = max(P_rate - I_rate, 0.0)

        # Total infiltration depth for this step
        F_i = I_rate * dt
        
        # --- Update Storage and Balance ---
        theta_new = _vg_theta(h, alpha, n_vg, m, theta_r, theta_s)
        
        # Storage capacity
        theta_avg = float(np.mean(theta_new))
        SC = max(soil.D * (theta_s - theta_avg), 0.0)

        F_cum += F_i

        # Store
        time_arr[step]    = t
        precip_int[step]  = P_rate
        theta_arr[step]   = theta_avg
        infil_rate[step]  = I_rate
        infil_depth[step] = F_i
        runoff_arr[step]  = Q_rate
        sc_arr[step]      = SC
        F_cum_arr[step]   = F_cum

    return {
        'time': time_arr, 'precip_intensity': precip_int,
        'theta': theta_arr, 'infil_rate': infil_rate,
        'infil_depth': infil_depth, 'runoff_intensity': runoff_arr,
        'storage_capacity': sc_arr, 'F_cumulative': F_cum_arr,
    }
