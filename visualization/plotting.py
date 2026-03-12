"""
Publication-quality 4-panel comparison plots.
Style: Arial 12pt, inward ticks, no grid lines.
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from infiltration_model.core.parameters import SoilParameters


def setup_plot_style():
    """Configure matplotlib for publication-quality output."""
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = ['Arial']
    mpl.rcParams['font.size'] = 12
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['xtick.bottom'] = True
    mpl.rcParams['ytick.left'] = True
    mpl.rcParams['ytick.right'] = True
    mpl.rcParams['xtick.major.size'] = 5
    mpl.rcParams['ytick.major.size'] = 5
    mpl.rcParams['axes.grid'] = False
    mpl.rcParams['axes.linewidth'] = 1.0


def plot_results(
    results_horton: dict,
    results_ga: dict,
    dt: float,
    soil: SoilParameters,
    save_path: str = "demo_output.png",
    results_richards: dict = None
) -> None:
    """Create 4-panel comparison: rainfall, infiltration, runoff, soil moisture."""
    setup_plot_style()
    time = results_horton['time']

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    title = 'Rainfall-Infiltration-Runoff Model: Horton vs Green-Ampt'
    if results_richards is not None:
        title += ' vs Richards'
    fig.suptitle(title, fontsize=14, fontweight='bold', y=0.98)

    # (a) Rainfall hyetograph
    ax1 = axes[0]
    ax1.bar(time - dt, results_horton['precip_intensity'], width=dt, align='edge',
            color='steelblue', alpha=0.7, edgecolor='navy', linewidth=0.5, label='Rainfall')
    ax1.set_ylabel('Rainfall Intensity\n[mm/h]')
    ax1.set_title('(a) Rainfall Hyetograph', fontsize=12, loc='left')
    max_hour = int(np.ceil(time[-1]))
    ax1.set_xticks(np.arange(0, max_hour + 1, 1))
    ax1.invert_yaxis()
    ax1.legend(loc='lower right', frameon=False)

    # (b) Infiltration rate
    ax2 = axes[1]
    ax2.plot(time, results_horton['infil_rate'], 'r-o', markersize=2.5, linewidth=1.2, label='Horton I(t)')
    ax2.plot(time, results_ga['infil_rate'], 'b-s', markersize=2.5, linewidth=1.2, label='Green-Ampt I(t)')
    if results_richards is not None:
        ax2.plot(time, results_richards['infil_rate'], 'g-^', markersize=2.5, linewidth=1.2, label='Richards I(t)')
    ax2.step(time - dt, results_horton['precip_intensity'], where='post',
             color='gray', linestyle='--', linewidth=0.8, alpha=0.6, label='Rainfall intensity')
    ax2.set_ylabel('Infiltration Rate\n[mm/h]')
    ax2.set_title('(b) Infiltration Rate Comparison', fontsize=12, loc='left')
    ax2.legend(loc='upper right', frameon=False)

    # (c) Soil moisture evolution
    ax3 = axes[2]
    ax3.plot(time, results_horton['theta'], 'r-o', markersize=2.5, linewidth=1.2, label=r'Horton $\theta$(t)')
    ax3.plot(time, results_ga['theta'], 'b-s', markersize=2.5, linewidth=1.2, label=r'Green-Ampt $\theta$(t)')
    if results_richards is not None:
        ax3.plot(time, results_richards['theta'], 'g-^', markersize=2.5, linewidth=1.2, label=r'Richards Avg $\theta$(t)')
    ax3.axhline(y=soil.theta_s, color='black', linestyle=':', linewidth=1, label=r'$\theta_s$ = ' + f'{soil.theta_s}')
    ax3.axhline(y=soil.theta_0, color='gray', linestyle=':', linewidth=1, label=r'$\theta_0$ = ' + f'{soil.theta_0}')
    ax3.set_ylabel(r'Soil Moisture $\theta$' + '\n' + r'[m$^3$/m$^3$]')
    ax3.set_xlabel('Time [h]')
    ax3.set_title(r'(c) Soil Moisture Evolution', fontsize=12, loc='left')
    ax3.legend(loc='lower right', frameon=False)

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\nPlot saved to: {save_path}")

def plot_moisture_profile(
    results_ga: dict,
    results_richards: dict,
    soil: SoilParameters,
    dt: float,
    target_time: float,
    save_path: str = "moisture_profile.png"
) -> None:
    """Create a depth vs moisture profile plot (z vs theta) matching the professor's diagram."""
    setup_plot_style()
    
    time_arr = results_richards['time']
    # Find the index closest to target time
    idx = np.argmin(np.abs(time_arr - target_time))
    actual_time = time_arr[idx]
    
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.set_title(f'Soil Moisture Profile at t = {actual_time:.1f} h', fontsize=14, fontweight='bold', pad=15)
    
    # Richards Profile (Continuous Curve)
    theta_z = results_richards['theta_z'][idx, :]
    z_nodes = results_richards['z_nodes']
    
    # We plot Z on the Y-axis inverted (0 at top, positive downwards)
    # and Theta on the X-axis
    ax.plot(theta_z, z_nodes, 'g-', linewidth=2.5, label='Richards Eqn')
    
    # Green-Ampt Profile (Step Function)
    # Z_f = F / Delta_theta
    F_ga = results_ga['F_cumulative'][idx]
    delta_theta = soil.theta_s - soil.theta_0
    Z_f = F_ga / delta_theta if delta_theta > 0 else 0.0
    
    # Create the rectangular shape for Green-Ampt
    ga_theta = [soil.theta_s, soil.theta_s, soil.theta_0, soil.theta_0]
    ga_z = [0.0, Z_f, Z_f, soil.D]
    ax.plot(ga_theta, ga_z, 'b--', linewidth=2.0, label='Green-Ampt (Step)')
    
    # Horton doesn't have a physical profile, so we omit it or just mention it.
    
    # Formatting to match the professor's diagram
    ax.set_ylim(soil.D, 0) # Invert Y axis
    ax.set_xlim(0, soil.theta_s + 0.05)
    
    ax.set_ylabel('Depth $z$ [mm]', fontsize=12)
    ax.set_xlabel(r'Soil Moisture $\theta$ [m$^3$/m$^3$]', fontsize=12)
    
    # Add vertical lines for theta_0 and theta_s
    ax.axvline(x=soil.theta_0, color='gray', linestyle=':', label=r'Initial $\theta_0$')
    ax.axvline(x=soil.theta_s, color='black', linestyle=':', label=r'Saturated $\theta_s$')
    
    # Fill under the curves to make it look like the diagram
    ax.fill_betweenx(z_nodes, soil.theta_0, theta_z, where=(theta_z > soil.theta_0), 
                     color='green', alpha=0.15, label='Infiltrated Water (Richards)')
    
    ax.legend(loc='lower right', frameon=True, facecolor='white', framealpha=0.9)
    plt.grid(True, linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Profile plot saved to: {save_path}")
