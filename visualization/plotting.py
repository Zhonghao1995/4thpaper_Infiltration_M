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
