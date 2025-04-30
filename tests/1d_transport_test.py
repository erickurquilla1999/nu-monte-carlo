import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from matplotlib.ticker import AutoLocator, AutoMinorLocator, LogLocator
import os
import glob
import yt
print("yt version:", yt.__version__)

# Font settings
mpl.rcParams['font.size'] = 22
mpl.rcParams['font.family'] = 'serif'
mpl.rc('text', usetex=True)

# Tick settings
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.major.pad'] = 8
mpl.rcParams['xtick.minor.size'] = 4

mpl.rcParams['xtick.minor.width'] = 2
mpl.rcParams['ytick.major.size'] = 7
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.size'] = 4
mpl.rcParams['ytick.minor.width'] = 2

# Axis linewidth
mpl.rcParams['axes.linewidth'] = 2

# Tick direction and enabling ticks on all sides
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True

# Function to apply custom tick locators and other settings to an Axes object
def apply_custom_settings(ax, leg, log_scale_y=False):

    if log_scale_y:
        # Use LogLocator for the y-axis if it's in log scale
        ax.set_yscale('log')
        ax.yaxis.set_major_locator(LogLocator(base=10.0))
        ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs='auto', numticks=100))
    else:
        # Use AutoLocator for regular scales
        ax.yaxis.set_major_locator(AutoLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
    
    # Apply the AutoLocator for the x-axis
    ax.xaxis.set_major_locator(AutoLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    
    # Legend settings
    leg.get_frame().set_edgecolor('w')
    leg.get_frame().set_linewidth(0.0)

# Load AMReX plotfile
ds = yt.load("plt300")

# List available fields in the dataset
print("Available fields:")
for field in ds.field_list:
    print(field)

# Access the full grid at the finest level
cg = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                      dims=ds.domain_dimensions)

# Read field data
n_invcm3 = cg["n_invcm3"]
print("n_data shape:", n_invcm3.shape)
print("n_data min:", n_invcm3.min())
print("n_data max:", n_invcm3.max())

n_xy_average_invcm3 = np.mean(n_invcm3, axis=(1, 2))
print("n_xy_average shape:", n_xy_average_invcm3.shape)

n_xy_average_invcm3 = n_xy_average_invcm3 / n_xy_average_invcm3[0]

cell_positions_cm = np.linspace(0, 9, 10)
k=1/5
theoryvalues = np.exp(-k*(cell_positions_cm))

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(cell_positions_cm, theoryvalues, color='red', label=r'$e^{-\kappa x}$', linestyle='solid', linewidth=3)
ax.plot(cell_positions_cm, n_xy_average_invcm3, marker='o', label='nu-monte-carlo', linestyle='dashed', linewidth=3)
ax.set_xlabel(r'$x$ (cm)')
ax.set_ylabel(r'$\left<n\right>/\left<n_{x_0}\right>$')
# ax.set_yscale('log')
leg = ax.legend(framealpha=0.0, ncol=1, fontsize=25)
apply_custom_settings(ax, leg, False)
plt.savefig("1dtransporttest.pdf", bbox_inches='tight')
plt.show()
plt.close(fig)

error = np.sum(np.abs((n_xy_average_invcm3 - theoryvalues)/theoryvalues))
print("Percentage Error:", error)
assert error <= 0.5, f"Error is too high: {error:.2f}%"