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
def apply_custom_settings(ax, log_scale_y=False):

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

# This function creates a pcolormesh plot with contour lines and saves it as a PNG file.
def plot_pcolormesh_with_contour(x, y, z, min_cb, max_cb, x_label, y_label, title, cbar_label, colormap, filename):

    fig, ax = plt.subplots(figsize=(12, 8))

    # Normalize Z data
    min_cb_scaled = min_cb
    max_cb_scaled = max_cb

    # Plot pcolormesh
    c = ax.pcolormesh(x, y, z, shading='auto', cmap=colormap, vmin=min_cb_scaled, vmax=max_cb_scaled)

    # Add contour lines
    # contour = ax.contour(x, y, z, colors='black', linewidths=0.5, levels=5)
    # ax.clabel(contour, inline=True, fontsize=15, fmt='%1.1f')

    # Plot settings
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)

    # Add color bar
    cbar = fig.colorbar(c, ax=ax, label=cbar_label)
    cbar.ax.yaxis.set_minor_locator(AutoMinorLocator())

    apply_custom_settings(ax, False)

    # Ensure equal aspect ratio
    ax.set_aspect('equal', 'box')

    # Save figure
    fig.savefig(filename, format='png', bbox_inches='tight')

    # Display figure
    # plt.show()
    # display(fig)
    
    # Close figure
    plt.close(fig)

# Load AMReX plotfile
ds = yt.load("pltmatter")
# ds = yt.load("plt0")

# List available fields in the dataset
print("Available fields:")
for field in ds.field_list:
    print(field)

# Access the full grid at the finest level
cg = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                      dims=ds.domain_dimensions)

# Read field data and convert to numpy arrays

rho_g_ccm = np.array(cg["rho_g_ccm"])
T_MeV = np.array(cg["T_MeV"])
chemical_potential_MeV = np.array(cg["chemical_potential_MeV"])
ye = np.array(cg["ye"])
IMFP_cm = np.array(cg["IMFP_cm"])

# Get physical coordinates of cell centers in the x-direction
x_centers = np.array(cg["x"])
y_centers = np.array(cg["y"])
z_centers = np.array(cg["z"])

# z slice
zslice = 15
zvalue = z_centers[0, 0, zslice]

# x and y coordinates
xplot = x_centers[:, :, zslice]
yplot = y_centers[:, :, zslice]

# TEMPERATURE

zplot = T_MeV[:, :, zslice]

plot_pcolormesh_with_contour(x=xplot, 
                            y=yplot, 
                            z=zplot, 
                            min_cb=None, 
                            max_cb=None, 
                            x_label=r'$x \, (\mathrm{km})$', 
                            y_label=r'$y \, (\mathrm{km})$', 
                            title=fr'$z = {zvalue:.2f} \, \mathrm{{km}}$',
                            cbar_label=fr'$T \, (\mathrm{{MeV}})$', 
                            colormap='viridis',
                            filename=f'temperature.png')

# IMFP

zplot = IMFP_cm[:, :, zslice]

plot_pcolormesh_with_contour(x=xplot, 
                            y=yplot, 
                            z=zplot, 
                            min_cb=None, 
                            max_cb=None, 
                            x_label=r'$x \, (\mathrm{km})$', 
                            y_label=r'$y \, (\mathrm{km})$', 
                            title=fr'$z = {zvalue:.2f} \, \mathrm{{km}}$',
                            cbar_label=fr'$\kappa \, (1/\mathrm{{cm}})$', 
                            colormap='viridis',
                            filename=f'imfp.png')

# CHEMICAL POTENTIAL

zplot = chemical_potential_MeV[:, :, zslice]

plot_pcolormesh_with_contour(x=xplot,
                            y=yplot, 
                            z=zplot, 
                            min_cb=None, 
                            max_cb=None, 
                            x_label=r'$x \, (\mathrm{km})$', 
                            y_label=r'$y \, (\mathrm{km})$', 
                            title=fr'$z = {zvalue:.2f} \, \mathrm{{km}}$',
                            cbar_label=fr'$\mu \, (\mathrm{{MeV}})$', 
                            colormap='viridis',
                            filename=f'chemical_potential.png')

# DENSITY

zplot = rho_g_ccm[:, :, zslice]

plot_pcolormesh_with_contour(x=xplot,
                            y=yplot, 
                            z=zplot, 
                            min_cb=None, 
                            max_cb=None, 
                            x_label=r'$x \, (\mathrm{km})$', 
                            y_label=r'$y \, (\mathrm{km})$', 
                            title=fr'$z = {zvalue:.2f} \, \mathrm{{km}}$',
                            cbar_label=fr'$\rho \, (\mathrm{{g/cm^3}})$', 
                            colormap='viridis',
                            filename=f'density.png')

# YE

zplot = ye[:, :, zslice]

plot_pcolormesh_with_contour(x=xplot,
                            y=yplot, 
                            z=zplot, 
                            min_cb=None, 
                            max_cb=None, 
                            x_label=r'$x \, (\mathrm{km})$', 
                            y_label=r'$y \, (\mathrm{km})$', 
                            title=fr'$z = {zvalue:.2f} \, \mathrm{{km}}$',
                            cbar_label=fr'$Y_e$', 
                            colormap='viridis',
                            filename=f'ye.png')

# NEUTRINO DATA

def plot_color_vector_map(x, y, z, vecx, vecy, ve, min_cb, max_cb, x_label, y_label, title, cbar_label, colormap, filename, doshow=True, dosave=True):

    fig, ax = plt.subplots(figsize=(12, 8))

    # Plot pcolormesh
    c = ax.pcolormesh(x, y, z, shading='auto', cmap=colormap, vmin=min_cb, vmax=max_cb)
    
    ax.quiver(x[::ve,::ve], y[::ve,::ve], ve*vecx[::ve,::ve], ve*vecy[::ve,::ve], color='white', scale=1, scale_units='xy', angles='xy')

    # Add contour lines
    # contour = ax.contour(x, y, z, colors='black', linewidths=0.5, levels=5)
    # ax.clabel(contour, inline=True, fontsize=10, fmt='%1.1e')

    # Plot settings
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)

    # Add color bar
    cbar = fig.colorbar(c, ax=ax, label=cbar_label)
    cbar.ax.yaxis.set_minor_locator(AutoMinorLocator())

    apply_custom_settings(ax, False)

    # Ensure equal aspect ratio
    ax.set_aspect('equal', 'box')

    # Save figure
    if dosave:
        fig.savefig(filename, format='png', bbox_inches='tight')

    # Display figure
    if doshow:
        plt.show()
    # display(fig)
    
    # Close figure
    plt.close(fig)

# Load AMReX plotfile
ds = yt.load("plt40")

# List available fields in the dataset
print("Available fields:")
for field in ds.field_list:
    print(field)

# Access the full grid at the finest level
cg = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                      dims=ds.domain_dimensions)

n_invcm3 = np.array(cg["n_invcm3"])
fx_invcm3 = np.array(cg["fx_invcm3"])
fy_invcm3 = np.array(cg["fy_invcm3"])
fz_invcm3 = np.array(cg["fz_invcm3"])

# Get physical coordinates of cell centers in the x-direction
x_centers = np.array(cg["x"])
y_centers = np.array(cg["y"])
z_centers = np.array(cg["z"])

# If you want directional resolution, you can use domain_width / domain_dimensions
domain_width = ds.domain_width
domain_dims = ds.domain_dimensions
cell_widths = domain_width / domain_dims
print("Cell widths per direction:", cell_widths.to("cm"))

# If your dataset is 1D or 2D, just access relevant axes:
if ds.dimensionality >= 1:
    print(f"dx = {cell_widths[0].to('cm')}")
if ds.dimensionality >= 2:
    print(f"dy = {cell_widths[1].to('cm')}")
if ds.dimensionality == 3:
    print(f"dz = {cell_widths[2].to('cm')}")

# time s
time = ds.current_time.to("s")
print(f"Time: {time:.2e}")

# z slice
zslice = 15
zvalue = z_centers[0, 0, zslice]

xplot = x_centers[:, :, zslice]
yplot = y_centers[:, :, zslice]

# PLOT NEUTRINO DENSITY AND FLUXES

zplot = n_invcm3[:, :, zslice]
fmag = np.sqrt(fx_invcm3[:, :, zslice]**2 + fy_invcm3[:, :, zslice]**2 + fz_invcm3[:, :, zslice]**2)
fxhat = np.where(fmag == 0, 0, fx_invcm3[:, :, zslice] / fmag) * cell_widths[0]
fyhat = np.where(fmag == 0, 0, fy_invcm3[:, :, zslice] / fmag) * cell_widths[1]
fzhat = np.where(fmag == 0, 0, fz_invcm3[:, :, zslice] / fmag) * cell_widths[2]

max_cb = np.max(zplot)
min_cb = np.min(zplot)

plot_color_vector_map(  x=xplot, 
                        y=yplot, 
                        z=zplot,
                        vecx=fxhat,
                        vecy=fyhat,
                        ve=2, 
                        min_cb=min_cb,
                        max_cb=max_cb*1.0e-1, 
                        x_label=r'$x \, (\mathrm{km})$', 
                        y_label=r'$y \, (\mathrm{km})$', 
                        title=rf'$z={zvalue:.1f}\,\mathrm{{km}}$, $t = {float(time):.2e}\,\mathrm{{s}}$', 
                        cbar_label=r'$n \, (\mathrm{cm^{-3}})$', 
                        colormap='viridis',
                        filename=f'neutrino_density_flux.png',
                        doshow=False,
                        dosave=True)