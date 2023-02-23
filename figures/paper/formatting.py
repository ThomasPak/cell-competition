import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def set_mpl_customisation(type='line'):

    # Plot customisation
    mpl.rcParams['lines.linewidth'] = 1
    mpl.rcParams['lines.markersize'] = 3
    mpl.rcParams['font.size'] = 7
    if type == 'line':
        mpl.rcParams['savefig.dpi'] = 1000
    elif type == 'combination':
        mpl.rcParams['savefig.dpi'] = 500
    else:
        mpl.rcParams['savefig.dpi'] = 300
    plt.rc('text', usetex=True)

def cbar_add_central_tick(cbar, central=None):

    # Get ticks
    ticks = cbar.get_ticks()

    # Define central
    if central == None:
        central = 0.5 * (ticks.min() + ticks.max())

    # Set central
    if not central in ticks:
        ticks = np.concatenate((ticks, [central]))
        cbar.set_ticks(ticks)
