import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from formatting import set_mpl_customisation
from sys import argv
from competition_regimes import plot_bounds

# Process arguments
if len(argv) != 3:
    print('Usage: {} CROSS_SECTION-CROSS_SECTION[-CROSS_SECTION-CROSS_SECTION] IMAGE_PNG'.format(argv[0]))
    exit(1)

cross_sections = argv[1].split('-')
image_png = argv[2]

eta_Bs = []
beta_Bs = []
for cross_section in cross_sections:
    if cross_section == 'I':
        eta_Bs.append(0.2)
        beta_Bs.append(0.2)
    elif cross_section == 'II':
        eta_Bs.append(0.2)
        beta_Bs.append(0.8)
    elif cross_section == 'IV':
        eta_Bs.append(0.175)
        beta_Bs.append(0.7)
    elif cross_section == 'V':
        eta_Bs.append(0.225)
        beta_Bs.append(0.9)
    elif cross_section == 'VI':
        #eta_Bs.append(0.24)
        eta_Bs.append(np.log(2))
        beta_Bs.append(0.7)
    elif cross_section == 'VII':
        #eta_Bs.append(0.09)
        eta_Bs.append(np.log(2))
        beta_Bs.append(0.9)
    else:
        raise ValueError('Cross Section {} not recognised'.format(cross_section))

# Customise matplotlib
set_mpl_customisation('combination')

# Create subplots
if len(cross_sections) == 2:
    fig, axes = plt.subplots(1, 2, figsize=(5.512, 2.6))
elif len(cross_sections) == 4:
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    axes = axes.flatten()
else:
    raise ValueError('Number of cross sections should be 2 or 4.  Given cross sections: {}'.format(
        cross_sections))

for cross_section, eta_B, beta_B, ax in zip(cross_sections, eta_Bs, beta_Bs, axes):

    plot_bounds(ax, eta_B, beta_B)

    if eta_B == np.log(2):
        eta_B = r'\ln(2)'

    ax.set_title(r'$\textrm{{Cross Section {}: }} \beta_B = {}, \eta_B = {}$'.format(cross_section,
        beta_B, eta_B))

fig.tight_layout()
fig.savefig(image_png)
