# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 10:36:44 2019

@author: Greg Smith

functions to work with CXRO atomic scattering files
"""

#%% load things
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys

# use fancy LaTeX plotting styles
plt.style.use(['science-fixed', 'high-vis'])
plt.rcParams['figure.dpi'] = 240  # fix high-dpi display scaling issues
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

# define useful file paths
script_dir = os.path.dirname(os.path.realpath(__file__))  # this script file
parent_dir = os.path.split(script_dir)[0]
funcs_dir = os.path.join(parent_dir, 'funcs')
CrossSection_dir = os.path.join(parent_dir, 'CrossSections')
index_dir = os.path.join(parent_dir, 'index')
trans_dir = os.path.join(parent_dir, 'Transmissions')
sys.path.append(funcs_dir)  # add funcs's parent directory to path

# load my modules
import CXRO_funcs as CXRO  # CXRO convenience functions

# parameters to save files
save_fig = True  # flag to save figs or not
diss_dir = r'C:\PhDPythonScripts\CXRO\figures'

#%% load data
os.chdir(CrossSection_dir)

myP = 760  # Torr
myL = 1e7  # nm
Ar = CXRO.load_CXRO_ASF('ar.nff.txt', P_Torr=myP, L_nm=myL)
N = CXRO.load_CXRO_ASF('n.nff.txt', P_Torr=myP, L_nm=myL)
O = CXRO.load_CXRO_ASF('o.nff.txt', P_Torr=myP, L_nm=myL)
Ne = CXRO.load_CXRO_ASF('ne.nff.txt', P_Torr=myP, L_nm=myL)
He = CXRO.load_CXRO_ASF('he.nff.txt', P_Torr=myP, L_nm=myL)

title_text = f'Pressure = {myP} Torr, Path = {myL*1e-7: 2.2f} cm'
df_list = [O, N]
CXRO.plot_trans(df_list, title_text, logy=True)
CXRO.plot_mu(df_list, title_text, logy=True)
CXRO.plot_index(Ar, title_text, logy=True)

#%% air mixture
os.chdir(CrossSection_dir)
myP = 760  # Torr
myL = 5e5  # nm
air = CXRO.load_air(myP, myL)

fig, ax = plt.subplots()
lns1 = ax.semilogy(air['E'], air['delta'], label=r'$\delta$')
lns2 = ax.semilogy(air['E'], air['beta'], label=r'$\beta$')
ax.set_xlabel('Energy [eV]')
ax.set_ylabel(r'$\delta, \beta$')

ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
lns3 = ax2.semilogy(air['E'], air['T'], label='T', c='g', ls='-.')
ax2.set_ylabel('T')
lns = lns1+lns2+lns3
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, frameon=True, loc=(0.75, 0.5))

ax.set_title(f'Air at STP (d ={myL*1e-6: 0.1f} mm)')
ax.text(0.3, 0.7, r'$n = 1-\delta - i \beta$', transform=ax.transAxes)
if save_fig:
    plt.savefig(os.path.join(diss_dir, 'AirAbs.pdf'), dpi=300)

#%% XUV transmission vs pressure

os.chdir(CrossSection_dir)

myL = 4*(1e9)  # nm
myP = 1e-7  # Torr

Ar = CXRO.load_CXRO_ASF('ar.nff.txt', P_Torr=myP, L_nm=myL)
He = CXRO.load_CXRO_ASF('he.nff.txt', P_Torr=myP, L_nm=myL)
N = CXRO.load_CXRO_ASF('n.nff.txt', P_Torr=2*myP, L_nm=myL)
H = CXRO.load_CXRO_ASF('h.nff.txt', P_Torr=2*myP, L_nm=myL)
O = CXRO.load_CXRO_ASF('o.nff.txt', P_Torr=2*myP, L_nm=myL)
N.name = 'N2'
H.name = 'H2'
O.name = 'O2'

# plot transmission vs. energy for each gas
df_list = [Ar, He, N, H, O]
fig, ax = CXRO.plot_trans(df_list, 'XUV propagation in partial vacuum\n')
ax.text(0.3,
        0.5,
        f'd ={myL*(1e-9): 0.0f} m\n P ={myP: 1.0e} Torr',
        transform=ax.transAxes)
fig.tight_layout()
if save_fig:
    plt.savefig(os.path.join(diss_dir, 'XUVinVacuum.pdf'), dpi=300)
#%% XUV transmission in HPC's outer shroud

os.chdir(CrossSection_dir)

myL = (6.88e-3)*(1e9)  # nm
myP = 5  # Torr

Ar = CXRO.load_CXRO_ASF('ar.nff.txt', P_Torr=myP, L_nm=myL)
He = CXRO.load_CXRO_ASF('he.nff.txt', P_Torr=myP, L_nm=myL)
N = CXRO.load_CXRO_ASF('n.nff.txt', P_Torr=2*myP, L_nm=myL)
H = CXRO.load_CXRO_ASF('h.nff.txt', P_Torr=2*myP, L_nm=myL)
O = CXRO.load_CXRO_ASF('o.nff.txt', P_Torr=2*myP, L_nm=myL)
N.name = 'N2'
H.name = 'H2'
O.name = 'O2'

# plot transmission vs. energy for each gas
df_list = [Ar, He, N, H, O]
fig, ax = CXRO.plot_trans(df_list, r'HPC: XUV reabsorption in $P_M$ region')
ax.text(0.4,
        0.5,
        f'd ={myL*(1e-9): 0.0f} m\n P ={myP: 1.0f} Torr',
        transform=ax.transAxes)
ax.legend(frameon=True, ncol=2)
fig.tight_layout()
if save_fig:
    plt.savefig(os.path.join(diss_dir, 'HPC_absorption.pdf'), dpi=300)


#%% sample transmission

os.chdir(trans_dir)

Ge100 = CXRO.load_CXRO_trans('Ge_100nm.dat')
GaAs50 = CXRO.load_CXRO_trans('GaAs_50nm.dat')
Si100 = CXRO.load_CXRO_trans('Si_100nm.dat')
SiN30 = CXRO.load_CXRO_trans('SiNx_30nm.dat')
Cr2O3 = CXRO.load_CXRO_trans('Cr2O3_10nm.dat')
WS2 = CXRO.load_CXRO_trans('WS2_20nm.dat')

df_list = [Ge100, Si100, SiN30, Cr2O3, WS2, GaAs50]


fig, ax = CXRO.plot_trans(df_list, 'Sample Transmission')
ax.set_xlim(10, 130)
ax.set_ylim(0, 1.0)

# resize figure
fig.set_size_inches(1.2*3.5, 2.625, forward=True)
legend_pos = 'bottom'

if legend_pos == 'right':
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=True)
elif legend_pos == 'bottom':
    # Shrink current axis's height by 10% on the bottom
    dx = 0.2
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * dx,
                     box.width, box.height * (1-dx)])

    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -dx),
              frameon=True, ncol=3)

if save_fig:
    plt.savefig(os.path.join(diss_dir, 'Sample_transmission_CXRO.pdf'), dpi=300)

#%% Filter transmission

os.chdir(trans_dir)

Al200 = CXRO.load_CXRO_trans('Al_200nm.dat')
Sn200 = CXRO.load_CXRO_trans('Sn_200nm.dat')
Zr200 = CXRO.load_CXRO_trans('Zr_200nm.dat')

df_list = [Al200, Sn200, Zr200]
fig, ax = CXRO.plot_trans(df_list, 'Filter Transmission')
ax.set_xlim(10, 250)
ax.set_ylim(0, 1.0)
# resize figure
fig.set_size_inches(1.2*3.5, 2.625, forward=True)
legend_pos = 'bottom'

if legend_pos == 'right':
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=True)
elif legend_pos == 'bottom':
    # Shrink current axis's height by 10% on the bottom
    dx = 0.2
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * dx,
                     box.width, box.height * (1-dx)])

    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -dx),
              frameon=True, ncol=3)

if save_fig:
    plt.savefig(os.path.join(diss_dir, 'Filter_transmission_CXRO.pdf'), dpi=300)

