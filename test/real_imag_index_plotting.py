# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 14:37:02 2019

@author: greg
"""

#%% load things
import matplotlib.pyplot as plt
plt.style.use(['science-fixed','high-vis'])
plt.rcParams['figure.dpi'] = 240  # fix high-dpi display scaling issues
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
import os
import sys

script_dir = os.path.dirname(os.path.realpath(__file__))  # location of this script file
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

#%% load refractive index info

os.chdir(index_dir)
my_h = 100  # sample thickness [nm]
Ge = CXRO.load_CXRO_index(fname='ge.dat', h=my_h)
Al = CXRO.load_CXRO_index(fname='al.dat', h=my_h)
Si = CXRO.load_CXRO_index(fname='si.dat', h=my_h)

#%% plot Ge data
fig, ax = plt.subplots(2, 1, sharex=True)

ax[0].plot(Ge['E'], Ge['delta'], label=r'$\delta$')
ax[0].plot(Ge['E'], Ge['beta'], label=r'$\beta$')
ax[0].legend(frameon=True)
ax[0].set_ylabel(r'$\tilde{n}$')
ax[0].set_ylim(0)
plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

ax[1].plot(label='bulk trans / total trans - 1')
ax[1].plot((Ge['trans_bulk']/Ge['trans_Infbounce'])-1, label='error')
ax[1].set_ylabel(r'rel. error, $\epsilon$')

plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

fig.suptitle(f'{my_h} nm Ge')
ax[-1].set_xlabel('Energy [eV]')
ax[1].set_xlim(min(Ge['E']), 50)
# if save_fig:
#    plt.savefig(os.path.join(diss_dir, 'Ge_transmission_Fresnel.pdf'), dpi=300)

#%% plot Si data
fig, ax = plt.subplots(2, 1, sharex=True)

ax[0].plot(Si['E'], Si['delta'], label=r'$\delta$')
ax[0].plot(Si['E'], Si['beta'], label=r'$\beta$')
ax[0].legend(frameon=True)
ax[0].set_ylabel(r'$\tilde{n}$')

ax[1].plot(label='bulk trans / total trans - 1')
ax[1].plot(Si['E'], Si['frac_err'], label='error')
ax[1].set_ylabel(r'rel. error, $\epsilon$')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

fig.suptitle(f'{my_h} nm Si')
ax[-1].set_xlabel('Energy [eV]')
ax[1].set_xlim(70, 130)
E_mask = (Si['E'] > 70) & (Si['E'] < 130)
ax[1].set_ylim(
    0.7*min(Si['frac_err'][E_mask]),
    1.1*max(Si['frac_err'][E_mask])
    )
ax[0].set_ylim(
    1.2*min([min(Si['delta'][E_mask]), min(Si['beta'][E_mask])]),
    1.1*max([max(Si['delta'][E_mask]), max(Si['beta'][E_mask])])
    )
# if save_fig:
#    plt.savefig(os.path.join(diss_dir, 'Si_transmission_Fresnel.pdf'), dpi=300)

#%% plot both Si and Ge
fig, ax = plt.subplots(2, 2, sharex='col', figsize=(1.5*3.5, 1.5*2.625))

# Ge data
ax[0, 0].plot(Ge['E'], Ge['delta'], label=r'$\delta$')
ax[0, 0].plot(Ge['E'], Ge['beta'], label=r'$\beta$')
ax[0, 0].legend(frameon=True)
ax[0, 0].set_ylabel(r'$\tilde{n} = (1-\delta) - i k$')
ax[0, 0].set_ylim(0)
ax[0, 0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1, 0].plot(Ge['E'], Ge['frac_err'], label='error', c='k')
ax[1, 0].set_ylabel(r'rel. error, $\epsilon$')
ax[1, 0].set_xlabel('Energy [eV]')
ax[1, 0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1, 0].set_ylim(2e-3, 3.2e-3)
ax[1, 0].set_ylim()
ax[0, 0].set_title(f'{my_h} nm Ge')

ax[1, 0].set_xlim(min(Ge['E']), 50)
E_mask = (Ge['E'] >= min(Ge['E'])) & (Ge['E'] <= 50)
ax[0, 0].set_ylim(
    0.5*min([min(Ge['delta'][E_mask]), min(Ge['beta'][E_mask])]),
    1.1*max([max(Ge['delta'][E_mask]), max(Ge['beta'][E_mask])])
    )
ax[1, 0].set_ylim(0, 3.4e-3)

# Si data
ax[0, 1].plot(Si['E'], Si['delta'], label=r'$\delta$')
ax[0, 1].plot(Si['E'], Si['beta'], label=r'$\beta$')
ax[0, 1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0, 1].set_title(f'{my_h} nm Si')
ax[1, 1].plot(Si['E'], Si['frac_err'], label='error', c='k')
ax[1, 1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1, 1].set_xlabel('Energy [eV]')
ax[1, 1].set_xlim(70, 130)
E_mask = (Si['E'] > 70) & (Si['E'] < 130)
ax[0, 1].set_ylim(
    1.2*min([min(Si['delta'][E_mask]), min(Si['beta'][E_mask])]),
    1.1*max([max(Si['delta'][E_mask]), max(Si['beta'][E_mask])])
    )
ax[1, 1].set_ylim(0, 2.5e-4)

if save_fig:
    plt.savefig(os.path.join(diss_dir, 'Ge_Si_transmission_Fresnel.pdf'),
                dpi=300)
