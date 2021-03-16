# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 13:56:52 2020

@author: greg

script file to calculate critical ionization fraction
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
from matplotlib.ticker import FormatStrFormatter

# define useful file paths
script_dir = os.path.dirname(os.path.realpath(__file__))  # this script file
parent_dir = os.path.split(script_dir)[0]
funcs_dir = os.path.join(parent_dir, 'funcs')
CrossSection_dir = os.path.join(parent_dir, 'CrossSections')
index_dir = os.path.join(parent_dir, 'index')
trans_dir = os.path.join(parent_dir, 'Transmissions')
sys.path.append(funcs_dir)  # add funcs's parent directory to path

# load my modules
import CXRO_funcs as CXRO  # CXRO file convenience functions
import ionization_phase_matching as IPM  # phase matching and ionization

# parameters to save files
save_fig = True  # flag to save figs or not
diss_dir = r'C:\PhDPythonScripts\CXRO\figures'


#%% heatmap: delta-n for Argon and Helium

# define wavelength ranges
lam1_range = np.linspace(0.8, 2.0, num=12)  # laser wavelength [um]
lamq_range = np.linspace(1240/30, 1240/300, num=50) * 1e-3  # harmonic wavelength [um]
lam11, lamqq = np.meshgrid(lam1_range, lamq_range)

# calculate Argon
z_Ar = np.vectorize(IPM.calc_deltan)('Ar', lam11, lamqq)

# plot Argon
fig, ax = plt.subplots()
pcm = ax.pcolormesh(lam1_range, 1.240/lamq_range, z_Ar, cmap='plasma',
                    rasterized=True)
ax.set_xlabel(r'fundamental wavelength [$\mu$m]')
ax.set_ylabel('HHG energy [eV]')
fig.colorbar(pcm, format='%.0e')
ax.set_title(r'$\Delta n$ for Argon')
if save_fig:
    plt.savefig(os.path.join(diss_dir, 'Ar_deltan.pdf'), dpi=300)


# calculate Helium
lam11, lamqq = np.meshgrid(lam1_range, lamq_range)
z_He = np.vectorize(IPM.calc_deltan)('He', lam11, lamqq)

# plot Helium
fig, ax = plt.subplots()
pcm = ax.pcolormesh(lam1_range, 1.240/lamq_range, z_He, cmap='plasma',
                    rasterized=True)
ax.set_xlabel(r'fundamental wavelength [$\mu$m]')
ax.set_ylabel('HHG energy [eV]')
fig.colorbar(pcm, format='%.0e')
ax.set_title(r'$\Delta n$ for Helium')
if save_fig:
    plt.savefig(os.path.join(diss_dir, 'He_deltan.pdf'), dpi=300)
#%% plot lineouts

# calculate delta-n for Ar & He at 0.8 micron
lamq_range = np.linspace(1240/30, 1240/300, num=50)*1e-3  # harmonic energy [um]
He_deltaN = np.vectorize(IPM.calc_deltan)('He', 0.8, lamq_range)
Ar_deltaN = np.vectorize(IPM.calc_deltan)('Ar', 0.8, lamq_range)

fig, ax = plt.subplots()
ax.plot(1.240/lamq_range, 1/He_deltaN, label='He', c='r', ls='solid')
ax.plot(1.240/lamq_range, 6*1/Ar_deltaN, label='Ar (x6)', c='b', ls='solid')
ax.set_xlim(30, 300)
ax.legend(frameon=True)
ax.set_xlabel('HH Photon Energy [eV]')
ax.set_ylabel(r'$(\Delta n)^{-1}$ for $\lambda_1 = 0.8 \ \mu$m')
ax.set_title('Optimal Phase Matching Density Scaling')
if save_fig:
    plt.savefig(os.path.join(diss_dir, 'recip_deltan_plot.pdf'), dpi=300)

#%% plot: optimum phase matching pressure at the focus

# set parameters
myw0 = 100  # spot size [um]
lamq_range = np.linspace(1240/30, 1240/300, num=50) * 1e-3  # harmonic wavelength [um]
E_range = 1.240 / lamq_range  # harmoinc energy [eV]

# calculate for 0.8, 1.3 and 1.8 micron for both He and Ar
Ar800 = IPM.calc_rhoopt_3('Ar', 0.8, lamq_range, myw0)
He800 = IPM.calc_rhoopt_3('He', 0.8, lamq_range, myw0)
Ar1300 = IPM.calc_rhoopt_3('Ar', 1.3, lamq_range, myw0)
He1300 = IPM.calc_rhoopt_3('He', 1.3, lamq_range, myw0)
Ar1800 = IPM.calc_rhoopt_3('Ar', 1.8, lamq_range, myw0)
He1800 = IPM.calc_rhoopt_3('He', 1.8, lamq_range, myw0)

# plot it
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(1*3.5, 1*2.625))

ax[0].plot(E_range, Ar800, label=r'0.8 $\mu$m', c='b', ls='solid')
ax[0].plot(E_range, Ar1300, label=r'1.3 $\mu$m', c='g', ls='solid')
ax[0].plot(E_range, Ar1800, label=r'1.8 $\mu$m', c='r', ls='solid')
ax[0].set_ylabel(r'$\rho_{opt}$/$\rho_0$, Ar')
ax[0].set_ylim(0)

ax[1].plot(E_range, He800, label=r'0.8 $\mu$m', c='b', ls='solid')
ax[1].plot(E_range, He1300, label=r'1.3 $\mu$m', c='g', ls='solid')
ax[1].plot(E_range, He1800, label=r'1.8 $\mu$m', c='r', ls='solid')
ax[1].set_ylabel(r'$\rho_{opt}$/$\rho_0$, He')
ax[1].set_ylim(0)

ax[1].legend(frameon=True, framealpha=1)
ax[0].set_xlim(30, 300)
ax[1].set_xlabel('HH Photon Energy [eV]')
plt.subplots_adjust(hspace=0)
ax[0].set_title('Optimal Phase Matching Density Scaling')

if save_fig:
    plt.savefig(os.path.join(diss_dir, 'Popt_Scaling.pdf'), dpi=300)

#%% spectrogram: critical ionization fraction for Ar, He

# define wavelength ranges
lam1_range = np.linspace(0.8, 2.0, num=50)  # [micron]
lamq_range = np.linspace(1240/30, 1240/300, num=50) * 1e-3  # micron
lam11, lamqq = np.meshgrid(lam1_range, lamq_range)

# do the calculations
z_Ar = np.vectorize(IPM.eta_c)('Ar', lam11, lamqq)
z_He = np.vectorize(IPM.eta_c)('He', lam11, lamqq)

# plot Argon
fig, ax = plt.subplots()
pcm = ax.pcolormesh(lam1_range, 1.240/lamq_range, 100*z_Ar, cmap='plasma',
                    rasterized=True)
ax.set_xlabel(r'fundamental wavelength [$\mu$m]')
ax.set_ylabel('HHG energy [eV]')
fig.colorbar(pcm)
ax.set_title(r'Argon: $\eta_c$ [\%]')
if save_fig:
    plt.savefig(os.path.join(diss_dir, 'Ar_critical_ionization.pdf'), dpi=300)
# plot Helium
fig, ax = plt.subplots()
pcm = ax.pcolormesh(lam1_range, 1.240/lamq_range, 100*z_He, cmap='plasma',
                    rasterized=True)
ax.set_xlabel(r'fundamental wavelength [$\mu$m]')
ax.set_ylabel('HHG energy [eV]')
fig.colorbar(pcm)
ax.set_title(r'Helium: $\eta_c$ [\%]')
if save_fig:
    plt.savefig(os.path.join(diss_dir, 'He_critical_ionization.pdf'), dpi=300)


#%% lineouts: critical ionization fraction for Ar and He

# define wavelengths
lam1_range = [0.8, 1.3, 1.8]
lamq_range = np.linspace(1240/30, 1240/300, num=50) * 1e-3  # micron
lam11, lamqq = np.meshgrid(lam1_range, lamq_range)

# do the calculations
z_Ar = np.vectorize(IPM.eta_c)('Ar', lam11, lamqq)
z_He = np.vectorize(IPM.eta_c)('He', lam11, lamqq)

# plot it
fig, ax = plt.subplots(2, 1, sharex=True)
for i, lami in enumerate(lam1_range):
    ax[0].plot(1.240/lamq_range, 100*z_Ar[:, i], label=f'{lami}' + r' $\mu$m')
    ax[1].plot(1.240/lamq_range, 100*z_He[:, i], label=f'{lami}' + r' $\mu$m')

# use a common legend for both plots
ax[0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
           frameon=True, ncol=3)
plt.subplots_adjust(hspace=0.5)
ax[1].set_xlabel('HHG energy [eV]')
ax[0].set_ylabel('Ar')
ax[1].set_ylabel('He')
ax[0].set_xlim([np.min(1.240/lamq_range), np.max(1.240/lamq_range)])
ax[0].set_ylim(0)
ax[1].set_ylim(0)
ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax[0].set_title(r'Critical Ionization Fraction, $\eta_c$ (\%)')
if save_fig:
    plt.savefig(os.path.join(diss_dir, 'crit_ion_frac.pdf'), dpi=300)


fig, ax = plt.subplots(figsize=(1.2*3.5, 1*2.625))
ax.semilogy(1.240/lamq_range, 100*z_Ar[:, 0], label=f'{lam1_range[0]}' + r' $\mu$m', ls='solid', c='b')
ax.semilogy(1.240/lamq_range, 100*z_Ar[:, 1], label=f'{lam1_range[1]}' + r' $\mu$m', ls='solid', c='g')
ax.semilogy(1.240/lamq_range, 100*z_Ar[:, 2], label=f'{lam1_range[2]}' + r' $\mu$m', ls='solid', c='r')
ax.semilogy(1.240/lamq_range, 100*z_He[:, 0], ls='dashed', c='b')
ax.semilogy(1.240/lamq_range, 100*z_He[:, 1], ls='dashed', c='g')
ax.semilogy(1.240/lamq_range, 100*z_He[:, 2], ls='dashed', c='r')
ax.set_xlabel('HHG energy [eV]')
ax.set_ylabel(r'$\eta_c$ (\%)')
ax.set_xlim([np.min(1.240/lamq_range), np.max(1.240/lamq_range)])
ax.grid(axis='y', which='both')
ax.grid(axis='x', which='major')
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.set_ylim(0.08, 10)

# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=True)
fig.suptitle(r'Critical Ionization Fraction, $\eta_c$ (\%)')

if save_fig:
    plt.savefig(os.path.join(diss_dir, 'crit_ion_frac_log.pdf'), dpi=300)
