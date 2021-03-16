# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 16:13:00 2019.

@author: Greg Smith
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

#  add CrossSections directory to PATH
sys.path.insert(0, CrossSection_dir)

# load modules
import Fresnel_funcs as FF

# should we save the figures?
save_fig = True
diss_dir = r'C:\PhDPythonScripts\CXRO\figures'

#%% definitions
re = 2.81794032e-6  # classical electron radius, nm
hc = 1239.84193  # eV*nm
NA = 6.0221409e+23  # Avogadro's number

# fancy LaTeX font for plotting
theta_t = r'$\theta_i$'
sigma_t = r'$\sigma$'

os.chdir(index_dir)
#%% Gold: refractive index and critical angle

# load data from CXRO file
Au = FF.calc_refl('au.dat', 85*(np.pi/180), 0)

# plot it
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(1*3.5, 1.5*2.625))

ax[0].semilogy(Au['E'], Au['delta'], label=r'$\delta$')
ax[0].semilogy(Au['E'], Au['beta'], label=r'$\beta$')
ax[0].set_ylim([1e-3, 1])
ax[0].legend(frameon=True)
ax[0].set_ylabel('Ref. Index')
ax[0].grid()

ax[1].plot(Au['E'], Au['th_c'])
ax[1].set_ylabel('Crit. Angle [deg]')
ax[1].grid()

ax[0].set_xlim([30, 500])
ax[1].set_xlabel('Energy [eV]')
ax[0].set_title('Au')

#%% Smooth mirror: reflectance vs material
th = 85*(np.pi/180)  # radian
sig = 0  # nm

Au = FF.calc_refl('au.dat', th, sig)
Ag = FF.calc_refl('ag.dat', th, sig)
Al = FF.calc_refl('al.dat', th, sig)
Si = FF.calc_refl('si.dat', th, sig)
Ir = FF.calc_refl('ir.dat', th, sig)
Cu = FF.calc_refl('cu.dat', th, sig)
W = FF.calc_refl('w.dat', th, sig)
Rh = FF.calc_refl('rh.dat', th, sig)
Pt = FF.calc_refl('pt.dat', th, sig)
Ru = FF.calc_refl('ru.dat', th, sig)
Mo = FF.calc_refl('mo.dat', th, sig)
Pd = FF.calc_refl('pd.dat', th, sig)
Ni = FF.calc_refl('ni.dat', th, sig)
C = FF.calc_refl('c.dat', th, sig)

fig, ax = plt.subplots()
ax.plot(Al['E'], Al['Rs'], label='Al')  # Z=13
# ax.plot(Si['E'], Si['Rs'], label='Si')  # Z=14
ax.plot(Cu['E'], Cu['Rs'], label='Cu')  # Z=29
ax.plot(Ag['E'], Ag['Rs'], label='Ag')  # Z=47
ax.plot(W['E'], W['Rs'], label='W')  # Z=74
ax.plot(Ir['E'], Ir['Rs'], label='Ir')  # Z=77
ax.plot(Au['E'], Au['Rs'], label='Au')  # Z=79
ax.legend(frameon=True, ncol=2)
ax.set_xlabel('Energy [eV]')
ax.set_ylabel('Reflectance')
ax.set_title(f's-pol, {theta_t}={th*(180/np.pi):.0f} degrees, {sigma_t}=0 nm')
if save_fig:
    plt.savefig(os.path.join(diss_dir, 'Fresnel_NoSigma.pdf'), dpi=300)

fig, ax = plt.subplots()
ax.plot(Au['E'], Au['Rs'], label='Au')
ax.plot(C['E'], C['Rs'], label='C')
ax.plot(Pt['E'], Pt['Rs'], label='Pt')
ax.plot(Rh['E'], Rh['Rs'], label='Rh')
ax.plot(Ni['E'], Ni['Rs'], label='Ni')
ax.plot(Mo['E'], Mo['Rs'], label='Mo')
ax.plot(Pd['E'], Pd['Rs'], label='Pd')
# ax.plot(Al['E'], Al['Rs'], label='Al')
ax.plot(Ru['E'], Ru['Rs'], label='Ru')
ax.plot(Si['E'], Si['Rs'], label='Si')
ax.legend(frameon=True, ncol=2)
ax.set_xlabel('Energy [eV]')
ax.set_ylabel('Reflectance')
ax.set_title(f's-pol, {theta_t}={th*(180/np.pi):.0f} degrees, {sigma_t}=0 nm')
# if save_fig:
#     plt.savefig('Fresnel_NoSigma.png', dpi=300)

#%% Smooth mirror: reflectance vs incident angle (1-bounce)
theta_list = [1, 3, 5, 10, 15]  # degrees
sig = 0  # nm
theta_list = list((90 - np.array(theta_list))*(np.pi/180))

fig, ax = plt.subplots()
for t in theta_list:
    Au = FF.calc_refl('au.dat', t, sig)
    ax.plot(Au['E'], Au['Rs'], label=f'{t*(180/np.pi):.0f} deg')
ax.legend(frameon=True, framealpha=0.9)
ax.set_xlabel('Energy [eV]')
ax.set_ylabel('Reflectance')
ax.set_title(f'Au, s-pol, {sigma_t}=0 nm, 1 bounce')
if save_fig:
    plt.savefig(os.path.join(diss_dir, 'Au_ReflvsAngle.pdf'), dpi=300)

#%% Smooth mirror: reflectance vs incident angle (2-bounces)
theta_list = [3, 5, 7.5, 10]
sig = 0  # nm
theta_list = list((90 - np.array(theta_list))*(np.pi/180))

fig, ax = plt.subplots()
for t in theta_list:
    Au = FF.calc_refl('au.dat', t, sig)
    if (90-(180/np.pi)*t).is_integer():
        ax.plot(Au['E'], Au['Rs']**2, label=f'{t*(180/np.pi):.0f} deg')
    else:
        ax.plot(Au['E'], Au['Rs']**2, label=f'{t*(180/np.pi):.1f} deg')
ax.legend(loc=1, frameon=True, framealpha=0.9)
ax.set_xlabel('Energy [eV]')
ax.set_ylabel('Reflectance')
ax.set_title(f'Au, s-pol, {sigma_t}=0 nm, 2 bounces')
if save_fig:
    plt.savefig(os.path.join(diss_dir, 'Au_ReflvsAngle_2bounce.pdf'), dpi=300)

#%% Gold mirror: reflectance vs surface roughness
th = 85*(np.pi/180)  # radian
sig = 0.3  # rms [nm]
Au = FF.calc_refl('au.dat', th, sig)

fig, ax = plt.subplots()
ax.plot(Au['E'], Au['Rs']*Au['DW'], label='Au, s-pol')
ax.legend(frameon=True)
ax.set_xlabel('Energy [eV]')
ax.set_ylabel('Reflectance')
ax.set_title(f'{theta_t}={th*(180/np.pi):.0f} degrees, {sigma_t}={sig} nm')

sig_list = [0, 0.3, 3, 10, 30]  # rms roughness [nm]

fig, ax = plt.subplots()
for s in sig_list:
    Au = FF.calc_refl('au.dat', th, s)
    # ax.plot(Au['E'], Au['Rs']*Au['DW'], label=f's-pol, {sigma_t}={s} nm')
    ax.plot(Au['E'], Au['Rs_DW'], label=f'{sigma_t}={s} nm')
    # ax.plot(Au['E'], Au['Rp']*Au['DW'], label=f'p-pol, {sigma_t}={s} nm')
ax.legend(frameon=True)
ax.set_xlabel('Energy [eV]')
ax.set_ylabel('Reflectance')
ax.set_title(f'Au, s-pol, {theta_t}={th*(180/np.pi):.0f} degrees')
if save_fig:
    plt.savefig(os.path.join(diss_dir, 'R_vs_roughness.pdf'), dpi=300)

fig, ax = plt.subplots()
for s in sig_list:
    Au = FF.calc_refl('au.dat', th, s)
    #ax.plot(Au['E'], Au['Rs']*Au['NC'], label=f's-pol, {sigma_t}={s} nm')
    ax.plot(Au['E'], Au['Rs_NC'], label=f's-pol, {sigma_t}={s} nm')
    # ax.plot(Au['E'], Au['Rp']*Au['NC'], label=f'p-pol, {sigma_t}={s} nm')
ax.legend(frameon=True)
ax.set_xlabel('Energy [eV]')
ax.set_ylabel('Reflectance')
ax.set_title(f'Au, NC-factor, {theta_t}={th*(180/np.pi):.0f} degrees')

#%% compare to CXRO's reflectance calculator
sig = 0.3  # rms [nm]
th = 85*(np.pi/180)  # radian

Au = FF.calc_refl('au.dat', th, sig)
Ag = FF.calc_refl('ag.dat', th, sig)

CXRO_0nm = pd.read_csv(
    'CXRO_Au_0nm.dat', delimiter=' ', skiprows=2, skipinitialspace=True,
    usecols=[0, 1], names=['E', 'R'], dtype={'E': np.float64, 'R': np.float64}
    )
CXRO_03nm = pd.read_csv(
    'CXRO_Au_0.3nm.dat', delimiter=' ', skiprows=2, skipinitialspace=True,
    usecols=[0, 1], names=['E', 'R'], dtype={'E': np.float64, 'R': np.float64}
    )
CXRO_3nm = pd.read_csv(
    'CXRO_Au_3nm.dat', delimiter=' ', skiprows=2, skipinitialspace=True,
    usecols=[0, 1], names=['E', 'R'], dtype={'E': np.float64, 'R': np.float64}
    )

fig, ax = plt.subplots()
ax.plot(CXRO_03nm['E'], CXRO_03nm['R'], label='CXRO')
# ax.plot(Au['E'], Au['Rs']*Au['NC'], label='NC')
ax.plot(Au['E'], Au['Rs_NC'], label='NC')
# ax.plot(Au['E'], Au['Rs']*Au['DW'], label='DW')
ax.plot(Au['E'], Au['Rs_DW'], label='DW')
ax.legend(frameon=True)
ax.set_xlabel('Energy [eV]')
ax.set_ylabel('Reflectance')
ax.set_title(f'Au, {sigma_t}={sig} nm, {theta_t}={th*(180/np.pi)} deg')

#%% plot extinction length

theta_list = [1, 3, 5]
theta_list = list((90 - np.array(theta_list))*(np.pi/180))
sig = 0  # nm

fig, ax = plt.subplots()
for th in theta_list:
    Au = FF.calc_refl('au.dat', th, sig)
    ax.plot(Au['E'], 1/Au['k2'].to_numpy().imag,
            label=f'{theta_t} = {th*(180/np.pi)} deg')
    # ax.plot(Ag['E'], 1/Ag['k2'].to_numpy().imag,
    # label=f'Ag, {theta_t}={th*(180/np.pi)} deg')
ax.legend(frameon=True)
ax.set_xlabel('Energy [eV]')
ax.set_ylabel(r'$1/k_2$ [nm]')
ax.set_title('Extinction Length in Gold')
if save_fig:
    plt.savefig(os.path.join(diss_dir, 'Au_ExtinctionLength.pdf'), dpi=300)
