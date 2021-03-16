# -*- coding: utf-8 -*-
"""
Created on Sat Jul 20 11:53:11 2019

@author: Greg Smith
"""

#%% load things
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
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
import geometry_funcs as geo

# save figures to disk?
save_fig = True
diss_dir = r'C:\PhDPythonScripts\CXRO\figures'

# create fancy text for plotting
l1_t = r'$l_1$'
l2_t = r'$l_2$'
alpha_t = r'$\alpha$'

#%% define mirror parameters
a = 1500  # ellipse width [mm]
b = 113.219   # ellipse height [mm]
x0 = 752.146  # impact point / x-position of center of mirror [mm]
L = 2 * 90  # length of mirror [mm]

#%% plot grazing angle over mirror surface
my_xL = geo.xL(a, b, x0, L)
my_xR = geo.xR(a, b, x0, L)
M = 3  # demagnification
y0 = geo.ellipse(x0, a, b)  # y-position of center of mirror
my_eps = geo.eps(a, b)
my_l1 = geo.l1(x0, y0, a, b, my_eps)
my_l2 = geo.l2(x0, y0, a, b, my_eps)

# grazing angle over mirror surface
xplot = np.linspace(my_xL, my_xR, 500)

fig, ax = plt.subplots()
ax.plot(xplot-x0, geo.alpha1(xplot, a, b)*(180/np.pi))
ax.plot(xplot-x0, [5]*len(xplot))
ax.set_xlabel('Horizontal Displacement from Center of Mirror [mm]')
ax.set_ylabel('Grazing Angle [deg]')
ax.set_title(f'Effect of Curvature on Grazing Angle\n{l1_t}={my_l1: 4.0f} mm, {l2_t}={my_l2: 3.0f} mm, {alpha_t} = 5 deg')
if save_fig:
    plt.savefig(os.path.join(diss_dir, 'EM_angle.pdf'), dpi=300)

#%% plot 2D representation of mirror
xplot = np.linspace(-a, a, 500)  # entire top half of ellipse
xplot2 = np.linspace(my_xL, my_xR, 500)  # optic surface only
f1 = np.array([-a*my_eps, 0])  # source point
f2 = np.array([a*my_eps, 0])   # focal point

fig, ax = plt.subplots()
ax.plot(xplot, geo.ellipse(xplot, a, b), c='k', ls=':')  # top of ellipse
ax.plot(xplot, -1.0*geo.ellipse(xplot, a, b), c='k', ls=':')  # bottom of ellipse
ax.plot(xplot2, geo.ellipse(xplot2, a, b), c='b', ls='-', linewidth=3)  # mirror surface
ax.scatter(x=-a*my_eps, y=0, c='k', s=6)  # source point
ax.scatter(x=a*my_eps, y=0, c='k', s=6)  # focal point
ax.plot([-a*my_eps, my_xL], [0, geo.ellipse(my_xL, a, b)], c='k', ls='-')  # line between source & xL
ax.plot([-a*my_eps, x0], [0, geo.ellipse(x0, a, b)], c='r', ls='-')  # line between source & x0
ax.plot([-a*my_eps, my_xR], [0, geo.ellipse(my_xR, a, b)], c='k', ls='-')  # line between source & xR
ax.plot([a*my_eps, my_xL], [0, geo.ellipse(my_xL, a, b)], c='k', ls='-')  # line between focus & xL
ax.plot([a*my_eps, x0], [0, geo.ellipse(x0, a, b)], c='r', ls='-')  # line between focus & x0
ax.plot([a*my_eps, my_xR], [0, geo.ellipse(my_xR, a, b)], c='k', ls='-')  # line between focus & xR
ax.set_ylim([-1.1*b, 1.1*b])
ax.text(x=-1375, y=-25, s=r'$f_1$')
ax.text(x=1275, y=-25, s=r'$f_2$')
ax.text(x=x0, y=105, s=r'$S$')
ax.set_xlabel('x [mm]')
ax.set_ylabel('y [mm]')
ax.set_title(f'Ellipse Geometry\n{l1_t}={my_l1: 4.0f} mm, {l2_t}={my_l2: 3.0f} mm, {alpha_t} = 5 deg')
if save_fig:
    plt.savefig(os.path.join(diss_dir, 'EM_2D.pdf'), dpi=300)

#%% plot mirror's deviation from ideal
os.chdir(parent_dir)
df = pd.read_csv('KP 15-0001.ASC', sep='\t', header=None,
                 names=['X', 'Y', 'Z'], skiprows=1)
df['Z'] = 1e6*df['Z']  # nm
x = np.unique(df['X'])
y = np.unique(df['Y'])
X, Y = np.meshgrid(x, y)
Z = df['Z'].values.reshape(len(y), len(x))

cmap = plt.get_cmap('seismic')
fig, ax = plt.subplots()
im = ax.pcolormesh(X, Y, Z, cmap=cmap,
                   norm=geo.MidpointNormalize(midpoint=0.0), rasterized=True)
ax.set_title('Deviation from Ideal Ellipsoid')
ax.set_xlabel('x [mm]')
ax.set_ylabel('y [mm]')
cbar = fig.colorbar(im, ax=ax)
cbar.ax.get_yaxis().labelpad = 15
cbar.ax.set_ylabel('height error [nm]', rotation=270)
if save_fig:
    plt.savefig(os.path.join(diss_dir, 'EM_error.pdf'), dpi=300)
