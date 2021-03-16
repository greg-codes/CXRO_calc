# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 15:07:18 2021

@author: smith
"""

#%% load modules
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


#%% function definitions
class MidpointNormalize(colors.Normalize):
    '''fixes colorbar'''

    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):  
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


def alpha1(x, a, b):
    """
    Calculate the grazing angle.

    Parameters
    ----------
    x : float
        x-position of impact point on ellipse.
    a : float
        semiaxis parameter for x-axis.
    b : float
        semiaxis parameter for y-axis.

    Returns
    -------
    α : float
        DESCRIPTION.
    """
    num = a**4 - 2*(a*b)**2 - (a*x)**2 + (b*x)**2
    arg1 = (-a * np.sqrt(1 - (b/a)**2) + x)**2 + b**2 * (1 - (x/a)**2)
    arg2 = (a * np.sqrt(1 - (b/a)**2) + x)**2 + b**2 * (1 - (x/a)**2)
    denom = a**2 * np.sqrt(arg1) * np.sqrt(arg2)
    alpha = 0.5 * np.arccos(num / denom)
    return alpha


def xL(a, b, x0, L):
    """
    Calculate the left (x) side of the mirror.

    Parameters
    ----------
    a : float
        semiaxis parameter for x-axis.
    b : float
        semiaxis parameter for y-axis.
    x0 : float
        offset amount in x-direction (center of mirror).
    L : float
        length of mirror.

    Returns
    -------
    xL : float
        DESCRIPTION.

    """
    num = (a*L)**2 * (x0**2 - a**2)
    denom = 2 * np.sqrt(
        (a*L)**2 * (a-x0) * (a+x0) * (a**4 - (a*x0)**2 + (b*x0)**2)
        )
    xL = x0 + num/denom
    return xL


def xR(a, b, x0, L):
    """
    Calculate the right (x) side of the mirror.

    Parameters
    ----------
    a : float
        semiaxis parameter for x-axis.
    b : float
        semiaxis parameter for y-axis.
    x0 : float
        offset amount in x-direction (center of mirror).
    L : float
        length of mirror.

    Returns
    -------
    xR : float
        DESCRIPTION.

    """
    num = (a*L)**2 * (a-x0) * (a+x0)
    denom = 2 * np.sqrt(
        (a*L)**2 * (a-x0) * (a+x0) * (a**4 - (a*x0)**2 + (b*x0)**2)
        )
    xR = x0 + num/denom
    return xR


def eps(a, b):
    """
    Calculate the eccentricity of the ellipse.

    Parameters
    ----------
    a : float
        semiaxis parameter for x-axis.
    b : float
        semiaxis parameter for y-axis.

    Returns
    -------
    eps : float
        eccentricity.

    """
    eps = np.sqrt(1 - (b/a)**2)
    return eps


def l1(x0, y0, a, b, eps):
    """
    Calculate the entrance arm length.

    Parameters
    ----------
    x0 : float
        DESCRIPTION.
    y0 : float
        DESCRIPTION.
    a : float
        semiaxis parameter for x-axis.
    b : float
        semiaxis parameter for y-axis.
    eps : float
        eccentricity.

    Returns
    -------
    l1 : float
        entrance arm length.

    """
    l1 = np.sqrt((x0 + a*eps)**2 + y0**2)
    return l1


def l2(x0, y0, a, b, eps):
    """
    Calculate the exit arm length.

    Parameters
    ----------
    x0 : float
        DESCRIPTION.
    y0 : float
        DESCRIPTION.
    a : float
        semiaxis parameter for x-axis.
    b : float
        semiaxis parameter for y-axis.
    eps : float
        eccentricity.

    Returns
    -------
    l2 : float
        exit arm length.

    """
    l2 = np.sqrt((x0 - a*eps)**2 + y0**2)
    return l2


def ellipse(x, a, b):
    """
    Compute y(x, a, b) for an ellipse.

    Parameters
    ----------
    x : float
        x-position in 2D plane.
    a : float
        semiaxis parameter for x-axis.
    b : float
        semiaxis parameter for y-axis.

    Returns
    -------
    y : float
        corresponding y-position of ellipse.

    """
    y = np.sqrt(b**2 * (1 - (x / a)**2))
    return y