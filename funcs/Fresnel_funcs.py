# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 14:24:19 2021

@author: smith
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

#%% function definitions

re = 2.81794032e-6  # classical electron radius, nm
hc = 1239.84193  # eV*nm
NA = 6.0221409e+23  # Avogadro's number


def calc_refl(fname, th, sig):
    """
    Master function that calculates Fresnel reflectivity.

    Parameters
    ----------
    fname : string
        filename of CXRO Energy, delta, beta values.
    th : TYPE
        incidence angle [radians]. normal = 0 degrees.
    sig : TYPE
        surface roughness [nm].

    Returns
    -------
    df : pandas dataframe
        columns are:
        index: energy [eV].
        delta: from n
        Beta: from n
        n: refractive index, n = 1 - delta - i*Beta
        WL: wavelength [nm]
        rs: Fresnel amplitude reflectivity (s-pol)
        rp: Fresnel amplitude reflectivity (p-pol)
        Rs: Fresnel intensity reflectivity (s-pol)
        Rp: Fresnel intensity reflectivity (p-pol)
        k1: normal component of wave vector in medium 1, Zangwill eqn 17.59
        k2: normal component of wave vector in medium 2, Zangwill eqn 17.59
        NC: Nevot-Croce factor (amplitude)
        DW: Debye-Waller factor (amplitude)
        Rs_NC: Fresnel (intensity) reflectivity, s-pol, NC factor
        Rp_NC: Fresnel (intensity) reflectivity, p-pol, NC factor
        Rs_DW: Fresnel (intensity) reflectivity, s-pol, DW factor
        Rp_DW: Fresnel (intensity) reflectivity, p-pol, DW factor
        th_c: Critical angle [degrees]. 0 degrees = glancing, 90 = normal.
        df.name: name the dataframe after the filename (element)
    """
    df = load_index(fname)
    df['rs'] = Fresnel_r(df, th, 's')
    df['rp'] = Fresnel_r(df, th, 'p')
    df['Rs'] = Fresnel_R(df, th, 's')
    df['Rp'] = Fresnel_R(df, th, 'p')
    df['k1'] = np.cos(th)*(2*np.pi/df['WL'])
    df['k2'] = 1j*(2*np.pi/df['WL'])*np.sqrt(
        np.sin(th)**2 - df['n'].to_numpy().real**2
        )
    df['NC'] = np.exp(-2*df['k1'] * df['k2'] * sig**2)
    df['DW'] = np.exp(-2*df['k1']**2 * sig**2)
    df['Rs_NC'] = np.abs(df['rs'] * df['NC'])**2
    df['Rp_NC'] = np.abs(df['rp'] * df['NC'])**2
    df['Rs_DW'] = np.abs(df['rs'] * df['DW'])**2
    df['Rp_DW'] = np.abs(df['rp'] * df['DW'])**2
    df['th_c'] = critical_angle(df)
    df.name = fname.split('.')[0].capitalize()
    return df


def critical_angle(df):
    """
    Compute the critical angle [degrees].

    Parameters
    ----------
    df : pandas dataframe
        DESCRIPTION.

    Returns
    -------
    float
        Critical angle [degrees] th_c = 0 is glancing, th_c = 90 is normal to
        the incident plane.
    """
    return 90 - (180/np.pi)*np.arcsin(1-df['delta'])


def load_CXRO_index(fname):
    """
    Load atomic scattering file from CXRO.

    Parameters
    ----------
    fname : CXRO file from online database
        columns:
        E: energy [eV]
        f1, f2: atomic scattering factors from f = f1 + i*f2
        WL: wavelength, from energy [nm]
        mu_a: photoatomic cross section [nm^2]

    Returns
    -------
    pandas dataframe
        index: energy [eV].
        WL: wavelength [nm]
        mu_a: photoatomic cross section
        n: refractive index, n = 1 - delta - i*Beta
        delta: from n
        Beta: from n
        df.name: name the dataframe after the element

    """
    def calc_N(rho, M):
        """
        Calculate the number density in atoms/nm^3.

        Parameters
        ----------
        rho : float
            massdensity [g/cm3].
        M : float
            atomic mass [g/mol].

        Returns
        -------
        N: float
            number density [atoms/nm3].
        """
        return (1e-21)*rho*NA/M
    # calculate N
    if fname.split('.')[0] == 'au':
        rho = 19.32  # mass density [gram/cm3]
        M = 196.9665687  # atomic weight [gram/mol]
    elif fname.split('.')[0] == 'ag':
        rho = 10.49
        M = 107.8682
    else:
        return 'unknown density!'
    N = calc_N(rho, M)  # number density [atoms/cm3]
    df = pd.read_csv(
        fname, delimiter='\t', skiprows=1, usecols=[0, 1, 2],
        names=['E', 'f1', 'f2'],
        dtype={'E': np.float64, 'f1': np.float64, 'f2': np.float64}
        )
    df['WL'] = hc / df['E']
    df['mu_a'] = 2 * re * df['WL'] * df['f2']
    df['n'] = 1 - (N*re*df['WL']**2*(df['f1']+1j*df['f2']))/(2*np.pi)
    df['delta'] = 1 - df['n'].to_numpy().real
    df['beta'] = -np.imag(df['n'])
    df.name = fname.split('.')[0].capitalize()
    return df


def load_index(fname):
    """
    Load a CXRO file of energy, delta and beta (real/imag indices).

    Parameters
    ----------
    fname : string
         CXRO file from online database.

    Returns
    -------
    df : pandas dataframe
        n: refractive index, n = 1 - delta - i*Beta
        WL: wavelength [nm]
        delta: from n
        beta: from n

    """
    df = pd.read_csv(
        fname, delimiter=' ', skiprows=2, skipinitialspace=True,
        usecols=[0, 1, 2], names=['E', 'delta', 'beta'],
        dtype={'E': np.float64, 'delta': np.float64, 'beta': np.float64}
        )
    df['n'] = 1 - df['delta'] - 1j*df['beta']
    df['WL'] = hc / df['E']
    return df


def Fresnel_r(df, th, type='s'):
    """
    Calculate the Fresnel reflectance factors for s and p polarized light.

    Assumes initial medium is vacuum (n1=1), and that the second medium is
    gold,silver,etc. (n2 from df)

    Parameters
    ----------
    df : pandas datafame
        Contains n vs E data, calculated from CXRO file.
    th : float
        incidence angle (th=0 is normal incidence) in radians.
    type : string, optional
        Polarization type: s or p. The default is 's'.

    Returns
    -------
    R : float
        Fresnel reflectance factor vs energy.

    """
    if type == 's':
        r = (
            (np.cos(th) - df['n']*np.sqrt(1-(np.sin(th)/df['n'])**2)) /
            (np.cos(th) + df['n']*np.sqrt(1-(np.sin(th)/df['n'])**2))
            )
    elif type == 'p':
        r = (
            (np.sqrt(1-(np.sin(th)/df['n'])**2) - df['n']*np.cos(th)) /
            (np.sqrt(1-(np.sin(th)/df['n'])**2) + df['n']*np.cos(th))
            )
    else:
        print('error: polarization must be s or p')
    return r


def Fresnel_R(df, th, type='s'):
    """
    Calculate the Fresnel reflectance factors for s and p polarized light.

    Assumes initial medium is vacuum (n1=1), and that the second medium is
    gold,silver,etc. (n2 from df)

    Parameters
    ----------
    df : pandas datafame
        Contains n vs E data, calculated from CXRO file.
    th : float
        incidence angle (th=0 is normal incidence) in radians.
    type : string, optional
        Polarization type: s or p. The default is 's'.

    Returns
    -------
    R : float
        Fresnel reflectance factor vs energy.

    """
    if type == 's':
        R = np.abs(
            (np.cos(th) - df['n']*np.sqrt(1-(np.sin(th)/df['n'])**2))
            / (np.cos(th) + df['n']*np.sqrt(1-(np.sin(th)/df['n'])**2))
            )**2
    elif type == 'p':
        R = np.abs(
            (np.sqrt(1-(np.sin(th)/df['n'])**2) - df['n']*np.cos(th))
            / (np.sqrt(1-(np.sin(th)/df['n'])**2) + df['n']*np.cos(th))
            )**2
    else:
        print('error: polarization must be s or p')
    return R


def plot_crit_angle(df):
    """
    Plot the critical angle.

    Parameters
    ----------
    df : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    fig, ax = plt.subplots()
    ax.plot(df['E'], np.arcsin(df['n'])*(180/np.pi))
    ax.axhline(85)
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('Critical Angle [deg]')
    return


def calc_roughness(df, th=85*(np.pi/180), sig=0.3):
    """
    Calculate roughness factors.

    Parameters
    ----------
    df : pandas dataframe
        from the CXRO database.
    th : FLOAT, optional
        incidence angle (0=normal) in radians. The default is 85*(np.pi/180).
    sig : FLOAT, optional
        surface roughness [nm]. The default is 0.3.

    Returns
    -------
    df : pandas dataframe
        Same as input, but with extra columns:
        k1: normal component of wave vector in medium 1, Zangwill eqn 17.59
        k2: normal component of wave vector in medium 2, Zangwill eqn 17.59
        NC: Nevot-Croce factor
        DW: Debye-Waller factor

    """
    df['k1'] = np.cos(th)*(2*np.pi/df['WL'])
    df['k2'] = 1j*(2*np.pi/df['WL'])*np.sqrt(
        np.sin(th)**2 - df['n'].to_numpy().real**2
        )
    df['NC'] = np.exp(-4*df['k1']*df['k2'] * sig**2)
    df['DW'] = np.exp(-4*df['k1']**2 * sig**2)
    return df