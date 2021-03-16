# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 10:37:24 2019

@author: greg
"""
#%% load things
import numpy as np
import matplotlib.pyplot as plt
plt.style.use(['science-fixed','high-vis'])
plt.rcParams['figure.dpi'] = 240  # fix high-dpi display scaling issues
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

import pandas as pd
import os, sys, itertools

# Add CrossSections directory to PATH
here = os.path.dirname(os.path.abspath(__file__))  # location of this file
sys.path.insert(0, os.path.normpath(os.path.join(here, '../CrossSections')))

#%% function definitions
re = 2.81794032e-6  # classical electron radius, nm
hc = 1239.84193  # eV*nm

def load_CXRO_ASF(fname, P_Torr, L_nm):
    """
    Load atomic scattering file (ASF) from CXRO database.
    
    Parameters
    ----------
    fname : string
        filename of ASF.
    P_Torr : float
        Pressure in Torr.
    L_nm : float
        Medium length in nm, used to calculate transmission.

    Returns
    -------
    df : pandas dataframe
    columns:
        E: energy [eV]
        f1, f2: atomic scattering factors from f = f1 + i*f2
        WL: wavelength, from energy [nm]
        mu_a: photoatomic cross section [nm^2]
        N: atomic number density [atoms/nm3^]
        n: refractive index
        T: transmission
    """
    # load file
    df = pd.read_csv(fname, delimiter='\t', skiprows=1, usecols=[0, 1, 2],
                     names=['E', 'f1', 'f2'],
                     dtype={'E':np.float64, 'f1':np.float64, 'f2':np.float64})
    # calculate parameters
    df['WL'] = hc / df['E']
    df['mu_a'] = 2 * re * df['WL'] * df['f2']
    df['N'] = (P_Torr/760) * 0.0269
    df['n'] = 1 - (df['N']*re*df['WL']**2*(df['f1']+1j*df['f2']))/(2*np.pi)
    df['T'] = np.exp( - df['N']*df['mu_a']*L_nm )
    
    # name the dataframe after fname
    df.name = fname.split('.')[0].capitalize()
    return df
    

def load_CXRO_trans(fname):
    """
    Load transmission file from CXRO database.

    Assumed filename: MAT_thicknessNM

    Parameters
    ----------
    fname : string
        Filename.

    Returns
    -------
    df : pandas dataframe
    columns:
        E: energy [eV]
        T: transmission [0 to 1]
    """

    df = pd.read_csv(fname, delimiter=' ', skiprows=[0, 1],
                     skipinitialspace=True, usecols=[0, 1],
                     names=['E', 'T'],
                     dtype={'E':np.float64, 'T':np.float64})
    
    # name the dataframe after fname
    df.name = fname.split('.')[0].replace('_', ' ').replace('nm', ' nm')
    return df


def load_CXRO_index(fname='ge.dat', h=100):
    """
    Calculate Frensel information from refractive index CXRO file.

    The complex refractive index is of the form:
        n = (1-delta) - i*beta

    Parameters
    ----------
    fname : string, optional
        Filename.
    h : float, optional
        sample thickness in nm. The default is 100.

    Returns
    -------
    df : pandas dataframe
    columns:
        E: energy [eV]
        delta: from n = (1-delta) - i*beta
        beta: from n = (1-delta) - i*beta
        wL: wavelength [nm]
        k: imaginary part of complex refractive index
        n: real part of complex refractive index
        alpha: base-e alpha [nm^-1]
        alpha10: base-10 alpha
        trans_bulk: bulk absorption
        T_F: Fresnel transmittance
        R_F: Fresnel reflectance
        trans_1bounce: transmission, assuming 1 bounce
        refl_1bounce: reflection, assuming 1 bounce
        trans_Infbounce: transmission, assuming infinite bounces
        refl_Infbounce: reflection, assuming infinite bounces
        frac_err: error introduced by ignoring edge effects
    """
    df = pd.read_csv(fname, delimiter=' ', skiprows=[0, 1],
                             skipinitialspace=True, usecols=[0, 1, 2],
                             names=['E', 'delta', 'beta'],
                             dtype={'E':np.float64, 'delta':np.float64,
                                    'beta':np.float64}
                             )
    df['wL'] = hc / df['E']
    df['k'] = df['beta']
    df['n'] = 1-df['delta']
    df['alpha'] = 4*np.pi*df['k'] / df['wL']
    df['alpha10'] = df['alpha'] / np.log(10)
    df['trans_bulk'] = np.exp(-1.0*df['alpha']*h)
    df['T_F'] = 4*df['n'] / np.abs(df['n'] - 1.j*df['k'] + 1)**2
    df['R_F'] = np.abs((df['n']-1.j*df['k']-1) / (df['n']-1.j*df['k']+1))**2
    df['trans_1bounce'] = df['T_F']**2 * df['trans_bulk']
    df['refl_1bounce'] = df['R_F'] + df['R_F']*df['T_F']**2 * df['trans_bulk']**2
    df['trans_Infbounce'] = (df['T_F']**2 * df['trans_bulk']) / (1 - df['R_F']**2 * df['trans_bulk']**2)
    df['refl_Infbounce'] = df['R_F'] + (df['R_F']*df['T_F']**2 * df['trans_bulk']**2) / (1 - df['R_F']**2 * df['trans_bulk']**2)
    df['frac_err'] = (df['trans_bulk']/df['trans_Infbounce']) - 1

    # name the dataframe using the thickness
    df.name = f'{h} nm'
    return df


def load_air(P_Torr, L_nm):
    """
    Load a gas mixture file.

    Note: factor of 0.0269 [atoms/nm3] corresponds to STP.

    air.txt comes from CXRO at STP.
    
    Parameters
    ----------
    P_Torr : float
        Pressure in Torr.
    L_nm : float
        Length of gas medium in nm.

    Returns
    -------
    df : pandas dataframe
    columns:
        n: complex refractive index
        WL: wavelength [nm]
        N: atomic number density at P_Torr [atoms/nm^3]
        f2: atomic scattering factors from f = f1 + i*f2
        mu_a: photoatomic cross section
        T: transmission

    """
    df = pd.read_csv('air.txt', delimiter=' ', skiprows=2,
                     skipinitialspace=True, usecols=[0, 1, 2],
                     names=['E', 'delta', 'beta'],
                     dtype={'E':np.float64, 'delta':np.float64,
                            'beta':np.float64})
    df['n'] = 1 - df['delta'] - 1j*df['beta']
    df['WL'] = hc / df['E']
    df['N'] = (P_Torr/760) * 0.0269
    df['f2'] = 2*np.pi*df['beta'] / (0.0269*re*df['WL']**2)
    df['mu_a'] = 2 * re * df['WL'] * df['f2']
    df['T'] = np.exp(-1.0*df['N']*df['mu_a']*L_nm)
    return df


def plot_trans(df_list, title_text, logy=False):
    """
    Plot the transmission of a list of materials vs Energy.

    Parameters
    ----------
    df_list : list
        list of dataframes.
    title_text : string
        Title of the Plot.
    logy : bool, optional
        Plot in logscale? The default is False.

    Comments refer to the output for an input:
        df.name = 'Cr2O3 10 nm'

    Returns
    -------
    ax : matplotlib axis

    """
    fig, ax = plt.subplots()
    for df in df_list:
        # create label text
        mat_name = df.name.split(' ')[0]  # 'Cr2O3'
        thickness = ' '.join(df.name.split(' ')[1:])  # '10 nm'
        if hasNumbers(mat_name):  # material name has numbers in it
            # make numbers subscripts, i.e. 'Cr2O3' --> r'Cr$_{2}$O$_{3}$'
            elem_list = ["".join(x) for _, x in itertools.groupby(mat_name, key=str.isdigit)]
            mat_name = r''.join(['$_{'+s+'}$' if s.isdigit() else s for s in elem_list])
        if logy:
            ax.semilogy(df['E'], df['T'], label=mat_name+' '+thickness)
        else:
            ax.plot(df['E'], df['T'], label=mat_name+' '+thickness)
    ax.set_xlim([0, 300])
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('Transmission')
    ax.legend(frameon=True)
    ax.set_title(title_text)
    plt.show()
    return fig, ax


def plot_scat_fact(df, title_text, logy=False):
    """
    Plot the atomic scattering factors (f1, f2) vs Energy.

    Parameters
    ----------
    df : pandas dataframe
        DESCRIPTION.
    title_text : string
        DESCRIPTION.
    logy : bool, optional
        Plot in logscale?. The default is False.

    Returns
    -------
    ax : matplotlib axis

    """
    fig, ax = plt.subplots()
    if logy:
        ax.semilogy(df['E'], df['f1'], label='f1')
        ax.semilogy(df['E'], df['f2'], label='f2')
    else:
        ax.plot(df['E'], df['f1'], label='f1')
        ax.plot(df['E'], df['f2'], label='f2')
    ax.set_xlim([0, 300])
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('Atomic Scattering Factor')
    ax.legend(frameon=True)
    ax.set_title(title_text)
    plt.show()
    return ax


def plot_mu(df_list, title_text, logy=False):
    """
    Plot the photoatomic cross section vs Energy.

    Parameters
    ----------
    df_list : list
        list of dataframes.
    title_text : string
        Plot title.
    logy : bool, optional
        Plot in logscale? The default is False.

    Returns
    -------
    ax : matplotlib axis

    """
    fig, ax = plt.subplots()
    for df in df_list:
        if logy:
            ax.semilogy(df['E'], df['mu_a'], label=df.name)
        else:
            ax.plot(df['E'], df['mu_a'], label=r'$\mu_a$')
    ax.set_xlim([0, 300])
    #ax.set_ylim([0, 2])
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel(r'$\mu_a \text{ [nm}^2$\text{]}')
    ax.legend(frameon=True)
    ax.set_title(title_text)
    plt.show()
    return ax


def plot_index(df, title_text, logy=False):
    """
    Plot the real & imag parts of the refractive index vs Energy.

    Parameters
    ----------
    df : pandas dataframe
        DESCRIPTION.
    title_text : string
        Plot title.
    logy : bool, optional
        Plot in logscale? The default is False.

    Returns
    -------
    fig : matplotlib figure
    ax : matplotlib axis

    """
    fig, ax = plt.subplots()
    if logy:
        ax.semilogy(df['E'], np.real(df['n']-1), label='n-1')
        ax.semilogy(df['E'], np.imag(df['n']), label='k')
    else:
        ax.plot(df['E'], np.real(df['n']-1), label='n-1')
        ax.plot(df['E'], np.imag(df['n']), label='k')
    ax.set_xlim([0, 300])
    #ax.set_ylim([0, 2])
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('Refractive Index')
    ax.legend(frameon=True)
    ax.set_title(title_text)
    plt.show()
    return fig, ax
    

def hasNumbers(inputString):
    """Convenience function to test if a string has numbers."""
    return any(char.isdigit() for char in inputString)