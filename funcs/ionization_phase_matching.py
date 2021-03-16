# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 12:47:33 2021

@author: smith
"""

#%% load things
import numpy as np
import os
import sys

# define useful file paths
funcs_dir = os.path.dirname(os.path.realpath(__file__))  # this script file
parent_dir = os.path.split(funcs_dir)[0]
CrossSection_dir = os.path.join(parent_dir, 'CrossSections')
index_dir = os.path.join(parent_dir, 'index')
trans_dir = os.path.join(parent_dir, 'Transmissions')
sys.path.append(funcs_dir)  # add funcs's parent directory to path

# load modules
import CXRO_funcs as CXRO

#%% function definitions
def eta_c(gas, lam1, lamq):
    """
    Calculate the critical ionization fraction, eta_c.

    Follows eqn (3) from Paul (2006).

    Parameters
    ----------
    gas : string
        which gas to use.
    lam1 : float
        fundamental wavelength (micron).
    lamq : float
        XUV wavelength (micron).

    Returns
    -------
    ans : float
        critical ionization fraction.
    """
    rho0 = 0.02504e27  # [1/m3] number density of ideal gas at STP
    re = 2.8179403227e-15  # [m] classical electron radius
    deltan = calc_deltan(gas, lam1, lamq)
    num = rho0*re*(lam1*1e-6)**2
    denom = 2*np.pi*deltan
    ans = (1 + num/denom)**(-1)
    return ans


def Sellmeier(gas, xum):
    """
    Calculate refractive index, n(xum).

    Parameters
    ----------
    gas : string
        which gas to use.
    xum : float
        wavelength [micron].

    Returns
    -------
    n : float
        real part of refractive index.
    """
    if gas == 'Ar':
        # https://refractiveindex.info/tmp/data/main/Ar/Peck-15C.html
        # xum: 0.47 - 2.06 micron
        n = 1 + 6.432135E-5 + 2.8606021E-2/(144-xum**-2)
    elif gas == 'He':
        # https://refractiveindex.info/tmp/data/main/He/Mansfield.html
        # xum: 0.47 - 2.06 micron
        n = 1 + 0.01470091/(423.98-xum**-2)
    elif gas == 'N2':
        # https://refractiveindex.info/tmp/data/main/N2/Peck-15C.html
        # xum: 0.47 - 2.06 micron
        n = 1 + 6.497378E-5 + 3.0738649E-2/(144-xum**-2)
    else:
        print('error! unrecognized gas!')
    return n


def calc_deltan(gas, lam1, lamq):
    """
    Calculate Delta n= n(lam1) - n(lamq).

    Parameters
    ----------
    gas : string, optional
        gas name. The default is 'Ar'.
    lam1 : float
        fundamental wavelength [micron].
    lamq : float
        harmonic wavelength [micron].

    Returns
    -------
    ans : float
        difference in refractive indices.

    """
    lamq = lamq * 1000.0  # convert lamq to nm
    lam1 = float(lam1)  # convert lam1 to float (int causes error)
    fname = os.path.join(CrossSection_dir, gas.lower() + '.nff.txt')
    nlam1 = Sellmeier(gas, lam1)
    df = CXRO.load_CXRO_ASF(fname=fname, P_Torr=760, L_nm=100)
    nlamq = df.iloc[(df['WL']-lamq).abs().argsort()[:1]]['n'].values[0].real
    ans = nlam1 - nlamq
    return ans


def calc_rhoopt(rho0, q, lam, deltan, w0, eta, etac):
    """
    Calculate rho_opt, assuming we are at the focus and on the optical axis.

    Parameters
    ----------
    rho0 : float
        density at STP.
    q : float
        harmonic order (can be fractional).
    lam : float
        fundamental wavelength [micron].
    deltan : float
        refractive index difference (n_XUV - n_IR).
    w0 : float
        beam waist [micron].
    eta : float
        current ionization fraction.
    etac : float
        critical ionization fraction.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    eta_frac = 1 - eta/etac
    num = rho0 * (q-1) * 2 * lam**2
    denom = q * deltan * w0**2 * eta_frac
    return num/denom


def calc_rhoopt_2(rho0, q, lam, deltan, w0, eta, etac, z, P0):
    """
    Calculate rho_opt, at arbitrary z, but on the optical axis.

    Parameters
    ----------
    rho0 : float
        density at STP.
    q : float
        harmonic order (can be fractional).
    lam : float
        fundamental wavelength [micron].
    deltan : float
        refractive index difference (n_XUV - n_IR).
    w0 : float
        beam waist [micron].
    eta : float
        current ionization fraction.
    etac : float
        critical ionization fraction.
    z : float
        position relative to focus. positive=downstream, negative=upstream
    P0 : float
        input power. units are [W/cm2]

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    aq = 2e-14  # long trajectory [cm^2/W]

    # calculate rayleigh range and confocal parameter
    zR = np.pi * w0**2 / lam  # [micrpn]
    b1 = 2*zR  # [micron]

    # calculate dI/dz [W/um^3]
    def calc_dIdz(P0, z, zR, lam):
        deriv = -4*P0*z*zR**2 / (lam * (z**2 + zR**2)**2)
        # print(f'deriv = {deriv}')
        return deriv

    dIdz = calc_dIdz(P0, z, zR, lam)

    # convert power and aq to microns from cm
    P0 = P0 * 1e-8  # convert to [W/micron2]
    aq = aq * 1e8  # convert to [micron^2/W]

    # eta_frac = 1 - eta/etac
    # fact1 = z*zR / (z**2+zR**2)**2 * 4*P0*aq / (np.pi*w0**2)
    # fact2 = -1.0*(q-1) / (z**2+zR**2)
    # prefact = rho0*zR*lam / (q*2*np.pi*deltan*eta_frac)
    # return prefact * (fact1 + fact2)

    num = etac*lam*rho0 * (2*b1*(q-1) + dIdz*aq*(b1**2 + 4*z**2))
    denom = 2*np.pi*q*(b1**2+4*z**2)*deltan*(eta-etac)
    return num / denom


def calc_rhoopt_3(gas, lam, lamq_range, w0):
    """
    Calculate rho_opt / rho_0 at the focus and on the optical axis.

    For a given gas, fundamental wavelength and spot size, calculate the 
    optimal density for each of the harmonic wavelengts in <lamq_range>.

    Parameters
    ----------
    gas : string
        name of gas (He, Ar, etc.).
    lam : float
        fundamental wavelength [micron].
    lamq_range : np.ndarray of floats
        array of harmonic wavelengths [micron].
    w0 : float
        beam wasit [micron].

    Returns
    -------
    popt : float
        optimum density.

    """
    q_range = lam / lamq_range  # convert to harmonic order [unitless]

    deltan = np.vectorize(calc_deltan)(gas, lam1=lam, lamq=lamq_range)
    etac = np.vectorize(eta_c)(gas, lam, lamq_range)
    eta = 0.5*etac.min()  # guarantee that we don't exceed eta_c
    popt = np.vectorize(calc_rhoopt)(
        rho0=1,
        q=q_range,
        lam=lam,
        deltan=deltan,
        w0=w0,
        eta=eta,
        etac=etac)
    return popt
