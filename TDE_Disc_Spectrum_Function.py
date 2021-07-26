"""
A Python function which returns the X-ray spectrum resulting from a cool (kT < 0.3 keV) accretion disc.  
Python implimentation of the model from the MNRAS Letter "Tidal disruption event discs are larger than they seem:
removing systematic biases in TDE X-ray spectral modelling", Andrew Mummery, 2021. 

This model will be particularly useful for interpreting TDE observations, 
as the derivation of this model does not make any assumptions about the 
accretion disc having evolved to the steady state. 

Standard disc models (e.g. the DISKBB model in XSPEC) assume that the radial 
profile of the accretion disc is that of a steady state configuration. 
This enforces radial structure onto the solutions which may well not be present 
in reality, particularly at the earliest (and brightest) times. 

This model therefore represents a more agnostic approach to modelling TDE X-ray spectra, 
assuming only that the each disc radius emits like a colour-corrected blackbody, 
and that there exists some temperature maximum within the accretion disc. 

Author: 
    Andrew Mummery, 
    Oxford University Astrophysics,
    andrew.mummery@wadham.ox.ac.uk

Paper:
    "Tidal disruption event discs are larger than they seem: removing systematic biases in TDE X-ray spectral modelling"
    Andrew Mummery, Oxford University Astrophysics.
    Submitted to MNRAS, June 2021. Accepted: july 2021. 
    arXiv:XXXX.XXXXX

git repo:
    https://github.com/andymummeryastro/TDEdiscXraySpectrum

"""

import numpy as np## For mathematical operators. 
from astropy import constants## For fundamental constants.

h = constants.h.value#Planck
c = constants.c.value#Speed of light
kb = constants.k_B.value#Boltzmann 
keV_to_Hz = 1000 * constants.e.value/h#Takes a photon energy in keV and converts it to a frequency in Hertz
Mpc_to_cm = 3.08567758128e24#Converts distances in Megaparsec to centimeters. 

def col_corr(temp):
    """
    The Done et al. (2012) model of the colour correction factor of an accretion disc.  
    Reference: Done C., Davis S.W., Jin C., Blaes O., Ward M., 2012, MNRAS, 420, 1848
    
    input:
        temp - Peak disc temperature, in units of Kelvin. 
    output:
        f_col(T) - disc colour correction factor, dimensionless. 
    """
    if temp < 3e4:
        ans = 1
    elif 3e4 < temp and temp < 1e5:
        ans = np.power(temp/3e4, 0.8333598980732597)# Power from math.log( 11598*72000/1e5, base = 1e5/3e4)
    else:
        ans = np.power(11598*72000/temp, 1/9)
    return ans

def TDEdiscXraySpectrum(E, Rp, Tp, gamma, D_Mpc=100):
    """
    Returns the observed flux from an accretion disc, according to the model of Mummery 2021c. 
    Reference: Submitted to MNRAS, June 2021. Accepted: XXXX., arXiv:XXXX.XXXXX 
    
    input:
        E - observed photon energy, in units of keV.
        Rp - Radius at which disc temperature peaks, in units of 10^12 cm. 
        Tp - Physical peak disc temperature, in units of 10^5 Kelvin. 
        gamma - dimensionless and restricted to the range ~ (1/2, 3/2). 
        D_mpc - the source-observer distance, in units of Mpc. 
    
    output:
        Fv - The observed accretion disc flux, in units of erg/s/cm^2/Hz. 
    
    Notes:
        The parameter Rp is related to the physical size of the disc (R0) approximately by Rp = R0 \sqrt{cos\theta},
        where \theta is the disc-observer inclination angle (See letter for more detail).
    
        Theoretical values for gamma: 
            gamma = 1/2 for a vanishing ISCO stess disc observed face on 
            gamma = 1 for a vanishing ISCO stess disc observed edge on 
            gamma = 1 for a finite ISCO stess disc observed face on 
            gamma = 3/2 for a finite ISCO stess disc observed edge on 
    
        This function makes no assumption about the properties of the disc's radial temperature profile, other than 
        that each disc radius emits like a colour-corrected blackbody, and that some temperature maximum Tp exists 
        within the disc. 
    """
    Rp *= 1e12## Corrects units. 
    Tp *= 1e5## Corrects units. 
    v = E * keV_to_Hz## Corrects units. 
    D = D_Mpc * Mpc_to_cm## Corrects units. 
    
    fcol = col_corr(Tp)## Gets colour correction factor, uses Done et al. 2012 model. 
    
    ## The peak effective temperature \tilde T_p requires the photons energy-shift, f_gamma. 
    ## This energy shift must in general be calculated numerically, which is beyond the scope of this 
    ## simple treatment. As an approximation we use the value of f_gamma which would occur if the peak
    ## disc temperature occured at the ISCO of a Schwarzschild black hole, which was then observed face on. 
    ## In this case f_gamma = 1/U^0(r_ISCO) = 1/sqrt{2}. Where U^0 is the time component of the discs 4-velocity.
    Tp = Tp/np.sqrt(2)
    
    xi1 = np.sqrt(2*np.pi/(3/4*7/4))## Determined analytically, see equation 87 of Mummery & Balbus 2020a.
    # Reference: Mummery A., Balbus S. A., 2020, MNRAS, 492, 5655
    xi2 = 3.501## Determined numerically, see Appendix B of Mummery & Balbus 2021a. 
    xi3 = 1.499## Determined numerically, see Appendix B of Mummery & Balbus 2021a. 
    # Reference: Mummery, A. & Balbus, S. A., arXiv:2104.06177
    
    higher_order = (1 + xi2 * kb*fcol*Tp/(h*v) + xi3 * (kb*fcol*Tp/(h*v)) ** 2)
    fv = 4 * np.pi * h * v**3 / (fcol**4 * c**2) * xi1 * (Rp/D)**2 * np.power(kb*fcol*Tp/(h*v), gamma) * np.exp(-h*v/(kb * fcol * Tp)) * higher_order * 1e3
    # The above equation is eq. 3 of the Mummery 2021c Letter. The factor 1000 corrects the units to erg/s/cm^2/Hz.  
    return fv



# End. 
