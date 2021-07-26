# README #
## A Fortran and Python implimentation of a tidal disruption event X-ray spectrum disc model
A Fortran90 and Python implimentation of a model for the X-ray spectrum resulting from a cool (kT < 0.3 keV) accretion disc.  This model will be particularly useful for interpreting TDE observations, 
as the derivation of this model does not make any assumptions about the accretion disc having evolved to the steady state. 

This model was presented in the MNRAS Letter "Tidal disruption event discs are larger than they seem:
removing systematic biases in TDE X-ray spectral modelling", Andrew Mummery, 2021. arxiv:
### Fortran implimentation 
* For use in the XSPEC software (https://heasarc.gsfc.nasa.gov/xanadu/xspec/)
* XSPEC requires the files lmodel.dat, load.xcm and diskmodel.f90
* diskmodel.f90 contains the actual fortran implementation of the model 
* diskmodel.f90 will need to be compiled. This can be done with e.g., the terminal command "gfortran diskmodel.f90"
* See https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixLocal.html for more information on how to load the TDE disc model into XSPEC
* For any questions about the model, please email andrew.mummery@wadham.ox.ac.uk 
* If you use this model in published work you must cite: Mummery, A., 2021, MNRAS, 
### Python implimentation
* Also included for ease of plotting and analysis
* Returns disc spectrum in units of erg/s/cm^2, which differs from the XSPEC Fortran implimentation which returns photons/s/cm^2
