!A Fortran90 subroutine which returns the X-ray spectrum resulting from a cool (kT < 0.3 keV) accretion disc.  
!Implimentation of the model from the MNRAS Letter "Tidal disruption event discs are larger than they seem:
!removing systematic biases in TDE X-ray spectral modelling", Andrew Mummery, 2021. 
!
! IF YOU USE THIS MODEL IN PUBLISHED WORK YOU MUST CITE: Mummery, A., 2021, MNRAS, 
! arXiv:XXXX.XXXXX
!
!
!The following parameters are required for model
!input:
!    Rp - Radius at which disc temperature peaks, in units of 10^12 cm. 
!    Tp - Physical peak disc temperature, in units of 10^5 Kelvin. 
!    gamma - dimensionless and restricted to the range ~ (1/2, 3/2). 
!    D_mpc - the source-observer distance, in units of Mpc. 
!
!output from model:
!    n(E)dE - The observed accretion disc spectrum, in units of photons/s/cm^2 
!
!Notes:
!    The parameter Rp is related to the physical size of the disc (R0) approximately by Rp = R0 \sqrt{cos\theta},
!    where \theta is the disc-observer inclination angle (See published letter for more detail).
!
!    Theoretical values for gamma: 
!        gamma = 1/2 for a vanishing ISCO stess disc observed face on 
!        gamma = 1 for a vanishing ISCO stess disc observed edge on 
!        gamma = 1 for a finite ISCO stess disc observed face on 
!        gamma = 3/2 for a finite ISCO stess disc observed edge on 
!
!    This function makes no assumption about the properties of the disc's radial temperature profile, other than 
!    that each disc radius emits like a colour-corrected blackbody, and that some temperature maximum Tp exists 
!    within the disc. 
!
!Author: 
!    Andrew Mummery, 
!    Oxford University Astrophysics,
!    andrew.mummery@wadham.ox.ac.uk
!
!Paper:
!    "Tidal disruption event discs are larger than they seem: removing systematic biases in TDE X-ray spectral modelling"
!    Andrew Mummery, Oxford University Astrophysics.
!    Submitted to MNRAS, June 2021. Accepted: July 2021. 
!    arXiv:XXXX.XXXXX
!
!git repo:
!    https://github.com/andymummeryastro/TDEdiscXraySpectrum
!
!
program wrapper  
! The wrapper program allows the user to run this model in the terminal. 
! This may be useful for outputting standard spectra and comparing models.
! The model used by XSPEC is the below subroutine 'diskmodel'.
  implicit none
  integer ne,i,ifl
  parameter (ne=6000)
  real Emax,Emin,ear(0:ne),param(4),photar(ne),E,dE

! Parameters (Default parameter values -- for running locally and testing).
  param(1) = 1.0      ! Radial size of disc in units of 10^12 cm
  param(2) = 1.0      ! Phyiscal peak disc temperature in units of 10^5 K 
  param(3) = 0.51     ! gamma (dimensionless), restricted to 0.5 ≤ gamma ≤ 1.5 (from theory).
  param(4) = 100.0    ! Source-observer distance, in Mpc. 

! Set energy grid  
  Emax  = 10.0
  Emin  = 1e-1
  do i = 0,ne
    ear(i) = Emin * (Emax/Emin)**(real(i)/real(ne))
  end do

! Call the model
  call diskmodel(ear,ne,param,ifl,photar)

! Write out model output
  do i = 1,ne
     E  = 0.5 * ( ear(i) + ear(i-1) )
     dE =         ear(i) - ear(i-1)
     write(99,*)E,E**2*photar(i)/dE
  end do
  write(99,*)"log"
  
end program wrapper

!=======================================================================
subroutine diskmodel(ear,ne,param,ifl,photar)
!A Fortran90 subroutine which returns the X-ray spectrum resulting from a cool (kT < 0.3 keV) accretion disc.  
!Implimentation of the model from the MNRAS Letter "Tidal disruption event discs are larger than they seem:
!removing systematic biases in TDE X-ray spectral modelling", Andrew Mummery, 2021. 
!
!This subroutine is written specifically for use with the XSPEC software. 
!It therefore returns the X-ray spectrum in units of photons/cm^2/s. 
!
!This model will be particularly useful for interpreting TDE observations, 
!as the derivation of this model does not make any assumptions about the 
!accretion disc having evolved to the steady state. 
!
!Standard disc models (e.g. the DISKBB model in XSPEC) assume that the radial 
!profile of the accretion disc is that of a steady state configuration. 
!This enforces radial structure onto the solutions which may well not be present 
!in reality, particularly at the earliest (and brightest) times. 
!
!This model therefore represents a more agnostic approach to modelling TDE X-ray spectra, 
!assuming only that the each disc radius emits like a colour-corrected blackbody, 
!and that there exists some temperature maximum within the accretion disc. 
  implicit none
  integer ne,ifl,i
  real ear(0:ne),param(4),photar(ne)
  real Tp,Tpp,Ep
  double precision pi, xi1, xi2, xi3, Rp, Dmpc, gamma, fc, ho, fo, amp
  real E,dE
  pi  = acos(-1.d0)
  ifl = 1
  
! Parameters
  Rp    = dble( param(1) )
  Tp  = dble( param(2) )
  gamma = dble( param(3) )
  Dmpc  = dble( param(4) )

! Get colour correction factor. 
  Tpp = Tp * 100000.0!Temperature in physical units, important for colour correction
  if (Tpp .lt. 30000.0) then
      fc = 1.0
  else
      if (Tpp .lt. 100000.0) then
          fc = (Tpp/30000.0) ** (0.833359898)! Power from math.log( 11598*72000/1e5, base = 1e5/3e4)
      else
          fc = (11598.0*72000.0/(Tpp)) ** (0.111111111)!Power is 1/9
      end if
  end if 
  
  Ep = Tp/1.41421356 * 0.00861732815!Final numerical factor takes 10^5 Kelvin to keV. 
! Factor of Root 2 is an approximation for the gravitational red shift. 
  
! Required unit conversions:
!  keV_to_Hz = 2.4179905e17!Takes a photon energy in keV to Hz.
!  Mpc_to_cm = 3.0856776e24!Takes a distance in Mpc to cm. 
!  c = 2.99792458e10!Units of cm/s
  amp = 1655969.941 * (Rp/Dmpc)**2! Numeric factor is (keV_to_Hz)^3 / c^2 * (10^12cm/1Mpc)^2, with c in cm/s. 

! Required non-fundamental constants
  xi1 = 2.188!sqrt(32*pi/21)
  xi2 = 3.501
  xi3 = 1.499
  
  ! Calculate n(E)dE for energy grid. 
  do i = 1,ne
    E        = 0.5 * ( ear(i) + ear(i-1) )
    dE       = ear(i) - ear(i-1)
    ho = (1.0 + xi2 * (fc*Ep)/E + xi3 * ((fc*Ep)/E)**2)!Higher order corrections. 
    fo = 4.0*pi/(fc**4) * xi1 * E**2 * (fc*Ep/E)**gamma * exp(-E/(fc*Ep))!Leading order, neglecting amplitude. 
    photar(i) = amp * fo * ho * dE!Spectrum is amplitude * first order * higher order corrections * dE
  end do

  return
end subroutine diskmodel
!=======================================================================


