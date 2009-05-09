!**************************************************************
! This module contains the conversion factors between normalised
! and internal code units
!**************************************************************

MODULE normalise
  USE constants
  USE shared_data
  IMPLICIT NONE

  SAVE

  ! These are the real SI physical constants
  ! Permiability of free space
  REAL(num), PARAMETER :: mu0_0 =  4.0e-7_num * pi

  ! Boltzmann's Constant
  REAL(num), PARAMETER :: kb_0 = 1.3806504e-23_num

  ! Mass of hydrogen ion
  REAL(num), PARAMETER :: mh_0 = 1.660538782e-27_num

  ! Mass of electron
  REAL(num), PARAMETER :: me_0 = 9.10938188e-31_num

  ! Planck's constant
  REAL(num), PARAMETER :: h_0 = 6.626068e-34_num

  ! Ionisation potential of hydrogen in J
  REAL(num), PARAMETER :: ionise_pot_0 = 2.17870364e-18_num

  ! These contain the correct value for the constant
  ! at the point where the code is running
  ! Permiability of free space
  REAL(num) :: mu0 =  mu0_0

  ! Boltzmann's Constant
  REAL(num) :: kb = kb_0

  ! Mass of hydrogen ion
  REAL(num) :: mh = mh_0

  ! Average mass of all particles in proton masses
  REAL(num) :: mf

  ! Average mass of a particle (assuming even mix of hydrogen ions & electrons)
  REAL(num) :: MBAR

  ! Magnetic field conversion factor (in T)
  REAL(num) :: B0 = 1.0_num

  ! Specific energy density conversion factor in K
  REAL(num) :: ENERGY0 = 1.0_num

  ! Velocity conversion factor in m / s
  REAL(num) :: VEL0 = 1.0_num

  ! Density conversion factor in kgm^{ - 3}
  REAL(num) :: RHO0 = 1.0_num

  ! Time conversion factor in units of s
  REAL(num) :: T0 = 1.0_num

  ! Length conversion factor in m
  REAL(num) :: L0 = 1.0_num

  ! Current conversion factor in A
  REAL(num) :: J0 = 1.0_num

  ! Temperature conversion factor
  REAL(num) :: TEMP0 = 1.0_num

  ! Pressure Conversion factor
  REAL(num) :: PRESSURE0 = 1.0_num

  ! These conversion factors are used to convert additional variables
  ! They are locked once the basic conversion variables are set
  ! Gravity conversion factor in ms^{ - 2}
  REAL(num) :: GRAV0 = 1.0_num

  ! Viscosity conversion factor in m^2 / s
  REAL(num) :: VISC0 = 1.0_num

  ! Reistivity conversion factor in m^2 / s
  REAL(num) :: RES0 = 1.0_num

  ! Thermal conductivity conversion factor in kg / m / s
  REAL(num) :: KAPPA0 = 1.0_num

CONTAINS

  SUBROUTINE set_normalisation

    ! Remember that the way to set the normalisations is to set
    ! B0 the normalising magnetic field (in Tesla)
    ! L0 the normalising length (in m)
    ! Rho0 the normalising density (in kgm^{ - 3})
    ! The code will generate the rest

    IF (.NOT. SI) THEN
      ! If not running as an SI code then force normalisation off
      ! Ignore any values read from the input deck
      B0 = 1.0_num
      RHO0 = 1.0_num
      L0 = 1.0_num

      ! Set the constants to one as well, so that
      ! you can use them in the main code
      mu0 = 1.0_num
      kb = 1.0_num
      mh = 1.0_num
      mf = 1.0_num
    END IF

    ! Calculate the derived quantities
    CALL derived_quantities

  END SUBROUTINE set_normalisation



  SUBROUTINE derived_quantities

    ! Average particle mass is mass of proton * mass of average particle
    ! in proton masses
    MBAR = mh * mf

    VEL0 = B0 / SQRT(MU0 * RHO0) ! Velocity
    ENERGY0 = VEL0**2            ! Specific energy density
    T0 = L0 / VEL0               ! Time

    ! Put code in here to normalise any derived quantities
    GRAV0  = VEL0**2 / L0      ! g in kgms^ - 2

    ! viscosity (Input as inverse Reynolds so L0 * VEL0 * RHO0 not needed)
    VISC0  = 1.0_num

    ! resistivity (Input as inverse Lundquist so L0 * VEL0 * MU0 not needed)
    RES0   = 1.0_num
    PRESSURE0 = B0**2 / MU0  ! Pressure
    TEMP0 = MBAR * PRESSURE0 / (KB * RHO0) ! Temperature in K
    J0 = B0 / (L0 * MU0)

    ! Thermal conductivity
    KAPPA0 = ENERGY0**(3.0_num / 2.0_num) * RHO0 * L0 &
        / (MBAR / KB * ENERGY0)**(7.0_num / 2.0_num)

  END SUBROUTINE derived_quantities

END MODULE normalise
