!******************************************************************************
! This module contains the conversion factors between normalised and
! internal code units
!******************************************************************************

MODULE normalise

  USE constants
  USE shared_data
  IMPLICIT NONE

CONTAINS

  SUBROUTINE normalise_transport

    REAL(num) :: eta_bar_0, pressure0, v0, energy0, mbar, kappa0

    pressure0 = B_norm**2 / mu0_si
    energy0 = pressure0 / rho_norm
    v0 = SQRT(energy0)
    mbar = mf * mh_si
    temp_norm = (mbar / kb_si) * energy0
    time_norm = L_norm / v0

    ! Normalise tbar, r_bar and eta_bar for including Cowling resistivity and
    ! neutrals
    t_bar = t_bar / temp_norm
    r_bar = r_bar * rho_norm / temp_norm**1.5_num
    eta_bar_0 = (mu0_si * L_norm * v0) * rho_norm**2 * SQRT(temp_norm) / B_norm**2
    eta_bar = eta_bar / eta_bar_0

    ! Normalise ionise_pot - inioisation potential of hydrogen
    ionise_pot = ionise_pot_si / (energy0 * mbar)
    IF (eos_number /= EOS_ION) ionise_pot = 0.0_num

    ! Normalise tr required for get_neutral etc.
    tr = tr / temp_norm

    ! Normalise parallel thermal conductivity
    kappa0 = energy0**1.5_num * rho_norm * L_norm / (mbar / kb_si * energy0)**3.5_num
    kappa_0 = 1.e-11_num / kappa0


    ! Constants used in radiative losses
    h_star = L_norm / (rho_norm * v0**3)
    lr_star = 1.148e-35_num * (rho_norm / mbar)**2

  END SUBROUTINE normalise_transport

END MODULE normalise
