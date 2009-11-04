MODULE EOS

  USE shared_data
  USE normalise

  IMPLICIT NONE

CONTAINS

  ! This module contains all the information about the equations of state
  ! used by LARE. These are expressed in normalised units except for
  ! get_energy which is only used in the initial conditions and MUST be
  ! in SI units using the define constants constants

  SUBROUTINE get_pressure(rho_in, en_in, m_in, ix, iy, iz, p)
    REAL(num), INTENT(IN) :: rho_in, en_in ! input energy & density
    INTEGER, INTENT(IN) :: m_in            ! EOS number
    INTEGER, INTENT(IN) :: ix, iy, iz
    REAL(num), INTENT(OUT) :: p            ! output pressure

    IF (m_in .EQ. EOS_IDEAL) THEN
      p = en_in * rho_in * (gamma - 1.0_num)
      RETURN
    END IF

    IF (m_in .EQ. EOS_PI) THEN
      p = en_in * rho_in * (gamma - 1.0_num)
      RETURN
    END IF

    IF (m_in .EQ. EOS_ION) THEN
      p = (en_in - (1.0_num - xi_n(ix, iy, iz)) * ionise_pot) &
          * (gamma - 1.0_num) * rho_in
      RETURN
    END IF

  END SUBROUTINE get_pressure



  SUBROUTINE get_temp(rho_in, energy_in, m_in, ix, iy, iz, temp_out)

    REAL(num), INTENT(IN) :: rho_in, energy_in
    INTEGER, INTENT(IN) :: m_in, ix, iy, iz
    REAL(num), INTENT(OUT) :: temp_out

    ! mbar and kb will be the correct form for the normalisation when
    ! the code is running, or when setting up initial conditions with
    ! SI_Code = F, mbar and kb are both 1.0, reducing to the normal case
    ! when setting up initial conditions with SI_Code = T, kb and mbar
    ! will have their normal SI values

    IF (m_in .EQ. EOS_IDEAL) THEN
      temp_out = energy_in * (gamma - 1.0_num) / 2.0_num
      RETURN
    END IF

    IF (m_in .EQ. EOS_PI) THEN
      temp_out = (gamma - 1.0_num) * energy_in / (2.0_num - xi_n(ix, iy, iz))
      RETURN
    END IF

    IF (m_in .EQ. EOS_ION) THEN
      temp_out = (gamma - 1.0_num) &
          * (energy_in - (1.0_num - xi_n(ix, iy, iz)) * ionise_pot) &
          / ((2.0_num - xi_n(ix, iy, iz)))
      RETURN
    END IF

  END SUBROUTINE get_temp



  SUBROUTINE get_cs(rho_in, energy_in, m_in, ix, iy, iz, cs_out)

    REAL(num), INTENT(IN) :: rho_in, energy_in
    INTEGER, INTENT(IN) :: m_in, ix, iy, iz
    REAL(num), INTENT(OUT) :: cs_out

    IF (m_in .EQ. EOS_IDEAL) THEN
      cs_out = SQRT(gamma * (gamma - 1.0_num) * energy_in)
      RETURN
    END IF

    IF (m_in .EQ. EOS_PI) THEN
      cs_out = SQRT(gamma * (gamma - 1.0_num) * energy_in)
      RETURN
    END IF

    IF (m_in .EQ. EOS_ION) THEN
      cs_out = SQRT(gamma * (gamma - 1.0_num) * energy_in)
      RETURN
    END IF

  END SUBROUTINE get_cs



  SUBROUTINE get_energy(rho_in, temp_in, m_in, ix, iy, iz, en_out)

    REAL(num), INTENT(IN) :: rho_in, temp_in
    INTEGER, INTENT(IN) :: m_in, ix, iy, iz
    REAL(num), INTENT(OUT) :: en_out
    REAL(num) :: xi_local, bof, r

    IF (m_in .EQ. EOS_IDEAL) THEN
      en_out = temp_in * kb / ((gamma - 1.0_num) * mbar / 2.0_num)
      RETURN
    END IF

    IF (m_in .EQ. EOS_PI) THEN
      ! Since we can't guarantee that the ionisation fraction already
      ! calculated is correct here, calculate it straight from the temperature
      bof = tr_bar / (f_bar * SQRT(temp_in)) &
          * EXP((0.25_num * (tr_bar * temp_in - 1.0_num) + 1.0_num) &
          * T_bar / temp_in)
      r = 0.5_num * (-1.0_num + SQRT(1.0_num + r_bar * rho_in * bof))
      xi_local = r / (1.0_num + r)

      en_out = (kb * temp_in * (2.0_num - xi_local)) &
          / (MBAR * (gamma - 1.0_num))
      RETURN
    END IF

    IF (m_in .EQ. EOS_ION) THEN
      ! Since we can't guarantee that the ionisation fraction already
      ! calculated is correct here, calculate it straight from the temperature
      bof = tr_bar / (f_bar * SQRT(temp_in)) &
          * EXP((0.25_num * (tr_bar * temp_in - 1.0_num) + 1.0_num) &
          * T_bar / temp_in)
      r = 0.5_num * (-1.0_num + SQRT(1.0_num + r_bar * rho_in * bof))
      xi_local = r / (1.0_num + r)

      en_out = (kb * temp_in * (2.0_num - xi_local) &
          + (1.0_num - xi_local) * ionise_pot * (gamma - 1.0_num)) &
          / (MBAR * (gamma - 1.0_num))
      RETURN
    END IF

  END SUBROUTINE get_energy

END MODULE EOS
