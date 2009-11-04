MODULE EOS

  USE shared_data
  USE normalise

  IMPLICIT NONE

CONTAINS

  ! This module contains all the information about the equations of state
  ! used by LARE. 

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


END MODULE EOS
