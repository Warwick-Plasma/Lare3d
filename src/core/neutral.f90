! All the subroutines in this module are for the partially ionised flux
! emergence simulations; see Leake & Arber, 2006

MODULE neutral

  USE shared_data
  USE boundary
  USE eos

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: perpendicular_resistivity, newton_relax, &
      neutral_fraction, setup_neutral, get_neutral, get_energy

CONTAINS


  SUBROUTINE setup_neutral
    ! Ion neutral collision cross section(m^2)
    REAL(num) :: sigma_in = 5.0e-19_num

    IF (include_neutrals) THEN
      ALLOCATE(xi_n(-1:nx+2, -1:ny+2, -1:nz+2))
      xi_n = 0.0_num
    END IF

    IF (cowling_resistivity) THEN
      ALLOCATE(eta_perp(-1:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(parallel_current(0:nx, 0:ny, 0:nz))
      ALLOCATE(perp_current(0:nx, 0:ny, 0:nz))
    END IF

    ! Ionisation potential of hydrogen(J)
    ionise_pot = ionise_pot_0

    ! Temperature of the photospheric radiation field
    tr = 6400.0_num 

    ! Calculate fbar^(2 / 3) in (k^-1 m^-2)
    f_bar = pi * (me_0 / h_0) * (kb_0 / h_0)
    f_bar = SQRT(2.0_num) * f_bar**(3.0_num / 2.0_num)

    ! Calculate tbar in (K)
    t_bar = ionise_pot / kb_0

    ! Calculate rbar in (kg^-1)
    mbar = mh * mf                     
    r_bar = 4.0_num / mbar

    ! Calculate eta_bar in (m^4 / (k s kg^2))
    eta_bar = 2.0_num * mbar &
        / (SQRT(16.0_num * kb_0 / (pi * mbar)) * sigma_in)

  END SUBROUTINE setup_neutral



  SUBROUTINE perpendicular_resistivity

    ! This subroutine calculates the cross field resistivity at the current
    ! temperature. 

    REAL(num) :: f, xi_v, bxv, byv, bzv, bfieldsq, rho_v, t_v, T
    INTEGER :: ixp, iyp, izp

    DO iz = 0, nz
      DO iy = 0, ny
        DO ix = 0, nx
          ixp = ix + 1
          iyp = iy + 1
          izp = iz + 1

          ! Get the vertex density
          rho_v = rho(ix, iy, iz) * cv(ix, iy, iz) &
              + rho(ixp, iy , iz ) * cv(ixp, iy , iz ) &
              + rho(ix , iyp, iz ) * cv(ix , iyp, iz ) &
              + rho(ixp, iyp, iz ) * cv(ixp, iyp, iz ) &
              + rho(ix , iy , izp) * cv(ix , iy , izp) &
              + rho(ixp, iy , izp) * cv(ixp, iy , izp) &
              + rho(ix , iyp, izp) * cv(ix , iyp, izp) &
              + rho(ixp, iyp, izp) * cv(ixp, iyp, izp)

          rho_v = rho_v / (cv(ix, iy, iz) + cv(ixp, iy, iz) &
              + cv(ix, iyp, iz ) + cv(ixp, iyp, iz ) &
              + cv(ix, iy , izp) + cv(ixp, iy , izp) &
              + cv(ix, iyp, izp) + cv(ixp, iyp, izp))

          ! Get the vertex magnetic field
          bxv = (bx(ix, iy, iz) + bx(ix, iyp, iz) + bx(ix, iy, izp) &
              + bx(ix, iyp, izp)) / 4.0_num
          byv = (by(ix, iy, iz) + by(ixp, iy, iz) + by(ix, iy, izp) &
              + by(ixp, iy, izp)) / 4.0_num
          bzv = (bz(ix, iy, iz) + bz(ixp, iy, iz) + bz(ix, iyp, iz) &
              + bz(ixp, iyp, iz)) / 4.0_num
          bfieldsq = bxv**2 + byv**2 + bzv**2

          t_v = 0.0_num
          ! Get the vertex temperature
          DO izp = iz, iz + 1
            DO iyp = iy, iy + 1
              DO ixp = ix, ix + 1
                CALL get_temp(rho(ixp, iyp, izp), energy(ixp, iyp, izp), &
                    eos_number, ixp, iyp, izp, T)
                t_v = t_v + T
              END DO
            END DO
          END DO
          t_v = t_v / 8.0_num

          xi_v = get_neutral(t_v, rho_v)

          f = MAX(1.0_num - xi_v, none_zero)
          IF (f .GT. 0) THEN
            eta_perp(ix, iy, iz) = eta_bar * xi_v / f * bfieldsq &
                / rho_v**2 / SQRT(t_v)
          ELSE
            eta_perp(ix, iy, iz) = 0.0_num
          END IF

!          eta_perp(ix, iy, iz) = MIN(eta_perp(ix, iy, iz), 0.0001_num)
        END DO
      END DO
    END DO

  END SUBROUTINE perpendicular_resistivity



  FUNCTION get_neutral(t_v, rho_v)

    REAL(num), INTENT(IN) :: t_v, rho_v
    REAL(num) :: get_neutral
    REAL(num) :: bof, r

    bof = 1.0_num / (f_bar * tr * SQRT(t_v)) &
        * EXP((0.25_num * (t_v / tr - 1.0_num) + 1.0_num) &
        * T_bar / t_v)
    r = 0.5_num * (-1.0_num + SQRT(1.0_num + r_bar * rho_v * bof))
    get_neutral = r / (1.0_num + r)

  END FUNCTION  get_neutral



  SUBROUTINE neutral_fraction(material)

    INTEGER, INTENT(IN) :: material
    REAL(num) :: t, rho0, e0, dx, x
    REAL(num), DIMENSION(2) :: ta, fa, xi_a
    REAL(num) :: ionise_pot_local
    INTEGER :: loop

    IF (material .EQ. EOS_ION) THEN
      ionise_pot_local = ionise_pot
    ELSE
      ionise_pot_local = 0.0_num
    END IF

    ! variable bof is b / f in the original version
    DO iz = -1, nz+2
      DO iy = -1, ny+2
        DO ix = -1, nx+2
          rho0 = rho(ix, iy, iz)
          e0 = energy(ix, iy, iz)
          ta = (gamma - 1.0_num) &
              * (/ MAX((e0 - ionise_pot_local) / 2.0_num, none_zero), e0 /)

          IF (ta(1) > ta(2)) THEN
            PRINT *, "Temperature bounds problem", ta
            STOP
          END IF

          dx = ta(2) - ta(1)
          T = ta(1)

          DO loop = 1, 100
            dx = dx / 2.0_num
            x = T  + dx
            xi_a(1) = get_neutral(x, rho0) 
            fa(1) = x - (gamma - 1.0_num) * (e0 &
                - (1.0_num - xi_a(1)) * ionise_pot_local) / (2.0_num - xi_a(1))
            IF (fa(1) <= 0.0_num) T = x
            IF (ABS(dx) < 1.e-8_num .OR. fa(1) == 0.0_num) EXIT
          END DO

          xi_n(ix, iy, iz) = get_neutral(t, rho0) 
        END DO
      END DO
    END DO

  END SUBROUTINE neutral_fraction



  SUBROUTINE get_energy(rho_in, temp_in, m_in, en_out)
    ! this routine is only used by initial conditions and must be in SI
    ! with variables as defined with kb etc. The routine needs to be here to call 
    ! get_neutral above
    REAL(num), INTENT(IN) :: rho_in, temp_in
    INTEGER, INTENT(IN) :: m_in
    REAL(num), INTENT(OUT) :: en_out
    REAL(num) :: xi_local, bof, r
  
    IF (m_in .EQ. EOS_IDEAL) THEN
      en_out = temp_in * kb / ((gamma - 1.0_num) * mbar * reduced_mass)
      RETURN
    END IF
  
    IF (m_in .EQ. EOS_PI) THEN
      ! Since we can't guarantee that the ionisation fraction already
      ! calculated is correct here, calculate it straight from the temperature
      xi_local = get_neutral(temp_in, rho_in)  
      en_out = (kb * temp_in * (2.0_num - xi_local)) &
          / (mbar * (gamma - 1.0_num))
      RETURN
    END IF
  
    IF (m_in .EQ. EOS_ION) THEN
      ! Since we can't guarantee that the ionisation fraction already
      ! calculated is correct here, calculate it straight from the temperature
      xi_local = get_neutral(temp_in, rho_in)    
      en_out = (kb * temp_in * (2.0_num - xi_local) &
          + (1.0_num - xi_local) * ionise_pot * (gamma - 1.0_num)) &
          / (mbar * (gamma - 1.0_num))
      RETURN
    END IF
  
  END SUBROUTINE get_energy



  SUBROUTINE newton_relax

!!$    INTEGER, DIMENSION(1) :: ref_index, z0(1) = 1
!!$    LOGICAL :: first_call = .TRUE., run_loop = .TRUE.
!!$
!!$    ! This should only be run above the photosphere so first call sets up
!!$    ! the lowest value of iz to use if at all, the -2 is due to zc starting
!!$    ! at -1
!!$    IF (first_call) THEN
!!$      z0 = MINLOC(ABS(zc - 0.0_num)) - 2
!!$      ! This process doesn't have any cells in the corona
!!$      IF (z0(1) > nz) run_loop = .FALSE.
!!$      IF (z0(1) < 1) z0(1) = 1 ! only need to run over the internal domain
!!$      first_call = .FALSE.
!!$    END IF
!!$
!!$    ! For every point need to find the reference density value and hence
!!$    ! the tau and temperature
!!$    IF (run_loop) THEN
!!$      DO iz = z0(1), nz
!!$        DO iy = 1, ny
!!$          DO ix = 1, nx
!!$            ! the 2 is subtracted due to rho_ref starting at -1
!!$            ref_index = MINLOC(ABS(rho(ix, iy) - rho_ref)) - 2
!!$            energy(ix, iy) = (energy(ix, iy) + dt / tau_ref(ref_index(1)) * &
!!$                T_ref(ref_index(1)) / (gamma - 1.0_num)) &
!!$                / (1.0_num + dt / tau_ref(ref_index(1)))
!!$          END DO
!!$        END DO
!!$      END DO
!!$    END IF
!!$
!!$    CALL energy_bcs

  END SUBROUTINE newton_relax

END MODULE neutral
