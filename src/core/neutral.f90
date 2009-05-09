! All the subroutines in this module are for the partially ionised flux
! emergence simulations; see Leake & Arber, 2006

MODULE neutral

  USE shared_data
  USE boundary
  USE eos

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: perpendicular_resistivity, newton_relax, &
      neutral_fraction, setup_neutral, get_neutral

CONTAINS

  SUBROUTINE setup_neutral
    ! Ion neutral collision cross section(m^2)
    REAL(num) :: sigma_in = 5.0e-19_num
    REAL(num) :: tr

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
    tr = 7230.85_num
    tr_bar = 1.0_num / tr

    ! Calculate fbar^(2 / 3) in (k^-1 m^-2)
    f_bar = (pi * me_0 * kb_0) / h_0**2
    f_bar = f_bar**(3.0_num / 2.0_num)

    ! Calculate tbar in (K)
    t_bar = ionise_pot / kb

    ! Calculate rbar in (kg^-1)
    r_bar = 4.0_num / MBAR

    ! Calculate eta_bar in (m^4 / (k s kg^2))
    eta_bar = (0.5_num / MBAR &
        * SQRT(16.0_num * kb_0 / (pi * MBAR)) * sigma_in)**(-1)

  END SUBROUTINE setup_neutral



  SUBROUTINE perpendicular_resistivity

    ! This subroutine calculates the cross field resistivity at the current
    ! temperature. If you're not using the "neutral_fraction" subroutine and
    ! the "Saha" equation of state, then this routine will give VERY strange
    ! results which are probably meaningless.

    ! normalising values are L0 = 150km, v0 = 6.5km / s, rho0 = 2.7e-4 kg / m3
    ! t0 = 23s, T0 = 6420K, P0 = 1.2e4 Pa, B0 = 1200G (0.12T)

    REAL(num) :: f, xi_v, bxv, byv, bzv, bfieldsq, rho_v, T_v, T
    INTEGER :: ixp, iyp, izp

    eta_perp = 0.0001
    RETURN

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

          T_v = 0.0_num
          ! Get the vertex temperature
          DO izp = iz, iz + 1
            DO iyp = iy, iy + 1
              DO ixp = ix, ix + 1
                CALL get_temp(rho(ixp, iyp, izp), energy(ixp, iyp, izp), &
                    eos_number, ixp, iyp, izp, T)
                T_v = T_v + T
              END DO
            END DO
          END DO
          T_v = T_v / 8.0_num

          xi_v = get_neutral(T_v, rho_v)

          f = MAX(1.0_num - xi_v, none_zero)
          IF (f .GT. 0) THEN
            eta_perp(ix, iy, iz) = eta_bar * xi_v / f * bfieldsq &
                / rho_v**2 / SQRT(T_v)
          ELSE
            eta_perp(ix, iy, iz) = 0.0_num
          END IF

          eta_perp(ix, iy, iz) = MIN(eta_perp(ix, iy, iz), 0.0001_num)
        END DO
      END DO
    END DO

  END SUBROUTINE perpendicular_resistivity



  FUNCTION get_neutral(T_v, rho_v)

    REAL(num), INTENT(IN) :: T_V, rho_v
    REAL(num) :: get_neutral
    REAL(num) :: f, b, r

    f = f_bar * SQRT(T_v) * EXP(-T_bar / T_v) ! T from b has been cancelled
    b = tr_bar * EXP(0.25_num * T_bar / T_v * (tr_bar * T_v - 1.0_num))
    r = 0.5_num * (-1.0_num + SQRT(1.0_num + r_bar * rho_v * b / f))
    get_neutral = r / (1.0_num + r)

  END FUNCTION  get_neutral



  SUBROUTINE neutral_fraction(material)

    INTEGER, INTENT(IN) :: material
    REAL(num) :: bof, r, T, rho0, e0, dx, x
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
              * (/ MAX((e0 - ionise_pot_local) / 2.0_num, none_zero), &
              e0 / 2.0_num /)

          IF (ta(1) > ta(2)) THEN
            PRINT *, "Temperature bounds problem", ta
            STOP
          END IF

          dx = ta(2) - ta(1)
          T = ta(1)

          DO loop = 1, 100
            dx = dx / 2.0_num
            x = T  + dx
            bof = tr_bar / (f_bar * SQRT(x)) &
                * EXP((0.25_num * (tr_bar * x - 1.0_num) + 1.0_num) * T_bar / x)
            r = 0.5_num * (-1.0_num + SQRT(1.0_num + r_bar * rho0 * bof))
            xi_a(1) = r / (1.0_num + r)
            fa(1) = x - (gamma - 1.0_num) * (e0 &
                - (1.0_num - xi_a(1)) * ionise_pot_local) / (2.0_num - xi_a(1))
            IF (fa(1) <= 0.0_num) T = x
            IF (ABS(dx) < 1.e-8_num .OR. fa(1) == 0.0_num) EXIT
          END DO

          bof = tr_bar / (f_bar * SQRT(T)) &
              * EXP((0.25_num * (tr_bar * T - 1.0_num) + 1.0_num) * T_bar / T)
          r = 0.5_num * (-1.0_num + SQRT(1.0_num + r_bar * rho0 * bof))
          xi_n(ix, iy, iz) = r / (1.0_num + r)
        END DO
      END DO
    END DO

  END SUBROUTINE neutral_fraction



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

    CALL energy_bcs

  END SUBROUTINE newton_relax

END MODULE neutral
