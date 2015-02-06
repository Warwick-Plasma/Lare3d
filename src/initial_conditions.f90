MODULE initial_conditions

  USE shared_data
  USE neutral

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_initial_conditions

CONTAINS

  !****************************************************************************
  ! This function sets up the initial condition for the code
  ! The variables which must be set are:
  !   rho - density
  !   v{x,y,z} - Velocities in x, y, z
  !   b{x,y,z} - Magnetic fields in x, y, z
  !   energy - Specific internal energy
  !
  ! You may also need the neutral fraction. This can be calculated by a
  ! function call to get_neutral(temperature, rho, z). This routine is in
  ! core/neutral.f90 and requires the local temperature and mass density.
  ! For example to set xi_n to the neutral fraction use:
  !   xi_n = get_neutral(temperature, rho, z)
  !****************************************************************************

  SUBROUTINE set_initial_conditions

    ! This is about the most complicated example for initial conditions
    ! used here as it covers including gravity and neutrals.
    ! The normalisation assumed is that from the defauls control.f90

    INTEGER :: loop
    INTEGER :: ix, iy, iz, iz1
    REAL(num) :: a1, a2, dg
    REAL(num) :: a = 1.0_num, Tph = 9.8_num, Tcor = 980.0_num, ycor = 11.0_num
    REAL(num) :: betafs = 0.25_num, yfsl = -10.0_num, yfsu = -1.0_num
    REAL(num) :: wtr = 0.6_num, wfsl = 0.5_num, wfsu = 0.5_num
    REAL(num) :: r1, maxerr, xi_v
    REAL(num) :: amp, wptb, fac, theta
    REAL(num), DIMENSION(:), ALLOCATABLE :: zc_global, dzb_global, dzc_global
    REAL(num), DIMENSION(:), ALLOCATABLE :: grav_ref, temp_ref, rho_ref
    REAL(num), DIMENSION(:), ALLOCATABLE :: beta_ref, mag_ref, mu_m

    ALLOCATE( zc_global(-1:nz_global+1))
    ALLOCATE(dzb_global(-1:nz_global+1))
    ALLOCATE(dzc_global(-1:nz_global+1))
    ALLOCATE(  grav_ref(-1:nz_global+2))
    ALLOCATE(  temp_ref(-1:nz_global+2))
    ALLOCATE(   rho_ref(-1:nz_global+2))
    ALLOCATE(   mag_ref(-1:nz_global+2))
    ALLOCATE(  beta_ref(-1:nz_global+2))
    ALLOCATE(      mu_m(-1:nz_global+2))

    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num
    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num

    ! Fill in zc_global with the positions central to the zb_global points
    DO iz = -1, nz_global + 1
      zc_global(iz) = 0.5_num * (zb_global(iz-1) + zb_global(iz))
    END DO

    ! Fill in dzb_global and dzc_global
    DO iz = -1, nz_global
      dzb_global(iz) = zb_global(iz) - zb_global(iz-1)
      dzc_global(iz) = zc_global(iz+1) - zc_global(iz)
    END DO

    ! Fill in the reference gravity array - lowering grav to zero at the top
    ! of the corona smoothly from a1 to grav = 0 at a2 and above
    grav_ref = 11.78_num
    a1 = zb_global(nz_global) - 20.0_num
    a2 = zb_global(nz_global) - 5.0_num
    DO iz = 0, nz_global + 2
      IF (zb_global(iz) > a1) THEN
        grav_ref(iz) = 11.78_num &
            * (1.0_num + COS(pi * (zb_global(iz) - a1) / (a2 - a1))) / 2.0_num
      END IF
      IF (zb_global(iz) > a2) THEN
        grav_ref(iz) = 0.0_num
      END IF
    END DO
    grav_ref(-1) = grav_ref(0)
    grav_ref(nz_global+1:nz_global+2) = grav_ref(nz_global)

    ! Beta profile from Archontis 2009 but in 2D
    ! Similar to that of Nozawa 1991
    ! NB : The variable beta used here is actually 1/beta
    beta_ref = 0.0_num
    DO iz = -1, nz_global + 1
      IF ((zc_global(iz) > yfsl) .AND. (zc_global(iz) < yfsu)) THEN
        beta_ref(iz) = betafs &
            * (0.5_num * (TANH((zc_global(iz) - yfsl) / wfsl) + 1.0_num)) &
            * (0.5_num * (1.0_num - TANH((zc_global(iz) - yfsu) / wfsu)))
      END IF
    END DO

    ! Calculate the density profile, starting from the refence density at the
    ! photosphere and calculating up and down from there including beta
    rho_ref = 1.0_num
    mu_m = 1.0_num
    IF (eos_number == EOS_IDEAL .AND. (.NOT. neutral_gas)) mu_m = 0.5_num

    DO loop = 1, 1
      maxerr = 0.0_num
      ! Go from photosphere down
      DO iz = -1, nz_global + 1
        IF (zc_global(iz) < 0.0_num) THEN
          temp_ref(iz) = Tph - a * (gamma - 1.0_num) &
              * zc_global(iz) * grav_ref(iz) * mu_m(iz) / gamma
        END IF
        IF (zc_global(iz) >= 0.0_num) THEN
          a1 = 0.5_num * (TANH((zc_global(iz) - ycor) / wtr) + 1.0_num)
          temp_ref(iz) = Tph - 1.0 + (Tcor - Tph)**a1
        END IF
      END DO
      temp_ref(nz_global+1:nz_global+2) = temp_ref(nz_global)

      DO iz = nz_global, 0, -1
        IF (zc_global(iz) < 0.0_num) THEN
          izm = iz - 1
          dg = 1.0_num / (dzb_global(iz) + dzb_global(izm))

          rho_ref(izm) = rho_ref(iz ) * (temp_ref(iz ) &
              * (1.0_num + beta_ref(iz )) / dzc_global(izm) / mu_m(iz ) &
              + grav_ref(izm) * dzb_global(iz ) * dg)

          rho_ref(izm) = rho_ref(izm) / (temp_ref(izm) &
              * (1.0_num + beta_ref(izm)) / dzc_global(izm) / mu_m(izm) &
              - grav_ref(izm) * dzb_global(izm) * dg)
        END IF
      END DO

      ! Now move from the photosphere up to the corona
      DO iz = 0, nz_global
        IF (zc_global(iz) >= 0.0_num) THEN
          izm = iz - 1
          dg = 1.0_num / (dzb_global(iz) + dzb_global(izm))

          rho_ref(iz)  = rho_ref(izm) * (temp_ref(izm) &
              * (1.0_num + beta_ref(izm)) / dzc_global(izm) / mu_m(izm) &
              - grav_ref(izm) * dzb_global(izm) * dg)

          rho_ref(iz)  = rho_ref(iz ) / (temp_ref(iz ) &
              * (1.0_num + beta_ref(iz )) / dzc_global(izm) / mu_m(iz ) &
              + grav_ref(izm) * dzb_global(iz ) * dg)
        END IF
      END DO

      IF (eos_number /= EOS_IDEAL) THEN
        DO iz = 0, nz_global
          xi_v = get_neutral(temp_ref(iz), rho_ref(iz), zb(iz))
          r1 = mu_m(iz)
          mu_m(iz) = 1.0_num / (2.0_num - xi_v)
          maxerr = MAX(maxerr, ABS(mu_m(iz) - r1))
        END DO
      END IF

      IF (maxerr < 1.e-16_num) EXIT
    END DO

    rho_ref(nz_global+1:nz_global+2) = rho_ref(nz_global)

    ! Magnetic flux sheet profile from Archontis2009
    ! Similar structure to the 2D version used in Nozawa1991 and Isobe2006
    DO iz= -1, nz_global + 2
      mag_ref(iz) = SQRT(2.0_num * beta_ref(iz) * temp_ref(iz) * rho_ref(iz) &
          / mu_m(iz))
    END DO

    ! Fill in all the final arrays from the ref arrays
    iz1 = coordinates(1) * nz - 1

    DO iz = -1, nz + 2
      grav(iz) = grav_ref(iz1)
      DO iy = -1, ny + 2
        DO ix = -1, nx + 2
          rho(ix,iy,iz) = rho_ref(iz1)
          energy(ix,iy,iz) = temp_ref(iz1)
          bx(ix,iy,iz) = mag_ref(iz1)

          IF (eos_number /= EOS_IDEAL) THEN
            xi_v = get_neutral(energy(ix,iy,iz), rho(ix,iy,iz), zb(iz))
          ELSE
            IF (neutral_gas) THEN
              xi_v = 1.0_num
            ELSE
              xi_v = 0.0_num
            END IF
          END IF

          energy(ix,iy,iz) = (energy(ix,iy,iz) * (2.0_num - xi_v) &
              + (1.0_num - xi_v) * ionise_pot * (gamma - 1.0_num)) &
              / (gamma - 1.0_num)
        END DO
      END DO
      iz1 = iz1 + 1
    END DO

    amp = 0.01_num
    wptb = 20.0_num

    DO iz = 1, nz
      IF (zc_global(iz) > -10.0_num .AND. zc_global(iz) < -1.0_num) THEN
        fac = 0.25_num * amp
        theta = 2.0_num * pi / wptb
        DO iy = 1, ny
          DO ix = 1,nx
            vz(ix,iy,iz) = fac * COS(theta * xb(ix)) * COS(theta * yb(iy)) &
              * (TANH(2.0_num * (zb(iz) - yfsl)) &
              -  TANH(2.0_num * (zb(iz) - yfsu)))
          END DO
        END DO
      END IF
    END DO

    DEALLOCATE(zc_global, dzb_global, dzc_global, mu_m)
    DEALLOCATE(grav_ref, temp_ref, rho_ref, beta_ref, mag_ref)

  END SUBROUTINE set_initial_conditions

END MODULE initial_conditions
