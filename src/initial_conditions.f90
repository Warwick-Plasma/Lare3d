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

    REAL(num) :: amp, wptb

    vx(:,:,:) = 0.0_num
    vy(:,:,:) = 0.0_num
    vz(:,:,:) = 0.0_num
    bx(:,:,:) = 0.0_num
    by(:,:,:) = 0.0_num
    bz(:,:,:) = 0.0_num

    rho(:,:,:) = 1.0_num
    energy(:,:,:) = 0.1_num
    grav(:) = 0.0_num

    amp = 0.01_num
    wptb = 20.0_num

    DO iz = -2, nz + 2
        DO iy = -2, ny + 2
          DO ix = -2,nx + 2
            vz(ix,iy,iz) = 10.0_num * EXP(-((xb(ix)-50.0_num)**2+(yb(iy)-50.0_num)**2 &
                +(zb(iz)-30.0_num)**2) / 20.0_num)
          END DO
        END DO
    END DO

  END SUBROUTINE set_initial_conditions

END MODULE initial_conditions
