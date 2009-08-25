MODULE initial_conditions

  USE shared_data
  USE eos
  USE neutral

  IMPLICIT NONE

  PRIVATE

  PUBLIC:: set_initial_conditions

CONTAINS

  !---------------------------------------------------------------------------
  ! This function sets up the initial condition for the code
  ! The variables which must be set are
  ! Rho - density
  ! V{x, y, z} - Velocities in x, y, z
  ! B{x, y, z} - Magnetic fields in x, y, z
  ! Energy - Specific internal energy
  ! Since temperature is a more intuitive quantity than specific internal energy
  ! There is a helper function get_energy which converts temperature to energy
  ! The syntax for this function is
  !
  ! CALL get_energy(density, temperature, equation_of_state, ix, iy, &
  !     output_energy)
  !
  ! REAL(num) :: density - The density at point (ix, iy) on the grid
  ! REAL(num) :: temperature - The temperature at point (ix, iy) on the grid
  ! INTEGER :: equation_of_state - The code for the equation of state to use.
  !            The global equation of state for the code is eos_number
  ! INTEGER :: ix - The current gridpoint in the x direction
  ! INTEGER :: iy - The current gridpoint in the y direction
  ! REAL(num) :: output_energy - The specific internal energy returned by
  !              the routine
  !---------------------------------------------------------------------------
  SUBROUTINE set_initial_conditions

    INTEGER:: ix, iy, iz
    REAL(num) :: T

    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num
    by = 0.0_num
    bx = 0.0_num
    bz = 0.0_num
    energy = 1.0_num
    rho = 1.0_num

    grav = 0.0_num

    rho = 1.0e-6_num

    DO iz = 0, nz+1
      DO iy = 0, ny+1
        DO ix = 0, nx+1
          T = 1.0e3_num + 1.0e6_num * EXP(-(xc(ix) / (length_x / 16.0_num))**2 &
              - (yc(iy) / (length_y / 16.0_num))**2 &
              - (zc(iz) / (length_z / 16.0_num))**2)
          CALL get_energy(rho(ix, iy, iz), T, EOS_IDEAL, ix, iy, iz, &
              energy(ix, iy, iz))
        END DO
      END DO
    END DO

  END SUBROUTINE set_initial_conditions

END MODULE initial_conditions
