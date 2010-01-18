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
      REAL(num) :: hh
      REAL(num) :: alpha1, alpha2, B1 
      REAL(num) :: rc, rb, rbx, rby, b_theta, b_z, delta, B2, C2
      REAL(num) :: k, amp, dx, dy, theta, v_perp, v_r, v_theta
      REAL(num) :: costh, sinth, coskz, sinkz, arg

      vx = 0.0_num
      vy = 0.0_num
      vz = 0.0_num

  ! Case 10
      alpha1 = 1.8_num
      alpha2 = alpha1*alpha1
      B1 = 1.0_num

  ! helicity
      hh = 0.0_num

      bx = 0.0_num
      by = 0.0_num
      bz = 0.0_num

      energy = 0.0_num
      rho = 1.0_num
      grav = 0.0_num

  ! Velocity perturbation
  !Case 3 perturbation
      k = 6.0_num*pi/20.0_num
      amp = 1.e-4_num
      dx = 0.1_num * dxb(nx/2)
      dy = 0.1_num * dyb(ny/2)
      DO ix = -1,nx+2
        DO iy = -1,ny+2
          DO iz = -1,nz+2

            rc  = SQRT(xc(ix)*xc(ix)+yc(iy)*yc(iy))
            rb  = SQRT(xb(ix)*xb(ix)+yb(iy)*yb(iy))
            rbx = SQRT(xb(ix)*xb(ix)+yc(iy)*yc(iy))
            rby = SQRT(xc(ix)*xc(ix)+yb(iy)*yb(iy))
            IF (rb .LE. sqrt(dx**2 + dy**2)) THEN
              costh = 1.0_num
              sinth = 0.0_num
            ELSE
              costh = xb(ix)/rb
              sinth = yb(iy)/rb
            ENDIF
            coskz = cos(k*zb(iz))
            sinkz = sin(k*zb(iz))
            v_r = EXP(-rb**2*(1+(rb/0.5_num)**6))*COS(pi*zb(iz)/length_z)*&
                  (costh*coskz + sinth*sinkz)
  !
  ! Define the velocity at the cell boundaries
  !
            IF (rb .LE. 1.0_num) THEN
              B2 = alpha2*((1.0_num-rb**2.0_num)**7.0_num-1.0_num)/7.0_num
              C2 = alpha2*rb*rb*(1.0_num-rb**2.0_num)**6.0_num
              b_z = sqrt(1.0_num + B2 - C2)
              b_theta = alpha1*rb*(1.0_num-rb*rb)**3.0_num
              v_perp = -((b_theta**2 + b_z**2)/&
               (b_z + k*rc*b_theta))*&
               (1.0_num - 2.0_num*rb**2 - 8.0_num*(rb/0.5_num)**6)*EXP(-4.0_num*rb**4)*&
               COS(pi*zb(iz)/length_z)*(sinth*coskz - costh*sinkz)

              v_theta = b_z*v_perp / (b_z**2 + b_theta**2)

              vz(ix,iy,iz) = amp*(-b_theta*v_perp / (b_z**2 +&
                    b_theta**2))
              vx(ix,iy,iz) = amp*(v_r*costh - v_theta*sinth)
              vy(ix,iy,iz) = amp*(v_r*sinth + v_theta*costh)
            ELSE
              b_z = sqrt(1.0_num - alpha2/7.0_num)
              b_theta = 0.0_num
              v_perp = -((b_theta**2 + b_z**2)/&
               (b_z + k*rb*b_theta))*&
               (1.0_num - 2.0_num*rb**2 - 8.0_num*(rb/0.5_num)**6)*EXP(-4.0_num*rb**4)*&
               COS(pi*zb(iz)/length_z)*(sinth*coskz - costh*sinkz)

              v_theta = b_z*v_perp / (b_z**2 + b_theta**2)

              vz(ix,iy,iz) = amp*(-b_theta*v_perp / (b_z**2 +&
                    b_theta**2))
              vx(ix,iy,iz) = amp*(v_r*costh - v_theta*sinth)
              vy(ix,iy,iz) = amp*(v_r*sinth + v_theta*costh)
            ENDIF
  !
  ! Define Bz on face centred at (xc,yc,zb)
  !
            IF (rc .LE. 1.0_num) THEN
              B2 = alpha2*((1.0_num-rc**2.0_num)**7.0_num-1.0_num)/7.0_num
              C2 = alpha2*rc*rc*(1.0_num-rc**2.0_num)**6.0_num
              bz(ix,iy,iz) = sqrt(1.0_num + B2 - C2)
            ELSE
              bz(ix,iy,iz) = sqrt(1.0_num - alpha2/7.0_num)
            ENDIF
  !
  ! Define Bx on face centred at (xb,yc,zc)
  !
            IF (rbx .LE. 1.0_num) THEN
              b_theta = alpha1*rbx*(1.0_num-rbx*rbx)**3.0_num
              bx(ix,iy,iz) = -b_theta * yc(iy) / rbx
            ELSE
              bx(ix,iy,iz) = 0.0_num
            ENDIF
  !
  ! Define By on face centred at (xc,yb,zc)
  !
            IF (rby .LE. 1.0_num) THEN
              b_theta = alpha1*rby*(1.0_num-rby*rby)**3.0_num
              by(ix,iy,iz) = b_theta * xc(ix) / rby
            ELSE
              by(ix,iy,iz) = 0.0_num
            ENDIF

          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE set_initial_conditions
  


END MODULE initial_conditions
