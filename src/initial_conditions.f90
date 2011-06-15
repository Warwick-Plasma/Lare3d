MODULE initial_conditions

  USE shared_data
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
  ! The syntax for this function (which is in core/neutral.f90) is
  !
  ! CALL get_energy(density, temperature, equation_of_state, &
  !     output_energy) 
  !
  ! REAL(num) :: density - The density at point (ix, iy) on the grid
  ! REAL(num) :: temperature - The temperature at point (ix, iy) on the grid
  ! INTEGER :: equation_of_state - The code for the equation of state to use.
  !            The global equation of state for the code is eos_number
  ! REAL(num) :: output_energy - The specific internal energy returned by
  !              the routine
  !
  ! You may also need the neutral fraction. This can be calculated by a function
  ! call to  get_neutral(temperature, rho, z). This routine is in core/neutral.f90
  ! and requires the local temperature and mass density. For example to set
  ! xi_n to the neutral fraction use
  ! xi_n = get_neutral(temperature, rho, z)
  !---------------------------------------------------------------------------
  
  SUBROUTINE set_initial_conditions

    INTEGER :: loop
    INTEGER :: ix, iy, iz
    REAL(num) :: a1, a2, dg
    REAL(num) :: a=1.0_num, Tph=9.8_num, Tcor=980.0_num, ycor=11.0_num, wtr=0.6_num
    REAL(num) :: betafs=0.25_num, yfsl=-10.0_num, yfsu=-1.0_num, wfsl=0.5_num, wfsu=0.5_num
    REAL(num) :: r1, maxerr, xi_v
    REAL(num) :: amp, wptb, yptb1, yptb2, yptb
    REAL(num), DIMENSION(:), ALLOCATABLE :: zc_global, dzb_global, dzc_global
    REAL(num), DIMENSION(:), ALLOCATABLE :: grav_ref, temp_ref, rho_ref
    REAL(num), DIMENSION(:), ALLOCATABLE :: beta_ref, mag_ref, mu_m
    
    ALLOCATE(zc_global(-1:nz_global+1))
    ALLOCATE(dzb_global(-1:nz_global+1), dzc_global(-1:nz_global))
    ALLOCATE(grav_ref(-1:nz_global+2), temp_ref(-1:nz_global+2))
    ALLOCATE(rho_ref(-1:nz_global+2), mag_ref(-1:nz_global+2))
    ALLOCATE(beta_ref(-1:nz_global+2), mu_m(-1:nz_global+2))
    
    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num
    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num
    
    !fill in zc_global with the positions central to the zb_global points
    DO iz = -1, nz_global+1
       zc_global(iz) = 0.5_num * (zb_global(iz-1) + zb_global(iz))
    END DO
    
    !fill in dzb_global and dzc_global
    DO iz = -1, nz_global
       dzb_global(iz) = zb_global(iz) - zb_global(iz-1)
       dzc_global(iz) = zc_global(iz+1) - zc_global(iz)
    END DO
    
    !fill in the reference gravity array - lowering grav to zero at the top 
    !of the corona smoothly from a1 to grav=0 at a2 and above
    grav_ref = 11.78_num
    a1 = zb_global(nz_global) - 20.0_num
    a2 = zb_global(nz_global) - 5.0_num
    DO iz = 0,nz_global+2
       IF (zb_global(iz) > a1) THEN
          grav_ref(iz) = 11.78_num * (1.0_num + COS(pi * (zb_global(iz) - a1) &
               / (a2-a1))) / 2.0_num
       END IF
       IF (zb_global(iz) > a2) THEN
          grav_ref(iz) = 0.0_num
       END IF
    END DO    
    grav_ref(-1) = grav_ref(0)
    grav_ref(nz_global+1:nz_global+2) = grav_ref(nz_global)
    
    !beta profile from Archontis 2009 but in 2D
    !similar to that of Nozawa 1991
    !NB : The variable beta used here is actually 1/beta
    beta_ref = 0.0_num
    DO iz = -1, nz_global+1
      IF ((zc_global(iz) .GT. yfsl) .AND. (zc_global(iz) .LT. yfsu)) THEN
        beta_ref(iz) = betafs * &
              (0.5_num * (TANH((zc_global(iz) - yfsl) / wfsl) + 1.0_num)) * &
              (0.5_num * (1.0_num - TANH((zc_global(iz) - yfsu) / wfsu)))
      END IF
    END DO
    
    !calculate the density profile, starting from the refence density at the
    !photosphere and calculating up and down from there including beta
    rho_ref = 1.0_num
    mu_m = 1.0_num
    IF (eos_number == EOS_IDEAL .AND. (.NOT. neutral_gas)) mu_m = 0.5_num
    DO loop = 1,1
       maxerr = 0.0_num
       !Go from photosphere down
       DO iz = -1, nz_global+1
          IF (zc_global(iz) < 0.0_num) THEN
             temp_ref(iz) = Tph - a * (gamma - 1.0_num) &
                  * zc_global(iz) * grav_ref(iz) * mu_m(iz) / gamma 
          END IF
          IF (zc_global(iz) >= 0.0_num) THEN
             temp_ref(iz) = Tph + (Tcor - Tph) * 0.5_num &
                  * (TANH((zc_global(iz) - ycor) / wtr) + 1.0_num)
          END IF
       END DO
       temp_ref(nz_global+1:nz_global+2) = temp_ref(nz_global)
       
       DO iz = nz_global,0,-1
          IF (zc_global(iz) < 0.0_num) THEN  
             dg = 1.0_num / (dzb_global(iz) + dzb_global(iz-1))
             rho_ref(iz-1) = rho_ref(iz) * (temp_ref(iz)*(1.0_num+beta_ref(iz)) &
                  /dzc_global(iz-1)/mu_m(iz)+grav_ref(iz-1)*dzb_global(iz)*dg)
             rho_ref(iz-1) = rho_ref(iz-1) / (temp_ref(iz-1)*(1.0_num+beta_ref(iz-1)) &
                  /dzc_global(iz-1)/mu_m(iz-1)-grav_ref(iz-1)*dzb_global(iz-1)*dg)
          END IF
       END DO
       !Now move from the photosphere up to the corona
       DO iz = 0, nz_global
          IF (zc_global(iz) >= 0.0_num) THEN
             dg = 1.0_num / (dzb_global(iz)+dzb_global(iz-1))
             rho_ref(iz) = rho_ref(iz-1) * (temp_ref(iz-1)*(1.0_num+beta_ref(iz-1)) &
                  /dzc_global(iz-1)/mu_m(iz-1)-grav_ref(iz-1)*dzb_global(iz-1)*dg)
             rho_ref(iz) = rho_ref(iz) / (temp_ref(iz)*(1.0_num+beta_ref(iz)) &
                  /dzc_global(iz-1)/mu_m(iz)+grav_ref(iz-1)*dzb_global(iz)*dg)
          END IF
       END DO
       IF (eos_number /= EOS_IDEAL) THEN
          DO iz = 0, nz_global
             xi_v = get_neutral(temp_ref(iz),rho_ref(iz),zb(iz))
             r1 = mu_m(iz)
             mu_m(iz) = 1.0_num / (2.0_num - xi_v)
             maxerr = MAX(maxerr, ABS(mu_m(iz) - r1))
          END DO
       END IF
       IF (maxerr < 1.e-16_num) EXIT
    END DO 
    
    rho_ref(nz_global+1:nz_global+2) = rho_ref(nz_global)
                                  
    !magnetic flux sheet profile from Archontis2009
    !similar structure to the 2D version used in Nozawa1991 and Isobe2006
    DO iz= -1,nz_global+2,1
       mag_ref(iz) = SQRT(2.0_num * beta_ref(iz) * temp_ref(iz) * rho_ref(iz) / mu_m(iz))
    END DO
    
    !fill in all the final arrays from the ref arrays
    grav(:) = grav_ref(coordinates(1)*nz-1:coordinates(1)*nz+nz+2)
    DO ix = -1, nx+2
    DO iy = -1, ny+2
       rho(ix,iy,:) = rho_ref(coordinates(1)*nz-1:coordinates(1)*nz+nz+2)
       energy(ix,iy,:) = temp_ref(coordinates(1)*nz-1:coordinates(1)*nz+nz+2)
       bx(ix,iy,:) = mag_ref(coordinates(1)*nz-1:coordinates(1)*nz+nz+2)
    END DO
    END DO
    DO ix = -1,nx+2
       DO iy = -1,ny+2 
         DO iz = -1, nz+2
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
    END DO
    
    amp = 0.01_num
    wptb = 20.0_num
    
    DO iz = 1, nz
      IF ((zc_global(iz) .GT. -10.0_num) .AND. (zc_global(iz) .LT. -1.0_num)) THEN
        DO iy = 1, ny
          DO ix = 1,nx
            vz(ix,iy,iz) = (amp / 4.0_num) * COS(2.0_num*pi*xb(ix)/wptb) & 
              * COS(2.0_num*pi*yb(iy)/wptb)  &
              * (TANH((zb(iz)-yfsl)/0.5_num)-TANH((zb(iz)-yfsu)/0.5_num))
          END DO 
        END DO
      END IF
    END DO  
    
    
    DEALLOCATE(zc_global, dzb_global, dzc_global, mu_m)
    DEALLOCATE(grav_ref, temp_ref, rho_ref, beta_ref, mag_ref)
    

  END SUBROUTINE set_initial_conditions


END MODULE initial_conditions
