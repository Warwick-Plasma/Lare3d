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
  
    INTEGER :: ix, iy, iz, loop
    REAL(num) :: rc, x1, y1, b_theta, amp, k, r0, a=1.1_num, r1, mu, m=1.5
    REAL(num):: b0, bz0, t_ph=1.0_num, t_cor=150.0_num, z_cor=25.0_num, wtr=5.0_num
    REAL(num) :: dg, w, q, lambda, r, bphi, b1, p0, p1, rho1, r_a, a1, a2, b 
    REAL(num) :: maxerr, xi_v
    REAL(num), DIMENSION(:), ALLOCATABLE :: dzb_global, dzc_global,zc_global,grav_global 
    
    REAL(num), DIMENSION(:), ALLOCATABLE :: rho_ref, energy_ref, t_ref, mu_m
    
    ALLOCATE(dzb_global(-1:nz_global+1), dzc_global(-1:nz_global), zc_global(-1:nz_global+1))
    ALLOCATE(grav_global(-1:nz_global+2), mu_m(-1:nz_global+2))
    ALLOCATE(rho_ref(-1:nz_global+2),energy_ref(-1:nz_global+2), t_ref(-1:nz_global+2))
    
    ! Changed from James's version to remove the need for a seperate routine setting up
    ! the newton cooling arrays and the MPI calls.
    
    ! Set up the initial 1D hydrostatic equilibrium
    
    grav_global = 0.9727_num
    
    ! example of lowering g to zero in corona
    a1 = 60.0_num
    a2 = 80.0_num
    WHERE (zb_global > a1) grav_global = grav_global(0)*(1.0_num+COS(pi*(zb_global-a1)/(a2-a1))) &
         /2.0_num
    WHERE (zb_global > a2) grav_global = 0.0_num
    
    !Y.Fan atmosphere temp profile
    DO iz = -1, nz_global+1 ! needs to be +1 for the dzc calculation
       zc_global(iz) = 0.5_num * (zb_global(iz) + zb_global(iz-1))
       IF (zc_global(iz) < 0.0_num) THEN
          t_ref(iz) = t_ph - (t_ph * a * zc_global(iz) * grav_global(iz) / (m+1.0_num))
       ELSE
          t_ref(iz) = t_ph + ((t_cor-t_ph) * 0.5_num * (TANH((zc_global(iz)-z_cor)/wtr)+1.0_num))
       ENDIF
    ENDDO
    t_ref(nz_global+2) = t_ref(nz_global+1)
  
    !solve HS eqn to get rho profile
    !density at bottom of domain
    rho_ref = 1.0_num
    
    DO iz = -1, nz_global
       dzb_global(iz) = zb_global(iz) - zb_global(iz-1)
       dzc_global(iz) = zc_global(iz+1) - zc_global(iz)
    ENDDO
                        
    !solve for density
    mu_m = 0.5_num    ! the reduced mass in units of proton mass
    IF (include_neutrals) xi_n = 0.0_num 
    DO loop = 1, 100
       maxerr = 0.0_num    
      DO iz = nz_global, 0, -1
        IF (zc_global(iz) < 0.0_num) THEN
          dg = 1.0_num / (dzb_global(iz)+dzb_global(iz-1))
          rho_ref(iz-1) = rho_ref(iz) * (T_ref(iz)/dzc_global(iz-1)/mu_m(iz)& 
               +grav_global(iz-1)*dzb_global(iz)*dg)
          rho_ref(iz-1) = rho_ref(iz-1) / (T_ref(iz-1)/dzc_global(iz-1)/mu_m(iz-1) & 
               -grav_global(iz-1)*dzb_global(iz-1)*dg)
        END IF
      END DO   
      !Now move from the photosphere up to the corona
      DO iz = 0, nz_global
        IF (zc_global(iz) > 0.0_num) THEN
          dg = 1.0_num / (dzb_global(iz)+dzb_global(iz-1))
          rho_ref(iz) = rho_ref(iz-1) * (T_ref(iz-1)/dzc_global(iz-1)/mu_m(iz-1)& 
               -grav_global(iz-1)*dzb_global(iz-1)*dg)
          rho_ref(iz) = rho_ref(iz) / (T_ref(iz)/dzc_global(iz-1)/mu_m(iz) & 
               +grav_global(iz-1)*dzb_global(iz)*dg)
        END IF
      END DO
      IF (include_neutrals) THEN     !note this always assumes EOS_PI
          DO iz = 0, nz_global 
             xi_v = get_neutral(t_ref(iz), rho_ref(iz)) 
             r1 = mu_m(iz)
             mu_m(iz) = 1.0_num / (2.0_num - xi_v)  
             maxerr = MAX(maxerr, ABS(mu_m(iz) - r1))
          END DO
       END IF 
       IF (maxerr < 1.e-10_num) EXIT
    END DO
    rho_ref(nz_global+1:nz_global+2) = rho_ref(nz_global)
    ! set the relaxation rate, James used an exponent of -1.67, here I use the value
    ! in the 2D paper of -1.67. The tau equation is scaled by 1/t0 * 0.1
    
    ! Convert into the local 3D arrays
    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num
    
    grav = grav_global(coordinates(1)*nz-1:coordinates(1)*nz+nz+2)
    
    DO iy = -1, ny + 2
       DO ix = -1, nx + 2
          !store temperature in energy array for a few lines    
          energy(ix,iy,:) = t_ref(coordinates(1)*nz-1:coordinates(1)*nz+nz+2) 
          rho(ix,iy,:) = rho_ref(coordinates(1)*nz-1:coordinates(1)*nz+nz+2)
       ENDDO
    ENDDO  
    
    DO iz = -1, nz + 2
      DO iy = -1, ny + 2
        DO ix = -1, nx + 2                
          r1 = energy(ix,iy,iz)
          CALL get_energy(rho(ix,iy,iz), r1, eos_number, ix, iy, iz, energy(ix,iy,iz))
        END DO
      END DO
    END DO   
    
    !add magnetic flux tube at (0,0,-10) and change pressure, dens, energy over it
    !grad(p1) matches lorentz force
    !MEQ at end of tube
    !pressure match at apex
    !seee Fan 2001 for details of initialisation
    !w = width of tube
    w = 2.0_num
    q = -(1.0_num/w)
    b0 = 5.0_num
    lambda = 20.0_num
    DO ix = -1,nx+2
       DO iy = -1,ny+2
          DO iz = -1,nz+2
             !define bx,by,bz at correct points
             r = SQRT(yc(iy)**2 + (zc(iz)+10.0_num)**2)
             bx(ix,iy,iz) =  b0 * EXP(-(r/w)**2)
             bphi =  bx(ix,iy,iz) * q * r
    
             r = SQRT(yb(iy)**2 + (zc(iz)+10.0_num)**2)
             b1 = b0 * EXP(-(r/w)**2)
             by(ix,iy,iz) = -b1 * q * (zc(iz)+10.0_num)
    
             r = SQRT(yc(iy)**2 + (zb(iz)+10.0_num)**2)
             b1 =  b0 * EXP(-(r/w)**2)
             bz(ix,iy,iz) =  b1 * q * yc(iy)
    
             !define gas pressure and magnetic pressure
             p0 =  rho(ix,iy,iz)*energy(ix,iy,iz)*(gamma-1.0_num)
             p1 =  -0.25_num * bx(ix,iy,iz)**2 - 0.5_num * bphi**2
             !change density and energy
             r1 = xc(ix)/lambda
             rho1 =  (p1/p0)*rho(ix,iy,iz)*EXP(-(r1**2))
             rho(ix,iy,iz)   =  rho(ix,iy,iz) + rho1
             energy(ix,iy,iz)= (p0 + p1) / (rho(ix,iy,iz) * (gamma - 1.0_num))
    
          END DO
       END DO
    END DO              
    
    DEALLOCATE(dzb_global, dzc_global, zc_global)
    DEALLOCATE(grav_global, mu_m)
    DEALLOCATE(rho_ref, energy_ref, t_ref)
  
  END SUBROUTINE set_initial_conditions

 
END MODULE initial_conditions
