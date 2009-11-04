! This file contains example initial conditions used in previous simulations

! should give the model atmosphere, with bouyant flux tube, as in Arber et al. ApJ
SUBROUTINE set_initial_conditions 

  INTEGER :: ix, iy, iz, loop
  REAL(num) :: rc, x1, y1, b_theta, amp, k, r0, a=1.1_num, r1, mu, m=1.5
  REAL(num):: b0, bz0, t_ph=1.0_num, t_cor=150.0_num, z_cor=25.0_num, wtr=5.0_num
  REAL(num) :: dg, w, q, lambda, r, bphi, b1, p0, p1, rho1, r_a, a1, a2, b
  REAL(num), DIMENSION(:), ALLOCATABLE :: dzb_global, dzc_global,zc_global,grav_global
  REAL(num), DIMENSION(:), ALLOCATABLE :: mu_m, t_ref, rho_ref
  REAL(num) :: f, xi_v, maxerr, change, T

  ALLOCATE(dzb_global(-1:nz_global+1), dzc_global(-1:nz_global), zc_global(-1:nz_global+1))
  ALLOCATE(grav_global(-1:nz_global+2), mu_m(-1:nz_global+2), t_ref(-1:nz_global+2))
  ALLOCATE(rho_ref(-1:nz_global+2))

  vx = 0.0_num
  vy = 0.0_num
  vz = 0.0_num
  bx = 0.0_num
  by = 0.0_num
  bz = 0.0_num
  rho = 1.0_num
  energy = 1.0_num
  rho_ref = 1.0_num
  
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
        t_ref(iz) = t_ph - zc_global(iz) * grav_global(iz) / (m + 1.0_num)
     ELSE
        t_ref(iz) = t_ph + ((t_cor-t_ph) * 0.5_num * (TANH((zc_global(iz)-z_cor)/wtr)+1.0_num))
     ENDIF
  ENDDO
  t_ref(nz_global+2) = t_ref(nz_global+1)
      
  DO iz = -1, nz_global
     dzb_global(iz) = zb_global(iz) - zb_global(iz-1)
     dzc_global(iz) = zc_global(iz+1) - zc_global(iz)
  ENDDO

  !solve for density
  mu_m = 1.0
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
            change = ABS(mu_m(iz) - r1) 
            maxerr = MAX(maxerr, change)
         END DO
      END IF 
      IF (maxerr < 1.e-10_num) EXIT    
  END DO

  rho_ref(nz_global+1:nz_global+2) = rho_ref(nz_global)
  DO iy = -1, ny + 2
     DO ix = -1, nx + 2
        rho(ix,iy,:) = rho_ref(coordinates(1)*nz-1:coordinates(1)*nz+nz+2) 
        !store temperature in energy array for a few lines
        energy(ix,iy,:) = t_ref(coordinates(1)*nz-1:coordinates(1)*nz+nz+2)
     ENDDO
  ENDDO 

  DO iz = -1, nz + 2
    DO iy = -1, ny + 2
      DO ix = -1, nx + 2                
        t = energy(ix,iy,iz)
        CALL get_energy(rho(ix,iy,iz), t, eos_number, ix, iy, iz, energy(ix,iy,iz))
      END DO
    END DO
  END DO          
       
   grav = grav_global(coordinates(1)*nz-1:coordinates(1)*nz+nz+2)

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
           p0 =  rho(ix,iy,iz)*(gamma-1.0_num)*energy(ix,iy,iz) 
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
  DEALLOCATE(grav_global, mu_m, t_ref, rho_ref)


  END SUBROUTINE set_initial_conditions



!Kink unstable loop from Arber et al, 1999
  SUBROUTINE equilibrium  

    INTEGER :: ix, iy, iz
    REAL(num) :: rc, x1, y1, b_theta, amp, k, r0, a, r1, mu
    REAL(num):: b0, bz0

    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num

    bx = 1.0_num
    by = 0.0_num
    bz = 0.0_num
    
    r0 = 1.0_num / SQRT(6.0_num)
    a = 5.0_num / (6.0_num * SQRT(6.0_num))

    bz0 = 1.0_num
    b0 = 4.3_num

    DO ix = -1, nx+2             ! setup static equilibrium values
       DO iy = -1, ny+2
          DO iz = -1, nz+2

             rc = SQRT(xc(ix)**2 + yc(iy)**2)
             IF (rc >= 1.0_num) rc = 1.0_num
             bz(ix,iy,iz) = (1.0/2.0)*rc**2 - (3.0/8.0)*(rc**4/r0**2)    &
                  + (7.0*a/25.0)*(rc**5/r0**3) + (1.0/12.0)*(rc**6/r0**4)     &
                  - (9.0*a/70.0)*(rc**7/r0**5) + (1.0*a**2/20.0)*(rc**8/r0**6)
             bz(ix,iy,iz) = SQRT(bz0**2 - b0**2*bz(ix,iy,iz))
             rho(ix,iy,iz) = 0.45_num*(1.0_num+COS(pi*rc))+0.1_num

             x1 = xb(ix) 
             y1 = yc(iy) 
             rc = SQRT(x1**2 + y1**2)
             IF (rc >= 1.0_num) rc = 1.0_num
             b_theta = rc/2.0 - rc**3/(4.0*r0**2) + a*rc**4/(5.0*r0**3)
             bx(ix,iy,iz) = - b0 * b_theta * y1 / rc  

             x1 = xc(ix) 
             y1 = yb(iy) 
             rc = SQRT(x1**2 + y1**2)
             IF (rc >= 1.0_num) rc = 1.0_num
             b_theta = rc/2.0 - rc**3/(4.0*r0**2) + a*rc**4/(5.0*r0**3)
             by(ix,iy,iz) = b0 * b_theta * x1 / rc

          END DO
       END DO
    END DO
    bx(-2,:,:) = bx(-1,:,:)
    by(:,-2,:) = by(:,-1,:)
    bz(:,:,-2) = bz(:,:,-1)


    WHERE (rho < 0.1_num) rho = 0.1_num
    energy = 0.01_num / (rho * (gamma-1.0_num))

    k = 2.0_num * pi / length_z     ! apply velocity perturbation 
    amp = 1.e-2_num
    r1 = 0.95_num
    mu = 0.2_num

    DO ix = -1, nx+2
       DO iy = -1, ny+2
          DO iz = -1, nz+2
             rc = SQRT(xb(ix)**2 + yb(iy)**2)
             IF (rc < r1) THEN
                vx(ix,iy,iz) = amp*COS(2.5_num*k*zc(iz))   &
                     * (1.0_num+COS(k*zc(iz)))*(1.0_num - (rc/r1)**2)**mu
                vy(ix,iy,iz) = amp*SIN(2.5_num*k*zc(iz))   &
                     * (1.0_num+COS(k*zc(iz)))*(1.0_num - (rc/r1)**2)**mu
             ELSE
                vx(ix,iy,iz) = 0.0_num
                vy(ix,iy,iz) = 0.0_num
             END IF
          END DO
       END DO
    END DO


  END SUBROUTINE equilibrium


!example of how to get B from A
  SUBROUTINE equilibrium  

    INTEGER :: ix, iy, iz
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: ax,ay,az
    REAL(num) :: mag_scale_height, b0, r, radius

    ALLOCATE(ax(-2:nx+2,-2:ny+2,-2:nz+2),ay(-2:nx+2,-2:ny+2,-2:nz+2),az(-2:nx+2,-2:ny+2,-2:nz+2))

    ax = 0.0_num
    ay = 0.0_num
    az = 0.0_num

    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num

    rho = 1.0_num
    energy = 0.01_num  / (gamma - 1.0_num)

    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num

    r = 10.0_num ! loop major radius, foot points are at +- r
    b0 = 1.0_num
    mag_scale_height = 2.0_num * r / pi
    grav = 0.0_num

    ! Define the vector potential
    DO iz = -2, nz + 2
       DO ix = -2, nx + 2
          ay(ix,:,iz) = b0 * mag_scale_height * COS(xb(ix) / mag_scale_height) * EXP(-zb(iz)/mag_scale_height)
       END DO
    END DO


    ! Take the curl of the vector potential to get B
    DO iz = -1, nz + 2
       DO iy = -1, ny + 2
          DO ix = -1, nx + 2
             ixm = ix - 1 
             iym = iy - 1
             izm = iz - 1
             bx(ix,iy,iz) = (az(ix,iy,iz) - az(ix,iym,iz)) / dyb(iy) - (ay(ix,iy,iz) - ay(ix,iy,izm)) / dzb(iz)
             by(ix,iy,iz) = (ax(ix,iy,iz) - ax(ix,iy,izm)) / dzb(iz) - (az(ix,iy,iz) - az(ixm,iy,iz)) / dxb(ix)
             bz(ix,iy,iz) = (ay(ix,iy,iz) - ay(ixm,iy,iz)) / dxb(ix) - (ax(ix,iy,iz) - ax(ix,iym,iz)) / dyb(iy)
          END DO
       END DO
    END DO

    DO iz = -1, nz + 2
       DO iy = -1, ny + 2
          iym = iy - 1
          izm = iz - 1
          bx(-2,iy,iz) = (az(-2,iy,iz) - az(-2,iym,iz)) / dyb(iy) - (ay(-2,iy,iz) - ay(-2,iy,izm)) / dzb(iz)
       END DO
    END DO

    DO iz = -1, nz + 2
       DO ix = -1, nx + 2
          by(ix,-2,iz) = (ax(ix,-2,iz) - ax(ix,-2,izm)) / dzb(iz) - (az(ix,-2,iz) - az(ixm,-2,iz)) / dxb(ix)
       END DO
    END DO

    DO iy = -1, ny + 2
       DO ix = -1, nx + 2
          bz(ix,iy,-2) = (ay(ix,iy,-2) - ay(ixm,iy,-2)) / dxb(ix) - (ax(ix,iy,-2) - ax(ix,iym,-2)) / dyb(iy)
       END DO
    END DO

    DEALLOCATE(ax,ay,az)


  END SUBROUTINE equilibrium
