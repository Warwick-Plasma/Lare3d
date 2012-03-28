MODULE conduct

  USE shared_data
  USE boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: conduct_heat

  REAL(num) :: heat0 = 0.0_num
  REAL(num) :: rho_cor = 0.0_num 
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE  :: temperature              
  
CONTAINS

  SUBROUTINE conduct_heat

    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: uxkx, uxky, uxkz  
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: uykx, uyky, uykz 
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: uzkx, uzky, uzkz  
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: energy0, limiter, temperature0
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: radiation, heat_in, alpha
    REAL(num) :: txb, tyb, tzb 
    REAL(num) :: b, bxc, byc, bzc, bpx, bpy, bpz 
    REAL(num) :: ux, uy, uz
    REAL(num) :: pow = 5.0_num / 2.0_num  
    REAL(num) :: a1, a2, a3, a4, error, errmax, rad, alf
    REAL(num) :: w, residual, q_shx, q_shy, q_shz, q_sh, q_f, q_nl

    INTEGER :: loop, redblack, x1, y1, z1 

    LOGICAL :: converged
    REAL(num), PARAMETER :: fractional_error = 1.e-4_num
    REAL(num), PARAMETER :: b_min = 1.e-3_num

    ALLOCATE(uxkx(-1:nx+1,-1:ny+1,-1:nz+1), uxky(-1:nx+1,-1:ny+1,-1:nz+1))
    ALLOCATE(uxkz(-1:nx+1,-1:ny+1,-1:nz+1))
    ALLOCATE(uykx(-1:nx+1,-1:ny+1,-1:nz+1), uyky(-1:nx+1,-1:ny+1,-1:nz+1))
    ALLOCATE(uykz(-1:nx+1,-1:ny+1,-1:nz+1))
    ALLOCATE(uzkx(-1:nx+1,-1:ny+1,-1:nz+1), uzky(-1:nx+1,-1:ny+1,-1:nz+1))
    ALLOCATE(uzkz(-1:nx+1,-1:ny+1,-1:nz+1))
    ALLOCATE(energy0(-1:nx+2,-1:ny+2,-1:nz+2))
    ALLOCATE(limiter(-1:nx+2,-1:ny+2,-1:nz+2))
    ALLOCATE(temperature0(-1:nx+2,-1:ny+2,-1:nz+2))
    ALLOCATE(temperature(-1:nx+2,-1:ny+2,-1:nz+2))
    ALLOCATE(radiation(-1:nx+2,-1:ny+2,-1:nz+2))
    ALLOCATE(heat_in(-1:nx+2,-1:ny+2,-1:nz+2))
    ALLOCATE(alpha(-1:nx+2,-1:ny+2,-1:nz+2))
            
    temperature = (gamma - 1.0_num) * (energy - (1.0_num - xi_n) * ionise_pot) &
                 / (2.0_num - xi_n)         

    IF (.NOT. conduction) kappa_0 = 0.0_num       

    DO iz = -1, nz + 1
      DO iy = -1, ny + 1
        DO ix = -1, nx + 1 
  
          ! x face centred B field
          bxc = bx(ix,iy,iz) 
          byc = (by(ix,iy,iz) + by(ix+1,iy,iz) + by(ix,iy-1,iz) + by(ix+1,iy-1,iz)) / 4.0_num
          bzc = (bz(ix,iy,iz) + bz(ix+1,iy,iz) + bz(ix,iy,iz-1) + bz(ix+1,iy,iz-1)) / 4.0_num
          bpx = SQRT(bxc**2 + byc**2 + bzc**2 + b_min**2)
          bpx = MAX(bpx, none_zero) 
  
          txb = (temperature(ix,iy,iz) + temperature(ix+1,iy,iz)) / 2.0_num
  
          ! Direction of magnetic field on x face 
          ux = bxc / bpx       
          uy = byc / bpx       
          uz = bzc / bpx       
          ! Kappa along magnetic field, now a vector  
          uxkx(ix,iy,iz) = ux * ux * kappa_0 * txb**pow 
          uxky(ix,iy,iz) = ux * uy * kappa_0 * txb**pow 
          uxkz(ix,iy,iz) = ux * uz * kappa_0 * txb**pow 
          ! add symmetic conduction near b=0 points 
          uxkx(ix,iy,iz) = uxkx(ix,iy,iz) + b_min**2 * kappa_0 * txb**pow &
              / (bpx**2 + b_min**2) 
  
          ! y face centred B field
          bxc = (bx(ix,iy,iz) + bx(ix,iy+1,iz) + bx(ix-1,iy,iz) + bx(ix-1,iy+1,iz)) / 4.0_num
          byc = by(ix,iy,iz)      
          bzc = (bz(ix,iy,iz) + bz(ix,iy+1,iz) + bz(ix,iy,iz-1) + bz(ix,iy+1,iz-1)) / 4.0_num
          bpy = SQRT(bxc**2 + byc**2 + bzc**2 + b_min**2)
          bpy = MAX(bpy, none_zero) 
  
          tyb = (temperature(ix,iy,iz) + temperature(ix,iy+1,iz)) / 2.0_num
  
          ! Direction of magnetic field on y face 
          ux = bxc / bpy      
          uy = byc / bpy       
          uz = bzc / bpy      
          ! Kappa along magnetic field, now a vector  
          uykx(ix,iy,iz) = uy * ux * kappa_0 * tyb**pow 
          uyky(ix,iy,iz) = uy * uy * kappa_0 * tyb**pow    
          uykz(ix,iy,iz) = uy * uz * kappa_0 * tyb**pow    
          ! add symmetic conduction near b=0 points 
          uyky(ix,iy,iz) = uyky(ix,iy,iz) + b_min**2 * kappa_0 * tyb**pow &
                / (bpy**2 + b_min**2)  
  
          ! z face centred B field
          bxc = (bx(ix,iy,iz) + bx(ix,iy,iz+1) + bx(ix-1,iy,iz) + bx(ix-1,iy,iz+1)) / 4.0_num
          byc = (by(ix,iy,iz) + by(ix,iy,iz+1) + by(ix,iy-1,iz) + by(ix,iy-1,iz+1)) / 4.0_num         
          bzc = bz(ix,iy,iz) 
          bpz = SQRT(bxc**2 + byc**2 + bzc**2 + b_min**2)
          bpz = MAX(bpz, none_zero) 
          
          tzb = (temperature(ix,iy,iz) + temperature(ix,iy,iz+1)) / 2.0_num
          
          ! Direction of magnetic field on z face 
          ux = bxc / bpz      
          uy = byc / bpz       
          uz = bzc / bpz      
          ! Kappa along magnetic field, now a vector  
          uzkx(ix,iy,iz) = uz * ux * kappa_0 * tzb**pow 
          uzky(ix,iy,iz) = uz * uy * kappa_0 * tzb**pow    
          uzkz(ix,iy,iz) = uz * uz * kappa_0 * tzb**pow    
          ! add symmetic conduction near b=0 points 
          uzkz(ix,iy,iz) = uzkz(ix,iy,iz) + b_min**2 * kappa_0 * tzb**pow &
              / (bpz**2 + b_min**2)   
  
        END DO
      END DO  
    END DO  

    IF (heat_flux_limiter) THEN
      DO iz = 0, nz + 1
        DO iy = 0, ny + 1
          DO ix = 0, nx + 1  
            ! estimate the parallel heat flux and the centre of a cell
            q_shx = &
                (uxkx(ix,iy,iz) + uxkx(ix-1,iy,iz))  &
                * (temperature(ix+1,iy,iz) - temperature(ix-1,iy,iz)) / dxc(ix) & 
              + (uxky(ix,iy,iz) + uxky(ix,iy-1,iz))  &
                * (temperature(ix,iy+1,iz) - temperature(ix,iy-1,iz)) / dyc(iy) &
              + (uxkz(ix,iy,iz) + uxkz(ix,iy,iz-1))  &
                * (temperature(ix,iy,iz+1) - temperature(ix,iy,iz-1)) / dzc(iz) 
            q_shy = &
                (uykx(ix,iy,iz) + uykx(ix-1,iy,iz))  &
                * (temperature(ix+1,iy,iz) - temperature(ix-1,iy,iz)) / dxc(ix) & 
              + (uyky(ix,iy,iz) + uyky(ix,iy-1,iz))  &
                * (temperature(ix,iy+1,iz) - temperature(ix,iy-1,iz)) / dyc(iy) &
              + (uykz(ix,iy,iz) + uykz(ix,iy,iz-1))  &
                * (temperature(ix,iy,iz+1) - temperature(ix,iy,iz-1)) / dzc(iz) 
            q_shz = &
                (uzkx(ix,iy,iz) + uzkx(ix-1,iy,iz))  &
                * (temperature(ix+1,iy,iz) - temperature(ix-1,iy,iz)) / dxc(ix) & 
              + (uzky(ix,iy,iz) + uzky(ix,iy-1,iz))  &
                * (temperature(ix,iy+1,iz) - temperature(ix,iy-1,iz)) / dyc(iy) &
              + (uzkz(ix,iy,iz) + uzkz(ix,iy,iz-1))  &
                * (temperature(ix,iy,iz+1) - temperature(ix,iy,iz-1)) / dzc(iz) 

            q_sh = SQRT(q_shx**2 + q_shy**2 + q_shz**2) / 16.0_num  
            ! estimate the free streaming limit
            ! 42.85 = SRQT(m_p/m_e)    
            q_f = 42.85_num * flux_limiter * rho(ix,iy,iz) &
              * MIN(temperature(ix,iy,iz), temperature_100mk)**1.5_num            
            q_nl = 1.0_num / (1.0_num / MAX(q_sh, none_zero) + 1.0_num / MAX(q_f, none_zero)) 
            limiter(ix,iy,iz) = q_nl / MAX(q_sh, none_zero) / 2.0_num
          END DO
        END DO  
      END DO  
      DO iz = 0, nz+1
        DO iy = 0, ny+1
          DO ix = 0, nx+1  
            uxkx(ix,iy,iz) = uxkx(ix,iy,iz) * (limiter(ix,iy,iz) + limiter(ix+1,iy,iz)) 
            uxky(ix,iy,iz) = uxky(ix,iy,iz) * (limiter(ix,iy,iz) + limiter(ix+1,iy,iz))
            uxkz(ix,iy,iz) = uxkz(ix,iy,iz) * (limiter(ix,iy,iz) + limiter(ix+1,iy,iz))
            uykx(ix,iy,iz) = uykx(ix,iy,iz) * (limiter(ix,iy,iz) + limiter(ix,iy+1,iz))
            uyky(ix,iy,iz) = uyky(ix,iy,iz) * (limiter(ix,iy,iz) + limiter(ix,iy+1,iz)) 
            uykz(ix,iy,iz) = uykz(ix,iy,iz) * (limiter(ix,iy,iz) + limiter(ix,iy+1,iz))
            uzkx(ix,iy,iz) = uzkx(ix,iy,iz) * (limiter(ix,iy,iz) + limiter(ix,iy,iz+1)) 
            uzky(ix,iy,iz) = uzky(ix,iy,iz) * (limiter(ix,iy,iz) + limiter(ix,iy,iz+1))
            uzkz(ix,iy,iz) = uzkz(ix,iy,iz) * (limiter(ix,iy,iz) + limiter(ix,iy,iz+1))
          END DO
        END DO
      END DO
    END IF  

    converged = .FALSE. 
    w = 1.6_num       ! initial over-relaxation parameter  
		! store energy^{n} 
		energy0 = energy   
		temperature0 = temperature  

    radiation = 0.0_num
    heat_in = 0.0_num 
    alpha = 0.0_num
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx  
          heat_in(ix,iy,iz) = heating(rho(ix,iy,iz), temperature0(ix,iy,iz))  
          CALL rad_losses(rho(ix,iy,iz), temperature0(ix,iy,iz), xi_n(ix,iy,iz), zc(iz), rad, alf)
          alpha(ix,iy,iz) = alf  
          radiation(ix,iy,iz) =  rad 
        END DO
      END DO    		
    END DO    		
		
		! interate to get energy^{n+1} by SOR Guass-Seidel 		
    iterate: DO loop = 1, 100
      errmax = 0.0_num 
      error = 0.0_num
      z1 = 1 
      DO redblack = 1, 2
        y1 = z1 
        DO iz = 1, nz
          x1 = z1
          DO iy = 1, ny               
             DO ix = x1, nx, 2
              
              ! terms containing energy(ix,iz) resulting from 
              ! d^2/dx^2 and d^2/dy^2 derivatives
              a1 = uxkx(ix,iy,iz)/(dxc(ix)*dxb(ix)) + uxkx(ix-1,iy,iz)/(dxc(ix-1)*dxb(ix)) &
                + uyky(ix,iy,iz)/(dyc(iy)*dyb(iy)) + uyky(ix,iy-1,iz)/(dyc(iy-1)*dyb(iy))  &
                + uzkz(ix,iy,iz)/(dzc(iz)*dzb(iz)) + uzkz(ix,iy,iz-1)/(dzc(iz-1)*dzb(iz))  
              
              ! terms not containing temperature(ix,iy,iz) resulting from 
              ! d^2/dx^2, d^2/dy^2 and d^2/dz^2 derivatives
              a2 = uxkx(ix,iy,iz)*temperature(ix+1,iy,iz)/(dxc(ix)*dxb(ix)) &
                + uxkx(ix-1,iy,iz)*temperature(ix-1,iy,iz)/(dxc(ix-1)*dxb(ix)) 
              a2 = a2 + uyky(ix,iy,iz)*temperature(ix,iy+1,iz)/(dyc(iy)*dyb(iy)) &
                + uyky(ix,iy-1,iz)*temperature(ix,iy-1,iz)/(dyc(iy-1)*dyb(iy))              
              a2 = a2 + uzkz(ix,iy,iz)*temperature(ix,iy,iz+1)/(dzc(iz)*dzb(iz)) &
                + uzkz(ix,iy,iz-1)*temperature(ix,iy,iz-1)/(dzc(iz-1)*dzb(iz))              
                
              ! terms not containing temperature(ix,iy,iz) resulting from 
              ! d^2/dxdy cross derivatives                  
              a2 = a2 + uxky(ix,iy,iz)  &
                * (temperature(ix+1,iy+1,iz) + temperature(ix,iy+1,iz) - temperature(ix+1,iy-1,iz) &
                - temperature(ix,iy-1,iz)) &
                / (2.0_num * dxb(ix) * (dyc(iy) + dyc(iy-1)))  
              a2 = a2 - uxky(ix-1,iy,iz)  &
                * (temperature(ix,iy+1,iz) + temperature(ix-1,iy+1,iz) - temperature(ix,iy-1,iz) &
                - temperature(ix-1,iy-1,iz)) &
                / (2.0_num * dxb(ix) * (dyc(iy) + dyc(iy-1)))  
              ! terms not containing temperature(ix,iy,iz) resulting from 
              ! d^2/dxdz cross derivatives                  
              a2 = a2 + uxkz(ix,iy,iz)  &
                * (temperature(ix+1,iy,iz+1) + temperature(ix,iy,iz+1) - temperature(ix+1,iy,iz-1) &
                - temperature(ix,iy,iz-1)) &
                / (2.0_num * dxb(ix) * (dzc(iz) + dzc(iz-1)))  
              a2 = a2 - uxkz(ix-1,iy,iz)  &
                * (temperature(ix,iy,iz+1) + temperature(ix-1,iy,iz+1) - temperature(ix,iy,iz-1) &
                - temperature(ix-1,iy,iz-1)) &
                / (2.0_num * dxb(ix) * (dzc(iz) + dzc(iz-1)))
        
              ! terms not containing temperature(ix,iy,iz) resulting from 
              ! d^2/dydx cross derivatives
              a2 = a2 + uykx(ix,iy,iz)  &
                * (temperature(ix+1,iy+1,iz) + temperature(ix+1,iy,iz) - temperature(ix-1,iy+1,iz) &
                - temperature(ix-1,iy,iz)) &
                / (2.0_num * dyb(iy) * (dxc(ix) + dxc(ix-1)))  
              a2 = a2 - uykx(ix,iy-1,iz)  &
                * (temperature(ix+1,iy,iz) + temperature(ix+1,iy-1,iz) - temperature(ix-1,iy,iz) &
                - temperature(ix-1,iy-1,iz)) &
                / (2.0_num * dyb(iy) * (dxc(ix) + dxc(ix-1)))                     
              ! terms not containing temperature(ix,iy,iz) resulting from 
              ! d^2/dydz cross derivatives
              a2 = a2 + uykz(ix,iy,iz)  &
                * (temperature(ix,iy+1,iz+1) + temperature(ix,iy,iz+1) - temperature(ix,iy+1,iz-1) &
                - temperature(ix,iy,iz-1)) &
                / (2.0_num * dyb(iy) * (dzc(iz) + dzc(iz-1)))  
              a2 = a2 - uykz(ix,iy-1,iz)  &
                * (temperature(ix,iy,iz+1) + temperature(ix,iy-1,iz+1) - temperature(ix,iy,iz-1) &
                - temperature(ix,iy-1,iz-1)) &
                / (2.0_num * dyb(iy) * (dzc(iz) + dzc(iz-1)))
              
              ! terms not containing temperature(ix,iy,iz) resulting from 
              ! d^2/dzdx cross derivatives
              a2 = a2 + uzkx(ix,iy,iz)  &
                * (temperature(ix+1,iy,iz+1) + temperature(ix+1,iy,iz) - temperature(ix-1,iy,iz+1) &
                - temperature(ix-1,iy,iz)) &
                / (2.0_num * dzb(iz) * (dxc(ix) + dxc(ix-1)))  
              a2 = a2 - uzkx(ix,iy,iz-1)  &
                * (temperature(ix+1,iy,iz) + temperature(ix+1,iy,iz-1) - temperature(ix-1,iy,iz) &
                - temperature(ix-1,iy,iz-1)) &
                / (2.0_num * dzb(iz) * (dxc(ix) + dxc(ix-1)))  
              ! terms not containing temperature(ix,iy,iz) resulting from 
              ! d^2/dzdy cross derivatives
              a2 = a2 + uzky(ix,iy,iz)  &
                * (temperature(ix,iy+1,iz+1) + temperature(ix,iy+1,iz) - temperature(ix,iy-1,iz+1) &
                - temperature(ix,iy-1,iz)) &
                / (2.0_num * dzb(iz) * (dyc(iy) + dyc(iy-1)))  
              a2 = a2 - uzky(ix,iy,iz-1)  &
                * (temperature(ix,iy+1,iz) + temperature(ix,iy+1,iz-1) - temperature(ix,iy-1,iz) &
                - temperature(ix,iy-1,iz-1)) &
                / (2.0_num * dzb(iz) * (dyc(iy) + dyc(iy-1)))     
              
              a3 = (a1 + radiation(ix,iy,iz) * alpha(ix,iy,iz) / temperature0(ix,iy,iz)) * temperature(ix,iy,iz) &
                  / energy(ix,iy,iz) 

              a4 = a2 + heat_in(ix,iy,iz)  - (1.0_num - alpha(ix,iy,iz)) * radiation(ix,iy,iz) 
  
              a3 = a3 * dt / rho(ix, iy, iz)     
              a4 = a4 * dt / rho(ix, iy, iz)  

              residual = energy(ix,iy,iz) &
                    - (energy0(ix,iy,iz)  + a4) / (1.0_num + a3) 
              energy(ix, iy, iz) = MAX(energy(ix, iy, iz) - w * residual, (1.0_num - xi_n(ix,iy,iz)) * ionise_pot)    
              error = ABS(residual) / energy0(ix,iy,iz)
              errmax = MAX(errmax, error)

            END DO 
            x1 = 3 - x1
          END DO 
          y1 = 3 - y1
        END DO 
        z1 = 3 - z1
        
        CALL energy_bcs

      END DO

      CALL MPI_ALLREDUCE(errmax, error, 1, mpireal, MPI_MAX, comm, errcode)
      errmax = error 
      
      IF (errmax .LT. fractional_error) THEN
        converged = .TRUE.  
        EXIT iterate
      END IF
    END DO iterate
    
    IF (rank == 0 .AND. .NOT. converged) PRINT * , "Conduction failed at t = ", time  

    DEALLOCATE(uxkx, uxky, uxkz)
    DEALLOCATE(uykx, uyky, uykz)
    DEALLOCATE(uzkx, uzky, uzkz)
    DEALLOCATE(energy0)
    DEALLOCATE(limiter)
    DEALLOCATE(radiation)
    DEALLOCATE(heat_in)
    DEALLOCATE(alpha)
    DEALLOCATE(temperature)
    DEALLOCATE(temperature0)     
    
  END SUBROUTINE conduct_heat
  


  SUBROUTINE rad_losses(density, temperature, xi, height, rad, alf)  
    ! returns the normalised RTV losses 
    REAL(num), INTENT(IN) :: density, temperature, xi, height  
    REAL(num), INTENT(OUT) :: rad, alf 

    REAL(num), DIMENSION(7) :: trange = (/0.02_num,0.0398_num,0.0794_num,0.251_num, 0.562_num,1.995_num,10.0_num /)
    REAL(num), DIMENSION(6) :: psi = (/1.2303_num, 870.96_num, 5.496_num, 0.3467_num, 1.0_num, 1.6218_num /)
    REAL(num), DIMENSION(6) :: alpha = (/0.0_num, 2.0_num, 0.0_num, -2.0_num, 0.0_num, -2.0_num/3.0_num /) 
    REAL(num) :: tmk, factor
    INTEGER :: i
                      
    rad = 0.0_num     
    alf = 0.0_num                           
    IF(.NOT. radiation) RETURN 
    IF (height < 1.0_num) RETURN

    tmk = temperature * t2tmk   
    IF (tmk < trange(1) .OR. tmk > trange(7)) RETURN
                                                      
    DO i = 1, 6
      IF (tmk >= trange(i) .AND. tmk <= trange(i+1)) EXIT
    END DO 
    
    ! account for reduced electron number density due to neutrals                                            
    factor = (1.0_num - xi)**2
    IF (eos_number == EOS_IDEAL) factor = 1.0_num

    rad = factor * density**2 * psi(i) * tmk**alpha(i) 
    rad = rad * h_star * lr_star 
    alf = alpha(i)

  END SUBROUTINE rad_losses  
       
  
  
  FUNCTION heating(density, t0)    
    ! for a given density and temperature returns a user specific heating function
    REAL(num), INTENT(IN) :: density, t0
    REAL(num) :: heating
    REAL(num) :: tmk, a1, a2, rad, alf, height
    LOGICAL :: first_call = .TRUE.
    REAL(num) :: heat0 = 0.0_num      
    REAL(num) :: rho_coronal = 0.0_num  
    LOGICAL :: store_state
    INTEGER :: loop
    
    IF (first_call) THEN
      a1 = 0.0_num     
      IF (proc_y_max == MPI_PROC_NULL) THEN  
         store_state = radiation    
         radiation = .TRUE.
         CALL rad_losses(rho(1,1,nz), temperature(1,1,nz), xi_n(1,1,nz), zc(nz), rad, alf) 
         radiation = store_state
         a1 = rad / rho(1,1,nz)**2 
      END IF               
      CALL MPI_ALLREDUCE(a1, heat0, 1, mpireal, MPI_MAX, comm, errcode)  
       
      ! choose a reference density based on height   
      height = 15.0_num  
      a2 = 0.0_num
      DO loop = 1, nz
         IF (zb(loop) >= height .AND. zb(loop-1) < height) a2 = rho(1,1,loop)  
      END DO 
      CALL MPI_ALLREDUCE(a2, rho_coronal, 1, mpireal, MPI_MAX, comm, errcode)    

      first_call = .FALSE.               
    END IF     

    heating = 0.0_num
    IF (.NOT. coronal_heating) RETURN

    tmk = t0 * t2tmk    
    ! For low density and high temeprature define a heating course term                   
    IF(density < rho_coronal .AND. tmk > 0.02_num) heating = 100.0_num * heat0 * density**2

  END FUNCTION heating
  
                
                          

END MODULE conduct
