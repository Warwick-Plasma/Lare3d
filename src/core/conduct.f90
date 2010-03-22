MODULE conduct

  USE shared_data
  USE boundary
  USE eos

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: conduct_heat

CONTAINS

  SUBROUTINE conduct_heat

    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: uxkx, uxky, uxkz  
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: uykx, uyky, uykz 
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: uzkx, uzky, uzkz  
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: energy0, limiter
    REAL(num) :: e2t, exb, eyb, ezb 
    REAL(num) :: b, bxc, byc, bzc, bpx, bpy, bpz 
    REAL(num) :: ux, uy, uz
    REAL(num) :: pow = 5.0_num / 2.0_num  
    REAL(num) :: a1, a2, a3, error, errmax, errmax_prev = 0.0_num
    REAL(num) :: w, residual, q_shx, q_shy, q_shz, q_sh, q_f, q_nl
    REAL(num) :: initial_energy, final_energy

    INTEGER :: loop, redblack, x1, y1, z1

    LOGICAL :: converged
    REAL(num), PARAMETER :: fractional_error = 1.0e-2_num
    REAL(num), PARAMETER :: b_min = 1.0e-5_num

    ALLOCATE(uxkx(-1:nx+1, -1:ny+1, -1:nz+1), uxky(-1:nx+1, -1:ny+1, -1:nz+1))
    ALLOCATE(uxkz(-1:nx+1, -1:ny+1, -1:nz+1))
    ALLOCATE(uykx(-1:nx+1, -1:ny+1, -1:nz+1), uyky(-1:nx+1, -1:ny+1, -1:nz+1))
    ALLOCATE(uykz(-1:nx+1, -1:ny+1, -1:nz+1))
    ALLOCATE(uzkx(-1:nx+1, -1:ny+1, -1:nz+1), uzky(-1:nx+1, -1:ny+1, -1:nz+1))
    ALLOCATE(uzkz(-1:nx+1, -1:ny+1, -1:nz+1))
    ALLOCATE(energy0(-1:nx+2, -1:ny+2, -1:nz+2))
    ALLOCATE(limiter(-1:nx+2, -1:ny+2, -1:nz+2))
            
    a1 = 0.0_num
    DO iz = 1, nz 
      DO iy = 1, ny 
        DO ix = 1, nx
          a1 = a1 + energy(ix,iy,iz) * dxb(ix) * dyb(iy) * dzb(iz)
        END DO
      END DO
    END DO
    CALL MPI_ALLREDUCE(a1, initial_energy, 1, mpireal, MPI_SUM, comm, errcode) 

		! find factor required to convert between energy and temperature
		! N.B. only works for simple EOS
    e2t = (gamma - 1.0_num) / 2.0_num  

    DO iz = -1, nz + 1
      DO iy = -1, ny + 1
        DO ix = -1, nx + 1 
  
          ! x face centred B field
          bxc = bx(ix,iy,iz) 
          byc = (by(ix,iy,iz) + by(ix+1,iy,iz) + by(ix,iy-1,iz) + by(ix+1,iy-1,iz)) / 4.0_num
          bzc = (bz(ix,iy,iz) + bz(ix+1,iy,iz) + bz(ix,iy,iz-1) + bz(ix+1,iy,iz-1)) / 4.0_num
          bpx = SQRT(bxc**2 + byc**2 + bzc**2 + b_min**2) 
  
          exb = (energy(ix,iy,iz) + energy(ix+1,iy,iz)) / 2.0_num
  
          ! Direction of magnetic field on x face 
          ux = bxc / bpx       
          uy = byc / bpx       
          uz = bzc / bpx       
          ! Kappa along magnetic field, now a vector  
          uxkx(ix,iy,iz) = ux * ux * kappa_0 * (e2t * exb)**pow 
          uxky(ix,iy,iz) = ux * uy * kappa_0 * (e2t * exb)**pow 
          uxkz(ix,iy,iz) = ux * uz * kappa_0 * (e2t * exb)**pow 
          ! add symmetic conduction near b=0 points 
          uxkx(ix,iy,iz) = uxkx(ix,iy,iz) + b_min**2 * kappa_0 * (e2t * exb)**pow / bpx**2 
  
          ! y face centred B field
          bxc = (bx(ix,iy,iz) + bx(ix,iy+1,iz) + bx(ix-1,iy,iz) + bx(ix-1,iy+1,iz)) / 4.0_num
          byc = by(ix,iy,iz)      
          bzc = (bz(ix,iy,iz) + bz(ix,iy+1,iz) + bz(ix,iy,iz-1) + bz(ix,iy+1,iz-1)) / 4.0_num
          bpy = SQRT(bxc**2 + byc**2 + bzc**2 + b_min**2)
  
          eyb = (energy(ix,iy,iz) + energy(ix,iy+1,iz)) / 2.0_num
  
          ! Direction of magnetic field on y face 
          ux = bxc / bpy      
          uy = byc / bpy       
          uz = bzc / bpy      
          ! Kappa along magnetic field, now a vector  
          uykx(ix,iy,iz) = uy * ux * kappa_0 * (e2t * eyb)**pow 
          uyky(ix,iy,iz) = uy * uy * kappa_0 * (e2t * eyb)**pow    
          uykz(ix,iy,iz) = uy * uz * kappa_0 * (e2t * eyb)**pow    
          ! add symmetic conduction near b=0 points 
          uyky(ix,iy,iz) = uyky(ix,iy,iz) + b_min**2 * kappa_0 * (e2t * eyb)**pow / bpy**2  
  
          ! z face centred B field
          bxc = (bx(ix,iy,iz) + bx(ix,iy,iz+1) + bx(ix-1,iy,iz) + bx(ix-1,iy,iz+1)) / 4.0_num
          byc = (by(ix,iy,iz) + by(ix,iy,iz+1) + by(ix,iy-1,iz) + by(ix,iy-1,iz+1)) / 4.0_num         
          bzc = bz(ix,iy,iz) 
          bpz = SQRT(bxc**2 + byc**2 + bzc**2 + b_min**2)
          
          ezb = (energy(ix,iy,iz) + energy(ix,iy,iz+1)) / 2.0_num
          
          ! Direction of magnetic field on z face 
          uy = byc / bpz       
          ux = bxc / bpz      
          uz = bzc / bpz      
          ! Kappa along magnetic field, now a vector  
          uzkx(ix,iy,iz) = uz * ux * kappa_0 * (e2t * ezb)**pow 
          uzky(ix,iy,iz) = uz * uy * kappa_0 * (e2t * ezb)**pow    
          uzkz(ix,iy,iz) = uz * uz * kappa_0 * (e2t * ezb)**pow    
          ! add symmetic conduction near b=0 points 
          uzkz(ix,iy,iz) = uzkz(ix,iy,iz) + b_min**2 * kappa_0 * (e2t * ezb)**pow / bpz**2   
  
        END DO
      END DO  
    END DO  

    IF (heat_flux_limiter) THEN
      DO iz = 0, nz + 1
        DO iy = 0, ny + 1
          DO ix = 0, nx + 1  
            ! estimate the parallel heat flux and the centre of a cell
            q_shx = (uxkx(ix,iy,iz) + uxkx(ix-1,iy,iz)) * e2t &
              * (energy(ix+1,iy,iz) - energy(ix-1,iy,iz)) / dxc(ix) 
            q_shy = (uyky(ix,iy,iz) + uyky(ix,iy-1,iz)) * e2t &
              * (energy(ix,iy+1,iz) - energy(ix,iy-1,iz)) / dyc(iy) 
            q_shz = (uzkz(ix,iy,iz) + uzkz(ix,iy,iz-1)) * e2t &
              * (energy(ix,iy,iz+1) - energy(ix,iy,iz-1)) / dzc(iz) 
            q_sh = SQRT(q_shx**2 + q_shy**2 + q_shz**2) / 16.0_num  
            ! estimate the free streaming limit
            ! 42.85 = SRQT(m_p/m_e)    
            q_f = 42.85_num * flux_limiter * rho(ix,iy,iz) &
              * (e2t * MIN(energy(ix,iy,iz), temperature_100mk))**1.5_num            
            q_nl = 1.0_num / (1.0_num / MAX(q_sh, none_zero) + 1.0_num / MAX(q_f, none_zero)) 
            limiter(ix,iy,iz) = q_nl / MAX(q_sh, none_zero) 
          END DO
        END DO  
      END DO  
      DO iz = 0, nz+1
        DO iy = 0, ny+1
          DO ix = 0, nx+1  
            uxkx(ix,iy,iz) = uxkx(ix,iy,iz) * (limiter(ix,iy,iz) + limiter(ix+1,iy,iz)) 
            uxky(ix,iy,iz) = uxky(ix,iy,iz) * (limiter(ix,iy,iz) + limiter(ix+1,iy,iz))
            uxkz(ix,iy,iz) = uxkz(ix,iy,iz) * (limiter(ix,iy,iz) + limiter(ix+1,iy,iz))
            uyky(ix,iy,iz) = uyky(ix,iy,iz) * (limiter(ix,iy,iz) + limiter(ix,iy+1,iz)) 
            uykx(ix,iy,iz) = uykx(ix,iy,iz) * (limiter(ix,iy,iz) + limiter(ix,iy+1,iz))
            uykz(ix,iy,iz) = uykz(ix,iy,iz) * (limiter(ix,iy,iz) + limiter(ix,iy+1,iz))
            uzkx(ix,iy,iz) = uzkx(ix,iy,iz) * (limiter(ix,iy,iz) + limiter(ix,iy,iz+1)) 
            uzky(ix,iy,iz) = uzky(ix,iy,iz) * (limiter(ix,iy,iz) + limiter(ix,iy,iz+1))
            uzkz(ix,iy,iz) = uzkz(ix,iy,iz) * (limiter(ix,iy,iz) + limiter(ix,iy,iz+1))
          END DO
        END DO
      END DO
    END IF  
    
    converged = .FALSE. 
    w = 1.5_num       ! initial over-relaxation parameter  
		! store energy^{n} 
		energy0 = energy  
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
              
              ! terms not containing energy(ix,iy,iz) resulting from 
              ! d^2/dx^2, d^2/dy^2 and d^2/dz^2 derivatives
              a2 = uxkx(ix,iy,iz)*e2t*energy(ix+1,iy,iz)/(dxc(ix)*dxb(ix)) &
                + uxkx(ix-1,iy,iz)*e2t*energy(ix-1,iy,iz)/(dxc(ix-1)*dxb(ix)) 
              a2 = a2 + uyky(ix,iy,iz)*e2t*energy(ix,iy+1,iz)/(dyc(iy)*dyb(iy)) &
                + uyky(ix,iy-1,iz)*e2t*energy(ix,iy-1,iz)/(dyc(iy-1)*dyb(iy))              
              a2 = a2 + uzkz(ix,iy,iz)*e2t*energy(ix,iy,iz+1)/(dzc(iz)*dzb(iz)) &
                + uzkz(ix,iy,iz-1)*e2t*energy(ix,iy,iz-1)/(dzc(iz-1)*dzb(iz))              
                
              ! terms not containing energy(ix,iy,iz) resulting from 
              ! d^2/dxdy cross derivatives                  
              a2 = a2 + uxky(ix,iy,iz) * e2t &
                * (energy(ix+1,iy+1,iz) + energy(ix,iy+1,iz) - energy(ix+1,iy-1,iz) - energy(ix,iy-1,iz)) &
                / (2.0_num * dxb(ix) * (dyc(iy) + dyc(iy-1)))  
              a2 = a2 - uxky(ix-1,iy,iz) * e2t &
                * (energy(ix,iy+1,iz) + energy(ix-1,iy+1,iz) - energy(ix,iy-1,iz) - energy(ix-1,iy-1,iz)) &
                / (2.0_num * dxb(ix) * (dyc(iy) + dyc(iy-1)))  
              ! terms not containing energy(ix,iy,iz) resulting from 
              ! d^2/dxdz cross derivatives                  
              a2 = a2 + uxkz(ix,iy,iz) * e2t &
                * (energy(ix+1,iy,iz+1) + energy(ix,iy,iz+1) - energy(ix+1,iy,iz-1) - energy(ix,iy,iz-1)) &
                / (2.0_num * dxb(ix) * (dzc(iz) + dzc(iz-1)))  
              a2 = a2 - uxkz(ix-1,iy,iz) * e2t &
                * (energy(ix,iy,iz+1) + energy(ix-1,iy,iz+1) - energy(ix,iy,iz-1) - energy(ix-1,iy,iz-1)) &
                / (2.0_num * dxb(ix) * (dzc(iz) + dzc(iz-1)))
        
              ! terms not containing energy(ix,iy,iz) resulting from 
              ! d^2/dydx cross derivatives
              a2 = a2 + uykx(ix,iy,iz) * e2t &
                * (energy(ix+1,iy+1,iz) + energy(ix+1,iy,iz) - energy(ix-1,iy+1,iz) - energy(ix-1,iy,iz)) &
                / (2.0_num * dyb(iy) * (dxc(ix) + dxc(ix-1)))  
              a2 = a2 - uykx(ix,iy-1,iz) * e2t &
                * (energy(ix+1,iy,iz) + energy(ix+1,iy-1,iz) - energy(ix-1,iy,iz) - energy(ix-1,iy-1,iz)) &
                / (2.0_num * dyb(iy) * (dxc(ix) + dxc(ix-1)))                     
              ! terms not containing energy(ix,iy,iz) resulting from 
              ! d^2/dydz cross derivatives
              a2 = a2 + uykz(ix,iy,iz) * e2t &
                * (energy(ix,iy+1,iz+1) + energy(ix,iy,iz+1) - energy(ix,iy+1,iz-1) - energy(ix,iy,iz-1)) &
                / (2.0_num * dyb(iy) * (dzc(iz) + dzc(iz-1)))  
              a2 = a2 - uykz(ix,iy-1,iz) * e2t &
                * (energy(ix+1,iy,iz) + energy(ix+1,iy-1,iz) - energy(ix-1,iy,iz) - energy(ix-1,iy-1,iz)) &
                / (2.0_num * dyb(iy) * (dzc(iz) + dzc(iz-1)))
              
              ! terms not containing energy(ix,iy,iz) resulting from 
              ! d^2/dzdx cross derivatives
              a2 = a2 + uzkx(ix,iy,iz) * e2t &
                * (energy(ix+1,iy,iz+1) + energy(ix+1,iy,iz) - energy(ix-1,iy,iz+1) - energy(ix-1,iy,iz)) &
                / (2.0_num * dzb(iz) * (dxc(ix) + dxc(ix-1)))  
              a2 = a2 - uzkx(ix,iy,iz-1) * e2t &
                * (energy(ix+1,iy,iz) + energy(ix+1,iy,iz-1) - energy(ix-1,iy,iz) - energy(ix-1,iy,iz-1)) &
                / (2.0_num * dzb(iz) * (dxc(ix) + dxc(ix-1)))  
              ! terms not containing energy(ix,iy,iz) resulting from 
              ! d^2/dzdy cross derivatives
              a2 = a2 + uzky(ix,iy,iz) * e2t &
                * (energy(ix,iy+1,iz+1) + energy(ix,iy+1,iz) - energy(ix,iy-1,iz+1) - energy(ix,iy-1,iz)) &
                / (2.0_num * dzb(iz) * (dyc(iy) + dyc(iy-1)))  
              a2 = a2 - uzky(ix,iy,iz-1) * e2t &
                * (energy(ix,iy+1,iz) + energy(ix,iy+1,iz-1) - energy(ix,iy-1,iz) - energy(ix,iy-1,iz-1)) &
                / (2.0_num * dzb(iz) * (dyc(iy) + dyc(iy-1)))     
              
              a1 = a1 * dt * e2t / rho(ix,iy,iz)     
              a2 = a2 * dt / rho(ix,iy,iz)         
  
              residual = energy(ix,iy,iz) &
                    - (energy0(ix,iy,iz)  + a2) / (1.0_num + a1) 
              energy(ix,iy,iz) = MAX(energy(ix,iy,iz) - w * residual, 0.0_num)
              error = ABS(residual / MAX(energy(ix,iy,iz), energy0(ix,iy,iz), none_zero))     
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

    a1 = 0.0_num
    DO iz = 1, nz 
      DO iy = 1, ny 
        DO ix = 1, nx 
          a1 = a1 + energy(ix,iy,iz) * dxb(ix) * dyb(iy) * dzb(iz)
        END DO
      END DO 
    END DO 
    CALL MPI_ALLREDUCE(a1, final_energy, 1, mpireal, MPI_SUM, comm, errcode)
    IF (initial_energy < final_energy) energy = energy * initial_energy / final_energy

    DEALLOCATE(uxkx, uxky, uxkz)
    DEALLOCATE(uykx, uyky, uykz)
    DEALLOCATE(uzkx, uzky, uzkz)
    DEALLOCATE(energy0)
    DEALLOCATE(limiter)

  END SUBROUTINE conduct_heat

END MODULE conduct
