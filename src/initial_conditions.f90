
MODULE initial_conditions

  USE shared_data
  USE neutral
  USE boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_initial_conditions

  REAL(num), DIMENSION(:), ALLOCATABLE :: axis, temperature_1d
  INTEGER :: table_count

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

    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: temperature
    REAL(num) :: xi_v, amp, centre, width

    ! Below are all the variables which must be defined and their sizes

    vx(-2:nx+2, -2:ny+2, -2:nz+2) = 0.0_num
    vy(-2:nx+2, -2:ny+2, -2:nz+2) = 0.0_num
    vz(-2:nx+2, -2:ny+2, -2:nz+2) = 0.0_num

    bx(-2:nx+2, -1:ny+2, -1:nz+2) = 0.0_num
    by(-1:nx+2, -2:ny+2, -1:nz+2) = 0.0_num
    bz(-1:nx+2, -1:ny+2, -2:nz+2) = 0.0_num

    rho(-1:nx+2, -1:ny+2, -2:nz+2) = 1.0_num
    energy(-1:nx+2, -1:ny+2, -2:nz+2) = 0.1_num

    grav(-1:nz+2) = 0.0_num

    ! If defining the initial conditions using temperature then use
    ALLOCATE(temperature(-1:nx+2, -1:ny+2, -1:nz+2))
    temperature(-1:nx+2, -1:ny+2, -1:nz+2) = 0.5_num

    ! If neutrals included xi_n is a function of temperature so iteration required
    ! Set the neutral fraction if needed
    DO iz = -1,nz+2
      DO iy = -1,ny+2
        DO ix = -1,nx+2
          IF (eos_number /= EOS_IDEAL) THEN         
            xi_v = get_neutral(temperature(ix,iy,iz), rho(ix,iy,iz))
          ELSE  
            IF (neutral_gas) THEN
              xi_v = 1.0_num
            ELSE
              xi_v = 0.0_num
            END IF
          END IF
          energy(ix,iy,iz) = (temperature(ix,iy,iz) * (2.0_num - xi_v) &
                + (1.0_num - xi_v) * ionise_pot * (gamma - 1.0_num)) &
                / (gamma - 1.0_num)
        END DO
      END DO
    END DO
    DEALLOCATE(temperature)

    ! set background, non-shock, viscosity
    visc3 = 0.0_num
    IF (use_viscous_damping) THEN
      width = length_z / 10.0_num
      centre = 0.4_num * length_z + width
      amp = 1.e2_num
      DO iz = -1, nz+1
        DO iy = -1, ny + 1
          DO ix = -1, nx + 1
            visc3(ix,iy,iz) = visc3(ix,iy,iz) + amp * (1.0_num + TANH((ABS(zb(iz)) - centre) / width)) 
          END DO
        END DO
      END DO
    END IF    
  
    ! An example for setting up simple potential field
    CALL potential_field

  END SUBROUTINE set_initial_conditions
  



  SUBROUTINE potential_field()

    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: phi
    REAL(num) :: w, errmax, error, residual, fractional_error
    REAL(num) :: bz_max, bz_max_local
    REAL(num) :: dx1, dy1, dz1
    INTEGER :: loop, x1, y1, z1, redblack
    LOGICAL :: converged

    ALLOCATE(phi(-1:nx+2,-1:ny+2, -1:nz+2))
    phi(:,:,:) = 0.0_num

    ! Only works for uniform x,y grids - no strecthed grids
    dx1 = (REAL(nx_global, num) / length_x)**2
    dy1 = (REAL(ny_global, num) / length_y)**2
    dz1 = (REAL(nz_global, num) / length_z)**2

    converged = .FALSE.
    w = 1.6_num   !2.0_num / (1.0_num + SIN(pi / REAL(nx_global,num)))
    fractional_error = 5.e-5_num !1.e-8_num

    ! Iterate to get phi^{n+1} by SOR Gauss-Seidel
    iterate: DO loop = 1, 1000000
      errmax = 0.0_num
      error = 0.0_num

      z1 = 1
      DO redblack = 1, 2
        y1 = z1
        DO iz = 1, nz
          x1 = z1
          izm = iz - 1
          izp = iz + 1
          DO iy = 1, ny
            iym = iy - 1
            iyp = iy + 1
            DO ix = x1, nx, 2
              ixm = ix - 1
              ixp = ix + 1
              residual = (((phi(ixp,iy,iz) - 2.0_num * phi(ix,iy,iz)) + phi(ixm,iy,iz)) * dx1 &
                       + ((phi(ix,iyp,iz) - 2.0_num * phi(ix,iy,iz)) + phi(ix,iym,iz)) * dy1 &
                       + ((phi(ix,iy,izp) - 2.0_num * phi(ix,iy,iz)) + phi(ix,iy,izm)) * dz1) &
                       / 2.0_num / (dx1 + dy1 + dz1)
              phi(ix,iy,iz) = phi(ix,iy,iz) + w * residual 
              error = ABS(residual) 
              errmax = MAX(errmax, error)
            END DO
            x1 = 3 - x1
          END DO
          y1 = 3 - y1
        END DO
        CALL phi_mpi
        z1 = 3 - z1
      END DO

      CALL MPI_ALLREDUCE(errmax, error, 1, mpireal, MPI_MAX, comm, errcode)
      IF (rank == 0 .AND. (MOD(loop,1000).EQ.0)) print *, 'loop, residual = ', loop, error
      IF (error < fractional_error) THEN
        converged = .TRUE.
        EXIT iterate
      END IF
    END DO iterate

    IF (rank == 0 .AND. .NOT. converged) PRINT*, 'potential_field failed'

    DO iz = 1, nz
     DO iy = 1, ny
       DO ix = 0, nx
         bx(ix,iy,iz) = -(phi(ix+1,iy,iz)-phi(ix,iy,iz))/dxc(ix)
       END DO
     END DO
    END DO

    DO iz = 1, nz
     DO iy = 0, ny
       DO ix = 1, nx
         by(ix,iy,iz) = -(phi(ix,iy+1,iz)-phi(ix,iy,iz))/dyc(iy)
       END DO
     END DO
    END DO

    DO iz = 0, nz
     DO iy = 1, ny
       DO ix = 1, nx
         bz(ix,iy,iz) = -(phi(ix,iy,iz+1)-phi(ix,iy,iz))/dzc(iz)
       END DO
     END DO
    END DO

    CALL bfield_bcs

    !Find maximum Bz on lower boundary
    bz_max_local = -largest_number
    IF (proc_z_min == MPI_PROC_NULL) THEN
      bz_max_local = MAXVAL(bz(1:nx,1:nz,0))
    END IF
    CALL MPI_ALLREDUCE(bz_max_local, bz_max, 1, mpireal, MPI_MAX, comm, errcode) 

    !Scale the field to one in normalised units
    bz(:,:,:) = bz(:,:,:) / bz_max 
    by(:,:,:) = by(:,:,:) / bz_max
    bx(:,:,:) = bx(:,:,:) / bz_max 

    CALL bfield_bcs

    DEALLOCATE(phi)

    CONTAINS

      SUBROUTINE phi_mpi

        REAL(num) :: local_flux
        REAL(num) :: centre, radius, amp, r1

        CALL MPI_SENDRECV( &
            phi(   1,-1,-1), 1, cell_xface, proc_x_min, tag, &
            phi(nx+1,-1,-1), 1, cell_xface, proc_x_max, tag, &
            comm, status, errcode)
        CALL MPI_SENDRECV( &
            phi(nx-1,-1,-1), 1, cell_xface, proc_x_max, tag, &
            phi(  -1,-1,-1), 1, cell_xface, proc_x_min, tag, &
            comm, status, errcode)

        CALL MPI_SENDRECV( &
            phi(-1,   1,-1), 1, cell_yface, proc_y_min, tag, &
            phi(-1,ny+1,-1), 1, cell_yface, proc_y_max, tag, &
            comm, status, errcode)
        CALL MPI_SENDRECV( &
            phi(-1,ny-1,-1), 1, cell_yface, proc_y_max, tag, &
            phi(-1,  -1,-1), 1, cell_yface, proc_y_min, tag, &
            comm, status, errcode)

        CALL MPI_SENDRECV( &
            phi(-1,-1,   1), 1, cell_zface, proc_z_min, tag, &
            phi(-1,-1,nz+1), 1, cell_zface, proc_z_max, tag, &
            comm, status, errcode)
        CALL MPI_SENDRECV( &
            phi(-1,-1,nz-1), 1, cell_zface, proc_z_max, tag, &
            phi(-1,-1,  -1), 1, cell_zface, proc_z_min, tag, &
            comm, status, errcode)

        !Unipolar flux
        local_flux = 0.0_num
        IF (proc_z_min == MPI_PROC_NULL) THEN
          DO iy = 1, ny
            DO ix = 1, nx

              centre = 10.0_num
              radius = 1.0_num
              amp = 1.0_num
              r1 = SQRT((xc(ix)-centre)**2 + yc(iy)**2)
              phi(ix,iy,0) = phi(ix,iy,1) &
                + amp * dzc(1) * (1.0_num - TANH((r1 - radius)/0.2_num))

              phi(ix,iy,-1) = phi(ix,iy,0)
            END DO
          END DO
        END IF

        IF (proc_z_max == MPI_PROC_NULL) THEN
          phi(:,:,nz+1) = 0.0_num
          phi(:,:,nz+2) = 0.0_num
        END IF        
        IF (proc_x_min == MPI_PROC_NULL) THEN
          phi(0,:,:) = phi(1,:,:) 
          phi(-1,:,:) = phi(1,:,:)
        END IF 
        IF (proc_x_max == MPI_PROC_NULL) THEN
          phi(nx+1,:,:) = phi(nx,:,:) 
          phi(nx+2,:,:) = phi(nx,:,:)
        END IF 
        IF (proc_y_min == MPI_PROC_NULL) THEN
          phi(:,0,:) = phi(:,1,:) 
          phi(:,-1,:) = phi(:,1,:)
        END IF 
        IF (proc_y_max == MPI_PROC_NULL) THEN
          phi(:,ny+1,:) = phi(:,ny,:) 
          phi(:,ny+2,:) = phi(:,ny,:)
        END IF 

      END SUBROUTINE phi_mpi

  END SUBROUTINE potential_field



END MODULE initial_conditions
