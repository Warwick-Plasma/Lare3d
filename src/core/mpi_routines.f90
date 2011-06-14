MODULE mpi_routines

  USE shared_data
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: mpi_initialise, mpi_close

  REAL(dbl) :: start_time, end_time

CONTAINS

  SUBROUTINE mpi_initialise

    INTEGER :: ndims, dims(3)
    LOGICAL :: periods(3), reorder
    INTEGER :: starts(3), sizes(3), subsizes(3)

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, errcode)

    ndims = 3
    dims = (/ nprocz, nprocy, nprocx /)

    IF (MAX(dims(1), 1)*MAX(dims(2), 1)*MAX(dims(3), 1) .GT. nproc) THEN
      dims = 0
      IF (rank .EQ. 0) THEN
        PRINT *, "Too many processors requested in override."
        PRINT *, "Reverting to automatic decomposition."
        PRINT *, "******************************************"
        PRINT *, ""
      END IF
    END IF

    CALL MPI_DIMS_CREATE(nproc, ndims, dims, errcode)

    nprocx = dims(3)
    nprocy = dims(2)
    nprocz = dims(1)
    nx = nx_global / nprocx
    ny = ny_global / nprocy
    nz = nz_global / nprocz

    periods = .TRUE.
    reorder = .TRUE.

    IF (xbc_left == BC_OTHER) periods(3) = .FALSE.
    IF (ybc_up == BC_OTHER) periods(2) = .FALSE.
    IF (zbc_front == BC_OTHER) periods(1) = .FALSE.
    IF (xbc_left == BC_OPEN) periods(3) = .FALSE.
    IF (ybc_up == BC_OPEN) periods(2) = .FALSE.
    IF (zbc_front == BC_OPEN) periods(1) = .FALSE. 
    
    IF (xbc_right == BC_OTHER) periods(3) = .FALSE.
    IF (ybc_down == BC_OTHER) periods(2) = .FALSE.
    IF (zbc_back == BC_OTHER) periods(1) = .FALSE.
    IF (xbc_right == BC_OPEN) periods(3) = .FALSE.
    IF (ybc_down == BC_OPEN) periods(2) = .FALSE.
    IF (zbc_back == BC_OPEN) periods(1) = .FALSE.

    CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods, &
        reorder, comm, errcode)

    CALL MPI_COMM_RANK(comm, rank, errcode)
    CALL MPI_CART_COORDS(comm, rank, 3, coordinates, errcode)
    CALL MPI_CART_SHIFT(comm, 2, 1, left, right, errcode)
    CALL MPI_CART_SHIFT(comm, 1, 1, down, up, errcode)
    CALL MPI_CART_SHIFT(comm, 0, 1, front, back, errcode)

    ! Create the subarray for this problem: subtype decribes where this
    ! process's data fits into the global picture.

    ! set up the starting point for my subgrid (assumes arrays start at 0)
    starts(1) = coordinates(3) * nx
    starts(2) = coordinates(2) * ny
    starts(3) = coordinates(1) * nz

    ! the grid sizes
    subsizes = (/ nx+1, ny+1, nz+1 /)
    sizes = (/ nx_global+1, ny_global+1, nz_global+1 /)

    ! set up and commit the subarray type
    CALL MPI_TYPE_CREATE_SUBARRAY(ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, subtype, errcode)

    CALL MPI_TYPE_COMMIT(subtype, errcode)

    ! Calculate initial displacement value:
    ! nx, ny, nz, (xb, yb, zb, time) * size of float
    initialdisp = 12 + (nx_global + ny_global + 3) * num

    ALLOCATE(rho(-1:nx+2, -1:ny+2, -1:nz+2))
    ALLOCATE(energy(-1:nx+2, -1:ny+2, -1:nz+2))
    ALLOCATE(vx(-2:nx+2, -2:ny+2, -2:nz+2))
    ALLOCATE(vy(-2:nx+2, -2:ny+2, -2:nz+2))
    ALLOCATE(vz(-2:nx+2, -2:ny+2, -2:nz+2))
    ALLOCATE(vx1(-2:nx+2, -2:ny+2, -2:nz+2))
    ALLOCATE(vy1(-2:nx+2, -2:ny+2, -2:nz+2))
    ALLOCATE(vz1(-2:nx+2, -2:ny+2, -2:nz+2))
    ALLOCATE(bx(-2:nx+2, -1:ny+2, -1:nz+2))
    ALLOCATE(by(-1:nx+2, -2:ny+2, -1:nz+2))
    ALLOCATE(bz(-1:nx+2, -1:ny+2, -2:nz+2))
    ALLOCATE(delta_ke(-1:nx+2, -1:ny+2, -1:nz+2))
    ALLOCATE(p_visc(-1:nx+2, -1:ny+2, -1:nz+2))
    ALLOCATE(eta(-1:nx+2, -1:ny+2, -1:nz+2))
    ALLOCATE(lambda_i(-1:nx+2, -1:ny+2, -1:nz+2))
    ! shocked and resistive need to be larger to allow offset = 4 in shock_test
    ALLOCATE(cv(-1:nx+2, -1:ny+2, -1:nz+2), cv1(-1:nx+2, -1:ny+2, -1:nz+2))
    ALLOCATE(xc(-1:nx+2), xb(-2:nx+2), dxb(-1:nx+2), dxc(-1:nx+2))
    ALLOCATE(yc(-1:ny+2), yb(-2:ny+2), dyb(-1:ny+2), dyc(-1:ny+2))
    ALLOCATE(zc(-1:nz+2), zb(-2:nz+2), dzb(-1:nz+2), dzc(-1:nz+2))
    ALLOCATE(grav(-1:nz+2))
    ALLOCATE(jx_r(0:nx+1, 0:ny+1, 0:nz+1))
    ALLOCATE(jy_r(0:nx+1, 0:ny+1, 0:nz+1))
    ALLOCATE(jz_r(0:nx+1, 0:ny+1, 0:nz+1))

    IF (rank == 0) start_time = MPI_WTIME()

    p_visc = 0.0_num
    eta = 0.0_num

  END SUBROUTINE mpi_initialise



  SUBROUTINE mpi_close

    INTEGER :: seconds, minutes, hours, total

    IF (rank == 0) THEN
      end_time = MPI_WTIME()
      total = INT(end_time - start_time)
      seconds = MOD(total, 60)
      minutes = MOD(total / 60, 60)
      hours = total / 3600
      WRITE(20, *)
      WRITE(20, '("runtime = ", i4, "h ", i2, "m ", i2, &
          & "s on ", i4, " process elements.")') hours, minutes, seconds, nproc
    END IF

    CALL MPI_BARRIER(comm, errcode)

    DEALLOCATE(rho, energy)
    DEALLOCATE(vx, vy, vz)
    DEALLOCATE(vx1, vy1, vz1)
    DEALLOCATE(bx, by, bz)
    DEALLOCATE(delta_ke, p_visc)
    DEALLOCATE(eta, lambda_i)
    DEALLOCATE(cv, cv1)
    DEALLOCATE(xc, xb, dxb, dxc)
    DEALLOCATE(yc, yb, dyb, dyc)
    DEALLOCATE(grav)
    DEALLOCATE(jx_r, jy_r, jz_r)  
    
    IF (ALLOCATED(xi_n)) DEALLOCATE(xi_n)

  END SUBROUTINE mpi_close

END MODULE mpi_routines
