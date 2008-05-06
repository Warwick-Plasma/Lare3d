MODULE mpi_routines
  USE shared_data
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: mpi_initialise, mpi_close,mpi_minimal_init

  REAL(dbl) :: start_time, end_time

CONTAINS

  SUBROUTINE mpi_minimal_init
    INTEGER s

    CALL MPI_INIT(errcode)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, errcode)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, errcode)

    IF (rank .EQ. 0) THEN
      open(unit=10,file="./input.deck",iostat=s,status='OLD')
      deckfile=(s .EQ. 0)
      close(10)
    ENDIF

  END SUBROUTINE mpi_minimal_init

  SUBROUTINE mpi_initialise

    INTEGER :: ndims, dims(3)
    LOGICAL :: periods(3), reorder
    INTEGER :: starts(3), sizes(3), subsizes(3)

    ndims = 3
    dims = (/nprocz,nprocy,nprocx/)
    IF (MAX(dims(1),1)*MAX(dims(2),1)*MAX(dims(3),1) .GT. nproc) THEN
       dims=0
       IF (rank .EQ. 0) THEN
          PRINT *,"Too many processors requested in override."
          PRINT *,"Reverting to automatic decomposition."
          PRINT *,"******************************************"
          PRINT *,""
       ENDIF
    ENDIF

    CALL MPI_DIMS_CREATE(nproc, ndims, dims, errcode)
    nprocx = dims(3)
    nprocy = dims(2)
    nprocz = dims(1)

    CALL check_mpi_decomposition

    periods = .TRUE.
    reorder = .TRUE.
    IF (xbc_left /= periodic) periods(3) = .FALSE.
    IF (ybc_up /= periodic) periods(2) = .FALSE.
    IF (zbc_front /= periodic) periods(1) = .FALSE.
    IF (xbc_left /= periodic) periods(3) = .FALSE.
    IF (ybc_up /= periodic) periods(2) = .FALSE.
    IF (zbc_front /= periodic) periods(1) = .FALSE.

    CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods,  &
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
    subsizes = (/nx+1,ny+1,nz+1/)
    sizes = (/nx_global+1, ny_global+1, nz_global+1/)

! set up and commit the subarray type
    CALL MPI_TYPE_CREATE_SUBARRAY(ndims, sizes, subsizes, starts, &
         MPI_ORDER_FORTRAN, mpireal, subtype, errcode)

    CALL MPI_TYPE_COMMIT(subtype,errcode)

! Calculate initial displacement value: nx,ny,nz,(xb,yb,zb,time)*size of float
    initialdisp = 16 + (nx_global + ny_global + nz_global + 4) * num


    ALLOCATE(rho(-1:nx+2,-1:ny+2,-1:nz+2), energy(-1:nx+2,-1:ny+2,-1:nz+2))
    ALLOCATE(vx(-2:nx+2,-2:ny+2,-2:nz+2), vy(-2:nx+2,-2:ny+2,-2:nz+2), &
         vz(-2:nx+2,-2:ny+2,-2:nz+2))
    ALLOCATE(vx1(-2:nx+2,-2:ny+2,-2:nz+2), vy1(-2:nx+2,-2:ny+2,-2:nz+2), &
         vz1(-2:nx+2,-2:ny+2,-2:nz+2))
    ALLOCATE(bx(-2:nx+2,-1:ny+2,-1:nz+2), by(-1:nx+2,-2:ny+2,-1:nz+2), &
         bz(-1:nx+2,-1:ny+2,-2:nz+2))
    ALLOCATE(delta_ke(-1:nx+2,-1:ny+2,-1:nz+2), p_visc(-1:nx+2,-1:ny+2,-1:nz+2))
    ALLOCATE(eta(-1:nx+2,-1:ny+2,-1:nz+2), curlb(-1:nx+2,-1:ny+2,-1:nz+2))
    ALLOCATE(bzone(-1:nx+2,-1:ny+2,-1:nz+2))

    ALLOCATE(dv_left(-2:ny+2,-2:nz+2), dv_right(-2:ny+2,-2:nz+2))
    ALLOCATE(dv_up(-2:nx+2,-2:nz+2), dv_down(-2:nx+2,-2:nz+2))
    ALLOCATE(dv_back(-2:nx+2,-2:ny+2), dv_front(-2:nx+2,-2:ny+2))

    ALLOCATE(cv(-1:nx+2,-1:ny+2,-1:nz+2), cv1(-1:nx+2,-1:ny+2,-1:nz+2))

    ALLOCATE(xc(-1:nx+2), xb(-2:nx+2), dxb(-1:nx+2), dxc(-1:nx+2), grav(-2:nz+2))
    ALLOCATE(yc(-1:ny+2), yb(-2:ny+2), dyb(-1:ny+2), dyc(-1:ny+2))
    ALLOCATE(zc(-1:nz+2), zb(-2:nz+2), dzb(-1:nz+2), dzc(-1:nz+2))

    IF(rank == 0) start_time = MPI_WTIME()

    ! The following are set to zero as they are used even when run as fully ionized
    ! This is memory/disk inefficient - might get sorted out in the future

  END SUBROUTINE mpi_initialise



  SUBROUTINE mpi_close
    INTEGER :: seconds, minutes, hours, total

    IF(rank == 0) THEN
       end_time = MPI_WTIME()
       total = INT(end_time - start_time)
       seconds = MOD(total,60)
       minutes = MOD(total / 60, 60)
       hours = total / 3600
       WRITE(20,*)
       WRITE(20,'("runtime = ",i4,"h ",i2,"m ",i2,"s on ",i4," process elements.")') &
            hours,minutes,seconds,nproc
    END IF

    CALL MPI_BARRIER(comm, errcode)
    
    DEALLOCATE(rho, energy)
    DEALLOCATE(vx, vy, vz)
    DEALLOCATE(vx1, vy1, vz1)
    DEALLOCATE(bx, by, bz)
    DEALLOCATE(delta_ke, p_visc)
    DEALLOCATE(eta, curlb, bzone)
    DEALLOCATE(dv_left, dv_right)
    DEALLOCATE(dv_up, dv_down)
    DEALLOCATE(dv_back, dv_front)
    DEALLOCATE(cv, cv1,grav)
    DEALLOCATE(xc, xb, dxb, dxc)
    DEALLOCATE(yc, yb, dyb, dyc)
    DEALLOCATE(zc, zb, dzb, dzc)
    DEALLOCATE(xb_global,yb_global,zb_global)

  END SUBROUTINE mpi_close



  SUBROUTINE check_mpi_decomposition

    nz = nz_global / nprocz
    IF (nz * nprocz /= nz_global) THEN
       WRITE(*,*) 'nz * nprocy /= nz_global'
       STOP
    END IF
    ny = ny_global / nprocy
    IF (ny * nprocy /= ny_global) THEN
       WRITE(*,*) 'ny * nprocy /= ny_global'
       STOP
    END IF
    nx = nx_global / nprocx
    IF (nx * nprocx /= nx_global) THEN
       WRITE(*,*) 'nx * nprocx /= nx_global'
       STOP
    END IF

    IF (nx<2 .OR. ny<2 .OR. nz<2) THEN
       WRITE(*,*) 'decomposition too fine'
       STOP
     END IF

  END SUBROUTINE check_mpi_decomposition


END MODULE mpi_routines


