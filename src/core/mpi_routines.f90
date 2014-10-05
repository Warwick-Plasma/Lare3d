!******************************************************************************
! Routines to set up the MPI routines and allocate dynamic arrays
!******************************************************************************

MODULE mpi_routines

  USE shared_data
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: mpi_initialise, mpi_close

  REAL(dbl) :: start_time, end_time

CONTAINS

  !****************************************************************************
  ! Start up the MPI layer, allocate arrays and set up MPI types
  !****************************************************************************

  SUBROUTINE mpi_initialise

    INTEGER, PARAMETER :: ng = 1
    LOGICAL, PARAMETER :: allow_cpu_reduce = .FALSE.
    INTEGER :: dims(c_ndims), icoord, old_comm, ierr
    LOGICAL :: periods(c_ndims), reorder, reset
    INTEGER :: ix, iy, iz
    INTEGER :: nx0, ny0, nz0
    INTEGER :: nxp, nyp, nzp
    INTEGER :: nxsplit, nysplit, nzsplit
    INTEGER :: x_coords, y_coords, z_coords
    INTEGER :: area, minarea, nprocyz
    INTEGER :: ranges(3,1), nproc_orig, oldgroup, newgroup

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, errcode)

    nproc_orig = nproc

    IF (nx_global < ng .OR. ny_global < ng .OR. nz_global < ng) THEN
      IF (rank == 0) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'Simulation domain is too small.'
        PRINT*,'There must be at least ', ng, ' cells in each direction.'
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
    ENDIF

    reset = .FALSE.
    IF (MAX(nprocx,1) * MAX(nprocy,1) * MAX(nprocz,1) > nproc) THEN
      reset = .TRUE.
    ELSE IF (nprocx * nprocy * nprocz > 0) THEN
      ! Sanity check
      nxsplit = nx_global / nprocx
      nysplit = ny_global / nprocy
      nzsplit = nz_global / nprocz
      IF (nxsplit < ng .OR. nysplit < ng .OR. nzsplit < ng) &
          reset = .TRUE.
    ENDIF

    IF (reset) THEN
      IF (rank == 0) THEN
        PRINT *, 'Unable to use requested processor subdivision. Using ' &
            // 'default division.'
      ENDIF
      nprocx = 0
      nprocy = 0
      nprocz = 0
    ENDIF

    IF (nprocx * nprocy * nprocz == 0) THEN
      DO WHILE (nproc > 1)
        ! Find the processor split which minimizes surface area of
        ! the resulting domain

        minarea = nx_global * ny_global + ny_global * nz_global &
            + nz_global * nx_global

        DO ix = 1, nproc
          nprocyz = nproc / ix
          IF (ix * nprocyz /= nproc) CYCLE

          nxsplit = nx_global / ix
          ! Actual domain must be bigger than the number of ghostcells
          IF (nxsplit < ng) CYCLE

          DO iy = 1, nprocyz
            iz = nprocyz / iy
            IF (iy * iz /= nprocyz) CYCLE

            nysplit = ny_global / iy
            nzsplit = nz_global / iz
            ! Actual domain must be bigger than the number of ghostcells
            IF (nysplit < ng .OR. nzsplit < ng) CYCLE

            area = nxsplit * nysplit + nysplit * nzsplit + nzsplit * nxsplit
            IF (area < minarea) THEN
              nprocx = ix
              nprocy = iy
              nprocz = iz
              minarea = area
            ENDIF
          ENDDO
        ENDDO

        IF (nprocx > 0) EXIT

        ! If we get here then no suitable split could be found. Decrease the
        ! number of processors and try again.

        nproc = nproc - 1
      ENDDO
    ENDIF

    IF (nproc_orig /= nproc) THEN
      IF (.NOT.allow_cpu_reduce) THEN
        IF (rank == 0) THEN
          PRINT*,'*** ERROR ***'
          PRINT*,'Cannot split the domain using the requested number of CPUs.'
          PRINT*,'Try reducing the number of CPUs to ', nproc
        ENDIF
        CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
        STOP
      ENDIF
      IF (rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Cannot split the domain using the requested number of CPUs.'
        PRINT*,'Reducing the number of CPUs to ', nproc
      ENDIF
      ranges(1,1) = nproc
      ranges(2,1) = nproc_orig - 1
      ranges(3,1) = 1
      old_comm = comm
      CALL MPI_COMM_GROUP(old_comm, oldgroup, errcode)
      CALL MPI_GROUP_RANGE_EXCL(oldgroup, 1, ranges, newgroup, errcode)
      CALL MPI_COMM_CREATE(old_comm, newgroup, comm, errcode)
      IF (comm == MPI_COMM_NULL) THEN
        CALL MPI_FINALIZE(errcode)
        STOP
      ENDIF
      CALL MPI_GROUP_FREE(oldgroup, errcode)
      CALL MPI_GROUP_FREE(newgroup, errcode)
      CALL MPI_COMM_FREE(old_comm, errcode)
    ENDIF

    dims = (/nprocz, nprocy, nprocx/)
    CALL MPI_DIMS_CREATE(nproc, c_ndims, dims, errcode)

    IF (PRODUCT(MAX(dims, 1)) > nproc) THEN
      dims = 0
      IF (rank == 0) THEN
        PRINT*, 'Too many processors requested in override.'
        PRINT*, 'Reverting to automatic decomposition.'
        PRINT*, '******************************************'
        PRINT*
      END IF
    END IF

    CALL MPI_DIMS_CREATE(nproc, c_ndims, dims, errcode)

    nprocx = dims(c_ndims  )
    nprocy = dims(c_ndims-1)
    nprocz = dims(c_ndims-2)

    ALLOCATE(cell_nx_mins(0:nprocx-1), cell_nx_maxs(0:nprocx-1))
    ALLOCATE(cell_ny_mins(0:nprocy-1), cell_ny_maxs(0:nprocy-1))
    ALLOCATE(cell_nz_mins(0:nprocz-1), cell_nz_maxs(0:nprocz-1))

    periods = .TRUE.
    reorder = .TRUE.

    IF (xbc_min == BC_OTHER) periods(c_ndims  ) = .FALSE.
    IF (xbc_max == BC_OTHER) periods(c_ndims  ) = .FALSE.
    IF (ybc_min == BC_OTHER) periods(c_ndims-1) = .FALSE.
    IF (ybc_max == BC_OTHER) periods(c_ndims-1) = .FALSE.
    IF (zbc_min == BC_OTHER) periods(c_ndims-2) = .FALSE.
    IF (zbc_max == BC_OTHER) periods(c_ndims-2) = .FALSE.

    IF (xbc_min == BC_OPEN) periods(c_ndims  ) = .FALSE.
    IF (xbc_max == BC_OPEN) periods(c_ndims  ) = .FALSE.
    IF (ybc_min == BC_OPEN) periods(c_ndims-1) = .FALSE.
    IF (ybc_max == BC_OPEN) periods(c_ndims-1) = .FALSE.
    IF (zbc_min == BC_OPEN) periods(c_ndims-2) = .FALSE.
    IF (zbc_max == BC_OPEN) periods(c_ndims-2) = .FALSE.

    CALL MPI_CART_CREATE(MPI_COMM_WORLD, c_ndims, dims, periods, reorder, &
        comm, errcode)

    CALL MPI_COMM_RANK(comm, rank, errcode)
    CALL MPI_CART_COORDS(comm, rank, c_ndims, coordinates, errcode)
    CALL MPI_CART_SHIFT(comm, c_ndims-1, 1, proc_x_min, proc_x_max, errcode)
    CALL MPI_CART_SHIFT(comm, c_ndims-2, 1, proc_y_min, proc_y_max, errcode)
    CALL MPI_CART_SHIFT(comm, c_ndims-3, 1, proc_z_min, proc_z_max, errcode)

    x_coords = coordinates(c_ndims  )
    y_coords = coordinates(c_ndims-1)
    z_coords = coordinates(c_ndims-2)

    ! Create the subarray for this problem: subtype decribes where this
    ! process's data fits into the global picture.

    nx0 = nx_global / nprocx
    ny0 = ny_global / nprocy
    nz0 = nz_global / nprocz

    ! If the number of gridpoints cannot be exactly subdivided then fix
    ! The first nxp processors have nx0 grid points
    ! The remaining processors have nx0+1 grid points
    IF (nx0 * nprocx /= nx_global) THEN
      nxp = (nx0 + 1) * nprocx - nx_global
    ELSE
      nxp = nprocx
    END IF

    IF (ny0 * nprocy /= ny_global) THEN
      nyp = (ny0 + 1) * nprocy - ny_global
    ELSE
      nyp = nprocy
    END IF

    IF (nz0 * nprocz /= nz_global) THEN
      nzp = (nz0 + 1) * nprocz - nz_global
    ELSE
      nzp = nprocz
    END IF

    ! Set up the starting point for my subgrid (assumes arrays start at 0)

    DO icoord = 0, nxp - 1
      cell_nx_mins(icoord) = icoord * nx0 + 1
      cell_nx_maxs(icoord) = (icoord + 1) * nx0
    END DO
    DO icoord = nxp, nprocx - 1
      cell_nx_mins(icoord) = nxp * nx0 + (icoord - nxp) * (nx0 + 1) + 1
      cell_nx_maxs(icoord) = nxp * nx0 + (icoord - nxp + 1) * (nx0 + 1)
    END DO

    DO icoord = 0, nyp - 1
      cell_ny_mins(icoord) = icoord * ny0 + 1
      cell_ny_maxs(icoord) = (icoord + 1) * ny0
    END DO
    DO icoord = nyp, nprocy - 1
      cell_ny_mins(icoord) = nyp * ny0 + (icoord - nyp) * (ny0 + 1) + 1
      cell_ny_maxs(icoord) = nyp * ny0 + (icoord - nyp + 1) * (ny0 + 1)
    END DO

    DO icoord = 0, nzp - 1
      cell_nz_mins(icoord) = icoord * nz0 + 1
      cell_nz_maxs(icoord) = (icoord + 1) * nz0
    END DO
    DO icoord = nzp, nprocz - 1
      cell_nz_mins(icoord) = nzp * nz0 + (icoord - nzp) * (nz0 + 1) + 1
      cell_nz_maxs(icoord) = nzp * nz0 + (icoord - nzp + 1) * (nz0 + 1)
    END DO

    n_global_min(1) = cell_nx_mins(x_coords) - 1
    n_global_max(1) = cell_nx_maxs(x_coords)

    n_global_min(2) = cell_ny_mins(y_coords) - 1
    n_global_max(2) = cell_ny_maxs(y_coords)

    n_global_min(3) = cell_nz_mins(z_coords) - 1
    n_global_max(3) = cell_nz_maxs(z_coords)

    nx = n_global_max(1) - n_global_min(1)
    ny = n_global_max(2) - n_global_min(2)
    nz = n_global_max(3) - n_global_min(3)

    ALLOCATE(energy(-1:nx+2, -1:ny+2, -1:nz+2))
    ALLOCATE(p_visc(-1:nx+2, -1:ny+2, -1:nz+2))
    ALLOCATE(rho(-1:nx+2, -1:ny+2, -1:nz+2))
    ALLOCATE(vx (-2:nx+2, -2:ny+2, -2:nz+2))
    ALLOCATE(vy (-2:nx+2, -2:ny+2, -2:nz+2))
    ALLOCATE(vz (-2:nx+2, -2:ny+2, -2:nz+2))
    ALLOCATE(vx1(-2:nx+2, -2:ny+2, -2:nz+2))
    ALLOCATE(vy1(-2:nx+2, -2:ny+2, -2:nz+2))
    ALLOCATE(vz1(-2:nx+2, -2:ny+2, -2:nz+2))
    ALLOCATE(bx (-2:nx+2, -1:ny+2, -1:nz+2))
    ALLOCATE(by (-1:nx+2, -2:ny+2, -1:nz+2))
    ALLOCATE(bz (-1:nx+2, -1:ny+2, -2:nz+2))
    ALLOCATE(eta(-1:nx+2, -1:ny+2, -1:nz+2))
    IF (rke) ALLOCATE(delta_ke(-1:nx+2, -1:ny+2, -1:nz+2))

    ! Shocked and resistive need to be larger to allow offset = 4 in shock_test
    ALLOCATE(cv(-1:nx+2, -1:ny+2, -1:nz+2), cv1(-1:nx+2, -1:ny+2, -1:nz+2))
    ALLOCATE(xc(-1:nx+2), xb(-2:nx+2), dxc(-1:nx+1), dxb(-1:nx+2))
    ALLOCATE(yc(-1:ny+2), yb(-2:ny+2), dyc(-1:ny+1), dyb(-1:ny+2))
    ALLOCATE(zc(-1:nz+2), zb(-2:nz+2), dzc(-1:nz+1), dzb(-1:nz+2))
    ALLOCATE(grav(-1:nz+2))
    ALLOCATE(jx_r(0:nx+1, 0:ny+1, 0:nz+1))
    ALLOCATE(jy_r(0:nx+1, 0:ny+1, 0:nz+1))
    ALLOCATE(jz_r(0:nx+1, 0:ny+1, 0:nz+1))

    ALLOCATE(xb_global(-2:nx_global+2))
    ALLOCATE(yb_global(-2:ny_global+2))
    ALLOCATE(zb_global(-2:nz_global+2))

    IF (rank == 0) start_time = MPI_WTIME()

    p_visc = 0.0_num
    eta = 0.0_num

    CALL mpi_create_types

  END SUBROUTINE mpi_initialise



  !****************************************************************************
  ! Shutdown the MPI layer, deallocate arrays and set up timing info
  !****************************************************************************

  SUBROUTINE mpi_close

    INTEGER :: seconds, minutes, hours, total

#ifndef NO_IO
    IF (rank == 0) THEN
      end_time = MPI_WTIME()
      total = INT(end_time - start_time)
      seconds = MOD(total, 60)
      minutes = MOD(total / 60, 60)
      hours = total / 3600
      WRITE(stat_unit,*)
      WRITE(stat_unit,'(''runtime = '', i4, ''h '', i2, ''m '', i2, ''s on '', &
          & i4, '' process elements.'')') hours, minutes, seconds, nproc
    END IF
#endif

    CALL MPI_BARRIER(comm, errcode)

    CALL mpi_destroy_types

    DEALLOCATE(rho, energy)
    DEALLOCATE(vx, vy, vz)
    DEALLOCATE(vx1, vy1, vz1)
    DEALLOCATE(bx, by, bz)
    DEALLOCATE(p_visc)
    DEALLOCATE(eta)
    DEALLOCATE(cv, cv1)
    DEALLOCATE(xc, xb, dxb, dxc)
    DEALLOCATE(yc, yb, dyb, dyc)
    DEALLOCATE(zc, zb, dzb, dzc)
    DEALLOCATE(grav)
    DEALLOCATE(jx_r, jy_r, jz_r)
    DEALLOCATE(xb_global, yb_global, zb_global)
    DEALLOCATE(cell_nx_mins, cell_nx_maxs)
    DEALLOCATE(cell_ny_mins, cell_ny_maxs)
    DEALLOCATE(cell_nz_mins, cell_nz_maxs)

    IF (ALLOCATED(xi_n)) DEALLOCATE(xi_n)
    IF (ALLOCATED(delta_ke)) DEALLOCATE(delta_ke)
    IF (ALLOCATED(eta_perp)) DEALLOCATE(eta_perp)
    IF (ALLOCATED(parallel_current)) DEALLOCATE(parallel_current)
    IF (ALLOCATED(perp_current)) DEALLOCATE(perp_current)

  END SUBROUTINE mpi_close



  SUBROUTINE mpi_create_types

    INTEGER :: sizes(c_ndims), subsizes(c_ndims), starts(c_ndims)
    INTEGER :: local_dims(c_ndims), global_dims(c_ndims)
    INTEGER :: local_size(c_ndims), global_size(c_ndims)
    INTEGER :: idir, vdir, mpitype
    INTEGER, PARAMETER :: ng = 2 ! Number of ghost cells

    local_dims = (/ nx, ny, nz /)
    global_dims = (/ nx_global, ny_global, nz_global /)

    ! MPI types for cell-centred variables
    local_size  = local_dims
    global_size = global_dims

    ! File view for cell-centred variables (excluding the ghost cells)
    sizes = global_size
    subsizes = local_size
    starts = n_global_min

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    cell_distribution = mpitype

    ! Subarray for cell-centred variable which has no ghost cells
    sizes = local_size
    starts = 0

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    cellng_subarray = mpitype

    ! Cell-centred array dimensions
    sizes = local_size + 2 * ng

    ! Subarray for cell-centred variable which excludes the ghost cells
    starts = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    cell_subarray = mpitype

    ! ng cells, 1d slice of cell-centred variable

    idir = 1
    subsizes = sizes
    subsizes(idir) = ng
    starts = 0

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    cell_xface = mpitype

    idir = 2
    subsizes = sizes
    subsizes(idir) = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    cell_yface = mpitype

    idir = 3
    subsizes = sizes
    subsizes(idir) = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    cell_zface = mpitype

    ! MPI types for node-centred variables
    local_size  = local_dims + 1
    global_size = global_dims + 1

    ! File view for node-centred variables (excluding the ghost cells)
    sizes = global_size
    subsizes = local_size
    starts = n_global_min

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    node_distribution = mpitype

    ! Subarray for node-centred variable which has no ghost cells
    sizes = local_size
    starts = 0

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    nodeng_subarray = mpitype

    ! Node-centred array dimensions
    sizes = local_size + 2 * ng

    ! Subarray for node-centred variable which excludes the ghost cells
    starts = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    node_subarray = mpitype

    ! ng cells, 1d slice of node-centred variable

    idir = 1
    subsizes = sizes
    subsizes(idir) = ng
    starts = 0

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    node_xface = mpitype

    idir = 2
    subsizes = sizes
    subsizes(idir) = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    node_yface = mpitype

    idir = 3
    subsizes = sizes
    subsizes(idir) = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    node_zface = mpitype

    ! ng+1 cells, 1d slice of node-centred variable

    idir = 1
    subsizes = sizes
    subsizes(idir) = ng + 1

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    node_xface1 = mpitype

    idir = 2
    subsizes = sizes
    subsizes(idir) = ng + 1

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    node_yface1 = mpitype

    idir = 3
    subsizes = sizes
    subsizes(idir) = ng + 1

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    node_zface1 = mpitype

    ! MPI types for Bx-sized variables
    vdir = 1
    local_size  = local_dims
    global_size = global_dims
    local_size (vdir) = local_size (vdir) + 1
    global_size(vdir) = global_size(vdir) + 1

    ! File view for Bx-sized variables (excluding the ghost cells)
    sizes = global_size
    subsizes = local_size
    starts = n_global_min

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    bx_distribution = mpitype

    ! Bx-sized array dimensions
    sizes = local_size + 2 * ng

    ! Subarray for Bx-sized variable which excludes the ghost cells
    starts = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    bx_subarray = mpitype

    ! ng cells, 1d slice of Bx-sized variable

    idir = 1
    subsizes = sizes
    subsizes(idir) = ng
    starts = 0

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    bx_xface = mpitype

    idir = 2
    subsizes = sizes
    subsizes(idir) = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    bx_yface = mpitype

    idir = 3
    subsizes = sizes
    subsizes(idir) = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    bx_zface = mpitype

    ! ng+1 cells, 1d slice of Bx-sized variable

    idir = vdir
    subsizes = sizes
    subsizes(idir) = ng + 1

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    bx_xface1 = mpitype

    ! MPI types for By-sized variables
    vdir = 2
    local_size  = local_dims
    global_size = global_dims
    local_size (vdir) = local_size (vdir) + 1
    global_size(vdir) = global_size(vdir) + 1

    ! File view for By-sized variables (excluding the ghost cells)
    sizes = global_size
    subsizes = local_size
    starts = n_global_min

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    by_distribution = mpitype

    ! By-sized array dimensions
    sizes = local_size + 2 * ng

    ! Subarray for By-sized variable which excludes the ghost cells
    starts = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    by_subarray = mpitype

    ! ng cells, 1d slice of By-sized variable

    idir = 1
    subsizes = sizes
    subsizes(idir) = ng
    starts = 0

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    by_xface = mpitype

    idir = 2
    subsizes = sizes
    subsizes(idir) = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    by_yface = mpitype

    idir = 3
    subsizes = sizes
    subsizes(idir) = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    by_zface = mpitype

    ! ng+1 cells, 1d slice of By-sized variable

    idir = vdir
    subsizes = sizes
    subsizes(idir) = ng + 1

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    by_yface1 = mpitype

    ! MPI types for Bz-sized variables
    vdir = 3
    local_size  = local_dims
    global_size = global_dims
    local_size (vdir) = local_size (vdir) + 1
    global_size(vdir) = global_size(vdir) + 1

    ! File view for Bz-sized variables (excluding the ghost cells)
    sizes = global_size
    subsizes = local_size
    starts = n_global_min

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    bz_distribution = mpitype

    ! Bz-sized array dimensions
    sizes = local_size + 2 * ng

    ! Subarray for Bz-sized variable which excludes the ghost cells
    starts = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    bz_subarray = mpitype

    ! ng cells, 1d slice of Bz-sized variable

    idir = 1
    subsizes = sizes
    subsizes(idir) = ng
    starts = 0

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    bz_xface = mpitype

    idir = 2
    subsizes = sizes
    subsizes(idir) = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    bz_yface = mpitype

    idir = 3
    subsizes = sizes
    subsizes(idir) = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    bz_zface = mpitype

    ! ng+1 cells, 1d slice of Bz-sized variable

    idir = vdir
    subsizes = sizes
    subsizes(idir) = ng + 1

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    bz_zface1 = mpitype

  END SUBROUTINE mpi_create_types



  SUBROUTINE mpi_destroy_types

    CALL MPI_TYPE_FREE(cell_xface, errcode)
    CALL MPI_TYPE_FREE(cell_yface, errcode)
    CALL MPI_TYPE_FREE(cell_zface, errcode)
    CALL MPI_TYPE_FREE(node_xface, errcode)
    CALL MPI_TYPE_FREE(node_yface, errcode)
    CALL MPI_TYPE_FREE(node_zface, errcode)
    CALL MPI_TYPE_FREE(node_xface1, errcode)
    CALL MPI_TYPE_FREE(node_yface1, errcode)
    CALL MPI_TYPE_FREE(node_zface1, errcode)
    CALL MPI_TYPE_FREE(bx_xface, errcode)
    CALL MPI_TYPE_FREE(bx_yface, errcode)
    CALL MPI_TYPE_FREE(bx_zface, errcode)
    CALL MPI_TYPE_FREE(by_xface, errcode)
    CALL MPI_TYPE_FREE(by_yface, errcode)
    CALL MPI_TYPE_FREE(by_zface, errcode)
    CALL MPI_TYPE_FREE(bz_xface, errcode)
    CALL MPI_TYPE_FREE(bz_yface, errcode)
    CALL MPI_TYPE_FREE(bz_zface, errcode)
    CALL MPI_TYPE_FREE(bx_xface1, errcode)
    CALL MPI_TYPE_FREE(by_yface1, errcode)
    CALL MPI_TYPE_FREE(bz_zface1, errcode)
    CALL MPI_TYPE_FREE(cell_subarray, errcode)
    CALL MPI_TYPE_FREE(node_subarray, errcode)
    CALL MPI_TYPE_FREE(cellng_subarray, errcode)
    CALL MPI_TYPE_FREE(nodeng_subarray, errcode)
    CALL MPI_TYPE_FREE(cell_distribution, errcode)
    CALL MPI_TYPE_FREE(node_distribution, errcode)
    CALL MPI_TYPE_FREE(bx_subarray, errcode)
    CALL MPI_TYPE_FREE(by_subarray, errcode)
    CALL MPI_TYPE_FREE(bz_subarray, errcode)
    CALL MPI_TYPE_FREE(bx_distribution, errcode)
    CALL MPI_TYPE_FREE(by_distribution, errcode)
    CALL MPI_TYPE_FREE(bz_distribution, errcode)

  END SUBROUTINE mpi_destroy_types

END MODULE mpi_routines
