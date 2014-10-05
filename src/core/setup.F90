MODULE setup

  USE shared_data
  USE iocommon
  USE iocontrol
  USE input
  USE input_cartesian

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: before_control, after_control
  PUBLIC :: grid
  PUBLIC :: open_files, close_files, restart_data

CONTAINS

  !****************************************************************************
  ! Any variables that need to have default values before the user specified
  ! ones in control.f90
  !****************************************************************************

  SUBROUTINE before_control

    ! Setup basic variables which have to have default values

    nprocx = 0
    nprocy = 0
    nprocz = 0

    time = 0.0_num
    gamma = 5.0_num / 3.0_num

  END SUBROUTINE before_control



  !****************************************************************************
  ! Variables which need to be specified after the control.f90 is run,
  ! MPI has been setup
  !****************************************************************************

  SUBROUTINE after_control

    ! Setup arrays and other variables which can only be set after
    ! user input

    IF (IAND(initial, IC_RESTART) == 0) restart_snapshot = 0

    p_visc = 0.0_num
    eta = 0.0_num
    grav = 0.0_num

    rho = 0.0_num
    energy = 0.0_num
    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num
    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num

  END SUBROUTINE after_control



  !****************************************************************************
  ! Stretched and staggered grid
  !****************************************************************************

  SUBROUTINE grid

    REAL(num) :: dx, dy, dz
    INTEGER :: ix, iy, iz

    length_x = x_max - x_min
    length_y = y_max - y_min
    length_z = z_max - z_min

    ! Initially assume uniform grid
    dx = length_x / REAL(nx_global, num)
    dy = length_y / REAL(ny_global, num)
    dz = length_z / REAL(nz_global, num)

    ! Grid cell boundary for x coordinates
    DO ix = -2, nx_global + 2
      xb_global(ix) = x_min + REAL(ix, num) * dx
    END DO

    IF (x_stretch) CALL stretch_x ! stretch grid ?

    ! Define position of ghost cells using sizes of adjacent cells
    IF (xbc_max == BC_PERIODIC) THEN
      xb_global(nx_global+1) = xb_global(nx_global) &
          + (xb_global(1) - xb_global(0))
      xb_global(nx_global+2) = xb_global(nx_global) &
          + (xb_global(2) - xb_global(0))
      xb_global(-1) = xb_global(0) &
          - (xb_global(nx_global) - xb_global(nx_global-1))
      xb_global(-2) = xb_global(0) &
          - (xb_global(nx_global) - xb_global(nx_global-2))
    ELSE
      xb_global(nx_global+1) = 2.0_num * xb_global(nx_global) &
          - xb_global(nx_global-1)
      xb_global(nx_global+2) = 2.0_num * xb_global(nx_global) &
          - xb_global(nx_global-2)
      xb_global(-1) = 2.0_num * xb_global(0) - xb_global(1)
      xb_global(-2) = 2.0_num * xb_global(0) - xb_global(2)
    END IF

    ! Setup local grid

    xb(-2) = xb_global(-2+n_global_min(1))

    DO ix = -1, nx + 2
      ixm = ix - 1

      xb(ix) = xb_global(ix+n_global_min(1))
      ! Cell centre
      xc(ix) = 0.5_num * (xb(ixm) + xb(ix))
      ! Cell width
      dxb(ix) = xb(ix) - xb(ixm)
    END DO

    DO ix = -1, nx + 1
      ixp = ix + 1
      ! Distance between centres
      dxc(ix) = xc(ixp) - xc(ix)
    END DO

    ! Repeat for y

    DO iy = -2, ny_global + 2
      yb_global(iy) = y_min + REAL(iy, num) * dy
    END DO

    IF (y_stretch) CALL stretch_y

    IF (ybc_max == BC_PERIODIC) THEN
      yb_global(ny_global+1) = yb_global(ny_global) &
          + (yb_global(1) - yb_global(0))
      yb_global(ny_global+2) = yb_global(ny_global) &
          + (yb_global(2) - yb_global(0))
      yb_global(-1) = yb_global(0) &
          - (yb_global(ny_global) - yb_global(ny_global-1))
      yb_global(-2) = yb_global(0) &
          - (yb_global(ny_global) - yb_global(ny_global-2))
    ELSE
      yb_global(ny_global+1) = 2.0_num * yb_global(ny_global) &
          - yb_global(ny_global-1)
      yb_global(ny_global+2) = 2.0_num * yb_global(ny_global) &
          - yb_global(ny_global-2)
      yb_global(-1) = 2.0_num * yb_global(0) - yb_global(1)
      yb_global(-2) = 2.0_num * yb_global(0) - yb_global(2)
    END IF

    yb(-2) = yb_global(-2+n_global_min(2))

    DO iy = -1, ny + 2
      iym = iy - 1

      yb(iy) = yb_global(iy+n_global_min(2))
      ! Cell centre
      yc(iy) = 0.5_num * (yb(iym) + yb(iy))
      ! Cell width
      dyb(iy) = yb(iy) - yb(iym)
    END DO

    DO iy = -1, ny + 1
      iyp = iy + 1
      dyc(iy) = yc(iyp) - yc(iy)
    END DO

    ! Repeat for z

    DO iz = -2, nz_global + 2
      zb_global(iz) = z_min + REAL(iz, num) * dz
    END DO

    IF (z_stretch) CALL stretch_z

    IF (zbc_max == BC_PERIODIC) THEN
      zb_global(nz_global+1) = zb_global(nz_global) &
          + (zb_global(1) - zb_global(0))
      zb_global(nz_global+2) = zb_global(nz_global) &
          + (zb_global(2) - zb_global(0))
      zb_global(-1) = zb_global(0) &
          - (zb_global(nz_global) - zb_global(nz_global-1))
      zb_global(-2) = zb_global(0) &
          - (zb_global(nz_global) - zb_global(nz_global-2))
    ELSE
      zb_global(nz_global+1) = 2.0_num * zb_global(nz_global) &
          - zb_global(nz_global-1)
      zb_global(nz_global+2) = 2.0_num * zb_global(nz_global) &
          - zb_global(nz_global-2)
      zb_global(-1) = 2.0_num * zb_global(0) - zb_global(1)
      zb_global(-2) = 2.0_num * zb_global(0) - zb_global(2)
    END IF

    zb(-2) = zb_global(-2+n_global_min(3))

    DO iz = -1, nz + 2
      izm = iz - 1

      zb(iz) = zb_global(iz+n_global_min(3))
      ! Cell centre
      zc(iz) = 0.5_num * (zb(izm) + zb(iz))
      ! Cell width
      dzb(iz) = zb(iz) - zb(izm)
    END DO

    DO iz = -1, nz + 1
      izp = iz + 1
      dzc(iz) = zc(izp) - zc(iz)
    END DO

    ! Define the cell area
    DO iz = -1, nz + 2
      DO iy = -1, ny + 2
        DO ix = -1, nx + 2
          cv(ix,iy,iz) = dxb(ix) * dyb(iy) * dzb(iz)
        END DO
      END DO
    END DO

  END SUBROUTINE grid



  !****************************************************************************
  ! Subroutine stretches the grid in the x direction
  ! Replace with any stretching algorithm as needed
  !****************************************************************************

  SUBROUTINE stretch_x

    REAL(num) :: width, dx, L, f, lx_new
    REAL(num), DIMENSION(:), ALLOCATABLE :: dxnew

    ALLOCATE(dxnew(-2:nx_global+2))

    ! New total length
    lx_new = 200.0_num

    ! Centre of tanh stretching in unstretched coordinates
    L = length_x / 1.5_num

    ! Width of tanh stretching in unstretched coordinates
    width = length_x / 10.0_num

    f = (lx_new - length_x) / (length_x - L) / 2.0_num

    dx = length_x / REAL(nx_global, num)
    dxnew = dx + f * (1.0_num + TANH((ABS(xb_global) - L) / width)) * dx

    DO ix = 1, nx_global + 2
      xb_global(ix) = xb_global(ix-1) + dxnew(ix)
    END DO

    length_x = lx_new

    DEALLOCATE(dxnew)

  END SUBROUTINE stretch_x



  !****************************************************************************
  ! Subroutine stretches the domain in the y direction
  ! Stretch domain upwards only
  !****************************************************************************

  SUBROUTINE stretch_y

    REAL(num) :: width, dy, L, f, ly_new
    REAL(num), DIMENSION(:), ALLOCATABLE :: dynew

    ALLOCATE(dynew(-2:ny_global+2))

    ! New total length
    ly_new = 100.0_num

    ! Centre of tanh stretching in unstretched coordinates
    L = length_y / 1.5_num

    ! Width of tanh stretching in unstretched coordinates
    width = length_y / 10.0_num

    f = (ly_new - length_y) / (length_y - L) / 2.0_num

    dy = length_y / REAL(ny_global, num)
    dynew = dy + f * (1.0_num + TANH((ABS(yb_global) - L) / width)) * dy

    DO iy = 1, ny_global + 2
      yb_global(iy) = yb_global(iy-1) + dynew(iy)
    END DO

    length_y = ly_new

    DEALLOCATE(dynew)

  END SUBROUTINE stretch_y



  !****************************************************************************
  ! Subroutine stretches the domain in the z direction
  ! Stretch domain upwards only
  !****************************************************************************

  SUBROUTINE stretch_z

    REAL(num) :: width, dz, L, f, lz_new
    REAL(num), DIMENSION(:), ALLOCATABLE :: dznew

    ALLOCATE(dznew(-2:nz_global+2))

    ! New total length
    lz_new = 100.0_num

    ! Centre of tanh stretching in unstretched coordinates
    L = length_z / 1.5_num

    ! Width of tanh stretching in unstretched coordinates
    width = length_z / 10.0_num

    f = (lz_new - length_z) / (length_z - L) / 2.0_num

    dz = length_z / REAL(nz_global, num)
    dznew = dz + f * (1.0_num + TANH((ABS(zb_global) - L) / width)) * dz

    DO iz = 1, nz_global + 2
      zb_global(iz) = zb_global(iz-1) + dznew(iz)
    END DO

    length_z = lz_new

    DEALLOCATE(dznew)

  END SUBROUTINE stretch_z



  !****************************************************************************
  ! Open the output diagnostic files.
  !****************************************************************************

  SUBROUTINE open_files

    CHARACTER(LEN=11+data_dir_max_length) :: file2
    CHARACTER(LEN=7+data_dir_max_length) :: file3
    INTEGER :: ios

    IF (rank == 0) THEN
      WRITE(file2, '(a, ''/lare3d.dat'')') TRIM(data_dir)
      OPEN(UNIT=stat_unit, STATUS='REPLACE', FILE=file2, iostat=ios)

      IF (ios /= 0) THEN
        PRINT*, 'Unable to open file lare3d.dat for writing. This is ', &
                'most commonly caused by the output directory not existing'
        PRINT*
        PRINT*
        CALL MPI_ABORT(comm, errcode)
      END IF

      WRITE(file3, '(a, ''/en.dat'')') TRIM(data_dir)
      OPEN(UNIT=en_unit, STATUS='REPLACE', FILE=file3, &
          FORM='UNFORMATTED', ACCESS='STREAM', iostat=ios)

      IF (ios /= 0) THEN
        PRINT*, 'Unable to open file en.dat for writing. This is ', &
                'most commonly caused by the output directory not existing'
        PRINT*
        PRINT*
        CALL MPI_ABORT(comm, errcode)
      END IF
    END IF

  END SUBROUTINE open_files



  !****************************************************************************
  ! Close the output diagnostic files
  !****************************************************************************

  SUBROUTINE close_files

    IF (rank == 0) THEN
      CLOSE(UNIT=stat_unit)
      CLOSE(UNIT=en_unit)
    END IF

  END SUBROUTINE close_files



  !****************************************************************************
  ! Function to compare two strings.
  !****************************************************************************

  FUNCTION str_cmp(str_in, str_test)

    CHARACTER(*), INTENT(IN) :: str_in, str_test
    CHARACTER(30) :: str_trim
    LOGICAL :: str_cmp

    str_trim = TRIM(ADJUSTL(str_in))

    IF (LEN(str_test) > LEN(str_in)) THEN
      str_cmp = .FALSE.
      RETURN
    END IF

    IF (str_trim(LEN(str_test)+1:LEN(str_test)+1) /= ' ') THEN
      str_cmp = .FALSE.
      RETURN
    END IF

    str_cmp = str_trim(1:LEN(str_test)) == str_test

  END FUNCTION str_cmp



  !****************************************************************************
  ! Restart the code from a previous output dump.
  !****************************************************************************

  SUBROUTINE restart_data

    CHARACTER(LEN=20+data_dir_max_length) :: filename
    CHARACTER(LEN=20) :: name, class, mesh_name, mesh_class
    INTEGER :: nblocks, type, nd, sof, snap
    INTEGER, DIMENSION(c_ndims) :: dims, global_dims
    REAL(dbl) :: time_d
    REAL(num), DIMENSION(c_ndims) :: extent
    REAL(num), DIMENSION(c_ndims) :: stagger
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: data

    ! Create the filename for the last snapshot
#ifdef MHDCLUSTER
    WRITE(filename, '(''nfs:'', a, ''/'', i4.4, ''.cfd'')') &
        TRIM(data_dir), restart_snapshot
#else
    WRITE(filename, '(a, ''/'', i4.4, ''.cfd'')') &
        TRIM(data_dir), restart_snapshot
#endif

    file_number = restart_snapshot
    global_dims = (/ nx_global+1, ny_global+1, nz_global+1 /)

    ALLOCATE(data(0:nx, 0:ny, 0:nz))
    CALL cfd_open(filename, rank, comm, MPI_MODE_RDONLY)
    ! Open the file
    nblocks = cfd_get_nblocks()

    DO ix = 1, nblocks
      CALL cfd_get_next_block_info_all(name, class, type)
      IF (rank == 0) PRINT*, ix, name, class, type

      IF (type == TYPE_SNAPSHOT) THEN
        CALL cfd_get_snapshot(time_d, snap)
        time = time_d
      END IF

      IF (type == TYPE_MESH) THEN
        ! Strangely, LARE doesn't actually read in the grid from a file
        ! This can be fixed, but for the moment, just go with the flow and
        ! replicate the old behaviour

        CALL cfd_skip_block
      ELSE IF (type == TYPE_MESH_VARIABLE) THEN
        CALL cfd_get_common_meshtype_metadata_all(type, nd, sof)

        IF (nd /= DIMENSION_3D) THEN
          IF (rank == 0) PRINT*, 'Non 3D Dataset found in input file, ', &
              'ignoring and continuting.'
          CALL cfd_skip_block
          CYCLE
        END IF

        IF (type /= VAR_CARTESIAN) THEN
          IF (rank == 0) PRINT*, 'Non - Cartesian variable block found ', &
              'in file, ignoring and continuing'
          CALL cfd_skip_block
          CYCLE
        END IF

        ! We now have a valid variable, let's load it up
        ! First error trapping
        CALL cfd_get_nd_cartesian_variable_metadata_all(nd, dims, extent, &
            stagger, mesh_name, mesh_class)

        IF (ANY(dims(1:c_ndims) /= global_dims(1:c_ndims))) THEN
          IF (rank == 0) PRINT*, 'Size of grid represented by one more ', &
              'variables invalid. Continuing'
          CALL cfd_skip_block
          CYCLE
        END IF

        IF (sof /= num) THEN
          IF (rank == 0) PRINT*, 'Precision of data does not match ', &
              'precision of code. Continuing.'
          CALL cfd_skip_block
        END IF

        ! We're not interested in the other parameters, so if we're here,
        ! load up the data

        CALL cfd_get_3d_cartesian_variable_parallel(data, subtype)

        ! Now have the data, just copy it to correct place

        IF (str_cmp(name(1:3), 'Rho')) THEN
          rho(0:nx, 0:ny, 0:nz) = data
        END IF

        IF (str_cmp(name(1:6), 'Energy')) THEN
          energy(0:nx, 0:ny, 0:nz) = data
        END IF

        IF (str_cmp(name(1:2), 'Vx')) THEN
          vx(0:nx, 0:ny, 0:nz) = data
        END IF

        IF (str_cmp(name(1:2), 'Vy')) THEN
          vy(0:nx, 0:ny, 0:nz) = data
        END IF

        IF (str_cmp(name(1:2), 'Vz')) THEN
          vz(0:nx, 0:ny, 0:nz) = data
        END IF

        IF (str_cmp(name(1:2), 'Bx')) THEN
          bx(0:nx, 0:ny, 0:nz) = data
        END IF

        IF (str_cmp(name(1:2), 'By')) THEN
          by(0:nx, 0:ny, 0:nz) = data
        END IF

        IF (str_cmp(name(1:2), 'Bz')) THEN
          bz(0:nx, 0:ny, 0:nz) = data
        END IF

        ! Should be at end of block, but force the point anyway
        CALL cfd_skip_block
      ELSE
        ! Unknown block, just skip it
        CALL cfd_skip_block
      END IF
    END DO

    DEALLOCATE(data)

    CALL cfd_close

    CALL MPI_BARRIER(comm, errcode)

  END SUBROUTINE restart_data

END MODULE setup
