MODULE setup

  USE shared_data
  USE sdf_job_info
  USE version_data
  USE welcome
  USE diagnostics

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

    CALL get_job_id(jobid)
    run_date = get_unix_time()

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

    CALL get_job_id(jobid)

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
    CHARACTER(LEN=3) :: magic
    REAL(num) :: time0, time1, dt_en
    INTEGER :: ios, num_sz_in, en_nvars_in, p1, p2, nrecs, nrec, recsz
    INTEGER :: version, revision, endianness, header_length, ierr
    LOGICAL :: exists

#ifdef NO_IO
    RETURN
#endif

    IF (rank /= 0) RETURN

    ! Setup lare3d.dat file

    WRITE(file2, '(a, ''/lare3d.dat'')') TRIM(data_dir)

    INQUIRE(FILE=file2, EXIST=exists)

    IF (.NOT.exists .OR. .NOT.restart) THEN
      OPEN(UNIT=stat_unit, STATUS='REPLACE', FILE=file2, &
          FORM='FORMATTED', iostat=ios)
    ELSE
      OPEN(UNIT=stat_unit, STATUS='OLD', POSITION='APPEND', FILE=file2, &
          FORM='FORMATTED', iostat=ios)
    END IF

    IF (ios /= 0) THEN
      PRINT*, 'Unable to open file lare3d.dat for writing. This is ', &
              'most commonly caused by the output directory not existing'
      PRINT*
      PRINT*
      CALL MPI_ABORT(comm, errcode, ierr)
    END IF

    ! Setup en.dat file

    WRITE(file3, '(a, ''/en.dat'')') TRIM(data_dir)

    INQUIRE(FILE=file3, EXIST=exists)

    IF (.NOT.exists .OR. .NOT.restart) THEN
      OPEN(UNIT=en_unit, STATUS='REPLACE', FILE=file3, &
          FORM='UNFORMATTED', ACCESS='STREAM', iostat=ios)
    ELSE
      OPEN(UNIT=en_unit, STATUS='OLD', FILE=file3, &
          FORM='UNFORMATTED', ACCESS='STREAM', iostat=ios)

      READ(en_unit) magic, version, revision
      READ(en_unit) endianness
      READ(en_unit) header_length
      READ(en_unit) num_sz_in, en_nvars_in

      IF (magic /= c_history_magic .OR. endianness /= c_endianness &
          .OR. num_sz_in /= num_sz .OR. en_nvars_in /= en_nvars) THEN
        PRINT*, 'WARNING: incompatible en.dat file found. ', &
                'File will be overwritten.'
        REWIND(en_unit)
      ELSE
        INQUIRE(en_unit, SIZE=p2)
        p1 = header_length

        recsz = num_sz * en_nvars
        nrecs = (p2 + 1 - p1) / recsz

        READ(en_unit) time0
        READ(en_unit,POS=p1+(nrecs-1)*recsz) time1

        ! Guess location of current time position and read it
        dt_en = (time1 - time0) / nrecs
        nrec = FLOOR((time - time0) / dt_en)

        READ(en_unit,POS=p1+nrec*recsz) time1

        ! Search forwards until we reach the first en.dat time after
        ! the current simulation time
        DO WHILE (time-time1 > 1.0e-20_num)
          nrec = nrec + 1
          READ(en_unit,POS=p1+nrec*recsz) time1
        END DO

        ! Search backwards until we reach the first en.dat time before
        ! the current simulation time
        DO WHILE (time1-time > 1.0e-20_num)
          nrec = nrec - 1
          READ(en_unit,POS=p1+nrec*recsz) time1
        END DO

        ! Set file position
        READ(en_unit,POS=p1+nrec*recsz-num_sz) time1
      END IF
    END IF

    IF (ios /= 0) THEN
      PRINT*, 'Unable to open file en.dat for writing. This is ', &
              'most commonly caused by the output directory not existing'
      PRINT*
      PRINT*
      CALL MPI_ABORT(comm, errcode, ierr)
    END IF

    CALL setup_files

  END SUBROUTINE open_files



  !****************************************************************************
  ! Close the output diagnostic files
  !****************************************************************************

  SUBROUTINE close_files

#ifdef NO_IO
    RETURN
#endif

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

    CHARACTER(LEN=c_id_length) :: code_name, block_id, mesh_id, str1
    CHARACTER(LEN=c_max_string_length) :: name
    CHARACTER(LEN=22) :: filename_fmt
    CHARACTER(LEN=5+n_zeros+c_id_length) :: filename
    CHARACTER(LEN=6+data_dir_max_length+n_zeros+c_id_length) :: full_filename
    INTEGER :: blocktype, datatype, code_io_version, string_len
    INTEGER :: ierr, iblock, nblocks, ndims, geometry
    INTEGER, DIMENSION(4) :: dims
    INTEGER, DIMENSION(c_ndims) :: global_dims
    REAL(num), DIMENSION(2*c_ndims) :: extents
    LOGICAL :: restart_flag
    TYPE(sdf_file_handle) :: sdf_handle

    step = -1
    file_number = restart_snapshot

    ! Set the filename. Allows a maximum of 10^999 output dumps.
    WRITE(filename_fmt, '(''(a, i'', i3.3, ''.'', i3.3, '', ".sdf")'')') &
        n_zeros, n_zeros
    WRITE(filename, filename_fmt) TRIM(file_prefix), file_number
    full_filename = TRIM(filesystem) // TRIM(data_dir) // '/' // TRIM(filename)

    IF (rank == 0) THEN
      PRINT*,'Attempting to restart from file: ',TRIM(full_filename)
    END IF

    CALL sdf_open(sdf_handle, full_filename, comm, c_sdf_read)

    CALL sdf_read_header(sdf_handle, step, time, code_name, code_io_version, &
        string_len, restart_flag)

    IF (.NOT. restart_flag) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF file is not a restart dump. Unable to continue.'
      END IF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
      STOP
    END IF

    IF (.NOT.str_cmp(code_name, TRIM(c_code_name))) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF restart file was not generated by ', &
            TRIM(c_code_name) // '.', 'Unable to ', 'continue.'
      END IF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
      STOP
    END IF

    IF (string_len > c_max_string_length) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF file string lengths are too large to read.'
        PRINT*, 'Please increase the size of "c_max_string_length" in ', &
            'shared_data.F90 to ','be at least ', string_len
      END IF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
      STOP
    END IF

    nblocks = sdf_read_nblocks(sdf_handle)
    jobid = sdf_read_jobid(sdf_handle)

    IF (rank == 0) THEN
      PRINT*, 'Loading snapshot for time', time
      CALL create_ascii_header
    END IF

    IF (rank == 0) PRINT*, 'Input file contains', nblocks, 'blocks'

    CALL sdf_read_blocklist(sdf_handle)

    CALL sdf_seek_start(sdf_handle)

    global_dims = (/ nx_global+1, ny_global+1, nz_global+1 /)

    DO iblock = 1, nblocks
      CALL sdf_read_next_block_header(sdf_handle, block_id, name, blocktype, &
          ndims, datatype)

      SELECT CASE(blocktype)
      CASE(c_blocktype_constant)
        IF (str_cmp(block_id, 'dt')) THEN
          CALL sdf_read_srl(sdf_handle, dt_from_restart)
        ELSE IF (str_cmp(block_id, 'time_prev')) THEN
          CALL sdf_read_srl(sdf_handle, time_prev)
        ELSE IF (str_cmp(block_id, 'visc_heating')) THEN
          CALL sdf_read_srl(sdf_handle, total_visc_heating)
          IF (rank /= 0) total_visc_heating = 0
        END IF
      CASE(c_blocktype_plain_mesh)
        IF (ndims /= c_ndims .OR. datatype /= sdf_num &
            .OR. .NOT.str_cmp(block_id, 'grid')) CYCLE

        CALL sdf_read_plain_mesh_info(sdf_handle, geometry, dims, extents)

        IF (geometry /= c_geometry_cartesian &
            .OR. ALL(dims(1:c_ndims) /= global_dims(1:c_ndims))) CYCLE

        ! Should read the grid from file at this point?
        x_min = extents(1)
        x_max = extents(c_ndims+1)
        y_min = extents(2)
        y_max = extents(c_ndims+2)
        z_min = extents(3)
        z_max = extents(c_ndims+3)

      CASE(c_blocktype_plain_variable)
        IF (ndims /= c_ndims .OR. datatype /= sdf_num) CYCLE

        CALL sdf_read_plain_variable_info(sdf_handle, dims, str1, mesh_id)

        IF (.NOT.str_cmp(mesh_id, 'grid')) CYCLE

        IF (str_cmp(block_id, 'Rho')) THEN
          CALL check_dims(dims)
          CALL sdf_read_plain_variable(sdf_handle, rho, &
              cell_distribution, cell_subarray)

        ELSE IF (str_cmp(block_id, 'Energy')) THEN
          CALL check_dims(dims)
          CALL sdf_read_plain_variable(sdf_handle, energy, &
              cell_distribution, cell_subarray)

        ELSE IF (str_cmp(block_id, 'Vx')) THEN
          dims = dims - 1
          CALL check_dims(dims)
          CALL sdf_read_plain_variable(sdf_handle, vx, &
              node_distribution, node_subarray)

        ELSE IF (str_cmp(block_id, 'Vy')) THEN
          dims = dims - 1
          CALL check_dims(dims)
          CALL sdf_read_plain_variable(sdf_handle, vy, &
              node_distribution, node_subarray)

        ELSE IF (str_cmp(block_id, 'Vz')) THEN
          dims = dims - 1
          CALL check_dims(dims)
          CALL sdf_read_plain_variable(sdf_handle, vz, &
              node_distribution, node_subarray)

        ELSE IF (str_cmp(block_id, 'Bx')) THEN
          dims(1) = dims(1) - 1
          CALL check_dims(dims)
          CALL sdf_read_plain_variable(sdf_handle, bx, &
              bx_distribution, bx_subarray)

        ELSE IF (str_cmp(block_id, 'By')) THEN
          IF (c_ndims >= 2) dims(2) = dims(2) - 1
          CALL check_dims(dims)
          CALL sdf_read_plain_variable(sdf_handle, by, &
              by_distribution, by_subarray)

        ELSE IF (str_cmp(block_id, 'Bz')) THEN
          IF (c_ndims >= 3) dims(3) = dims(3) - 1
          CALL check_dims(dims)
          CALL sdf_read_plain_variable(sdf_handle, bz, &
              bz_distribution, bz_subarray)

        END IF

      END SELECT
    END DO

    CALL sdf_close(sdf_handle)

  END SUBROUTINE restart_data



  SUBROUTINE check_dims(dims)

    INTEGER, INTENT(IN) :: dims(:)
    INTEGER, DIMENSION(c_ndims) :: global_dims
    INTEGER :: ierr

    global_dims = (/ nx_global, ny_global, nz_global /)

    IF (ALL(dims(1:c_ndims) == global_dims(1:c_ndims))) RETURN

    IF (rank == 0) THEN
      PRINT*, '*** ERROR ***'
      PRINT*, 'Number of gridpoints in restart dump does not match', &
          ' the control parameters.'
      PRINT*, 'Control grid: ', nx_global, ',', ny_global, ',', nz_global
      PRINT*, 'Restart dump grid: ', dims(1), ',', dims(2), ',', dims(3)
    END IF

    CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
    STOP

  END SUBROUTINE check_dims

END MODULE setup
