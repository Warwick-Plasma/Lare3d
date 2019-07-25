!******************************************************************************
! Controls all I/O and diagnostics. Output files are 'lare3d.dat',
! 'control.dat', 'en.dat' and a series of snapshots in 'nnnn.sdf'
!******************************************************************************

MODULE diagnostics

  USE shared_data
  USE boundary
  USE conduct
  USE sdf
  USE version_data

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: output_routines, energy_correction, write_file, setup_files

  REAL(dbl) :: visc_heating
  LOGICAL, SAVE :: visc_heating_updated = .FALSE.

CONTAINS

  !****************************************************************************
  ! Call the output routines
  !****************************************************************************

  SUBROUTINE output_routines(step)

    INTEGER, INTENT(IN) :: step
    INTEGER, PARAMETER :: history_frequency = 1
    INTEGER, PARAMETER :: dump_frequency = 100
    INTEGER, SAVE :: ndump = 1
    INTEGER :: i
    LOGICAL :: print_arrays, last_call
    REAL(num) :: en_ke = 0.0_num, en_int = 0.0_num, en_b = 0.0_num
    REAL(num), ALLOCATABLE, SAVE :: t_out(:)
    REAL(dbl), ALLOCATABLE, SAVE :: var_local(:,:), var_sum(:,:)
    LOGICAL, SAVE :: first = .TRUE.

#ifdef NO_IO
    RETURN
#endif

    visc_heating_updated = .FALSE.

    ! Check if snapshot is needed
    CALL io_test(step, print_arrays, last_call)

    ! Do every (history_frequency) steps
    IF (MOD(step, history_frequency) == 0 .OR. last_call) THEN
      IF (first) THEN
        ALLOCATE(var_local(en_nvars-1,dump_frequency))
        IF (rank == 0) THEN
          ALLOCATE(var_sum(en_nvars-1,dump_frequency))
          ALLOCATE(t_out(dump_frequency))
        END IF
        first = .FALSE.
      END IF

      CALL energy_account(en_b, en_ke, en_int, .FALSE.)

      IF (rank == 0) t_out(ndump) = time
      var_local(1,ndump) = en_b
      var_local(2,ndump) = en_ke
      var_local(3,ndump) = en_int
      var_local(4,ndump) = total_visc_heating
      var_local(5,ndump) = total_ohmic_heating

      IF (ndump == dump_frequency .OR. last_call) THEN
        CALL MPI_REDUCE(var_local, var_sum, (en_nvars-1) * ndump, &
            MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, errcode)

        visc_heating_updated = .TRUE.

        IF (rank == 0) THEN
          visc_heating = var_sum(4,ndump)
          DO i = 1, ndump
            WRITE(en_unit) t_out(i), REAL(var_sum(:,i), num)
          END DO
        END IF

        ndump = 0
      END IF

      ndump = ndump + 1
    END IF

    IF (print_arrays) CALL write_file(step)

    ! Output energy diagnostics etc
    IF (last_call .AND. rank == 0) THEN
      WRITE(stat_unit,*) 'final nsteps / time = ', step, time

      IF (ALLOCATED(var_local)) DEALLOCATE(var_local)
      IF (ALLOCATED(var_sum)) DEALLOCATE(var_sum)
      IF (ALLOCATED(t_out)) DEALLOCATE(t_out)
    END IF

  END SUBROUTINE output_routines



  !****************************************************************************
  ! Write SDF file
  !****************************************************************************

  SUBROUTINE write_file(step)

    INTEGER, INTENT(IN) :: step
    REAL(num), DIMENSION(:), ALLOCATABLE :: work
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: array
    LOGICAL :: restart_flag, convert
    INTEGER, DIMENSION(c_ndims) :: global_dims, dims
    CHARACTER(LEN=22) :: filename_fmt
    CHARACTER(LEN=5+n_zeros+c_id_length) :: filename
    CHARACTER(LEN=6+data_dir_max_length+n_zeros+c_id_length) :: full_filename
    CHARACTER(LEN=c_id_length) :: varname, units
    TYPE(sdf_file_handle) :: sdf_handle
    LOGICAL, SAVE :: first = .TRUE.

    global_dims = (/ nx_global, ny_global, nz_global /)

    IF (first) THEN
      ! Resize the {x,y,z}b_global to be the correct size for output
      ALLOCATE(work(-2:MAX(nx_global,ny_global,nz_global)+2))

      work(-2:nx_global+2) = xb_global
      DEALLOCATE(xb_global)
      ALLOCATE(xb_global(0:nx_global))
      xb_global(:) = work(0:nx_global)

      work(-2:ny_global+2) = yb_global
      DEALLOCATE(yb_global)
      ALLOCATE(yb_global(0:ny_global))
      yb_global(:) = work(0:ny_global)

      work(-2:nz_global+2) = zb_global
      DEALLOCATE(zb_global)
      ALLOCATE(zb_global(0:nz_global))
      zb_global(:) = work(0:nz_global)

      DEALLOCATE(work)
      first = .FALSE.
    END IF

    ! Output a snapshot of arrays
    IF (rank == 0) THEN
      WRITE(stat_unit,*) 'Dumping ', file_number, ' at time', time
      IF (conduction) WRITE(stat_unit,*) 'Number of super-steps = ', n_s_stages
      CALL FLUSH(stat_unit)
    END IF

    ! Set the filename. Allows a maximum of 10^999 output dumps.
    WRITE(filename_fmt, '(''(a, i'', i3.3, ''.'', i3.3, '', ".sdf")'')') &
        n_zeros, n_zeros
    WRITE(filename, filename_fmt) TRIM(file_prefix), file_number
    full_filename = TRIM(filesystem) // TRIM(data_dir) // '/' // TRIM(filename)

    ! If dump_mask(1:8) are true then this file can be used for restarting
    restart_flag = ALL(dump_mask(1:8))

    convert = .FALSE.

    IF (.NOT.visc_heating_updated) THEN
      CALL MPI_REDUCE(total_visc_heating, visc_heating, 1, &
          MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, errcode)
      visc_heating_updated = .TRUE.
    END IF

    CALL sdf_open(sdf_handle, full_filename, comm, c_sdf_write)
    CALL sdf_set_string_length(sdf_handle, c_max_string_length)
    CALL sdf_write_header(sdf_handle, TRIM(c_code_name), 1, step, time, &
        restart_flag, jobid)
    CALL sdf_write_run_info(sdf_handle, c_version, c_revision, c_minor_rev, &
        c_commit_id, '', c_compile_machine, c_compile_flags, 0_8, &
        c_compile_date, run_date)
    CALL sdf_write_cpu_split(sdf_handle, 'cpu_rank', 'CPUs/Original rank', &
        cell_nx_maxs, cell_ny_maxs, cell_nz_maxs)
    CALL sdf_write_srl(sdf_handle, 'dt', 'Time increment', dt)
    CALL sdf_write_srl(sdf_handle, 'time_prev', 'Last dump time requested', &
        time_prev)
    CALL sdf_write_srl(sdf_handle, 'visc_heating', 'Viscous heating total', &
        visc_heating)

    CALL sdf_write_srl_plain_mesh(sdf_handle, 'grid', 'Grid/Grid', &
        xb_global, yb_global, zb_global, convert)

    IF (dump_mask(1)) THEN
      varname = 'Rho'
      units = 'kg/m^3'
      dims = global_dims

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Fluid/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_cell_centre, 'grid', rho, &
          cell_distribution, cell_subarray, convert)
    END IF

    IF (dump_mask(2)) THEN
      varname = 'Energy'
      units = 'J'
      dims = global_dims

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Fluid/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_cell_centre, 'grid', energy, &
          cell_distribution, cell_subarray, convert)
    END IF

    IF (dump_mask(3)) THEN
      varname = 'Vx'
      units = 'm/s'
      dims = global_dims + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Velocity/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_vertex, 'grid', vx, &
          node_distribution, node_subarray, convert)
    END IF

    IF (dump_mask(4)) THEN
      varname = 'Vy'
      units = 'm/s'
      dims = global_dims + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Velocity/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_vertex, 'grid', vy, &
          node_distribution, node_subarray, convert)
    END IF

    IF (dump_mask(5)) THEN
      varname = 'Vz'
      units = 'm/s'
      dims = global_dims + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Velocity/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_vertex, 'grid', vz, &
          node_distribution, node_subarray, convert)
    END IF

    IF (dump_mask(6)) THEN
      varname = 'Bx'
      units = 'T'
      dims = global_dims
      dims(1) = dims(1) + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Magnetic_Field/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_bx, 'grid', bx, &
          bx_distribution, bx_subarray, convert)
    END IF

    IF (dump_mask(7)) THEN
      varname = 'By'
      units = 'T'
      dims = global_dims
      dims(2) = dims(2) + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Magnetic_Field/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_by, 'grid', by, &
          by_distribution, by_subarray, convert)
    END IF

    IF (dump_mask(8)) THEN
      varname = 'Bz'
      units = 'T'
      dims = global_dims
      dims(3) = dims(3) + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Magnetic_Field/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_bz, 'grid', bz, &
          bz_distribution, bz_subarray, convert)
    END IF

    IF (dump_mask(9)) THEN
      varname = 'Temperature'
      units = 'K'
      dims = global_dims

      IF (.NOT.ALLOCATED(array)) ALLOCATE(array(nx,ny,nz))

      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            array(ix,iy,iz) = (gamma - 1.0_num) / (2.0_num - xi_n(ix,iy,iz)) &
                * (energy(ix,iy,iz) - (1.0_num - xi_n(ix,iy,iz)) * ionise_pot)
          END DO
        END DO
      END DO

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Fluid/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_cell_centre, 'grid', array, &
          cell_distribution, cellng_subarray, convert)
    END IF

    IF (dump_mask(10)) THEN
      varname = 'Pressure'
      units = 'Pa'
      dims = global_dims

      IF (.NOT.ALLOCATED(array)) ALLOCATE(array(nx,ny,nz))

      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            array(ix,iy,iz) = (gamma - 1.0_num) * rho(ix,iy,iz) &
                * (energy(ix,iy,iz) - (1.0_num - xi_n(ix,iy,iz)) * ionise_pot)
          END DO
        END DO
      END DO

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Fluid/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_cell_centre, 'grid', array, &
          cell_distribution, cellng_subarray, convert)
    END IF

    IF (dump_mask(11)) THEN
      varname = 'Cs'
      units = 'm/s'
      dims = global_dims

      IF (.NOT.ALLOCATED(array)) ALLOCATE(array(nx,ny,nz))

      array(:,:,:) = SQRT(gamma * (gamma - 1.0_num) * energy(1:nx,1:ny,1:nz))

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Fluid/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_cell_centre, 'grid', array, &
          cell_distribution, cellng_subarray, convert)
    END IF

    IF (dump_mask(12)) THEN
      varname = 'j_par'
      units = 'A/m^2'
      dims = global_dims + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'PIP/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_vertex, 'grid', parallel_current, &
          node_distribution, nodeng_subarray, convert)
    END IF

    IF (dump_mask(13)) THEN
      varname = 'j_perp'
      units = 'A/m^2'
      dims = global_dims + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'PIP/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_vertex, 'grid', perp_current, &
          node_distribution, nodeng_subarray, convert)
    END IF

    IF (dump_mask(14)) THEN
      varname = 'neutral_fraction'
      units = '%'
      dims = global_dims

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'PIP/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_cell_centre, 'grid', xi_n, &
          cell_distribution, cell_subarray, convert)
    END IF

    IF (dump_mask(15)) THEN
      varname = 'eta_perp'
      units = ''
      dims = global_dims

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'PIP/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_cell_centre, 'grid', eta_perp, &
          cell_distribution, cell_subarray, convert)
    END IF

    IF (dump_mask(16)) THEN
      varname = 'eta'
      units = ''
      dims = global_dims

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'PIP/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_cell_centre, 'grid', eta, &
          cell_distribution, cell_subarray, convert)
    END IF

    IF (dump_mask(17)) THEN
      varname = 'Jx'
      units = 'A/m^2'
      dims = global_dims + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Current/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_vertex, 'grid', jx_r, &
          node_distribution, nodeng_subarray, convert)
    END IF

    IF (dump_mask(18)) THEN
      varname = 'Jy'
      units = 'A/m^2'
      dims = global_dims + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Current/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_vertex, 'grid', jy_r, &
          node_distribution, nodeng_subarray, convert)
    END IF

    IF (dump_mask(19)) THEN
      varname = 'Jz'
      units = 'A/m^2'
      dims = global_dims + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Current/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_vertex, 'grid', jz_r, &
          node_distribution, nodeng_subarray, convert)
    END IF

    IF (dump_mask(20)) THEN
      varname = 'Resistive_Heat'
      units = 'J/m^3'
      dims = global_dims

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Fluid/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_cell_centre, 'grid', ohmic_dep, &
          cell_distribution, cell_subarray, convert)

      varname = 'Viscous_Heat'
      units = 'J/m^3'
      dims = global_dims

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Fluid/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_cell_centre, 'grid', visc_dep, &
          cell_distribution, cell_subarray, convert)

      ohmic_dep = 0.0_num
      visc_dep = 0.0_num
    END IF

    IF (ALLOCATED(array)) DEALLOCATE(array)

    ! Close the file
    CALL sdf_close(sdf_handle)

    file_number = file_number + 1

  END SUBROUTINE write_file



  !****************************************************************************
  ! Test whether any of the conditions for doing output on the current
  ! iteration are met
  !****************************************************************************

  SUBROUTINE io_test(step, print_arrays, last_call)

    INTEGER, INTENT(IN) :: step
    LOGICAL, INTENT(OUT) :: print_arrays, last_call

    REAL(num), SAVE :: t1 = 0.0_num

    IF (restart) THEN
      t1 = time_prev + dt_snapshots
    END IF

    print_arrays = .FALSE.
    last_call = .FALSE.

    IF (time >= t1) THEN
      print_arrays = .TRUE.
      time_prev = t1
      t1 = t1 + dt_snapshots
    END IF

    IF (time >= t_end .OR. step == nsteps) THEN
      last_call = .TRUE.
      print_arrays = .TRUE.
    END IF

    IF (restart) THEN
      print_arrays = .FALSE.
      file_number = file_number + 1
      restart = .FALSE.
    END IF

  END SUBROUTINE io_test




  SUBROUTINE energy_account(energy_b, energy_ke, energy_int, do_sum)

    REAL(dbl), INTENT(OUT) :: energy_b, energy_ke, energy_int
    LOGICAL, INTENT(IN) :: do_sum
    REAL(dbl) :: energy_b_local, energy_ke_local, energy_int_local
    REAL(dbl) :: energy_local(3), energy_sum(3)
    REAL(dbl) :: cv_v, rho_v, w1, w2, w3

    energy_b_local   = 0.0_dbl
    energy_ke_local  = 0.0_dbl
    energy_int_local = 0.0_dbl

    DO iz = 1, nz
      izm = iz - 1
      DO iy = 1, ny
        iym = iy - 1
        DO ix = 1, nx
          ixm = ix - 1

          w1 = (bx(ix,iy,iz)**2 + bx(ixm,iy ,iz )**2) * 0.5_num
          w2 = (by(ix,iy,iz)**2 + by(ix ,iym,iz )**2) * 0.5_num
          w3 = (bz(ix,iy,iz)**2 + bz(ix ,iy ,izm)**2) * 0.5_num
          w1 = (w1 + w2 + w3) * 0.5_dbl
          energy_b_local = energy_b_local + w1 * cv(ix,iy,iz)

          energy_int_local = energy_int_local &
              + energy(ix,iy,iz) * rho(ix,iy,iz) * cv(ix,iy,iz)
        END DO
      END DO
    END DO

    DO iz = 0, nz
      izp = iz + 1
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = 0, nx
          ixp = ix + 1

          ! WARNING the KE is summed on the vertices
          rho_v = rho(ix ,iy ,iz ) * cv(ix ,iy ,iz ) &
                + rho(ixp,iy ,iz ) * cv(ixp,iy ,iz ) &
                + rho(ix ,iyp,iz ) * cv(ix ,iyp,iz ) &
                + rho(ixp,iyp,iz ) * cv(ixp,iyp,iz ) &
                + rho(ix ,iy ,izp) * cv(ix ,iy ,izp) &
                + rho(ixp,iy ,izp) * cv(ixp,iy ,izp) &
                + rho(ix ,iyp,izp) * cv(ix ,iyp,izp) &
                + rho(ixp,iyp,izp) * cv(ixp,iyp,izp)

          cv_v = cv(ix,iy ,iz ) + cv(ixp,iy ,iz ) &
               + cv(ix,iyp,iz ) + cv(ixp,iyp,iz ) &
               + cv(ix,iy ,izp) + cv(ixp,iy ,izp) &
               + cv(ix,iyp,izp) + cv(ixp,iyp,izp)

          rho_v = rho_v / cv_v
          cv_v = cv_v * 0.125_dbl
          w1 = rho_v * cv_v &
              * (vx(ix,iy,iz)**2 + vy(ix,iy,iz)**2 + vz(ix,iy,iz)**2)

          IF (ix == 0 .OR. ix == nx) THEN
            w1 = w1 * 0.5_dbl
          END IF

          IF (iy == 0 .OR. iy == ny) THEN
            w1 = w1 * 0.5_dbl
          END IF

          IF (iz == 0 .OR. iz == nz) THEN
            w1 = w1 * 0.5_dbl
          END IF

          energy_ke_local = energy_ke_local + w1 * 0.5_dbl
        END DO
      END DO
    END DO

    IF (do_sum) THEN
      energy_local(1) = energy_b_local
      energy_local(2) = energy_ke_local
      energy_local(3) = energy_int_local

      CALL MPI_ALLREDUCE(energy_local, energy_sum, 3, MPI_DOUBLE_PRECISION, &
          MPI_SUM, comm, errcode)

      energy_b   = energy_sum(1)
      energy_ke  = energy_sum(2)
      energy_int = energy_sum(3)
    ELSE
      energy_b   = energy_b_local
      energy_ke  = energy_ke_local
      energy_int = energy_int_local
    END IF

  END SUBROUTINE energy_account



  SUBROUTINE energy_correction

    delta_ke = -delta_ke
    WHERE (delta_ke < 0.0_num) delta_ke = 0.0_num
    delta_ke(:,:,:) = delta_ke(:,:,:) / (rho(:,:,:) * cv(:,:,:))

    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          energy(ix,iy,iz) = energy(ix,iy,iz) + delta_ke(ix,iy,iz)
        END DO
      END DO
    END DO

    CALL energy_bcs

  END SUBROUTINE energy_correction



  SUBROUTINE setup_files

    INTEGER :: p1, header_length
    CHARACTER(LEN=c_id_length) :: varnames(en_nvars)

    CALL output_log

    INQUIRE(en_unit, POS=p1)
    IF (p1 /= 1) RETURN

    varnames(1) = 'time'
    varnames(2) = 'en_b'
    varnames(3) = 'en_ke'
    varnames(4) = 'en_int'
    varnames(5) = 'heating_visc'
    varnames(6) = 'heating_ohmic'

    header_length = 3 + 7 * 4 + en_nvars * c_id_length
    ! Write history file header if not appending to file
    WRITE(en_unit) c_history_magic, c_history_version, c_history_revision
    WRITE(en_unit) c_endianness
    WRITE(en_unit) header_length
    WRITE(en_unit) num_sz, en_nvars
    WRITE(en_unit) c_id_length
    WRITE(en_unit) varnames

  END SUBROUTINE setup_files



  SUBROUTINE output_log

    ! Writes basic data to 'lare3d.dat'

    IF (restart) THEN
      WRITE(stat_unit,*)
      WRITE(stat_unit,*) '#####################################################'
      WRITE(stat_unit,*)
      WRITE(stat_unit,*) 'Restarting from step ', step, ' and time ', time
      WRITE(stat_unit,*)
    END IF

    WRITE(stat_unit,*) ascii_header
    WRITE(stat_unit,*)
    WRITE(stat_unit,*) 'nprocx, nprocy, nprocz = ', nprocx, nprocy, nprocz
    WRITE(stat_unit,*) 'nx, ny, nz = ', nx_global, ny_global, nz_global
    WRITE(stat_unit,*)
    WRITE(stat_unit,*) 'length_x = ', length_x
    WRITE(stat_unit,*) 'length_y = ', length_y
    WRITE(stat_unit,*) 'length_z = ', length_z
    WRITE(stat_unit,*)
#ifdef QMONO
    WRITE(stat_unit,*) 'q_mono viscosity (-DQMONO)'
#else
    WRITE(stat_unit,*) 'tensor shock viscosity'
#endif
#ifdef FOURTHORDER
    WRITE(stat_unit,*) '4th-order resistive update (-DFOURTHORDER)'
#endif
#ifdef SINGLE
    WRITE(stat_unit,*) 'single precision (-DSINGLE)'
#endif
    WRITE(stat_unit,*) 'linear viscosity coeff = ', visc1
    WRITE(stat_unit,*) 'quadratic viscosity coeff = ', visc2
    WRITE(stat_unit,*) 'j_max = ', j_max
    WRITE(stat_unit,*) 'eta0 = ', eta0
    WRITE(stat_unit,*) 'eta_background = ', eta_background
    WRITE(stat_unit,*) 'kappa = ', kappa_0
    WRITE(stat_unit,*)
    WRITE(stat_unit,*) 'mass_fraction = ', mf
    WRITE(stat_unit,*) 'normalising B = ', B_norm
    WRITE(stat_unit,*) 'normalising L = ', L_norm
    WRITE(stat_unit,*) 'normalising density = ', rho_norm
    WRITE(stat_unit,*) 'normalising speed = ', B_norm / SQRT(mu0_SI * rho_norm)
    WRITE(stat_unit,*) 'normalising time = ', L_norm / (B_norm / SQRT(mu0_SI * rho_norm))
    WRITE(stat_unit,*) 'normalising temperature = ', temp_norm
    WRITE(stat_unit,*)
    WRITE(stat_unit,*) 't_start, t_end = ', time, t_end
    WRITE(stat_unit,*) 'nsteps =', nsteps
    WRITE(stat_unit,*)

  END SUBROUTINE output_log

END MODULE diagnostics
