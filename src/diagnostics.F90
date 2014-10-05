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

  PUBLIC :: set_dt, output_routines, energy_correction

CONTAINS

  !****************************************************************************
  ! Call the output routines
  !****************************************************************************

  SUBROUTINE output_routines(i) ! i = step index

    INTEGER, INTENT(IN) :: i
    INTEGER, PARAMETER :: outstep = 1
    REAL(num), DIMENSION(:), ALLOCATABLE :: work
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: array
    LOGICAL :: print_arrays, last_call, restart_flag, convert
    INTEGER, DIMENSION(c_ndims) :: global_dims, dims

    CHARACTER(LEN=22) :: filename_fmt
    CHARACTER(LEN=5+n_zeros+c_id_length) :: filename
    CHARACTER(LEN=6+data_dir_max_length+n_zeros+c_id_length) :: full_filename
    CHARACTER(LEN=c_id_length) :: varname, units
    TYPE(sdf_file_handle) :: sdf_handle

    REAL(num) :: t_out = 0.0_num
    REAL(num) :: en_ke = 0.0_num, en_int = 0.0_num
    REAL(num) :: en_b = 0.0_num, heating_visc = 0.0_num
    REAL(num) :: heating_ohmic = 0.0_num
    REAL(num) :: total
    LOGICAL, SAVE :: first = .TRUE.

    global_dims = (/ nx_global, ny_global, nz_global /)

#ifdef NO_IO
    RETURN
#endif

    ! Done just once at the start
    IF (i == 0 .AND. rank == 0) THEN
      CALL output_log
      IF (.NOT. restart) WRITE(en_unit) num, 6
    END IF

    ! Do every (outstep) steps
    IF (MOD(i, outstep) == 0 .OR. last_call) THEN
      t_out = time
      CALL energy_account(en_b, en_ke, en_int)

      CALL MPI_ALLREDUCE(total_visc_heating, total, 1, mpireal, MPI_SUM, &
          comm, errcode)

      heating_visc = total

      CALL MPI_ALLREDUCE(total_ohmic_heating, total, 1, mpireal, MPI_SUM, &
          comm, errcode)

      heating_ohmic = total

      IF (rank == 0) THEN
        WRITE(en_unit) t_out, en_b, en_ke, en_int
        WRITE(en_unit) heating_visc, heating_ohmic
      END IF
    END IF

    ! Check if snapshot is needed
    CALL io_test(i, print_arrays, last_call)

    ! Output energy diagnostics etc
    IF (last_call .AND. rank == 0) THEN
      WRITE(stat_unit,*) 'final nsteps / time = ', i, time
    END IF

    IF (.NOT.print_arrays) RETURN

    IF (first) THEN
      ! Resize the {x,y,z}b_global to be the correct size for output
      ALLOCATE(work(-2:MAX(nx_global,ny_global,nz_global)+2))

      work(-2:nx_global+2) = xb_global
      DEALLOCATE(xb_global)
      ALLOCATE(xb_global(0:nx_global))
      xb_global = work(0:nx_global)

      work(-2:ny_global+2) = yb_global
      DEALLOCATE(yb_global)
      ALLOCATE(yb_global(0:ny_global))
      yb_global = work(0:ny_global)

      work(-2:nz_global+2) = zb_global
      DEALLOCATE(zb_global)
      ALLOCATE(zb_global(0:nz_global))
      zb_global = work(0:nz_global)

      DEALLOCATE(work)
      first = .FALSE.
    END IF

    ! Output a snapshot of arrays
    IF (rank == 0) THEN
      WRITE(stat_unit,*) 'Dumping ', file_number, ' at time', time
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

    CALL sdf_open(sdf_handle, full_filename, comm, c_sdf_write)
    CALL sdf_write_header(sdf_handle, TRIM(c_code_name), 1, i, time, &
        restart_flag, jobid)
    CALL sdf_write_run_info(sdf_handle, c_version, c_revision, c_minor_rev, &
        c_commit_id, '', c_compile_machine, c_compile_flags, 0_8, &
        c_compile_date, run_date)
    CALL sdf_write_cpu_split(sdf_handle, 'cpu_rank', 'CPUs/Original rank', &
        cell_nx_maxs, cell_ny_maxs, cell_nz_maxs)

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

      array = SQRT(gamma * (gamma - 1.0_num) * energy(1:nx,1:ny,1:nz))

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
      dims = global_dims + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'PIP/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_vertex, 'grid', xi_n, &
          node_distribution, nodeng_subarray, convert)
    END IF

    IF (dump_mask(15)) THEN
      varname = 'eta_perp'
      units = ''
      dims = global_dims + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'PIP/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_vertex, 'grid', eta_perp, &
          node_distribution, nodeng_subarray, convert)
    END IF

    IF (dump_mask(16)) THEN
      varname = 'eta'
      units = ''
      dims = global_dims + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'PIP/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_vertex, 'grid', eta, &
          node_distribution, nodeng_subarray, convert)
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

    IF (ALLOCATED(array)) DEALLOCATE(array)

    ! Close the file
    CALL sdf_close(sdf_handle)

    file_number = file_number + 1

  END SUBROUTINE output_routines



  !****************************************************************************
  ! Test whether any of the conditions for doing output on the current
  ! iteration are met
  !****************************************************************************

  SUBROUTINE io_test(i, print_arrays, last_call)

    INTEGER, INTENT(IN) :: i
    LOGICAL, INTENT(OUT) :: print_arrays, last_call

    REAL(num), SAVE :: t1 = 0.0_num

    IF (restart) THEN
      t1 = time
      restart = .FALSE.
    END IF

    print_arrays = .FALSE.
    last_call = .FALSE.

    IF (time >= t1) THEN
      print_arrays = .TRUE.
      t1 = t1 + dt_snapshots
    END IF

    IF (time >= t_end .OR. i == nsteps) THEN
      last_call = .TRUE.
      print_arrays = .TRUE.
    END IF

  END SUBROUTINE io_test



  SUBROUTINE set_dt ! sets CFL limited step

    ! Assumes all variables are defined at the same point. Be careful with
    ! setting 'dt_multiplier' if you expect massive changes across cells.

    REAL(num) :: vxbm, vxbp, avxm, avxp, dvx, ax
    REAL(num) :: vybm, vybp, avym, avyp, dvy, ay
    REAL(num) :: vzbm, vzbp, avzm, avzp, dvz, az
    REAL(num) :: cons, cs, volume
    REAL(num) :: dxlocal, dt_local, dtr_local, dt1, dt2, dt3
    REAL(num) :: dt_locals(2), dt_min(2)

    dt_local = largest_number
    dtr_local = largest_number
    cons = gamma * (gamma - 1.0_num)

    DO iz = 0, nz
      izm = iz - 1
      DO iy = 0, ny
        iym = iy - 1
        DO ix = 0, nx
          ixm = ix - 1

          ! Fix dt for Lagrangian step
          w1 = bx(ix,iy,iz)**2 + by(ix,iy,iz)**2 + bz(ix,iy,iz)**2
          ! Sound speed squared
          cs = cons * energy(ix,iy,iz)

          w2 = SQRT(cs + w1 / MAX(rho(ix,iy,iz), none_zero) &
              + 2.0_num * p_visc(ix,iy,iz) / MAX(rho(ix,iy,iz), none_zero))

          ! Find ideal MHD CFL limit for Lagrangian step
          dt1 = MIN(dxb(ix), dyb(iy), dzb(iz)) / w2
          dt_local = MIN(dt_local, dt1)

          ! Now find dt for remap step
          ax = 0.25_num * dyb(iy) * dzb(iz)
          ay = 0.25_num * dxb(ix) * dzb(iz)
          az = 0.25_num * dxb(ix) * dyb(iy)
          vxbm = (vx(ixm,iy ,iz ) + vx(ixm,iym,iz ) &
              +   vx(ixm,iy ,izm) + vx(ixm,iym,izm)) * ax
          vxbp = (vx(ix ,iy ,iz ) + vx(ix ,iym,iz ) &
              +   vx(ix ,iy ,izm) + vx(ix ,iym,izm)) * ax
          vybm = (vy(ix ,iym,iz ) + vy(ixm,iym,iz ) &
              +   vy(ix ,iym,izm) + vy(ixm,iym,izm)) * ay
          vybp = (vy(ix ,iy ,iz ) + vy(ixm,iy ,iz ) &
              +   vy(ix ,iy ,izm) + vy(ixm,iy ,izm)) * ay
          vzbm = (vz(ix ,iy ,izm) + vz(ixm,iy ,izm) &
              +   vz(ix ,iym,izm) + vz(ixm,iym,izm)) * az
          vzbp = (vz(ix ,iy ,iz ) + vz(ixm,iy ,iz ) &
              +   vz(ix ,iym,iz ) + vz(ixm,iym,iz )) * az

          dvx = ABS(vxbp - vxbm)
          dvy = ABS(vybp - vybm)
          dvz = ABS(vzbp - vzbm)
          avxm = ABS(vxbm)
          avxp = ABS(vxbp)
          avym = ABS(vybm)
          avyp = ABS(vybp)
          avzm = ABS(vzbm)
          avzp = ABS(vzbp)

          volume = ax * dxb(ix)
          dt1 = volume / MAX(avxm, avxp, dvx, 1.e-10_num * volume)
          dt2 = volume / MAX(avym, avyp, dvy, 1.e-10_num * volume)
          dt3 = volume / MAX(avzm, avzp, dvz, 1.e-10_num * volume)

          ! Fix dt for remap step
          dt_local = MIN(dt_local, dt1, dt2, dt3)

          ! Note resistive limits assumes uniform resistivity hence cautious
          ! factor 0.2
          dxlocal = 1.0_num / (1.0_num / dxb(ix)**2 &
              + 1.0_num / dyb(iy)**2 + 1.0_num / dzb(iz)**2)

          IF (cowling_resistivity) THEN
            dt1 = 0.2_num * dxlocal &
                / MAX(MAX(eta(ix,iy,iz), eta_perp(ix,iy,iz)), none_zero)
          ELSE
            dt1 = 0.2_num * dxlocal / MAX(eta(ix,iy,iz), none_zero)
          END IF

          ! Adjust to accomodate resistive effects
          dtr_local = MIN(dtr_local, dt1)
        END DO
      END DO
    END DO

    dt_locals(1) = dt_local
    dt_locals(2) = dtr_local

    CALL MPI_ALLREDUCE(dt_locals, dt_min, 2, mpireal, MPI_MIN, comm, errcode)

    dt  = dt_multiplier * dt_min(1)
    dtr = dt_multiplier * dt_min(2)

    time = time + dt

  END SUBROUTINE set_dt



  SUBROUTINE energy_account(energy_b, energy_ke, energy_int)

    REAL(num), INTENT(OUT) :: energy_b, energy_ke, energy_int
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

    energy_local(1) = energy_b_local
    energy_local(2) = energy_ke_local
    energy_local(3) = energy_int_local

    CALL MPI_ALLREDUCE(energy_local, energy_sum, 3, MPI_DOUBLE_PRECISION, &
        MPI_SUM, comm, errcode)

    energy_b   = REAL(energy_sum(1), num)
    energy_ke  = REAL(energy_sum(2), num)
    energy_int = REAL(energy_sum(3), num)

  END SUBROUTINE energy_account



  SUBROUTINE energy_correction

    delta_ke = -delta_ke
    WHERE (delta_ke < 0.0_num) delta_ke = 0.0_num
    delta_ke = delta_ke / (rho * cv)

    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          energy(ix,iy,iz) = energy(ix,iy,iz) + delta_ke(ix,iy,iz)
        END DO
      END DO
    END DO

    CALL energy_bcs

  END SUBROUTINE energy_correction



  SUBROUTINE output_log

    ! Writes basic data to 'lare3d.dat'

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
    WRITE(stat_unit,*) 'uniform tensor viscosity coeff = ', visc3
    WRITE(stat_unit,*) 'j_max = ', j_max
    WRITE(stat_unit,*) 'eta0 = ', eta0
    WRITE(stat_unit,*) 'eta_background = ', eta_background
    WRITE(stat_unit,*) 'kappa = ', kappa_0
    WRITE(stat_unit,*)
    WRITE(stat_unit,*) 't_start, t_end = ', time, t_end
    WRITE(stat_unit,*) 'nsteps =', nsteps
    WRITE(stat_unit,*)

  END SUBROUTINE output_log

END MODULE diagnostics
