!******************************************************************************
! Controls all I/O and diagnostics. Output files are 'lare3d.dat',
! 'control.dat', 'en.dat' and a series of snapshots in 'nnnn.sdf'
!******************************************************************************

MODULE diagnostics

  USE shared_data
  USE boundary
  USE conduct
  USE output_cartesian
  USE output
  USE iocontrol

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_dt, output_routines, energy_correction

CONTAINS

  !****************************************************************************
  ! Call the output routines
  !****************************************************************************

  SUBROUTINE output_routines(i) ! i = step index

    INTEGER, INTENT(IN) :: i

    INTEGER, PARAMETER :: out = 1000
    INTEGER, SAVE :: index = 1, step = 1
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: array
    LOGICAL :: print_arrays, last_call
    REAL(num), DIMENSION(c_ndims) :: stagger = 0.0_num
    INTEGER, DIMENSION(c_ndims) :: dims

    ! This output routine uses the same structure as needed for MPI output.
    ! This is more complicated than need for the serial code
    ! rank always equals zero in this serial code
    CHARACTER(LEN=9+data_dir_max_length+n_zeros) :: filename
    CHARACTER(LEN=35) :: filename_desc

    REAL(num) :: t_out = 0.0_num
    REAL(num) :: en_ke = 0.0_num, en_int = 0.0_num
    REAL(num) :: en_b = 0.0_num, heating_visc = 0.0_num
    REAL(num) :: heating_ohmic = 0.0_num
    REAL(num) :: total

    dims = (/nx_global+1, ny_global+1, nz_global+1/)

    ! Make sure output fits arrays
    IF (nsteps >= out) step = nsteps / out + 1

    ! Done just once at the start
    IF (i == 0 .AND. rank == 0) THEN
      CALL output_log
      IF (.NOT. restart) WRITE(en_unit) num, 6
    END IF

    ! Do every (step) steps
    IF (MOD(i, step) == 0 .OR. last_call) THEN
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

      index = index + 1
    END IF

    ! Check if snapshot is needed
    CALL io_test(i, print_arrays, last_call)

    ! Output a snapshot of arrays
    IF (print_arrays) THEN
      IF (rank == 0) THEN
        WRITE(stat_unit,*) 'Dumping ', file_number, ' at time', time
        CALL FLUSH(stat_unit)
      END IF

      ! Set the filename
      WRITE(filename_desc, '("(''nfs:'', a, ''/'', i", i3.3, ".", i3.3, &
          & ", ''.cfd'')")'), n_zeros, n_zeros
      WRITE(filename, filename_desc) TRIM(data_dir), file_number

      CALL cfd_open(filename, rank, comm, MPI_MODE_CREATE + MPI_MODE_WRONLY)
      CALL cfd_write_snapshot_data(REAL(time, dbl), i, 0)

      ALLOCATE(array(0:nx,0:ny,0:nz))

      CALL cfd_write_3d_cartesian_grid('Grid', 'Grid', &
          xb_global(0:nx_global), yb_global(0:ny_global), &
          zb_global(0:nz_global), 0)

      IF (dump_mask(1)) THEN
        array = rho(0:nx,0:ny,0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel('Rho', 'Fluid', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(2)) THEN
        array = energy(0:nx,0:ny,0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel('Energy', 'Fluid', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(3)) THEN
        array = vx(0:nx,0:ny,0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel('Vx', 'Velocity', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(4)) THEN
        array = vy(0:nx,0:ny,0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel('Vy', 'Velocity', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(5)) THEN
        array = vz(0:nx,0:ny,0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel('Vz', 'Velocity', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(6)) THEN
        array = bx(0:nx,0:ny,0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel('Bx', 'Magnetic_Field', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(7)) THEN
        array = by(0:nx,0:ny,0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel('By', 'Magnetic_Field', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(8)) THEN
        array = bz(0:nx,0:ny,0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel('Bz', 'Magnetic_Field', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(9)) THEN
        DO iz = 0, nz
          DO iy = 0, ny
            DO ix = 0, nx
              array(ix,iy,iz) = (gamma - 1.0_num) / (2.0_num - xi_n(ix,iy,iz)) &
                  * (energy(ix,iy,iz) - (1.0_num - xi_n(ix,iy,iz)) * ionise_pot)
            END DO
          END DO
        END DO
        CALL cfd_write_3d_cartesian_variable_parallel('Temperature', 'Fluid', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(10)) THEN
        DO iz = 0, nz
          DO iy = 0, ny
            DO ix = 0, nx
              array(ix,iy,iz) = (gamma - 1.0_num) * rho(ix,iy,iz) &
                  * (energy(ix,iy,iz) - (1.0_num - xi_n(ix,iy,iz)) * ionise_pot)
            END DO
          END DO
        END DO
        CALL cfd_write_3d_cartesian_variable_parallel('Pressure', 'Fluid', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(11)) THEN
        array = SQRT(gamma * (gamma-1.0_num) * energy(0:nx,0:ny,0:nz))
        CALL cfd_write_3d_cartesian_variable_parallel('cs', 'Fluid', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(12)) THEN
        array = parallel_current(0:nx,0:ny,0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel('j_par', 'PIP', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(13)) THEN
        array = perp_current(0:nx,0:ny,0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel('j_perp', 'PIP', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(14)) THEN
        array = xi_n(0:nx,0:ny,0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel('neutral_fraction', &
            'PIP', dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(15)) THEN
        array = eta_perp(0:nx,0:ny,0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel('eta_perp', 'PIP', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(16)) THEN
        array = eta(0:nx,0:ny,0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel('eta', 'PIP', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(17)) THEN
        array = jx_r(0:nx,0:ny,0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel('jx', 'current', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(18)) THEN
        array = jy_r(0:nx,0:ny,0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel('jy', 'current', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(19)) THEN
        array = jz_r(0:nx,0:ny,0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel('jz', 'current', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      DEALLOCATE(array)

      ! Close the file
      CALL cfd_close

      file_number = file_number + 1
    END IF

    ! Output energy diagnostics etc
    IF (last_call .AND. rank == 0) THEN
      WRITE(stat_unit,*) 'final nsteps / time = ', i, time
    END IF

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

    WRITE(stat_unit,*) 'nprocx, nprocy, nprocz = ', nprocx, nprocy, nprocz
    WRITE(stat_unit,*) 'nx, ny, nz = ', nx, ny, nz
    WRITE(stat_unit,*)
    WRITE(stat_unit,*) 'length_x = ', length_x
    WRITE(stat_unit,*) 'length_y = ', length_y
    WRITE(stat_unit,*) 'length_z = ', length_z
    WRITE(stat_unit,*)
#ifdef QMONO
    WRITE(stat_unit,*) 'q_mono viscosity'
#else
    WRITE(stat_unit,*) 'tensor shock viscosity'
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
