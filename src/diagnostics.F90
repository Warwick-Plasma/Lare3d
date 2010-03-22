!*************************************************************************
! Controls all I/O and diagnostics. Output files are 'lare2d.dat',
! 'control.dat', 'en.dat' and a series of snapshots in 'fort.5x'
! The idl package in 'plot.pro' gives simple loading and surface
! plotting based on these files. This isn't documented but is very simple!
!*************************************************************************
MODULE diagnostics

  USE shared_data; USE boundary; USE normalise
  USE eos
  USE output_cartesian
  USE output
  USE iocontrol

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_dt, output_routines, energy_correction

CONTAINS

  SUBROUTINE output_routines(i) ! i = step index

    ! if halt set to false then code stops
    INTEGER, INTENT(IN) :: i

    INTEGER, PARAMETER :: out = 1000
    INTEGER, SAVE :: index = 1, step = 1
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: data
    LOGICAL :: print_arrays, last_call
    REAL(num), DIMENSION(3) :: stagger = 0.0_num
    INTEGER, DIMENSION(3) :: dims

    ! this output routine uses the same sturcture as needed for mpi output
    ! this is more complicated than need for the serial code
    ! rank always equals zero in this serial code
    CHARACTER(LEN = 9+data_dir_max_length+n_zeros) :: filename
    CHARACTER(LEN = 35) :: filename_desc

    REAL(num) :: t_out = 0.0_num
    REAL(num) :: en_ke = 0.0_num, en_int = 0.0_num
    REAL(num) :: en_b = 0.0_num, heating_visc = 0.0_num
    REAL(num) :: heating_ohmic = 0.0_num
    REAL(num) :: total

    dims = (/ nx_global+1, ny_global+1, nz_global+1 /)

    IF (nsteps >= out) step = nsteps / out + 1 ! make sure output fits arrays
    IF (i == 0 .AND. rank == 0) THEN ! done just once at the start
      CALL output_log
      IF (.NOT. restart) WRITE(30) num, 6
    END IF

    IF (MOD(i, step) .EQ. 0 .OR. last_call) THEN ! do every (step) steps
      t_out = time
      CALL energy_account(en_b, en_ke, en_int)

      CALL MPI_ALLREDUCE(total_visc_heating, total, 1, mpireal, MPI_SUM, &
          comm, errcode)

      heating_visc = total

      CALL MPI_ALLREDUCE(total_ohmic_heating, total, 1, mpireal, MPI_SUM, &
          comm, errcode)

      heating_ohmic = total

      IF (rank .EQ. 0) THEN
        WRITE(30) t_out, en_b, en_ke, en_int
        WRITE(30) heating_visc, heating_ohmic
      END IF

      index = index + 1
    END IF

    CALL io_test(i, print_arrays, last_call) ! check if snapshot is needed

    IF (print_arrays) THEN ! output a snapshot of arrays
      IF (rank .EQ. 0) THEN
        WRITE(20, *) "Dumping ", output_file, " at time", time * t0
        CALL FLUSH(20)
      END IF

      ! Set the filename
      WRITE(filename_desc, '("(''nfs:'', a, ''/'', i", i3.3, ".", i3.3, &
          & ", ''.cfd'')")'), n_zeros, n_zeros
      WRITE(filename, filename_desc) TRIM(data_dir), output_file

      CALL cfd_open(filename, rank, comm, MPI_MODE_CREATE + MPI_MODE_WRONLY)
      CALL cfd_write_snapshot_data(REAL(time * t0,dbl), i, 0)

      ALLOCATE(data(0:nx, 0:ny, 0:nz))
      CALL cfd_write_3d_cartesian_grid("Grid", "Grid", &
          xb_global(0:nx_global) * L0, yb_global(0:ny_global) * L0, &
          zb_global(0:nz_global) * L0, 0)

      IF (dump_mask(1))  THEN
        data = rho(0:nx, 0:ny, 0:nz) * rho0
        CALL cfd_write_3d_cartesian_variable_parallel("Rho", "Fluid", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(2))  THEN
        data = energy(0:nx, 0:ny, 0:nz) * energy0
        CALL cfd_write_3d_cartesian_variable_parallel("Energy", "Fluid", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(3))  THEN
        data = vx(0:nx, 0:ny, 0:nz) * vel0
        CALL cfd_write_3d_cartesian_variable_parallel("Vx", "Velocity", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(4))  THEN
        data = vy(0:nx, 0:ny, 0:nz) * vel0
        CALL cfd_write_3d_cartesian_variable_parallel("Vy", "Velocity", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(5))  THEN
        data = vz(0:nx, 0:ny, 0:nz) * vel0
        CALL cfd_write_3d_cartesian_variable_parallel("Vz", "Velocity", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(6))  THEN
        data = bx(0:nx, 0:ny, 0:nz) * B0
        CALL cfd_write_3d_cartesian_variable_parallel("Bx", "Magnetic_Field", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(7))  THEN
        data = by(0:nx, 0:ny, 0:nz) * B0
        CALL cfd_write_3d_cartesian_variable_parallel("By", "Magnetic_Field", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(8))  THEN
        data = bz(0:nx, 0:ny, 0:nz) * B0
        CALL cfd_write_3d_cartesian_variable_parallel("Bz", "Magnetic_Field", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(9)) THEN
        DO iz = 0, nz
          DO iy = 0, ny
            DO ix = 0, nx
              CALL get_temp(rho(ix, iy, iz), energy(ix, iy, iz), &
                  eos_number, ix, iy, iz, data(ix, iy, iz))
            END DO
          END DO
        END DO
        data = data * temp0
        CALL cfd_write_3d_cartesian_variable_parallel("Temperature", "Fluid", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(10))  THEN
        DO iz = 0, nz
          DO iy = 0, ny
            DO ix = 0, nx
              CALL get_pressure(rho(ix, iy, iz), energy(ix, iy, iz), &
                  eos_number, ix, iy, iz, data(ix, iy, iz))
            END DO
          END DO
        END DO
        data = data * pressure0
        CALL cfd_write_3d_cartesian_variable_parallel("Pressure", "Fluid", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(11)) THEN
        data = SQRT(gamma*(gamma - 1.0_num) * energy(0:nx,0:ny,0:nz)) * vel0
        CALL cfd_write_3d_cartesian_variable_parallel("cs", "Fluid", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(12)) THEN
        data = parallel_current(0:nx, 0:ny, 0:nz) * j0
        CALL cfd_write_3d_cartesian_variable_parallel("j_par", "PIP", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(13)) THEN
        data = perp_current(0:nx, 0:ny, 0:nz) * j0
        CALL cfd_write_3d_cartesian_variable_parallel("j_perp", "PIP", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(14)) THEN
        data = xi_n(0:nx, 0:ny, 0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel("neutral_fraction", &
            "PIP", dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(15)) THEN
        data = eta_perp(0:nx, 0:ny, 0:nz) * res0
        CALL cfd_write_3d_cartesian_variable_parallel("eta_perp", "PIP", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(16)) THEN
        data = eta(0:nx, 0:ny, 0:nz) * res0
        CALL cfd_write_3d_cartesian_variable_parallel("eta", "PIP", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(17)) THEN
        data = jx_r(0:nx, 0:ny, 0:nz) * j0
        CALL cfd_write_3d_cartesian_variable_parallel("jx", "current", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(18)) THEN
        data = jy_r(0:nx, 0:ny, 0:nz) * j0
        CALL cfd_write_3d_cartesian_variable_parallel("jy", "current", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(19)) THEN
        data = jz_r(0:nx, 0:ny, 0:nz) * j0
        CALL cfd_write_3d_cartesian_variable_parallel("jz", "current", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF
#ifdef PARALLEL_DEBUG
      data = rank
      CALL cfd_write_3d_cartesian_variable_parallel("rank", "parallel", &
          dims, stagger, "Grid", "Grid", data, subtype)
#endif

      ! Close the file
      CALL cfd_close()

      output_file = output_file + 1

    END IF

    IF (last_call .AND. rank == 0) THEN ! output energy diagnostics etc
      WRITE(20, *) 'final nsteps / time =', i, time * t0
    END IF

  END SUBROUTINE output_routines



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

    ! Assumes all variables are defined at the same point. Be careful
    ! with setting 'dt_multiplier' if you expect massive changes across
    ! cells.

    REAL(num) :: cons, dt1, dt3, dt4, dt5, dt_local, dxlocal
    REAL(num) :: dtr_local, dth_local, cs

    dt_local = largest_number
    dtr_local = largest_number
    dth_local = largest_number
    cons = gamma * (gamma - 1.0_num)

    DO iz = -1, nz+2
      izm = iz - 1
      DO iy = -1, ny+2
        iym = iy - 1
        DO ix = -1, nx+2
          ixm = ix - 1

          w1 = bx(ix, iy, iz)**2 + by(ix, iy, iz)**2 + bz(ix, iy, iz)**2 
          cs = cons * energy(ix,iy,iz)      ! sound speed squared

          w2 = SQRT(cs + w1 / MAX(rho(ix, iy, iz), none_zero)) &
              + 2.0_num * SQRT(p_visc(ix, iy, iz) &
              / MAX(rho(ix, iy, iz), none_zero))

          w2 = w2 * (1.0_num + visc1)

          ! find ideal MHD CFL limit
          dt1 = dxb(ix) / (w2 + ABS(vx(ix, iy, iz)))
          dt_local = MIN(dt_local, dt1)
          dt1 = dyb(iy) / (w2 + ABS(vy(ix, iy, iz)))
          dt_local = MIN(dt_local, dt1)
          dt1 = dzb(iz) / (w2 + ABS(vz(ix, iy, iz)))
          dt_local = MIN(dt_local, dt1)

          ! note resistive limits assumes uniform resistivity hence cautious
          ! factor 0.2
          dxlocal = 1.0_num / (1.0_num / dxb(ix)**2 &
              + 1.0_num / dyb(iy)**2 + 1.0_num / dzb(iz)**2)

          IF (cowling_resistivity) THEN
            dt3 = 0.2_num * dxlocal &
                / MAX(MAX(eta(ix, iy, iz), eta_perp(ix, iy, iz)), none_zero)
          ELSE
            dt3 = 0.2_num * dxlocal / MAX(eta(ix, iy, iz), none_zero)
          END IF

          ! Hall MHD CFL limit
          dt4 = 0.75_num * rho(ix, iy, iz) * MIN(dxb(ix), dyb(iy), dzb(iz))**2 &
              / MAX(lambda_i(ix, iy, iz) * SQRT(w1), none_zero)

          ! adjust to accomodate resistive effects
          dtr_local = MIN(dtr_local, dt3)
          dth_local = MIN(dth_local, dt4)

          ! correct to stop overlapping of Lagrangian cells

          w1 = ABS(vx(ix, iy, iz) - vx(ixm, iy, iz)) / dxb(ix) &
              + ABS(vy(ix, iy, iz) - vy(ix, iym, iz)) / dyb(iy) &
              + ABS(vz(ix, iy, iz) - vz(ix, iy, izm)) / dzb(iz)

          dt5 = 1.0_num / MAX(w1, none_zero)
          dt_local = MIN(dt_local, dt5)
        END DO
      END DO
    END DO

    CALL MPI_ALLREDUCE(dt_local, dt, 1, mpireal, MPI_MIN, comm, errcode)
    CALL MPI_ALLREDUCE(dtr_local, dtr, 1, mpireal, MPI_MIN, comm, errcode)
    CALL MPI_ALLREDUCE(dth_local, dth, 1, mpireal, MPI_MIN, comm, errcode)

    dtr = dt_multiplier * dtr
    dth = dt_multiplier * dth
    dt = dt_multiplier * dt

    time = time + dt

  END SUBROUTINE set_dt



  SUBROUTINE energy_account(energy_b, energy_ke, energy_int)

    REAL(num), INTENT(OUT) :: energy_b, energy_ke, energy_int
    REAL(dbl) :: energy_b_local, energy_ke_local, energy_int_local
    REAL(num) :: cv_v, rho_v, a, b, c

    energy_b_local = 0.0_dbl
    energy_ke_local = 0.0_dbl
    energy_int_local = 0.0_dbl

    DO iz = 1, nz
      izm = iz - 1
      DO iy = 1, ny
        iym = iy - 1
        DO ix = 1, nx
          ixm = ix - 1

          w2 = (bx(ix, iy, iz)**2 + bx(ixm, iy, iz)**2) / 2.0_dbl
          w3 = (by(ix, iy, iz)**2 + by(ix, iym, iz)**2) / 2.0_dbl
          w4 = (bz(ix, iy, iz)**2 + bz(ix, iy, izm)**2) / 2.0_dbl
          w1 = (w2 + w3 + w4) / 2.0_dbl
          energy_b_local = energy_b_local + w1 * cv(ix, iy, iz)

          energy_int_local = energy_int_local &
              + energy(ix, iy, iz) * rho(ix, iy, iz) * cv(ix, iy, iz)
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
          rho_v = (rho(ix, iy, iz) * cv(ix, iy, iz) &
              + rho(ixp, iy , iz ) * cv(ixp, iy , iz ) &
              + rho(ix , iyp, iz ) * cv(ix , iyp, iz ) &
              + rho(ixp, iyp, iz ) * cv(ixp, iyp, iz ) &
              + rho(ix , iy , izp) * cv(ix , iy , izp) &
              + rho(ixp, iy , izp) * cv(ixp, iy , izp) &
              + rho(ix , iyp, izp) * cv(ix , iyp, izp) &
              + rho(ixp, iyp, izp) * cv(ixp, iyp, izp))

          cv_v = (cv(ix, iy, iz) + cv(ixp, iy, iz) &
              + cv(ix, iyp, iz ) + cv(ixp, iyp, iz ) &
              + cv(ix, iy , izp) + cv(ixp, iy , izp) &
              + cv(ix, iyp, izp) + cv(ixp, iyp, izp))

          rho_v = rho_v / cv_v
          cv_v = cv_v / 8.0_dbl
          w1 = rho_v * cv_v &
              * (vx(ix, iy, iz)**2 + vy(ix, iy, iz)**2 + vz(ix, iy, iz)**2)

          IF ((ix == 0) .OR. (ix == nx)) THEN
            w1 = w1 / 2.0_dbl
          END IF

          IF ((iy == 0) .OR. (iy == ny)) THEN
            w1 = w1 / 2.0_dbl
          END IF

          IF ((iz == 0) .OR. (iz == nz)) THEN
            w1 = w1 / 2.0_dbl
          END IF

          energy_ke_local = energy_ke_local + w1 / 2.0_dbl
        END DO
      END DO
    END DO

    a = REAL(energy_ke_local, num)
    b = REAL(energy_b_local, num)
    c = REAL(energy_int_local, num)

    CALL MPI_ALLREDUCE(a, energy_ke, 1, mpireal, MPI_SUM, &
        comm, errcode)
    CALL MPI_ALLREDUCE(b, energy_b, 1, mpireal, MPI_SUM, &
        comm, errcode)
    CALL MPI_ALLREDUCE(c, energy_int, 1, mpireal, MPI_SUM, &
        comm, errcode)

  END SUBROUTINE energy_account



  SUBROUTINE energy_correction

    delta_ke = -delta_ke
    WHERE (delta_ke < 0.0_num) delta_ke = 0.0_num
    delta_ke = delta_ke / (rho * cv)

    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          energy(ix, iy, iz) = energy(ix, iy, iz) + delta_ke(ix, iy, iz)
        END DO
      END DO
    END DO
    CALL energy_bcs

  END SUBROUTINE energy_correction



  SUBROUTINE output_log ! writes basic data to 'lare3d.dat'

    WRITE(20, *) 'Density normalisation = ', RHO0, ' kg m^(-3)'
    WRITE(20, *) 'Specific energy density normalisation = ', ENERGY0, ' K'
    WRITE(20, *) 'Velocity normalisation = ', VEL0, ' m s^(-1)'
    WRITE(20, *) 'Magnetic field normalisation = ', B0, ' T'
    WRITE(20, *) 'Length normalisation = ', L0, ' m'
    WRITE(20, *) 'Time normalisation = ', T0, ' s'
    WRITE(20, *) 'Viscosity normalisation = ', VISC0, ' m^2s^(-1)'
    WRITE(20, *) 'Gravity normalisation = ', GRAV0, ' ms^(-2)'
    WRITE(20, *) 'Resistivity normalisation = ', RES0, ' m^2s^(-1)'
    WRITE(20, *) 'Thermal conductivity normalisation = ', KAPPA0, &
        ' kgm^(-1)s^(-1) '
    WRITE(20, *) 'Temperature normalisation =', TEMP0, 'K'
    WRITE(20, *) 'Pressure normalisation =', PRESSURE0, 'Pa'

    WRITE(20, *) 'nprocx, nprocy, nproca = ', nprocx, nprocy, nprocz
    WRITE(20, *) 'nx, ny, nz = ', nx, ny, nz
    WRITE(20, *)
    WRITE(20, *) 'length_x = ', length_x * L0
    WRITE(20, *) 'length_y = ', length_y * L0
    WRITE(20, *) 'length_z = ', length_z * L0
    WRITE(20, *)
#ifndef Q_MONO
    WRITE(20, *) 'tensor shock viscosity'
#else
    WRITE(20, *) 'q_mono viscosity'
#endif
    WRITE(20, *) 'linear viscosity coeff = ', visc1
    WRITE(20, *) 'quadratic viscosity coeff = ', visc2
    WRITE(20, *) 'uniform tensor viscosity coeff = ', visc3 * VISC0
    WRITE(20, *) 'j_max = ', j_max * B0 / L0
    WRITE(20, *) 'vc = ', vc * B0 / L0
    WRITE(20, *) 'eta0 = ', eta0 * RES0
    WRITE(20, *) 'eta_background = ', eta_background * RES0
    WRITE(20, *) 'kappa = ', kappa_0 * KAPPA0
    WRITE(20, *)
    WRITE(20, *) 't_start, t_end = ', time * T0, t_end * T0
    WRITE(20, *) 'nsteps =', nsteps
    WRITE(20, *)

  END SUBROUTINE output_log

END MODULE diagnostics
