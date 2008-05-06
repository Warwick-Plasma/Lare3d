
MODULE diagnostics

  USE shared_data ; USE boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_dt, output_routines, energy_correction

CONTAINS


  SUBROUTINE output_routines(i)   ! i=step index
    ! if halt set to false then code stops
    INTEGER, INTENT(IN) :: i

    INTEGER, SAVE :: output_file = 0
    INTEGER, SAVE :: index = 0, step = 1
    INTEGER :: filehandle, localcellcount
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: data
    LOGICAL :: print_arrays, last_call

    !this output routine uses the same sturcture as needed for mpi output
    !this is more complicated than need for the serial code 
    !rank always equals zero in this serial code
    CHARACTER(LEN=21+Data_Dir_Max_Length) :: filename

    REAL(num) :: j_max, en_ke, en_int
    REAL(num) :: en_b, heating_visc
    REAL(num) :: heating_ohmic, eta_perp_max
    REAL(num) :: j_max_local, total

    IF (i == 0) THEN         ! do just once at the start
       IF (rank == 0) CALL output_log
       IF ((rank == 0) .AND. (.NOT. restart)) WRITE(30) num, 8 ! Change if number of columns in energy diag change
       IF (restart) output_file = restart_snapshot + 1
    END IF

    CALL io_test(i, print_arrays, last_call) !check if snapshot is needed
    IF ((MOD(i,step) .EQ. 0) .OR. last_call) THEN       ! do every (step) steps
       IF (ALLOCATED(curlb)) THEN
         j_max_local = MAXVAL(curlb(1:nx,1:ny,1:nz))
       ELSE
         j_max_local = 0.0_num
       ENDIF
       CALL MPI_ALLREDUCE(j_max_local, j_max, 1, mpireal, MPI_MAX, &
            comm, errcode)

       CALL energy_account(en_b, en_ke, en_int)

       CALL MPI_ALLREDUCE(total_visc_heating, total, 1, mpireal, MPI_SUM, &
            comm, errcode)
       heating_visc = total

       CALL MPI_ALLREDUCE(total_ohmic_heating, total, 1, mpireal, MPI_SUM, &
            comm, errcode)
       heating_ohmic = total

       eta_perp_max = MAXVAL(eta_perp)
       CALL MPI_ALLREDUCE(eta_perp_max, total, 1, mpireal,MPI_MAX,comm,errcode)
       eta_perp_max = total

       IF (rank == 0) WRITE(30) time, j_max, en_b, en_ke, en_int, heating_visc, heating_ohmic,&
            eta_perp_max
       CALL FLUSH(30)
    END IF


    IF (print_arrays) THEN   ! output a snapshot of arrays   
       IF (rank == 0) THEN
          WRITE(20,*) 'output at nstep/time/sbsteps = ',i,time,peak_substeps
          peak_substeps=0
          CALL flush(20)
       END IF

#ifdef NONMPIIO
       ! The old output routines
       WRITE(filename, '(a,"/",i3.3,i4.4,".dat")') TRIM(data_dir), rank, output_file
       !after all of this charcter manipulation the filename for output
       !is of the form rrrnnnn.dat where rrr is the rank of the process
       !and nnnn is the snapshot number

       OPEN(unit=50, FORM = 'UNFORMATTED',STATUS = 'REPLACE', &
            FILE = TRIM(filename))
       WRITE(50) REAL(nx_global,num)-0.1_num,REAL(ny_global,num)-0.1_num, &
            REAL(nz_global,num)-0.1_num
       WRITE(50)  REAL(nx,num)-0.1_num,REAL(ny,num)-0.1_num, REAL(nz,num)-0.1_num
       WRITE(50) time
       WRITE(50) rho(1:nx,1:ny,1:nz),(gamma - 1.0_num) * energy(1:nx,1:ny,1:nz)
       WRITE(50) vx(1:nx,1:ny,1:nz), vy(1:nx,1:ny,1:nz), &
            vz(1:nx,1:ny,1:nz)
       WRITE(50) bx(1:nx,1:ny,1:nz), by(1:nx,1:ny,1:nz), &
            bz(1:nx,1:ny,1:nz)
       WRITE(50) xb(1:nx), yb(1:ny), zb(1:nz)
       WRITE(50) xb(1:nx+1), yb(1:ny+1), zb(1:nz+1)
       CLOSE(50)
#else       
       ! MPI file output routines  --------------------------------------

       ! Create the filename
#ifdef MHDCLUSTER
       WRITE(filename, '("nfs:",a,"/",i4.4,".llld")') TRIM(data_dir), output_file
#else
       WRITE(filename, '(a,"/",i4.4,".llld")') TRIM(data_dir), output_file
#endif
       ! Rank 0 deletes the file as there is no replace option on the file open,
       ! not deleting causes literal overwriting and potential confusion.
       IF (rank == 0) CALL MPI_FILE_DELETE(TRIM(filename), MPI_INFO_NULL, errcode)

       ! Wait for the deletetion before we open the file...
       CALL MPI_BARRIER(comm,errcode)

       ! Open the file, filehandle is the MPI file unit number.
       CALL MPI_FILE_OPEN(comm, TRIM(filename), MPI_MODE_CREATE + MPI_MODE_WRONLY, &
            MPI_INFO_NULL, filehandle, errcode)

       ! If rank 0 then write the file header.
       IF (rank == 0) THEN
          CALL MPI_FILE_WRITE(filehandle,(/nx_global, ny_global, nz_global /), 3, MPI_INTEGER, &
               status, errcode)
          CALL MPI_FILE_WRITE(filehandle,num,1,MPI_INTEGER,status,errcode)
          CALL MPI_FILE_WRITE(filehandle,time,1,mpireal,status,errcode)
          CALL MPI_FILE_WRITE(filehandle,xb_global(0:nx_global),nx_global+1,mpireal,status,errcode)
          CALL MPI_FILE_WRITE(filehandle,yb_global(0:ny_global),ny_global+1,mpireal,status,errcode)
          CALL MPI_FILE_WRITE(filehandle,zb_global(0:nz_global),nz_global+1,mpireal,status,errcode)
       ENDIF

       ! Set my view of the file
       CALL MPI_FILE_SET_VIEW(filehandle, initialdisp, mpireal, subtype,&
            "native", MPI_INFO_NULL, errcode)

       localcellcount = (nx+1) * (ny+1) * (nz+1)
       ALLOCATE(data(0:nx,0:ny,0:nz))

       data = rho(0:nx,0:ny,0:nz)
       CALL MPI_FILE_WRITE_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
       data =(gamma - 1.0_num) * energy(0:nx,0:ny,0:nz) 
       CALL MPI_FILE_WRITE_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
       data = vx(0:nx,0:ny,0:nz)
       CALL MPI_FILE_WRITE_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
       data = vy(0:nx,0:ny,0:nz)
       CALL MPI_FILE_WRITE_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
       data = vz(0:nx,0:ny,0:nz)
       CALL MPI_FILE_WRITE_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
       data = bx(0:nx,0:ny,0:nz)
       CALL MPI_FILE_WRITE_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
       data = by(0:nx,0:ny,0:nz)
       CALL MPI_FILE_WRITE_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
       data = bz(0:nx,0:ny,0:nz)
       CALL MPI_FILE_WRITE_ALL(filehandle,data, localcellcount, mpireal, status, errcode)

       CALL MPI_FILE_CLOSE(filehandle, errcode)
       DEALLOCATE(data)
#endif

       output_file = output_file + 1
    END IF

    IF (last_call .AND. rank == 0) THEN   ! output energy diagnostics etc
       WRITE(20,*) 'final nsteps / time =', i, time
    END IF

  END SUBROUTINE output_routines



  SUBROUTINE io_test(i, print_arrays, last_call)

    INTEGER, INTENT(IN) :: i
    LOGICAL, INTENT(OUT) :: print_arrays, last_call

    REAL(num), SAVE :: t1=0.0_num

    print_arrays = .FALSE.
    last_call = .FALSE.

    IF ((i == 0) .AND. restart) THEN
       t1 = REAL(restart_snapshot + 1,num) *  dt_snapshots
    ENDIF

    IF (time >= t1) THEN
       print_arrays = .TRUE.
       t1 = t1 + dt_snapshots
    END IF

    IF ((time >= t_end) .OR. (i == nsteps)) THEN
       last_call = .TRUE.
       print_arrays = .TRUE.
    END IF

  END SUBROUTINE io_test



  SUBROUTINE set_dt        ! sets CFL limited step
    ! Assumes all variables are defined at the same point. Be careful
    ! with setting 'dt_multiplier' if you expect massive changes across
    ! cells. 

    REAL(num) :: cons, dt1, dt2, dt3, dt4, dt5, dt6, dt7, dt_local, dtr_local
    REAL(num) :: p



    dt_local = largest_number
    dtr_local = largest_number
    cons = gamma * (gamma - 1.0_num)

    DO iz = -1, nz+1
       DO iy = -1, ny+1
          DO ix = -1, nx+1
             ixm = ix - 1
             iym = iy - 1
             izm = iz - 1
             p = (energy(ix,iy,iz)) * rho(ix,iy,iz)
             w1 = bx(ix,iy,iz)**2 + by(ix,iy,iz)**2 + bz(ix,iy,iz)**2
             w2 = SQRT((w1 + p * cons)  &
                  / MAX(rho(ix,iy,iz), none_zero))  &
                  +2.0_num*SQRT(p_visc(ix,iy,iz)/MAX(rho(ix,iy,iz), none_zero))
             w2 = w2 * (1.0_num + visc1)
             dt1 = dxb(ix) / (w2 + ABS(vx(ix,iy,iz)))
             dt2 = dyb(iy) / (w2 + ABS(vy(ix,iy,iz))) 
             dt6 = dzb(iz) / (w2 + ABS(vz(ix,iy,iz))) 


             ! find ideal MHD CFL limit    
             dt_local = MIN(dt_local, dt1, dt2, dt6) 
             ! assumes uniform resistivity hence cautious factor 0.2
             dt3 = 0.2_num * dxb(ix)**2 / MAX(eta(ix,iy,iz),none_zero)
             dt4 = 0.2_num * dyb(iy)**2 / MAX(eta(ix,iy,iz),none_zero)
             dt7 = 0.2_num * dzb(iz)**2 / MAX(eta(ix,iy,iz),none_zero)
             ! adjust to accomodate resistive effects 
             dtr_local = MIN(dtr_local, dt3, dt4, dt7) 
             ! correct to stop overlapping of Lagrangian cells	
             w1 = ABS(vx(ix,iy,iz) - vx(ixm,iy,iz)) / dxb(ix)   &
                  + ABS(vy(ix,iy,iz) - vy(ix,iym,iz)) / dyb(iy)   &
                  + ABS(vz(ix,iy,iz) - vz(ix,iy,izm)) / dzb(iz)   
             dt5 = 1.0_num / MAX(w1, none_zero)
             dt_local = MIN(dt_local, dt5)
          END DO
       END DO
    END DO


    CALL MPI_ALLREDUCE(dt_local, dt, 1, mpireal, MPI_MIN, comm, errcode)
    CALL MPI_ALLREDUCE(dtr_local, dtr, 1, mpireal, MPI_MIN, comm, errcode)

    dt = dt_multiplier * dt
    dtr = dt_multiplier * dtr

  END SUBROUTINE set_dt



  SUBROUTINE energy_account(energy_b, energy_ke, energy_int)

    REAL(num), INTENT(OUT) :: energy_b, energy_ke, energy_int
    REAL(dbl) :: energy_b_local, energy_ke_local, energy_int_local
    REAL(num) :: cv_v, rho_v, a, b, c

    energy_b_local = 0.0_dbl
    energy_ke_local= 0.0_dbl
    energy_int_local = 0.0_dbl

    DO iz = 1, nz
       DO iy = 1, ny
          DO ix = 1, nx
             ixm = ix - 1
             iym = iy - 1
             izm = iz - 1
             w2 = (bx(ix,iy,iz)**2 + bx(ixm,iy,iz)**2)/2.0_dbl   
             w3 = (by(ix,iy,iz)**2 + by(ix,iym,iz)**2)/2.0_dbl   
             w4 = (bz(ix,iy,iz)**2 + bz(ix,iy,izm)**2)/2.0_dbl 
             w1 = (w2 + w3 + w4) / 2.0_dbl  
             energy_b_local = energy_b_local + w1*cv(ix,iy,iz)

             energy_int_local = energy_int_local &
                  + energy(ix,iy,iz)*rho(ix,iy,iz)*cv(ix,iy,iz)
          END DO
       END DO
    END DO

    DO iz = 0, nz
       DO iy = 0, ny
          DO ix = 0, nx
             ixp = ix + 1
             iyp = iy + 1
             izp = iz + 1
             !WARNING the KE is summed on the vertices
             rho_v = (rho(ix,iy,iz)*cv(ix,iy,iz)+rho(ixp,iy,iz)*cv(ixp,iy,iz) &
                  +rho(ix,iyp,iz)*cv(ix,iyp,iz)+rho(ixp,iyp,iz)*cv(ixp,iyp,iz)&
                  +rho(ix,iy,izp)*cv(ix,iy,izp)+rho(ixp,iy,izp)*cv(ixp,iy,izp)&
                  +rho(ix,iyp,izp)*cv(ix,iyp,izp)  &
                  +rho(ixp,iyp,izp)*cv(ixp,iyp,izp))
             cv_v = (cv(ix,iy,iz)+cv(ixp,iy,iz)+cv(ix,iyp,iz)+cv(ixp,iyp,iz)  &
                  +cv(ix,iy,izp)+cv(ixp,iy,izp)+cv(ix,iyp,izp)+cv(ixp,iyp,izp))
             rho_v = rho_v / cv_v
             cv_v = cv_v / 8.0_dbl
             w1 = rho_v*(vx(ix,iy,iz)**2+vy(ix,iy,iz)**2+vz(ix,iy,iz)**2)*cv_v
             IF((ix == 0) .OR. (ix == nx)) THEN
                w1 = w1 / 2.0_dbl
             END IF
             IF((iy == 0) .OR. (iy == ny)) THEN
                w1 = w1 / 2.0_dbl
             END IF
             IF((iz == 0) .OR. (iz == nz)) THEN
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

    delta_ke = - delta_ke
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




  SUBROUTINE output_log  !writes basic data to 'lare2d.dat'

    REAL(num) :: temp
    IF (.NOT. restart) THEN
       WRITE(20,*) 'nprocx, nprocy, nprocz = ',nprocx, nprocy, nprocz
       WRITE(20,*) 'nx, ny, nz = ', nx, ny, nz
       WRITE(20,*)
       WRITE(20,*) 'length_x = ', length_x
       WRITE(20,*) 'length_y = ', length_y
       WRITE(20,*) 'length_z = ', length_z
       WRITE(20,*)
       WRITE(20,*) 'x boundary condition = ', xbc_right
       WRITE(20,*) 'y boundary condition = ', ybc_up
       WRITE(20,*) 'z boundary condition = ', zbc_back
       WRITE(20,*)
       WRITE(20,*) 'linear viscosity coeff = ', visc1
       WRITE(20,*) 'quadratic viscosity coeff = ', visc2
       WRITE(20,*) 'background viscosity coeff = ', visc3
       WRITE(20,*) 'eta0 = ', eta0
       WRITE(20,*)
       WRITE(20,*) 'Resistive MHD =', resistiveMHD
       WRITE(20,*)
       WRITE(20,*) 't_start, t_end = ', time, t_end
       WRITE(20,*) 'nsteps =',nsteps
       WRITE(20,*)
    ELSE
       WRITE(20,*) '*********************'
       WRITE(20,*) 'Restarted from output', restart_snapshot,' on ', nproc, 'process elements.'
       WRITE(20,*) 'nprocx, nprocy, nprocz = ',nprocx, nprocy, nprocz
    ENDIF

    CALL flush (20)

  END SUBROUTINE output_log


END MODULE diagnostics







