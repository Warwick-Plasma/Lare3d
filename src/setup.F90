MODULE setup
  USE shared_data
  IMPLICIT NONE

  PRIVATE :: initialise_arrays
  PUBLIC :: minimal_init,after_control, grid
  PUBLIC :: open_files, close_files,restart_data

  REAL(num), DIMENSION(:), ALLOCATABLE :: dxnew, dynew, dznew

CONTAINS

  SUBROUTINE minimal_init

    nprocx=0
    nprocy=0
    nprocz=0

    length_x=1.0_num
    length_y=1.0_num
    length_z=1.0_num

    IF (KIND(1.0_num) .EQ. 4) mpireal=MPI_REAL

  END SUBROUTINE minimal_init

  SUBROUTINE after_control

    IF (.NOT. restart)  restart_snapshot = 0

    CALL initialise_arrays 

  END SUBROUTINE after_control

  SUBROUTINE grid                 ! stretched and staggered grid

    REAL(num) :: dx, dy, dz, dxmin, dzmin, dymin, xcstar, ycstar, zcstar
    INTEGER :: ix, iy

    ALLOCATE(xb_global(-2:nx_global+2), dxnew(-2:nx_global+2))
    ALLOCATE(yb_global(-2:ny_global+2), dynew(-2:ny_global+2))
    ALLOCATE(zb_global(-2:nz_global+2), dznew(-2:nz_global+2))

    dx = 1.0_num / REAL(nx_global)       ! initially assume uniform grid
    dy = 1.0_num / REAL(ny_global)
    dz = 1.0_num / REAL(nz_global)

    xb_global(0) = 0.0_num               !grid cell boundary for x coordinates
    DO ix = -2, nx_global+2
       xb_global(ix) =  xb_global(0) + REAL(ix, num)*dx
    END DO
    xb_global = xb_global * (x_end-x_start) + x_start
    IF (x_stretch) CALL stretch_x     ! stretch grid ?
    !define position of ghost cells using sizes of adjacent cells
    xb_global(nx_global+1) = xb_global(nx_global) + &
         (xb_global(nx_global) - xb_global(nx_global-1)) 
    ! needed for ghost cell
    xb_global(nx_global+2) = xb_global(nx_global+1) &
         + (xb_global(nx_global+1) - xb_global(nx_global))
    xb_global(-1) = xb_global(0) &
         - (xb_global(1) - xb_global(0))
    xb_global(-2) = xb_global(-1) &
         - (xb_global(0) - xb_global(-1))
    xb = xb_global(coordinates(3)*nx-2:coordinates(3)*nx+nx+2)

    DO ix = -1, nx+2
       ixm = ix - 1
       xc(ix) = 0.5_num*(xb(ixm) + xb(ix))     ! cell centre
    END DO
    DO ix = -1, nx+1
       ixp = ix + 1
       dxc(ix) = xc(ixp) - xc(ix)    ! distance between centres
    END DO
    IF (coordinates(3)==nprocx-1) THEN
       dxc(nx+2) = dxc(nx+1)
    ELSE
       xcstar = 0.5_num*(xb(nx+2) + xb_global(coordinates(3)*nx+nx+3))
       dxc(nx+2) = xcstar - xc(nx+2)
    END IF
    DO ix = -1, nx+2
       ixm = ix - 1
       dxb(ix) = xb(ix) - xb(ixm)    ! cell width
    END DO

    yb_global(0) = 0.0_num               !grid cell boundary for y coordinates
    DO iy = -2, ny_global+2
       yb_global(iy) =  yb_global(0) + REAL(iy, num)*dy
    END DO
    yb_global = yb_global * (y_end-y_start) + y_start
    IF (y_stretch) CALL stretch_y     ! stretch grid ?
    !define position of ghost cells using sizes of adjacent cells
    yb_global(ny_global+1) = yb_global(ny_global) + &
         (yb_global(ny_global) - yb_global(ny_global-1)) 
    ! needed for ghost cell
    yb_global(ny_global+2) = yb_global(ny_global+1) &
         + (yb_global(ny_global+1) - yb_global(ny_global))
    yb_global(-1) = yb_global(0) &
         - (yb_global(1) - yb_global(0))
    yb_global(-2) = yb_global(-1) &
         - (yb_global(0) - yb_global(-1))
    yb = yb_global(coordinates(2)*ny-2:coordinates(2)*ny+ny+2)
    DO iy = -1, ny+2
       iym = iy - 1
       yc(iy) = 0.5_num*(yb(iym) + yb(iy))     ! cell centre
    END DO
    DO iy = -1, ny+1
       iyp = iy + 1
       dyc(iy) = yc(iyp) - yc(iy)    ! distance between centres
    END DO
    IF (coordinates(2)==nprocy-1) THEN
       dyc(ny+2) = dyc(ny+1)
    ELSE
       ycstar = 0.5_num*(yb(ny+2) + yb_global(coordinates(2)*ny+ny+3))
       dyc(ny+2) = ycstar - yc(ny+2)
    END IF
    DO iy = -1, ny+2
       iym = iy - 1
       dyb(iy) = yb(iy) - yb(iym)    ! cell width
    END DO

    zb_global(0) = 0.0_num               !grid cell boundary for z coordinates
    DO iz = -2, nz_global+2
       zb_global(iz) =  zb_global(0) + REAL(iz, num)*dz
    END DO
    zb_global = zb_global * (z_end-z_start) + z_start
    IF (z_stretch) CALL stretch_z     ! stretch grid ?
    !define position of ghost cells using sizes of adjacent cells
    zb_global(nz_global+1) = zb_global(nz_global) + &
         (zb_global(nz_global) - zb_global(nz_global-1)) 
    ! needed for ghost cell
    zb_global(nz_global+2) = zb_global(nz_global+1) &
         + (zb_global(nz_global+1) - zb_global(nz_global))
    zb_global(-1) = zb_global(0) &
         - (zb_global(1) - zb_global(0))
    zb_global(-2) = zb_global(-1) &
         - (zb_global(0) - zb_global(-1))
    zb = zb_global(coordinates(1)*nz-2:coordinates(1)*nz+nz+2)
    DO iz = -1, nz+2
       izm = iz - 1
       zc(iz) = 0.5_num*(zb(izm) + zb(iz))     ! cell centre
    END DO
    DO iz = -1, nz+1
       izp = iz + 1
       dzc(iz) = zc(izp) - zc(iz)    ! distance between centres
    END DO
    IF (coordinates(1)==nprocz-1) THEN
       dzc(nz+2) = dzc(nz+1)
    ELSE
       zcstar = 0.5_num*(zb(nz+2) + zb_global(coordinates(1)*nz+nz+3))
       dzc(nz+2) = zcstar - zc(nz+2)
    END IF
    DO iz = -1, nz+2
       izm = iz - 1
       dzb(iz) = zb(iz) - zb(izm)    ! cell width
    END DO

    DO ix = -1, nx+2
       DO iy = -1, ny+2
          DO iz = -1, nz+2
             cv(ix,iy,iz) = dxb(ix) * dyb(iy) * dzb(iz)    ! define the cell area
          END DO
       END DO
    END DO

    CALL MPI_ALLREDUCE(MINVAL(dxb(1:nx)), dxmin, 1, mpireal, MPI_MIN, &
         comm, errcode)
    CALL MPI_ALLREDUCE(MINVAL(dyb(1:ny)), dymin, 1, mpireal, MPI_MIN, &
         comm, errcode)
    CALL MPI_ALLREDUCE(MINVAL(dzb(1:nz)), dzmin, 1, mpireal, MPI_MIN, &
         comm, errcode)
    min_grid_spacing = MIN(dxmin, dymin, dzmin)

    DEALLOCATE(dxnew, dynew, dznew)


  END SUBROUTINE grid



  SUBROUTINE stretch_x   ! replace with any stretching algorithm as needed

    REAL(num) :: width, dx, L, f, lx_new

    lx_new = 200.0_num                ! new tolal length
    L = length_x / 2.0_num       ! centre of tanh stretching in unstretched coordinates
    width = length_x / 10.0_num       ! width of tanh stretching in unstretched coordinates

    f = (lx_new - length_x)/(length_x - L)/2.0_num

    dx = length_x / REAL(nx_global,num)  
    dxnew = dx + f*(1.0_num+TANH((ABS(xb_global)-L)/width))*dx

    DO ix = 1, nx_global+2
       xb_global(ix) = xb_global(ix-1) + dxnew(ix)
    ENDDO

  END SUBROUTINE stretch_x



  SUBROUTINE stretch_y 

    REAL(num) :: width, dy, L, f, ly_new

    ly_new = 150.0_num                ! new tolal length
    L = length_y / 2.0_num       ! centre of tanh stretching in unstretched coordinates
    width = length_y / 10.0_num       ! width of tanh stretching in unstretched coordinates

    f = (ly_new - length_y)/(length_y - L)/2.0_num

    dy = length_y / REAL(ny_global,num)  
    dynew = dy + f*(1.0_num+TANH((ABS(yb_global)-L)/width))*dy

    DO iy = 1, ny_global+2
       yb_global(iy) = yb_global(iy-1) + dynew(iy)
    ENDDO

  END SUBROUTINE stretch_y




  SUBROUTINE stretch_z 

    REAL(num) :: width, dz, L, f, lz_new

    lz_new = 33.0_num                ! new tolal length
    L = 2.0_num * length_z / 3.0_num       ! centre of tanh stretching in unstretched coordinates
    width = length_z / 10.0_num       ! width of tanh stretching in unstretched coordinates

    f = (lz_new - length_z)/(length_z - L)/2.0_num

    dz = length_z / REAL(nz_global,num)  
    dznew = dz + f*(1.0_num+TANH((ABS(zb_global)-L)/width))*dz

    DO iz = 1, nz_global+2
       zb_global(iz) = zb_global(iz-1) + dznew(iz)
    ENDDO

  END SUBROUTINE stretch_z


  SUBROUTINE open_files

    CHARACTER(LEN=11+Data_Dir_Max_Length) :: file2
    CHARACTER(LEN=7+Data_Dir_Max_Length) :: file3

    IF (rank == 0) THEN
       IF (.NOT. restart) THEN
          WRITE(file2, '(a,"/lare3d.dat")') TRIM(data_dir)
          OPEN(unit=20, STATUS='REPLACE',FILE=file2)
          WRITE(file3, '(a,"/en.dat")') TRIM(data_dir)
          OPEN(unit=30, STATUS = 'REPLACE', FORM='BINARY', FILE = file3)
       ELSE
          WRITE(file2, '(a,"/lare3d.dat")') TRIM(data_dir)
          OPEN(unit=20, STATUS='OLD',POSITION='APPEND',FILE=file2)
          WRITE(file3, '(a,"/en.dat")') TRIM(data_dir)
          OPEN(unit=30, STATUS = 'OLD',POSITION='APPEND',FORM='BINARY', FILE = file3)
       ENDIF
    END IF

  END SUBROUTINE open_files



  SUBROUTINE close_files

    IF (rank == 0) THEN
       CLOSE(unit=20)
       CLOSE(unit=30)
    END IF

  END SUBROUTINE close_files



  SUBROUTINE restart_data

    CHARACTER(LEN=21+Data_Dir_Max_Length) :: filename
    INTEGER :: filehandle, localcellcount, tn(3),sof
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: data

#ifdef NONMPIIO
    print*, "Cannot restart from non-MPI io yet"
    STOP
#else

    ! Create the filename for the last snapshot
#ifdef MHDCLUSTER
    WRITE(filename, '("nfs:",a,"/",i4.4,".llld")') TRIM(data_dir), restart_snapshot
#else
    WRITE(filename, '(a,"/",i4.4,".llld")') TRIM(data_dir), restart_snapshot
#endif
    ! Open the file
    CALL MPI_FILE_OPEN(comm, TRIM(filename), MPI_MODE_RDONLY, &
         MPI_INFO_NULL, filehandle, errcode)

    CALL MPI_FILE_READ_ALL(filehandle,tn,3,MPI_INTEGER,status,errcode)
    CALL MPI_FILE_READ_ALL(filehandle,sof,1,MPI_INTEGER,status,errcode)
    CALL MPI_FILE_READ_ALL(filehandle,time,1,mpireal,status,errcode)

    IF ((tn(1) /= nx_global).OR.(tn(2)/=ny_global).OR.(tn(3) /= nz_global)) THEN
       IF (rank ==0) THEN
          WRITE(20,*) 'Error: global dimensions do not match restart file.'
       ENDIF
       CALL MPI_FILE_CLOSE(filehandle, errcode)
       CALL close_files
       CALL MPI_FINALIZE(errcode)
       STOP
    ENDIF

    IF (sof /= num) THEN
       IF (rank == 0) THEN
        WRITE(20,*) "Error precision of restart file doesn't match precision of code"
       ENDIF
       CALL MPI_FILE_CLOSE(filehandle,errcode)
       CALL close_files
       CALL MPI_FINALIZE(errcode)
       STOP
    ENDIF


    ! Set my view of the file
    CALL MPI_FILE_SET_VIEW(filehandle, initialdisp, mpireal, subtype,&
         "native", MPI_INFO_NULL, errcode)

    localcellcount = (nx+1) * (ny+1) * (nz+1)
    ALLOCATE(data(0:nx,0:ny,0:nz))
    CALL MPI_FILE_READ_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
    rho(0:nx,0:ny,0:nz) = data
    CALL MPI_FILE_READ_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
    energy(0:nx,0:ny,0:nz) = data / (gamma - 1.0_num)
    CALL MPI_FILE_READ_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
    vx(0:nx,0:ny,0:nz) = data
    CALL MPI_FILE_READ_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
    vy(0:nx,0:ny,0:nz) = data
    CALL MPI_FILE_READ_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
    vz(0:nx,0:ny,0:nz) = data
    CALL MPI_FILE_READ_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
    bx(0:nx,0:ny,0:nz) = data
    CALL MPI_FILE_READ_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
    by(0:nx,0:ny,0:nz) = data
    CALL MPI_FILE_READ_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
    bz(0:nx,0:ny,0:nz) = data

    DEALLOCATE(data)

#endif

  CALL MPI_BARRIER(comm,errcode)

  END SUBROUTINE restart_data



  ! This routine initialises arrays to sensible starting values

  SUBROUTINE initialise_arrays

    curlb = 0.0_num
    p_visc = 0.0_num
    eta = 0.0_num
    eta_perp = 0.0_num

    grav = 0.0_num

  END SUBROUTINE initialise_arrays





END MODULE setup

