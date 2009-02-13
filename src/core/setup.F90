MODULE setup

  USE shared_data
  USE normalise
  USE iocommon
  USE iocontrol
  USE input
  USE input_cartesian

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: before_control, after_control
  PUBLIC :: grid
  PUBLIC :: open_files, close_files, restart_data
  PUBLIC :: normalise_code

  REAL(num), DIMENSION(:), ALLOCATABLE:: dxnew, dynew, dznew

CONTAINS


  SUBROUTINE before_control
    !Setup basic variables which have to have default values

    nprocx = 0
    nprocy = 0
    nprocz = 0

    time = 0.0_num
    gamma = 5.0_num/3.0_num

    IF (num .EQ. 4) mpireal = MPI_REAL

  END SUBROUTINE before_control



  SUBROUTINE after_control
    !Setup arrays and other variables which can only be set after
    !user input

    IF (IAND(initial,IC_RESTART) .EQ. 0)  restart_snapshot = 0

    p_visc = 0.0_num
    eta = 0.0_num
    grav = 0.0_num
    lambda_i = 0.0_num

    rho = 0.0_num
    energy = 0.0_num
    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num
    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num

  END SUBROUTINE after_control



  SUBROUTINE normalise_code

    !Normalise physical parameters (eta, grav, visc3 etc.)
    CALL normalise_constants    ! setup.f90
    !Normalise the grid
    CALL normalise_grid         ! setup.f90
    !Normalise the actual initial conditions
    CALL normalise_eqm          ! setup.f90
    !Normalise the constants needed for the
    !partially ionised plasma routines
    IF (include_neutrals) CALL normalise_neutral      ! setup.f90

  END SUBROUTINE normalise_code



  SUBROUTINE normalise_grid

    !Normalise the grid, including control volumes

    xb = xb / L0
    xc = xc / L0
    xb_global = xb_global / L0
    dxc = dxc / L0
    dxb = dxb / L0

    yb = yb / L0
    yc = yc / L0
    yb_global = yb_global / L0
    dyc = dyc / L0
    dyb = dyb / L0

    zb = zb / L0
    zc = zc / L0
    zb_global = zb_global / L0
    dzc = dzc / L0
    dzb = dzb / L0

    cv = cv / L0**3

  END SUBROUTINE normalise_grid



  SUBROUTINE normalise_eqm

    !Normalise the initial conditions

    rho = rho / RHO0
    energy = energy / ENERGY0
    vx = vx / VEL0
    vy = vy / VEL0
    vz = vz / VEL0
    bx = bx / B0
    by = by / B0
    bz = bz / B0

  END SUBROUTINE normalise_eqm



  SUBROUTINE normalise_constants

    !Normalise gravity etc.

    grav = grav / GRAV0
    visc3 = visc3 / VISC0
    eta0 = eta0 / RES0
    eta_background = eta_background / RES0
    
    kappa = 1.0e-11_num  ! the SI value for the constant in the conductivity assuming ln(Lambda)=18.4
    kappa = kappa / KAPPA0

    time = time / T0
    t_end = t_end / T0
    dt_snapshots = dt_snapshots / T0
    lambda_i = lambda_i / L0

  END SUBROUTINE normalise_constants



  SUBROUTINE Normalise_Neutral

    !Normalise constants used in the calculation of properties
    !Of partially ionised plasmas

    !Normalised mass
    REAL(num) :: MASS0

    REAL(num) :: Tr !Temperature of photospheric radiation field
    REAL(num) :: eta_bar_0

    MASS0 = RHO0 * L0**2

    f_bar=f_bar * L0**2 * TEMP0**(3.0_num/2.0_num)

    !Normalise tbar
    t_bar=t_bar / TEMP0

    !Normalise rbar
    r_bar=r_bar * MASS0

    !Normalise eta_bar
    eta_bar_0 = 1.0_num*(RHO0**2 * SQRT(TEMP0) * RES0/B0**2)
    eta_bar=eta_bar / eta_bar_0

    !Finally normalise ion_mass and ionise_pot which are needed in the code
    ionise_pot=ionise_pot/(ENERGY0 * MASS0)

    Tr=7230.85_num/TEMP0
    Tr_bar=1.0_num/Tr

  END SUBROUTINE Normalise_Neutral


 SUBROUTINE grid                 ! stretched and staggered grid

    REAL(num) :: dx, dy, dz, dxmin, dzmin, dymin, xcstar, ycstar, zcstar
    INTEGER :: ix, iy

    ALLOCATE(xb_global(-2:nx_global+2), dxnew(-2:nx_global+2))
    ALLOCATE(yb_global(-2:ny_global+2), dynew(-2:ny_global+2))
    ALLOCATE(zb_global(-2:nz_global+2), dznew(-2:nz_global+2))

    dx = 1.0_num / REAL(nx_global,num)       ! initially assume uniform grid
    dy = 1.0_num / REAL(ny_global,num)
    dz = 1.0_num / REAL(nz_global,num)

    length_x=x_end-x_start
    length_y=y_end-y_start
    length_z=z_end-z_start

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

    DEALLOCATE(dxnew, dynew, dznew)

  END SUBROUTINE grid




  !Subroutine stretches the grid in the x direction
  SUBROUTINE stretch_x   ! replace with any stretching algorithm as needed

    REAL(num) :: width, dx, L, f, lx_new

    lx_new = 200.0_num                ! new tolal length
    L = length_x / 1.5_num       ! centre of tanh stretching in unstretched coordinates
    width = length_x / 10.0_num       ! width of tanh stretching in unstretched coordinates

    f = (lx_new - length_x)/(length_x - L)/2.0_num

    dx = length_x / REAL(nx_global,num)  
    dxnew = dx + f*(1.0_num+TANH((ABS(xb_global)-L)/width))*dx

!!$    DO ix = nx_global/2+1, nx_global+2
!!$       xb_global(ix) = xb_global(ix-1) + dxnew(ix)
!!$    ENDDO
!!$    DO ix = nx_global/2-1, -2, -1
!!$       xb_global(ix) = xb_global(ix+1) - dxnew(ix)
!!$    ENDDO

    DO ix = 1, nx_global+2
       xb_global(ix) = xb_global(ix-1) + dxnew(ix)
    ENDDO

    length_x = lx_new

  END SUBROUTINE stretch_x


  !Subroutine stretches the domain in the y direction
  SUBROUTINE stretch_y !stretch domain upwards only

    REAL(num) :: width, dy, L, f, ly_new

    ly_new = 100.0_num                ! new tolal length
    L = length_y / 1.5_num       ! centre of tanh stretching in unstretched coordinates
    width = length_y / 10.0_num       ! width of tanh stretching in unstretched coordinates

    f = (ly_new - length_y)/(length_y - L)/2.0_num

    dy = length_y / REAL(ny_global,num)  
    dynew = dy + f*(1.0_num+TANH((ABS(yb_global)-L)/width))*dy

!!$    DO iy = ny_global/2+1, ny_global+2
!!$       yb_global(iy) = yb_global(iy-1) + dynew(iy)
!!$    ENDDO
!!$    DO iy = ny_global/2-1, -2, -1
!!$       yb_global(iy) = yb_global(iy+1) - dynew(iy)
!!$    ENDDO

    DO iy = 1, ny_global+2
       yb_global(iy) = yb_global(iy-1) + dynew(iy)
    ENDDO

    length_y = ly_new

  END SUBROUTINE stretch_y

  !Subroutine stretches the domain in the z direction
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


  !Open the output diagnostic files
  SUBROUTINE open_files

    CHARACTER(LEN=11+Data_Dir_Max_Length) :: file2
    CHARACTER(LEN=7+Data_Dir_Max_Length) :: file3
    INTEGER :: ios

    IF (rank == 0) THEN
       WRITE(file2, '(a,"/lare3d.dat")') TRIM(data_dir)
       OPEN(unit=20, STATUS = 'REPLACE',FILE = file2,iostat=ios)
       IF (ios .NE. 0) THEN
          PRINT *,"Unable to open file lare3d.dat for writing. This is most commonly caused by the output directory not existing"
          PRINT *," "
          PRINT *," "
          CALL MPI_ABORT(comm,errcode)
       ENDIF
       WRITE(file3, '(a,"/en.dat")') TRIM(data_dir)
       OPEN(unit=30, STATUS = 'REPLACE',FILE = file3,FORM="binary",iostat=ios)
       IF (ios .NE. 0) THEN
          PRINT *,"Unable to open file en.dat for writing. This is most commonly caused by the output directory not existing"
          PRINT *," "
          PRINT *," "
          CALL MPI_ABORT(comm,errcode)
       ENDIF
    END IF

  END SUBROUTINE open_files


  !Close the output diagnostic files
  SUBROUTINE close_files

    IF (rank == 0) THEN
       CLOSE(unit=20)
       CLOSE(unit=30)
    END IF

  END SUBROUTINE close_files


  !Subroutine to perform string comparisons
  FUNCTION StrCmp(StrIn,StrTest)

    CHARACTER(*),INTENT(IN) ::  StrIn,StrTest
    CHARACTER(30) :: StrTrim
    LOGICAL :: StrCmp

    StrTrim=TRIM(ADJUSTL(StrIn))

    IF (LEN(StrTest) .GT. LEN(StrIn)) THEN
       StrCmp=.FALSE.
       return
    ENDIF

    IF (StrTrim(LEN(StrTest)+1:LEN(StrTest)+1) .NE. " ") THEN
       StrCmp=.FALSE.
       RETURN
    ENDIF

    StrCmp=StrTrim(1:Len(StrTest)) == StrTest

  END FUNCTION StrCmp


  !Restart from previous output dumps
  SUBROUTINE restart_data
    CHARACTER(LEN=20+Data_Dir_Max_Length) :: filename
    CHARACTER(LEN=20) :: Name,Class,MeshName,MeshClass
    INTEGER :: nBlocks, Type, nd, sof, snap
    INTEGER,DIMENSION(3) :: dims
    REAL(num), DIMENSION(3) :: extent
    REAL(num), DIMENSION(3) :: stagger
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: data

    ! Create the filename for the last snapshot
#ifdef MHDCLUSTER
    WRITE(filename, '("nfs:",a,"/",i4.4,".cfd")') TRIM(data_dir), restart_snapshot
#else
    WRITE(filename, '(a,"/",i4.4,".cfd")') TRIM(data_dir), restart_snapshot
#endif

    output_file=restart_snapshot

    ALLOCATE(Data(0:nx,0:ny,0:nz))
    CALL cfd_Open(filename,rank,comm,MPI_MODE_RDONLY)
    ! Open the file
    nBlocks=cfd_Get_nBlocks()

    DO ix=1,nBlocks
       CALL cfd_Get_Next_Block_Info_All(Name,Class,Type)
       IF (rank == 0) PRINT *,ix,Name,Class,Type
       IF (Type == TYPE_SNAPSHOT) THEN
          CALL cfd_Get_Snapshot(time,snap)
       ENDIF
       IF (Type == TYPE_MESH) THEN
          !Strangely, LARE doesn't actually read in the grid from a file
          !This can be fixed, but for the moment, just go with the flow and
          !Replicate the old behaviour

          CALL cfd_Skip_Block()
       ELSE IF (Type == TYPE_MESH_VARIABLE) THEN
	  CALL cfd_Get_Common_MeshType_MetaData_All(Type,nd,sof)
	  IF (nd /= DIMENSION_3D) THEN
	     IF (rank == 0) PRINT *,"Non 3D Dataset found in input file, ignoring and continuting."
             CALL cfd_Skip_Block()
             CYCLE
	  ENDIF
	  IF (Type /= VAR_CARTESIAN) THEN
	     IF (rank == 0) PRINT *,"Non-Cartesian variable block found in file, ignoring and continuing"
             CALL cfd_Skip_Block()
             CYCLE
          ENDIF
	  !We now have a valid variable, let's load it up
	  !First error trapping
	  CALL  cfd_Get_nD_Cartesian_Variable_MetaData_All(nd,dims,extent,stagger,MeshName,MeshClass)

	  IF (dims(1) /= nx_global + 1 .OR. dims(2) /= ny_global + 1 .OR. dims(3) /= nz_global) THEN
	     IF (rank == 0) PRINT *,"Size of grid represented by one more variables invalid. Continuing"
	     CALL cfd_Skip_Block
	     CYCLE
	  ENDIF

	  IF (sof /= num) THEN
	     IF (rank == 0) PRINT *,"Precision of data does not match precision of code. Continuing."
	     CALL cfd_Skip_Block
	  ENDIF

   !We're not interested in the other parameters, so if we're here,
   !load up the data

	  CALL cfd_Get_3D_Cartesian_Variable_Parallel(Data,subtype)

   !Now have the data, just copy it to correct place

	  IF (StrCmp(Name(1:3),"Rho")) THEN
	     rho(0:nx,0:ny,0:nz)=data
	  ENDIF

	  IF (StrCmp(Name(1:6),"Energy")) THEN
             energy(0:nx,0:ny,0:nz)=data
	  ENDIF

	  IF (StrCmp(Name(1:2),"Vx")) THEN
             Vx(0:nx,0:ny,0:nz)=data
	  ENDIF

	  IF (StrCmp(Name(1:2),"Vy")) THEN
             Vy(0:nx,0:ny,0:nz)=data
	  ENDIF

	  IF (StrCmp(Name(1:2),"Vz")) THEN
             Vz(0:nx,0:ny,0:nz)=data
	  ENDIF

	  IF (StrCmp(Name(1:2),"Bx")) THEN
             Bx(0:nx,0:ny,0:nz)=data
	  ENDIF

	  IF (StrCmp(Name(1:2),"By")) THEN
             By(0:nx,0:ny,0:nz)=data
	  ENDIF

	  IF (StrCmp(Name(1:2),"Bz")) THEN
             Bz(0:nx,0:ny,0:nz)=data
	  ENDIF

   !Should be at end of block, but force the point anyway
          CALL cfd_Skip_Block()
       ELSE
          !Unknown block, just skip it
          CALL cfd_Skip_Block()
       ENDIF
    ENDDO

    DEALLOCATE(Data)

    CALL cfd_Close()

    CALL MPI_BARRIER(comm,errcode)
  END SUBROUTINE restart_data

END MODULE setup
