!****************************************************************
! All global variables defined here (cf F77 COMMON block).
! All the names in here are public provided the MODULE is USE'd
!**************************************************************** 
MODULE constants
  IMPLICIT NONE
  INTEGER, PARAMETER :: num = KIND(1.D0)
  INTEGER, PARAMETER :: dbl = KIND(1.D0)
  REAL(num), PARAMETER :: pi = 3.14159265358979323_num
  REAL(num), PARAMETER :: none_zero = TINY(1.0_num)
  REAL(num), PARAMETER :: largest_number = HUGE(1.0_num)
  INTEGER, PARAMETER :: periodic = 1, other = 2
  INTEGER, PARAMETER :: open = 3
  INTEGER, PARAMETER :: Version = 2, Revision = 0

  INTEGER,PARAMETER ::ERR_NONE=0,ERR_UNKNOWN_BLOCK=1,ERR_UNKNOWN_ELEMENT=2
  INTEGER,PARAMETER ::ERR_PRESET_ELEMENT=4,ERR_PRESET_ELEMENT_USE_LATER=8
  INTEGER,PARAMETER ::ERR_BAD_VALUE=16,ERR_OTHER=1024
END MODULE constants 


MODULE shared_data
  USE constants


!  USE mpi
  IMPLICIT NONE
  INCLUDE 'mpif.h'

TYPE :: Entry
   character(30) :: Value
END TYPE Entry


  INTEGER :: mpireal = MPI_DOUBLE_PRECISION

  INTEGER :: nx_global,ny_global,nz_global
! NB: as there are now 2 ghost celss so indexing will fail if (nx,ny,nz)<MAX(2,offset)
  INTEGER :: nx, ny, nz

  INTEGER  :: nsteps 
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: rho, energy
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: bx, vx, vx1
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: by, vy, vy1  
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: bz, vz, vz1 

  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: delta_ke, p_visc
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: eta, eta_perp, curlb
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: cv, cv1, bzone

  REAL(num), DIMENSION(:,:), ALLOCATABLE :: dv_right, dv_left, dv_up, dv_down
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: dv_back, dv_front

  INTEGER,PARAMETER :: Data_Dir_Max_Length = 64
  CHARACTER(LEN=Data_Dir_Max_Length) :: data_dir

  REAL(num), DIMENSION(:), ALLOCATABLE :: xc, xb, grav
  REAL(num), DIMENSION(:), ALLOCATABLE :: dxb, dxc
  REAL(num), DIMENSION(:), ALLOCATABLE :: yc, yb
  REAL(num), DIMENSION(:), ALLOCATABLE :: dyb, dyc
  REAL(num), DIMENSION(:), ALLOCATABLE :: zc, zb
  REAL(num), DIMENSION(:), ALLOCATABLE :: dzb, dzc
  REAL(num), DIMENSION(:), ALLOCATABLE :: xb_global, yb_global, zb_global

  REAL(num) :: w1, w2, w3, w4, w5, w6, w7, w8
  REAL(num)  :: dt, dt2, dt_multiplier, t_end, time, min_grid_spacing
  REAL(num) :: length_x, length_y, length_z, visc1, visc2, visc3, visc4
  REAL(num) :: x_start,x_end,y_start,y_end,z_start,z_end
  REAL(num) :: gamma, eta0, j_max, dt_snapshots
  REAL(num) :: dtr, eta_background, ion_mass
  REAL(num) :: total_visc_heating =0.0_num, total_ohmic_heating=0.0_num

  INTEGER :: xbc_left, xbc_right, ybc_down, ybc_up, zbc_front, zbc_back
  INTEGER :: ix, iy, iz, ixp, iyp, izp, ixm, iym, izm, xpass, ypass, zpass
  INTEGER :: restart_snapshot
  INTEGER :: peak_substeps = 0
  LOGICAL :: x_stretch, y_stretch, z_stretch, jxB, rke, tensor_shock_visc
  LOGICAL :: resistiveMHD, any_open, restart
  LOGICAL :: include_neutrals,farfield,damping,deckfile

! MPI data

  INTEGER :: rank, left, right, up, down, back, front, coordinates(3)
  INTEGER :: errcode, comm, tag, nproc, nprocx, nprocy, nprocz
  INTEGER :: status(MPI_STATUS_SIZE)
  REAL(num) :: walltime_end,walltime_start

! file handling
  INTEGER :: subtype, obstype
  INTEGER(KIND=MPI_OFFSET_KIND) :: initialdisp 

!Problem Dependant input deck demo
  INTEGER :: Problem
  REAL(num) :: Bx_g,By_g,Bz_g
  REAL(num) :: p_factor,p_dropoff

END MODULE shared_data

! The pre-processor removes the following line so it compiles without error
! unless the pre-processor hasn't been run over it

#ifdef PRECHECK
This line deliberately breaks the compile if the preprocessor has not worked.
#endif
