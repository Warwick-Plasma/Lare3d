!****************************************************************
! All global variables defined here (cf F77 COMMON block).
! All the names in here are public provided the MODULE is USE'd
!****************************************************************
MODULE constants

  IMPLICIT NONE

#ifdef Q_SINGLE
  INTEGER, PARAMETER :: num = KIND(1.0) 
#else
  INTEGER, PARAMETER :: num = KIND(1.D0) 
#endif  
  INTEGER, PARAMETER :: dbl = KIND(1.D0)
  REAL(num), PARAMETER :: pi = 3.14159265358979323_num
  REAL(num), PARAMETER :: none_zero = TINY(1.0_num) 
  REAL(num), PARAMETER :: largest_number = HUGE(1.0_num)
  INTEGER, PARAMETER :: BC_PERIODIC = 1, BC_OTHER = 2
  INTEGER, PARAMETER :: BC_OPEN = 3

  INTEGER, PARAMETER :: version = 2, revision = 3

  ! IC codes
  ! This is a bitmask, remember that
  INTEGER, PARAMETER :: IC_NEW = 1, IC_RESTART = 2

  ! Equation of state codes
  INTEGER, PARAMETER :: EOS_IDEAL = 1, EOS_ION = 2, EOS_PI = 4

END MODULE constants



MODULE shared_data

  USE constants
  IMPLICIT NONE
  INCLUDE 'mpif.h'

#ifdef Q_SINGLE
  INTEGER :: mpireal = MPI_REAL
#else
  INTEGER :: mpireal = MPI_DOUBLE_PRECISION
#endif  

  INTEGER :: nx_global, ny_global, nz_global
  ! NB: as there are now 2 ghost celss so indexing will fail if (nx, ny, nz)<2
  INTEGER :: nx, ny, nz
  INTEGER  :: nsteps
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: rho, energy
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: bx, vx, vx1
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: by, vy, vy1
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: bz, vz, vz1
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: bx1, by1, bz1
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: jx_r, jy_r, jz_r

  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: delta_ke, p_visc
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: eta, lambda_i
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: cv, cv1

  REAL(num), DIMENSION(:), ALLOCATABLE :: xc, xb
  REAL(num), DIMENSION(:), ALLOCATABLE :: dxb, dxc
  REAL(num), DIMENSION(:), ALLOCATABLE :: xb_global, yb_global, zb_global
  REAL(num), DIMENSION(:), ALLOCATABLE :: yc, yb
  REAL(num), DIMENSION(:), ALLOCATABLE :: dyb, dyc, grav
  REAL(num), DIMENSION(:), ALLOCATABLE :: zc, zb
  REAL(num), DIMENSION(:), ALLOCATABLE :: dzb, dzc

  INTEGER, PARAMETER :: data_dir_max_length = 64
  CHARACTER(LEN = data_dir_max_length) :: data_dir

  REAL(num) :: w1, w2, w3, w4, w5, w6, w7, w8
  REAL(num) :: dt, dt2, dtr, dth, t_end, time
  REAL(num) :: dt_multiplier = 0.8_num
  REAL(num) :: length_x, length_y, length_z, visc1, visc2, visc3
  REAL(num) :: x_start, x_end, y_start, y_end, z_start, z_end
  REAL(num) :: gamma, eta0, j_max, dt_snapshots, lambda0, eta_background
  REAL(num) :: total_visc_heating = 0.0_num, total_ohmic_heating = 0.0_num

  INTEGER :: xbc_right, ybc_up, xbc_left, ybc_down, zbc_front, zbc_back
  INTEGER :: ix, iy, iz, ixp, iyp, izp, ixm, iym, izm, xpass, ypass, zpass
  INTEGER :: restart_snapshot
  INTEGER :: peak_substeps = 0
  LOGICAL :: x_stretch, y_stretch, z_stretch, rke
  LOGICAL :: resistive_mhd, any_open, hall_mhd
  LOGICAL :: restart

  ! Heat conduction
  LOGICAL :: conduction
  LOGICAL :: heat_flux_limiter
  REAL(num) :: kappa_0, flux_limiter, temperature_100mk  

  ! RTV radiation
  LOGICAL :: radiation

  ! Code normalisation
  LOGICAL :: SI

  ! Driving
  REAL(num) :: drv_amp = 0.0_num, omega = 2.0_num * pi / 1.0_num

  ! Equation of state
  INTEGER :: eos_number = EOS_IDEAL

  ! Damping boundary variables
  LOGICAL :: damping

  ! Partially ionised plasma
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: eta_perp, xi_n, eta_perp0
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: parallel_current, perp_current
  LOGICAL :: include_neutrals, cowling_resistivity
  REAL(num) :: f_bar, t_bar, tr, ionise_pot, r_bar
  REAL(num) :: eta_bar

  ! MPI data
  INTEGER :: rank, left, right, up, down, front, back, coordinates(3)
  INTEGER :: errcode, comm, tag, nproc, nprocx, nprocy, nprocz
  INTEGER :: status(MPI_STATUS_SIZE)

  ! file handling
  INTEGER :: subtype, obstype
  INTEGER(KIND = MPI_OFFSET_KIND) :: initialdisp
  INTEGER :: initial
  INTEGER :: n_zeros = 4
  INTEGER :: output_file = 0

  ! Number of variables to dump
  LOGICAL, DIMENSION(19) :: dump_mask

END MODULE shared_data



! The pre-processor removes the following line so it compiles without error
! unless the pre-processor hasn't been run over it

#ifdef PRECHECK
This line deliberately breaks the compile IF the preprocessor has not worked.
#endif
