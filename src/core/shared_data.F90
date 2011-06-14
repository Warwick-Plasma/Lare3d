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

  ! These are the real SI physical constants
  ! Permiability of free space
  REAL(num), PARAMETER :: mu0_si =  4.0e-7_num * pi
  
  ! Boltzmann's Constant
  REAL(num), PARAMETER :: kb_si = 1.3806504e-23_num
  
  ! Mass of hydrogen ion
  REAL(num), PARAMETER :: mh_si = 1.67262158e-27_num
  
  ! Mass of electron
  REAL(num), PARAMETER :: me_si = 9.10938188e-31_num
  
  ! Planck's constant
  REAL(num), PARAMETER :: hp_si = 6.626068e-34_num
  
  ! Ionisation potential of hydrogen in J
  REAL(num), PARAMETER :: ionise_pot_si = 2.17870364e-18_num
  

  REAL(num), PARAMETER :: dt_multiplier = 0.9_num

  REAL(num), PARAMETER :: none_zero = TINY(1.0_num) 
  REAL(num), PARAMETER :: largest_number = HUGE(1.0_num)  
  REAL(num), PARAMETER :: third = 1.0_num / 3.0_num, sixth = 1.0_num / 6.0_num
  INTEGER, PARAMETER :: BC_PERIODIC = 1, BC_OTHER = 2
  INTEGER, PARAMETER :: BC_OPEN = 3

  INTEGER, PARAMETER :: version = 2, revision = 8

  ! IC codes
  ! This is a bitmask, remember that
  INTEGER, PARAMETER :: IC_NEW = 1, IC_RESTART = 2

  ! Equation of state codes
  INTEGER, PARAMETER :: EOS_IDEAL = 1, EOS_ION = 2, EOS_PI = 3


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

  REAL(num) :: w1, w2, w3, w4, w5, w6, w7, w8, w9
  REAL(num) :: dt, dt2, dtr, dth, t_end, time
  REAL(num) :: length_x, length_y, length_z, visc1, visc2, visc3
  REAL(num) :: x_start, x_end, y_start, y_end, z_start, z_end
  REAL(num) :: gamma, eta0, j_max, dt_snapshots, lambda0, eta_background
  REAL(num) :: total_visc_heating = 0.0_num, total_ohmic_heating = 0.0_num

  INTEGER :: xbc_max, ybc_max, xbc_min, ybc_min, zbc_min, zbc_max
  INTEGER :: ix, iy, iz, ixp, iyp, izp, ixm, iym, izm, xpass, ypass, zpass
  INTEGER :: restart_snapshot
  INTEGER :: peak_substeps = 0
  LOGICAL :: x_stretch, y_stretch, z_stretch, rke
  LOGICAL :: resistive_mhd, any_open, hall_mhd
  LOGICAL :: restart

  ! normalising constants                 
  REAL(num) :: B0, L0, rho0
  ! mass fraction - mass of ions in units of proton mass
  REAL(num) :: mf  
  !convertion factor to get normalised temperature from normalised energy  
  REAL(num) :: e2t
  !convertion factor to get temperature in MK from normalised energy  
  REAL(num) :: e2tmk
  ! normalisation used for radiative losses
  REAL(num) :: lr_star
  ! normalisation used for coronal heating
  REAL(num) :: h_star
  
  ! Heat conduction
  LOGICAL :: conduction
  LOGICAL :: heat_flux_limiter
  REAL(num) :: kappa_0, flux_limiter, temperature_100mk  

  ! Equation of state
  INTEGER :: eos_number = EOS_IDEAL

  ! Damping boundary variables
  LOGICAL :: damping

  ! Partially ionised plasma
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: eta_perp, xi_n, eta_perp0
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: parallel_current, perp_current
  LOGICAL :: cowling_resistivity, neutral_gas
  REAL(num) :: f_bar, t_bar, tr, ionise_pot, r_bar
  REAL(num) :: eta_bar

  ! MPI data
  INTEGER :: rank, proc_x_min, proc_x_max, proc_y_min, proc_y_max, proc_z_min, proc_z_max, coordinates(3)
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
