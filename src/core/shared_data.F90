!******************************************************************************
! All global variables defined here (cf F77 COMMON block).
! All the names in here are public provided the MODULE is USE'd
!******************************************************************************

MODULE constants

  IMPLICIT NONE

#ifdef SINGLE
  INTEGER, PARAMETER :: num = KIND(1.0)
  INTEGER, PARAMETER :: num_sz = 4
#else
  INTEGER, PARAMETER :: num = KIND(1.D0)
  INTEGER, PARAMETER :: num_sz = 8
#endif
  INTEGER, PARAMETER :: dbl = KIND(1.D0)

  ! Code dimensions
  INTEGER, PARAMETER :: c_ndims = 3

  REAL(num), PARAMETER :: pi = 3.141592653589793238462643383279503_num

  ! These are the real SI physical constants
  ! Taken from http://physics.nist.gov/cuu/Constants (05/07/2012)

  ! Permeability of free space
  REAL(num), PARAMETER :: mu0_si =  4.0e-7_num * pi ! N/A^2 (exact)

  ! Boltzmann's Constant
  REAL(num), PARAMETER :: kb_si = 1.3806488e-23_num  ! J/K (+/- 1.3e-29)

  ! Mass of hydrogen ion
  REAL(num), PARAMETER :: mh_si = 1.672621777e-27_num ! kg (+/- 7.4e-35)

  ! Mass of electron
  REAL(num), PARAMETER :: me_si = 9.10938291e-31_num ! kg (+/- 4e-38)

  ! Planck's constant
  REAL(num), PARAMETER :: hp_si = 6.62606957e-34_num ! J s (+/- 2.9e-41)

  ! Ionisation potential of hydrogen in J
  REAL(num), PARAMETER :: ionise_pot_si = 2.17870364e-18_num


  REAL(num), PARAMETER :: dt_multiplier = 0.9_num

  REAL(num), PARAMETER :: none_zero = TINY(1.0_num)
  REAL(num), PARAMETER :: largest_number = HUGE(1.0_num)
  REAL(num), PARAMETER :: third = 1.0_num / 3.0_num, sixth = 1.0_num / 6.0_num
  INTEGER, PARAMETER :: BC_PERIODIC = 1, BC_OTHER = 2
  INTEGER, PARAMETER :: BC_OPEN = 3

  ! IC codes
  ! This is a bitmask, remember that
  INTEGER, PARAMETER :: IC_NEW = 1, IC_RESTART = 2

  ! Equation of state codes
  INTEGER, PARAMETER :: EOS_IDEAL = 1, EOS_ION = 2, EOS_PI = 3

END MODULE constants



MODULE shared_data

  USE constants
  USE sdf_job_info
  USE sdf

  IMPLICIT NONE

  INCLUDE 'mpif.h'

#ifdef SINGLE
  INTEGER, PARAMETER :: mpireal = MPI_REAL
  INTEGER, PARAMETER :: sdf_num = c_datatype_real4
#else
  INTEGER, PARAMETER :: mpireal = MPI_DOUBLE_PRECISION
  INTEGER, PARAMETER :: sdf_num = c_datatype_real8
#endif

  INTEGER :: nx_global, ny_global, nz_global

  ! NB: as there are now 2 ghost cells, indexing will fail if (nx,ny,nz) < 2
  INTEGER :: nx, ny, nz
  INTEGER :: nsteps, step
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: rho, energy
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: bx, vx, vx1
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: by, vy, vy1
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: bz, vz, vz1
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: bx1, by1, bz1
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: jx_r, jy_r, jz_r

  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: delta_ke, p_visc
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: eta, cv, cv1

  REAL(num), DIMENSION(:), ALLOCATABLE :: xc, xb, dxb, dxc, xb_global
  REAL(num), DIMENSION(:), ALLOCATABLE :: yc, yb, dyb, dyc, yb_global
  REAL(num), DIMENSION(:), ALLOCATABLE :: zc, zb, dzb, dzc, zb_global
  REAL(num), DIMENSION(:), ALLOCATABLE :: grav

  INTEGER, PARAMETER :: data_dir_max_length = 64
  CHARACTER(LEN = data_dir_max_length) :: data_dir

  INTEGER, PARAMETER :: c_max_string_length = 128

  REAL(num) :: w1, w2, w3, w4, w5, w6, w7, w8, w9
  REAL(num) :: dt, dt2, dtr, dth, t_end, time
  REAL(num) :: dt_from_restart, time_prev
  REAL(num) :: visc1, visc2
  REAL(num) :: x_min, x_max, length_x
  REAL(num) :: y_min, y_max, length_y
  REAL(num) :: z_min, z_max, length_z
  REAL(num) :: gamma, eta0, j_max, dt_snapshots, eta_background
  REAL(num) :: total_visc_heating = 0.0_num, total_ohmic_heating = 0.0_num

  INTEGER :: xbc_min, xbc_max, ix, ixm, ixp, xpass
  INTEGER :: ybc_min, ybc_max, iy, iym, iyp, ypass
  INTEGER :: zbc_min, zbc_max, iz, izm, izp, zpass
  INTEGER :: restart_snapshot
  INTEGER :: peak_substeps = 0
  LOGICAL :: x_stretch, y_stretch, z_stretch
  LOGICAL :: resistive_mhd, any_open, rke
  LOGICAL :: restart

  ! Normalising constants
  REAL(num) :: B0, L0, rho0
  ! Mass fraction - mass of ions in units of proton mass
  REAL(num) :: mf
  ! Conversion factor to get temperature in MK from normalised energy
  REAL(num) :: t2tmk
  ! Normalisation used for radiative losses
  REAL(num) :: lr_star
  ! Normalisation used for coronal heating
  REAL(num) :: h_star

  ! Heat conduction
  LOGICAL :: conduction, heat_flux_limiter, radiation, coronal_heating
  REAL(num) :: kappa_0, flux_limiter, temperature_100mk

  ! Equation of state
  INTEGER :: eos_number = EOS_IDEAL

  ! Damping boundary variables
  LOGICAL :: damping

  ! Partially ionised plasma
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: eta_perp, xi_n, eta_perp0
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: parallel_current, perp_current
  LOGICAL :: cowling_resistivity, neutral_gas
  REAL(num) :: f_bar, t_bar, tr, ionise_pot, r_bar
  REAL(num) :: eta_bar

  ! MPI data
  INTEGER :: coordinates(c_ndims), n_global_min(c_ndims), n_global_max(c_ndims)
  INTEGER, DIMENSION(:), ALLOCATABLE :: cell_nx_mins, cell_nx_maxs
  INTEGER, DIMENSION(:), ALLOCATABLE :: cell_ny_mins, cell_ny_maxs
  INTEGER, DIMENSION(:), ALLOCATABLE :: cell_nz_mins, cell_nz_maxs
  INTEGER :: nprocx, proc_x_min, proc_x_max
  INTEGER :: nprocy, proc_y_min, proc_y_max
  INTEGER :: nprocz, proc_z_min, proc_z_max
  INTEGER :: rank, errcode, comm, tag, nproc
  INTEGER :: status(MPI_STATUS_SIZE)
  INTEGER :: cell_subarray, cellng_subarray, cell_distribution
  INTEGER :: node_subarray, nodeng_subarray, node_distribution
  INTEGER :: bx_subarray, bx_distribution
  INTEGER :: by_subarray, by_distribution
  INTEGER :: bz_subarray, bz_distribution
  INTEGER :: cell_xface, node_xface, node_xface1
  INTEGER :: cell_yface, node_yface, node_yface1
  INTEGER :: cell_zface, node_zface, node_zface1
  INTEGER :: bx_xface, by_xface, bz_xface, bx_xface1
  INTEGER :: bx_yface, by_yface, bz_yface, by_yface1
  INTEGER :: bx_zface, by_zface, bz_zface, bz_zface1

  ! File handling
  INTEGER :: subtype, obstype
  INTEGER :: initial
  INTEGER, PARAMETER :: n_zeros = 4
  INTEGER, PARAMETER :: en_nvars = 6
  INTEGER :: file_number = 0
#ifdef FILEPREFIX
  CHARACTER(LEN=4), PARAMETER :: filesystem = 'nfs:'
#else
  CHARACTER(LEN=1), PARAMETER :: filesystem = ''
#endif
  CHARACTER(LEN=6) :: file_prefix = ''
  TYPE(jobid_type) :: jobid
  INTEGER :: run_date = 0

  ! History file header
  CHARACTER(LEN=3) :: c_history_magic = 'HIS'
  INTEGER, PARAMETER :: c_history_version = 1
  INTEGER, PARAMETER :: c_history_revision = 0
  INTEGER, PARAMETER :: c_endianness = 16911887

  INTEGER, PARAMETER :: c_stagger_bx = c_stagger_face_x
  INTEGER, PARAMETER :: c_stagger_by = c_stagger_face_y
  INTEGER, PARAMETER :: c_stagger_bz = c_stagger_face_z

  ! Number of variables to dump
  LOGICAL, DIMENSION(20) :: dump_mask

  INTEGER, PARAMETER :: stat_unit = 20
  INTEGER, PARAMETER :: en_unit = 30

END MODULE shared_data



! The pre-processor removes the following line so it compiles without error
! unless the pre-processor hasn't been run over it

#ifdef PRECHECK
This line deliberately breaks the compile IF the preprocessor has not worked.
#endif
