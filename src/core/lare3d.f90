PROGRAM lare3d

  USE shared_data
  USE initial_conditions
  USE setup
  USE boundary
  USE openboundary
  USE diagnostics
  USE lagran
  USE remap
  USE mpi_routines
  USE welcome
  USE normalise
  USE eos
  USE neutral
  USE control

  IMPLICIT NONE

  INTEGER :: i = 0

  CALL MPI_INIT(errcode)

  CALL before_control      ! setup.F90
  CALL user_normalisation  ! control.f90
  CALL control_variables   ! control.f90
  CALL set_output_dumps    ! control.f90
  CALL mpi_initialise      ! mpi_routines.f90
  CALL after_control       ! setup.f90

  CALL welcome_message     ! welcome.f90

  CALL set_normalisation         ! normalise.f90
  CALL set_boundary_conditions   ! boundary.f90
  CALL open_files                ! setup.f90
  CALL grid                      ! setup.f90
  IF (include_neutrals .OR. &
      cowling_resistivity) CALL setup_neutral ! neutral.f90

  IF (IAND(initial, IC_RESTART) .NE. 0) THEN
    CALL restart_data            ! setup.f90
    restart = .TRUE.
  END IF

  IF (IAND(initial, IC_NEW) .NE. 0) THEN
    CALL set_initial_conditions  ! initial_conditions.f90
  END IF

  ! Initial conditions, parameters etc. specified in SI
  IF (SI) THEN
    ! Normalise everything
    CALL normalise_code          ! setup.f90
  END IF

  CALL set_boundary_conditions   ! boundary.f90
  CALL boundary_conditions       ! boundary.f90
  CALL eta_calc                  ! lagran.f90

  IF (include_neutrals) CALL neutral_fraction(eos_number) ! neutral.f90
  IF (cowling_resistivity) CALL perpendicular_resistivity ! neutral.f90

  IF (rank .EQ. 0) PRINT *, "Initial conditions setup OK. Running Code"

  CALL output_routines(i)        ! diagnostics.f90

  DO
    IF ((i >= nsteps .AND. nsteps >= 0) .OR. (time >= t_end)) EXIT
    i = i + 1
    CALL eta_calc                    ! lagran.f90
    CALL set_dt                      ! diagnostics.f90
    CALL lagrangian_step             ! lagran.f90
    CALL eulerian_remap(i)           ! remap.f90
    IF (rke) CALL energy_correction  ! diagnostics.f90
    IF (any_open) THEN
      CALL open_bcs                  ! openboundary.f90
    END IF
    CALL boundary_conditions         ! boundary.f90
    CALL output_routines(i)          ! diagnostics.f90
  END DO

  IF (rank .EQ. 0) PRINT *, "Code Terminated normally"
  CALL mpi_close                     ! mpi_routines.f90
  CALL close_files                   ! setup.f90
  CALL MPI_FINALIZE(errcode)

END PROGRAM lare3d
