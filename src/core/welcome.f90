MODULE welcome

  USE version_data
  USE shared_data

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: welcome_message, create_ascii_header

CONTAINS

  SUBROUTINE welcome_message

    CHARACTER, DIMENSION(4) :: clrstr = (/ ' ', '[', '2', 'J' /)

    IF (rank /= 0) RETURN

    clrstr(1) = CHAR(27)
    WRITE(*,'(1x, 4a1)') clrstr

    PRINT*,'    @@          @@@@    @@@@@@@@  @@@@@@@@    @@@@@@@@  @@@@@@   '
    PRINT*,'    @@          @@@@    @@@@@@@@  @@@@@@@@    @@@@@@@@  @@@@@@   '
    PRINT*,'    @@        @@@@@@@@  @@    @@  @@                @@  @@  @@   '
    PRINT*,'    @@        @@@@@@@@  @@    @@  @@                @@  @@  @@   '
    PRINT*,'    @@        @@    @@  @@    @@  @@                @@  @@    @@ '
    PRINT*,'    @@        @@    @@  @@    @@  @@                @@  @@    @@ '
    PRINT*,'    @@        @@@@@@@@  @@@@@@@@  @@@@@@@@    @@@@@@@@  @@    @@ '
    PRINT*,'    @@        @@@@@@@@  @@@@@@@@  @@@@@@@@    @@@@@@@@  @@    @@ '
    PRINT*,'    @@        @@    @@  @@@@      @@                @@  @@    @@ '
    PRINT*,'    @@        @@    @@  @@@@      @@                @@  @@    @@ '
    PRINT*,'    @@        @@    @@  @@  @@    @@                @@  @@  @@   '
    PRINT*,'    @@        @@    @@  @@  @@    @@                @@  @@  @@   '
    PRINT*,'    @@@@@@@@  @@    @@  @@    @@  @@@@@@@@    @@@@@@@@  @@@@@@   '
    PRINT*,'    @@@@@@@@  @@    @@  @@    @@  @@@@@@@@    @@@@@@@@  @@@@@@   '
    PRINT*
    PRINT*

    CALL create_ascii_header

    WRITE(*,*)
    WRITE(*,*) 'Welcome to ', TRIM(c_code_name), ' version ', &
        TRIM(version_string) // '   (commit ' // TRIM(c_commit_id) // ')'
    WRITE(*,*)

    CALL mpi_status_message

  END SUBROUTINE welcome_message



  SUBROUTINE mpi_status_message

    CHARACTER(LEN=8) :: string

    CALL integer_as_string(nproc, string)

    WRITE(*,*) 'Code is running on ', TRIM(string), ' processing elements'
    WRITE(*,*)

  END SUBROUTINE mpi_status_message



  SUBROUTINE create_ascii_header

    CHARACTER(LEN=11) :: ver, rev, minor_rev

    CALL integer_as_string(c_version, ver)
    CALL integer_as_string(c_revision, rev)
    CALL integer_as_string(c_minor_rev, minor_rev)
    version_string = TRIM(ver) // '.' // TRIM(ADJUSTL(rev)) // '.' &
      // TRIM(ADJUSTL(minor_rev))
    !CALL integer_as_string(jobid%start_seconds, ver)
    !CALL integer_as_string(jobid%start_milliseconds, rev)
    ascii_header = c_code_name // ' v' // TRIM(version_string) // ' ' &
        // c_commit_id // ' ' // TRIM(ver) // '.' // TRIM(ADJUSTL(rev))

  END SUBROUTINE create_ascii_header



  SUBROUTINE integer_as_string(int_in, string)

    INTEGER, INTENT(IN) :: int_in
    CHARACTER(LEN=*), INTENT(OUT) :: string

    INTEGER :: n_nums
    CHARACTER(LEN=9) :: numfmt

    IF (int_in == 0) THEN
      n_nums = 1
    ELSE
      n_nums = 1 + INT(LOG10(REAL(ABS(int_in), num)))
    END IF
    WRITE(numfmt, '(''(I'', I6.6, '')'')') n_nums
    WRITE(string, numfmt) int_in

  END SUBROUTINE integer_as_string

END MODULE welcome
