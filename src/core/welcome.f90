!  Copyright 2020 University of Warwick

!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at

!      http://www.apache.org/licenses/LICENSE-2.0

!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.

!******************************************************************************
! Welcome message routines
!******************************************************************************

MODULE welcome

  USE version_data
  USE shared_data

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: welcome_message, create_ascii_header

CONTAINS

  !****************************************************************************
  ! This routine prints the welcome message, MPI status
  !****************************************************************************

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

    CHARACTER(LEN=16) :: job1, job2
    CHARACTER(LEN=4) :: str
    INTEGER :: i, strmin, strmax, strlen

    ! Parse commit string to get version number
    ! Commit ID begins with the string v[0-9].[0-9].[0-9]-
    strlen = LEN_TRIM(c_commit_id)
    strmin = 1
    strmax = strmin + 12
    c_revision = -1

    IF (c_commit_id(1:1) == 'v') THEN
      strmin = 2
    END IF

    ! Version
    DO i = strmin, MIN(strmax,strlen)
      IF (c_commit_id(i:i) == '.') THEN
        str = c_commit_id(strmin:i-1)
        READ(str, '(i9)') c_version
        strmin = i + 1
        strmax = strmin + 4
        EXIT
      END IF
    END DO

    ! Revision
    DO i = strmin, MIN(strmax,strlen)
      IF (c_commit_id(i:i) == '.') THEN
        str = c_commit_id(strmin:i-1)
        READ(str, '(i9)') c_revision
        strmin = i + 1
        strmax = strmin + 4
        EXIT
      END IF
    END DO

    ! Minor revision
    DO i = strmin, MIN(strmax,strlen)
      IF (c_commit_id(i:i) == '-') THEN
        str = c_commit_id(strmin:i-1)
        READ(str, '(i9)') c_minor_rev
        strmax = i - 1
        EXIT
      END IF
    END DO

    IF (c_revision < 0) THEN
      c_revision = c_minor_rev
      c_minor_rev = 0
    END IF

    version_string = c_commit_id(1:strmax)

    CALL integer_as_string(jobid%start_seconds, job1)
    CALL integer_as_string(jobid%start_milliseconds, job2)
    ascii_header = c_code_name // ' v' // TRIM(version_string) // '   ' &
        // c_commit_id // ' ' // TRIM(job1) // '.' // TRIM(ADJUSTL(job2))

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
