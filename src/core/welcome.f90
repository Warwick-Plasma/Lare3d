MODULE welcome

  USE shared_data

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: welcome_message

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

    WRITE(*,*)
    WRITE(*,'(''Welcome to Lare2D Version '', I1, ''.'', I1)') version, revision
    WRITE(*,*)

  END SUBROUTINE welcome_message

END MODULE welcome
