!*******************************************************************
! All the ghost cell values are controlled by these routines.
! To speed things up it may be worth having this routine hard coded
! for each particular run, i.e. remove all if statements.
!*******************************************************************

MODULE boundary

  USE shared_data
  USE mpiboundary

  IMPLICIT NONE

  INTEGER :: ndx, ndy, ndz     

CONTAINS

  SUBROUTINE set_boundary_conditions

    REAL(num) :: a
    LOGICAL :: first_call = .TRUE.

    IF (first_call) THEN
      any_open = .FALSE.
      first_call = .FALSE.
      IF ((xbc_right == BC_OPEN) .OR. (xbc_left == BC_OPEN) &
          .OR. (ybc_up == BC_OPEN) .OR. (ybc_down == BC_OPEN) &
          .OR. (zbc_front == BC_OPEN) .OR. (zbc_back == BC_OPEN)) &
          any_open = .TRUE.
    ELSE
      ! when bzone = 0 uses first order remap scheme so that farfield not
      ! used in remap
      bzone = 1.0_num

      IF (xbc_right == BC_OPEN .AND. right == MPI_PROC_NULL) &
          bzone(nx-4:nx+2, :, :) = 0.0_num
      IF (xbc_left == BC_OPEN .AND. left == MPI_PROC_NULL) &
          bzone(-1:4, :, :) = 0.0_num
      IF (ybc_up == BC_OPEN .AND. up == MPI_PROC_NULL) &
          bzone(:, ny-4:ny+2, :) = 0.0_num
      IF (ybc_down == BC_OPEN .AND. down == MPI_PROC_NULL) &
          bzone(:, -1:4, :) = 0.0_num
      IF (zbc_back == BC_OPEN .AND. back == MPI_PROC_NULL) &
          bzone(:, :, nz-4:nz+1) = 0.0_num
      IF (zbc_front == BC_OPEN .AND. front == MPI_PROC_NULL) &
          bzone(:, :, -1:4) = 0.0_num

      ! define region for linear damping of solution at boundaries
      ! note if ndx > nx damping only in outer process
      a = length_x / 10.0_num  
      DO ix = 0, nx
        ndx = ix
        IF ((xb(ix) - xb(0)) > a) EXIT
      END DO

      a = length_y / 10.0_num  
      DO iy = 0, ny
        ndy = iy
        IF ((yb(iy) - yb(0)) > a) EXIT
      END DO

      a = length_z / 10.0_num
      DO iz = 0, nz
        ndz = iz
        IF ((zb(iz) - zb(0)) > a) EXIT
      END DO
   
    END IF

  END SUBROUTINE set_boundary_conditions



  SUBROUTINE boundary_conditions

    CALL damp_boundaries
    CALL bfield_bcs
    CALL energy_bcs
    CALL density_bcs
    CALL velocity_bcs

  END SUBROUTINE boundary_conditions



  SUBROUTINE damp_boundaries
  
    REAL(num) :: a
  
    IF (damping) THEN
  
      IF (right == MPI_PROC_NULL) THEN
        DO iz = -2, nz+2
          DO iy = -2, ny+2
            DO ix = nx-ndx, nx+2
              a = dt * REAL(ix - (nx-ndx), num) / REAL(ndx, num)
              vx(ix, iy, iz) = vx(ix, iy, iz) / (1.0_num + a)
              vy(ix, iy, iz) = vy(ix, iy, iz) / (1.0_num + a)
              vz(ix, iy, iz) = vz(ix, iy, iz) / (1.0_num + a)
            END DO
          END DO
        END DO
      END IF
        
      IF (left == MPI_PROC_NULL) THEN
        DO iz = -2, nz+2
          DO iy = -2, ny+2
            DO ix = -2, ndx
              a = dt * REAL((ndx - ix), num) / REAL(ndx, num)
              vx(ix, iy, iz) = vx(ix, iy, iz) / (1.0_num + a)
              vy(ix, iy, iz) = vy(ix, iy, iz) / (1.0_num + a)
              vz(ix, iy, iz) = vz(ix, iy, iz) / (1.0_num + a)
            END DO
          END DO
        END DO
      END IF
  
      IF (up == MPI_PROC_NULL) THEN
        DO iz = -2, nz+2
          DO iy = ny-ndy, ny+2
            DO ix = -2, nx+2
              a = dt * REAL(iy - (ny-ndy), num) / REAL(ndy, num)
              vx(ix, iy, iz) = vx(ix, iy, iz) / (1.0_num + a)
              vy(ix, iy, iz) = vy(ix, iy, iz) / (1.0_num + a)
              vz(ix, iy, iz) = vz(ix, iy, iz) / (1.0_num + a)
            END DO
          END DO
        END DO
      END IF
  
      IF (down == MPI_PROC_NULL) THEN
        DO iz = -2, nz+2
          DO iy = -2, ndy
            DO ix = -2, nx+2
              a = dt * REAL((ndy - iy), num) / REAL(ndy, num)
              vx(ix, iy, iz) = vx(ix, iy, iz) / (1.0_num + a)
              vy(ix, iy, iz) = vy(ix, iy, iz) / (1.0_num + a)
              vz(ix, iy, iz) = vz(ix, iy, iz) / (1.0_num + a)
            END DO
          END DO
        END DO
      END IF
  
      IF (back == MPI_PROC_NULL) THEN
        DO iz = nz-ndz, nz+2
          DO iy = -2, ny+2
            DO ix = -2, nx+2
              a = dt * REAL(iz - (nz-ndz), num) / REAL(ndz, num)
              vx(ix, iy, iz) = vx(ix, iy, iz) / (1.0_num + a)
              vy(ix, iy, iz) = vy(ix, iy, iz) / (1.0_num + a)
              vz(ix, iy, iz) = vz(ix, iy, iz) / (1.0_num + a)
            END DO
          END DO
        END DO
      END IF
   
!       IF (front == MPI_PROC_NULL) THEN
!         DO iz = -2, ndz
!           DO iy = -2, ny+2
!             DO ix = -2, nx+2
!               a = dt * REAL((ndz - iz), num) / REAL(ndz, num)
!               vx(ix, iy, iz) = vx(ix, iy, iz) / (1.0_num + a)
!               vy(ix, iy, iz) = vy(ix, iy, iz) / (1.0_num + a)
!               vz(ix, iy, iz) = vz(ix, iy, iz) / (1.0_num + a)
!             END DO
!           END DO
!         END DO
!       END IF        
  
     END IF      
  
  END SUBROUTINE damp_boundaries



  SUBROUTINE bfield_bcs

    CALL bfield_MPI

    IF (front == MPI_PROC_NULL .AND. zbc_front == BC_OTHER) THEN
      bx(:, :, -1) = bx(:, :, 2)
      bx(:, :,  0) = bx(:, :, 1)
      by(:, :, -1) = by(:, :, 2)
      by(:, :,  0) = by(:, :, 1)
      bz(:, :, -1) = bz(:, :, 1)
      bz(:, :, -2) = bz(:, :, 2)
    END IF
    IF (back == MPI_PROC_NULL .AND. zbc_back == BC_OTHER) THEN
      bx(:, :, nz+1) = bx(:, :, nz  )
      bx(:, :, nz+2) = bx(:, :, nz-1)
      by(:, :, nz+1) = by(:, :, nz  )
      by(:, :, nz+2) = by(:, :, nz-1)
      bz(:, :, nz+1) = bz(:, :, nz-1)
      bz(:, :, nz+2) = bz(:, :, nz-2)
    END IF

    IF (right == MPI_PROC_NULL .AND. xbc_right == BC_OTHER) THEN
      bx(nx+1, :, :) = bx(nx-1, :, :)
      bx(nx+2, :, :) = bx(nx-2, :, :)
      by(nx+1, :, :) = by(nx  , :, :)
      by(nx+2, :, :) = by(nx-1, :, :)
      bz(nx+1, :, :) = bz(nx  , :, :)
      bz(nx+2, :, :) = bz(nx-1, :, :)
    END IF
    IF (left == MPI_PROC_NULL .AND. xbc_left == BC_OTHER) THEN
      bx(-1, :, :) = bx(1, :, :)
      bx(-2, :, :) = bx(2, :, :)
      by( 0, :, :) = by(1, :, :)
      by(-1, :, :) = by(2, :, :)
      bz( 0, :, :) = bz(1, :, :)
      bz(-1, :, :) = bz(2, :, :)
    END IF

    IF (down == MPI_PROC_NULL .AND. ybc_down == BC_OTHER) THEN
      bx(:,  0, :) = bx(:, 1, :)
      bx(:, -1, :) = bx(:, 2, :)
      by(:, -1, :) = by(:, 1, :)
      by(:, -2, :) = by(:, 2, :)
      bz(:,  0, :) = bz(:, 1, :)
      bz(:, -1, :) = bz(:, 2, :)
    END IF
    IF (up == MPI_PROC_NULL .AND. ybc_up == BC_OTHER) THEN
      bx(:, ny+1, :) = bx(:, ny  , :)
      bx(:, ny+2, :) = bx(:, ny-1, :)
      by(:, ny+1, :) = by(:, ny-1, :)
      by(:, ny+2, :) = by(:, ny-2, :)
      bz(:, ny+1, :) = bz(:, ny  , :)
      bz(:, ny+2, :) = bz(:, ny-1, :)
    END IF

  END SUBROUTINE bfield_bcs



  SUBROUTINE energy_bcs

    CALL energy_MPI

    IF (front == MPI_PROC_NULL .AND. zbc_front == BC_OTHER) THEN
      energy(:, :,  0) = energy(:, :, 1)
      energy(:, :, -1) = energy(:, :, 2)
    END IF
    IF (back == MPI_PROC_NULL .AND. zbc_back == BC_OTHER) THEN
      energy(:, :, nz+1) = energy(:, :, nz  )
      energy(:, :, nz+2) = energy(:, :, nz-1)
    END IF

    IF (right == MPI_PROC_NULL .AND. xbc_right == BC_OTHER) THEN
      energy(nx+1, :, :) = energy(nx  , :, :)
      energy(nx+2, :, :) = energy(nx-1, :, :)
    END IF
    IF (left == MPI_PROC_NULL .AND. xbc_left == BC_OTHER) THEN
      energy( 0, :, :) = energy(1, :, :)
      energy(-1, :, :) = energy(2, :, :)
    END IF

    IF (down == MPI_PROC_NULL .AND. ybc_down == BC_OTHER) THEN
      energy(:,  0, :) = energy(:, 1, :)
      energy(:, -1, :) = energy(:, 2, :)
    END IF
    IF (up == MPI_PROC_NULL .AND. ybc_up == BC_OTHER) THEN
      energy(:, ny+1, :) = energy(:, ny  , :)
      energy(:, ny+2, :) = energy(:, ny-1, :)
    END IF

  END SUBROUTINE energy_bcs



  SUBROUTINE velocity_bcs

    CALL velocity_MPI

    IF (front == MPI_PROC_NULL .AND. zbc_front == BC_OTHER) THEN
      vx(:, :, -2:0) = 0.0_num
      vy(:, :, -2:0) = 0.0_num
      vz(:, :, -2:0) = 0.0_num
    END IF
    IF (back == MPI_PROC_NULL .AND. zbc_back == BC_OTHER) THEN
      vx(:, :, nz:nz+2) = 0.0_num
      vy(:, :, nz:nz+2) = 0.0_num
      vz(:, :, nz:nz+2) = 0.0_num
    END IF

    IF (right == MPI_PROC_NULL .AND. xbc_right == BC_OTHER) THEN
      vx(nx:nx+2, :, :) = 0.0_num
      vy(nx:nx+2, :, :) = 0.0_num
      vz(nx:nx+2, :, :) = 0.0_num
    END IF
    IF (left == MPI_PROC_NULL .AND. xbc_left == BC_OTHER) THEN
      vx(-2:0, :, :) = 0.0_num
      vy(-2:0, :, :) = 0.0_num
      vz(-2:0, :, :) = 0.0_num
    END IF

    IF (up == MPI_PROC_NULL .AND. ybc_up == BC_OTHER) THEN
      vx(:, ny:ny+2, :) = 0.0_num
      vy(:, ny:ny+2, :) = 0.0_num
      vz(:, ny:ny+2, :) = 0.0_num
    END IF
    IF (down == MPI_PROC_NULL .AND. ybc_down == BC_OTHER) THEN
      vx(:, -2:0, :) = 0.0_num
      vy(:, -2:0, :) = 0.0_num
      vz(:, -2:0, :) = 0.0_num
    END IF

  END SUBROUTINE velocity_bcs



  SUBROUTINE remap_v_bcs

    CALL remap_v_MPI

    IF (front == MPI_PROC_NULL .AND. zbc_front == BC_OTHER) THEN
      vx1(:, :, -2:0) = 0.0_num
      vy1(:, :, -2:0) = 0.0_num
      vz1(:, :, -2:0) = 0.0_num
    END IF
    IF (back == MPI_PROC_NULL .AND. zbc_back == BC_OTHER) THEN
      vx1(:, :, nz:nz+2) = 0.0_num
      vy1(:, :, nz:nz+2) = 0.0_num
      vz1(:, :, nz:nz+2) = 0.0_num
    END IF

    IF (right == MPI_PROC_NULL .AND. xbc_right == BC_OTHER) THEN
      vx1(nx:nx+2, :, :) = 0.0_num
      vy1(nx:nx+2, :, :) = 0.0_num
      vz1(nx:nx+2, :, :) = 0.0_num
    END IF
    IF (left == MPI_PROC_NULL .AND. xbc_left == BC_OTHER) THEN
      vx1(-2:0, :, :) = 0.0_num
      vy1(-2:0, :, :) = 0.0_num
      vz1(-2:0, :, :) = 0.0_num
    END IF

    IF (up == MPI_PROC_NULL .AND. ybc_up == BC_OTHER) THEN
      vx1(:, ny:ny+2, :) = 0.0_num
      vy1(:, ny:ny+2, :) = 0.0_num
      vz1(:, ny:ny+2, :) = 0.0_num
    END IF
    IF (down == MPI_PROC_NULL .AND. ybc_down == BC_OTHER) THEN
      vx1(:, -2:0, :) = 0.0_num
      vy1(:, -2:0, :) = 0.0_num
      vz1(:, -2:0, :) = 0.0_num
    END IF

  END SUBROUTINE remap_v_bcs



  SUBROUTINE density_bcs

    CALL density_MPI

    IF (front == MPI_PROC_NULL .AND. zbc_front == BC_OTHER) THEN
      rho(:, :, -1) = rho(:, :, 2)
      rho(:, :,  0) = rho(:, :, 1)
    END IF
    IF (back == MPI_PROC_NULL .AND. zbc_back == BC_OTHER) THEN
      rho(:, :, nz+1) = rho(:, :, nz  )
      rho(:, :, nz+2) = rho(:, :, nz-1)
    END IF

    IF (right == MPI_PROC_NULL .AND. xbc_right == BC_OTHER) THEN
      rho(nx+1, :, :) = rho(nx  , :, :)
      rho(nx+2, :, :) = rho(nx-1, :, :)
    END IF
    IF (left == MPI_PROC_NULL .AND. xbc_left == BC_OTHER) THEN
      rho( 0, :, :) = rho(1, :, :)
      rho(-1, :, :) = rho(2, :, :)
    END IF

    IF (down == MPI_PROC_NULL .AND. ybc_down == BC_OTHER) THEN
      rho(:,  0, :) = rho(:, 1, :)
      rho(:, -1, :) = rho(:, 2, :)
    END IF
    IF (up == MPI_PROC_NULL .AND. ybc_up == BC_OTHER) THEN
      rho(:, ny+1, :) = rho(:, ny  , :)
      rho(:, ny+2, :) = rho(:, ny-1, :)
    END IF

  END SUBROUTINE density_bcs

END MODULE boundary
