MODULE remap

  USE shared_data; USE xremap ;  USE yremap ; USE zremap

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eulerian_remap

CONTAINS


  SUBROUTINE eulerian_remap(i)  ! Strang splitting

    INTEGER, INTENT(IN) :: i
    INTEGER :: case_test


    delta_ke = 0.0_num
    xpass = 1
    ypass = 1
    zpass = 1

    DO iz = -1, nz+2
       DO iy = -1, ny+2
          DO ix = -1, nx+2
             bx(ix,iy,iz) = bx(ix,iy,iz) * (dyb(iy) * dzb(iz))  
             by(ix,iy,iz) = by(ix,iy,iz) * (dxb(ix) * dzb(iz))
             bz(ix,iy,iz) = bz(ix,iy,iz) * (dxb(ix) * dyb(iy))
          END DO
       END DO
    END DO
    DO iz = -1, nz+2
       DO iy = -1, ny+2
          bx(-2,iy,iz) = bx(-2,iy,iz) * (dyb(iy) * dzb(iz))  
       END DO
    END DO
    DO iz = -1, nz+2
       DO ix = -1, nx+2
          by(ix,-2,iz) = by(ix,-2,iz) * (dxb(ix) * dzb(iz))
       END DO
    END DO
    DO iy = -1, ny+2
       DO ix = -1, nx+2
          bz(ix,iy,-2) = bz(ix,iy,-2) * (dxb(ix) * dyb(iy))
       END DO
    END DO

    case_test = MODULO(i,6) 
    SELECT CASE(case_test)  ! Strang ordering
    CASE (0)
       CALL remap_x
       CALL remap_y
       CALL remap_z
    CASE (1)
       CALL remap_y
       CALL remap_z
       CALL remap_x
    CASE (2)
       CALL remap_z
       CALL remap_x
       CALL remap_y
    CASE (3)
       CALL remap_x
       CALL remap_z
       CALL remap_y
    CASE (4)
       CALL remap_z
       CALL remap_y
       CALL remap_x
    CASE (5)
       CALL remap_y
       CALL remap_x
       CALL remap_z
    END SELECT

    DO iz = -1, nz+2
       DO iy = -1, ny+2
          DO ix = -1, nx+2
             bx(ix,iy,iz) = bx(ix,iy,iz) / (dyb(iy) * dzb(iz))  
             by(ix,iy,iz) = by(ix,iy,iz) / (dxb(ix) * dzb(iz))
             bz(ix,iy,iz) = bz(ix,iy,iz) / (dxb(ix) * dyb(iy))
          END DO
       END DO
    END DO
    DO iz = -1, nz+2
       DO iy = -1, ny+2
          bx(-2,iy,iz) = bx(-2,iy,iz) / (dyb(iy) * dzb(iz))  
       END DO
    END DO
    DO iz = -1, nz+2
       DO ix = -1, nx+2
          by(ix,-2,iz) = by(ix,-2,iz) / (dxb(ix) * dzb(iz))
       END DO
    END DO
    DO iy = -1, ny+2
       DO ix = -1, nx+2
          bz(ix,iy,-2) = bz(ix,iy,-2) / (dxb(ix) * dyb(iy))
       END DO
    END DO

  END SUBROUTINE eulerian_remap

END MODULE remap













