MODULE openboundary

  USE shared_data

  IMPLICIT NONE

  REAL(num) :: direction,  pfar, rhofar, uxfar, uyfar, uzfar
  REAL(num) :: bxfar, byfar, bzfar, pbc, vnorm

  REAL(num), DIMENSION(0:1) :: vxbc, vybc, vzbc
  REAL(num), DIMENSION(0:1) :: bxbc, bybc, bzbc
  REAL(num), DIMENSION(0:1) :: rbc, ebc

  REAL(num) :: v0, v1
  INTEGER :: ndx, ndy, ndz

CONTAINS

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
   
      IF (front == MPI_PROC_NULL) THEN
        DO iz = -2, ndz
          DO iy = -2, ny+2
            DO ix = -2, nx+2
              a = dt * REAL((ndz - iz), num) / REAL(ndz, num)
              vx(ix, iy, iz) = vx(ix, iy, iz) / (1.0_num + a)
              vy(ix, iy, iz) = vy(ix, iy, iz) / (1.0_num + a)
              vz(ix, iy, iz) = vz(ix, iy, iz) / (1.0_num + a)
            END DO
          END DO
        END DO
      END IF
 
     END IF      

  END SUBROUTINE damp_boundaries



  SUBROUTINE open_bcs

    REAL(num) :: bperp

    ! update ghost cells based on Riemann problem with farfield
    ! only expected to work perfectly for prblems with straight B field
    ! through boundaries which do not drastically change shape during the
    ! simulation.

    ! right boundary
    IF (xbc_right == BC_OPEN .AND. right == MPI_PROC_NULL) THEN
      DO iz = -1, nz+1
        DO iy = -1, ny+1
          ! variables carried out of domain by Riemann invariants
          vxbc(1) = vx(nx, iy, iz)
          vybc(1) = vy(nx, iy, iz)
          vzbc(1) = vz(nx, iy, iz)
          bxbc(1) = bx(nx, iy, iz)
          bybc(1) = by(nx, iy, iz)
          bzbc(1) = bz(nx, iy, iz)
          rbc(1) = rho(nx, iy, iz)
          rbc(0) = rho(nx+1, iy, iz)
          ebc(1) = energy(nx, iy, iz)

          pbc = (gamma - 1.0_num) * energy(nx, iy, iz) * rho(nx, iy, iz)

          ! farfield values carried into domain
          pfar = (gamma - 1.0_num) * energy(nx+2, iy, iz) * rho(nx+2, iy, iz)

          rhofar = rho(nx+2, iy, iz)
          uxfar = vx(nx+2, iy, iz)
          uyfar = vy(nx+2, iy, iz)
          uzfar = vz(nx+2, iy, iz)
          bxfar = bx(nx+2, iy, iz)
          byfar = by(nx+2, iy, iz)
          bzfar = bz(nx+2, iy, iz)

          ! direction of boundary (1 = right ; -1 = left)
          direction = 1.0_num
          vnorm = vx(nx, iy, iz)

          ! select correct open bc solver
          bperp = SQRT(byfar**2 + bzfar**2)

          IF (ABS(bxfar) <= 0.01_num * bperp) THEN
            CALL open_bcs_1
          ELSE IF (bperp <= 0.01_num * ABS(bxfar)) THEN
            CALL open_bcs_2
          ELSE
            CALL open_bcs_3
          END IF

          v0 = SQRT(vx(nx+1, iy, iz)**2 &
              + vy(nx+1, iy, iz)**2 + vz(nx+1, iy, iz)**2)

          v1 = SQRT(vxbc(0)**2 + vybc(0)**2 + vzbc(0)**2)

          IF (ABS(v0 - v1) < ABS(dv_right(iy, iz))) THEN
            rho(nx+1, iy, iz) = rbc(0)
            energy(nx+1, iy, iz) = ebc(0)
            bz(nx+1, iy, iz) = bzbc(0)
            by(nx+1, iy, iz) = bybc(0)
            bx(nx+1, iy, iz) = bxbc(0)
            vx(nx+1, iy, iz) = vxbc(0)
            vy(nx+1, iy, iz) = vybc(0)
            vz(nx+1, iy, iz) = vzbc(0)
          ELSE
            vx(nx+1, iy, iz) = vxbc(1)
            vy(nx+1, iy, iz) = vybc(1)
            vz(nx+1, iy, iz) = vzbc(1)
          END IF
        END DO
      END DO
    END IF

    ! left bounday
    IF (xbc_left == BC_OPEN .AND. left == MPI_PROC_NULL) THEN
      DO iz = -1, nz+1
        DO iy = -1, ny+1
          vxbc(1) = vx(0, iy, iz)
          vybc(1) = vy(0, iy, iz)
          vzbc(1) = vz(0, iy, iz)
          bxbc(1) = bx(0, iy, iz)
          bybc(1) = by(1, iy, iz)
          bzbc(1) = bz(1, iy, iz)
          rbc(1) = rho(1, iy, iz)
          rbc(0) = rho(0, iy, iz)
          ebc(1) = energy(1, iy, iz)

          pbc = (gamma - 1.0_num) * energy(1, iy, iz) * rho(1, iy, iz)

          pfar = (gamma - 1.0_num) * energy(-1, iy, iz) * rho(-1, iy, iz)

          rhofar = rho(-1, iy, iz)
          uxfar = vx(-2, iy, iz)
          uyfar = vy(-2, iy, iz)
          uzfar = vz(-2, iy, iz)
          bxfar = bx(-2, iy, iz)
          byfar = by(-1, iy, iz)
          bzfar = bz(-1, iy, iz)

          direction = -1.0_num
          vnorm = vx(0, iy, iz)
          bperp = SQRT(byfar**2 + bzfar**2)

          IF (ABS(bxfar) <= 0.01_num * bperp) THEN
            CALL open_bcs_1
          ELSE IF (bperp <= 0.01_num * ABS(bxfar)) THEN
            CALL open_bcs_2
          ELSE
            CALL open_bcs_3
          END IF

          v0 = SQRT(vx(-1, iy, iz)**2 + vy(-1, iy, iz)**2 + vz(-1, iy, iz)**2)
          v1 = SQRT(vxbc(0)**2 + vybc(0)**2 + vzbc(0)**2)

          IF (ABS(v0 - v1) < ABS(dv_left(iy, iz))) THEN
            rho(0, iy, iz) = rbc(0)
            energy(0, iy, iz) = ebc(0)
            bz(0, iy, iz) = bzbc(0)
            bx(-1, iy, iz) = bxbc(0)
            vx(-1, iy, iz) = vxbc(0)
            vy(-1, iy, iz) = vybc(0)
            vz(-1, iy, iz) = vzbc(0)
            by(0, iy, iz) = bybc(0)
          ELSE
            vx(-1, iy, iz) = vxbc(1)
            vy(-1, iy, iz) = vybc(1)
            vz(-1, iy, iz) = vzbc(1)
          END IF
        END DO
      END DO
    END IF

    ! top boundary
    IF (ybc_up == BC_OPEN .AND. up == MPI_PROC_NULL) THEN
      DO iz = -1, nz+1
        DO ix = -1, nx+1
          vxbc(1) = vy(ix, ny, iz)
          vybc(1) = vx(ix, ny, iz)
          vzbc(1) = vz(ix, ny, iz)
          bxbc(1) = by(ix, ny, iz)
          bybc(1) = bx(ix, ny, iz)
          bzbc(1) = bz(ix, ny, iz)
          rbc(1) = rho(ix, ny, iz)
          rbc(0) = rho(ix, ny+1, iz)
          ebc(1) = energy(ix, ny, iz)

          pbc = (gamma - 1.0_num) * energy(ix, ny, iz) * rho(ix, ny, iz)

          pfar = (gamma - 1.0_num) * energy(ix, ny+2, iz) * rho(ix, ny+2, iz)

          rhofar = rho(ix, ny+2, iz)
          uxfar = vy(ix, ny+2, iz)
          uyfar = vx(ix, ny+2, iz)
          uzfar = vz(ix, ny+2, iz)
          bxfar = by(ix, ny+2, iz)
          byfar = bx(ix, ny+2, iz)
          bzfar = bz(ix, ny+2, iz)

          direction = 1.0_num
          vnorm = vy(ix, ny, iz)
          bperp = SQRT(byfar**2 + bzfar**2)

          IF (ABS(bxfar) <= 0.01_num * bperp) THEN
            CALL open_bcs_1
          ELSE IF (bperp <= 0.01_num * ABS(bxfar)) THEN
            CALL open_bcs_2
          ELSE
            CALL open_bcs_3
          END IF

          v0 = SQRT(vx(ix, ny+1, iz)**2 &
              + vy(ix, ny+1, iz)**2 + vz(ix, ny+1, iz)**2)
          v1 = SQRT(vxbc(0)**2 + vybc(0)**2 + vzbc(0)**2)

          IF (ABS(v0 - v1) < ABS(dv_up(ix, iz))) THEN
            rho(ix, ny+1, iz) = rbc(0)
            energy(ix, ny+1, iz) = ebc(0)
            bz(ix, ny+1, iz) = bzbc(0)
            by(ix, ny+1, iz) = bxbc(0)
            vx(ix, ny+1, iz) = vybc(0)
            vy(ix, ny+1, iz) = vxbc(0)
            vz(ix, ny+1, iz) = vzbc(0)
            bx(ix, ny+1, iz) = bybc(0)
          ELSE
            vx(ix, ny+1, iz) = vybc(1)
            vy(ix, ny+1, iz) = vxbc(1)
            vz(ix, ny+1, iz) = vzbc(1)
          END IF
        END DO
      END DO
    END IF

    ! bottom boundary
    IF (ybc_down == BC_OPEN .AND. down == MPI_PROC_NULL) THEN
      DO iz = -1, nz+1
        DO ix = -1, nx+1
          vxbc(1) = vy(ix, 0, iz)
          vybc(1) = vx(ix, 0, iz)
          vzbc(1) = vz(ix, 0, iz)
          bxbc(1) = by(ix, 0, iz)
          bybc(1) = bx(ix, 1, iz)
          bzbc(1) = bz(ix, 1, iz)
          rbc(1) = rho(ix, 1, iz)
          rbc(0) = rho(ix, 0, iz)
          ebc(1) = energy(ix, 1, iz)

          pbc = (gamma - 1.0_num) * energy(ix, 1, iz) * rho(ix, 1, iz)

          pfar = (gamma - 1.0_num) * energy(ix, -1, iz) * rho(ix, -1, iz)

          rhofar = rho(ix, -1, iz)
          uxfar = vy(ix, -2, iz)
          uyfar = vx(ix, -2, iz)
          uzfar = vz(ix, -2, iz)
          bxfar = by(ix, -2, iz)
          byfar = bx(ix, -1, iz)
          bzfar = bz(ix, -1, iz)

          direction = -1.0_num
          vnorm = vy(ix, 0, iz)
          bperp = SQRT(byfar**2 + bzfar**2)

          IF (ABS(bxfar) <= 0.01_num * bperp) THEN
            CALL open_bcs_1
          ELSE IF (bperp <= 0.01_num * ABS(bxfar)) THEN
            CALL open_bcs_2
          ELSE
            CALL open_bcs_3
          END IF

          v0 = SQRT(vx(ix, -1, iz)**2 + vy(ix, -1, iz)**2 + vz(ix, -1, iz)**2)
          v1 = SQRT(vxbc(0)**2 + vybc(0)**2 + vzbc(0)**2)

          IF (ABS(v0 - v1) < ABS(dv_down(ix, iz))) THEN
            rho(ix, 0, iz) = rbc(0)
            energy(ix, 0, iz) = ebc(0)
            bz(ix, 0, iz) = bzbc(0)
            by(ix, -1, iz) = bxbc(0)
            vx(ix, -1, iz) = vybc(0)
            vy(ix, -1, iz) = vxbc(0)
            vz(ix, -1, iz) = vzbc(0)
            bx(ix, 0, iz) = bybc(0)
          ELSE
            vx(ix, -1, iz) = vybc(1)
            vy(ix, -1, iz) = vxbc(1)
            vz(ix, -1, iz) = vzbc(1)
          END IF
        END DO
      END DO
    END IF

    ! back boundary
    IF (zbc_back == BC_OPEN .AND. back == MPI_PROC_NULL) THEN
      DO iy = -1, ny+1
        DO ix = -1, nx+1
          vxbc(1) = vz(ix, iy, nz)
          vybc(1) = vy(ix, iy, nz)
          vzbc(1) = vx(ix, iy, nz)
          bxbc(1) = bz(ix, iy, nz)
          bybc(1) = by(ix, iy, nz)
          bzbc(1) = bx(ix, iy, nz)
          rbc(1) = rho(ix, iy, nz)
          rbc(0) = rho(ix, iy, nz+1)
          ebc(1) = energy(ix, iy, nz)

          pbc = (gamma - 1.0_num) * energy(ix, iy, nz) * rho(ix, iy, nz)

          pfar = (gamma - 1.0_num) * energy(ix, iy, nz+2) * rho(ix, iy, nz+2)

          rhofar = rho(ix, iy, nz+2)
          uxfar = vz(ix, iy, nz+2)
          uyfar = vy(ix, iy, nz+2)
          uzfar = vx(ix, iy, nz+2)
          bxfar = bz(ix, iy, nz+2)
          byfar = by(ix, iy, nz+2)
          bzfar = bx(ix, iy, nz+2)

          direction = 1.0_num
          vnorm = vz(ix, iy, nz)
          bperp = SQRT(byfar**2 + bzfar**2)

          IF (ABS(bxfar) <= 0.01_num * bperp) THEN
            CALL open_bcs_1
          ELSE IF (bperp <= 0.01_num * ABS(bxfar)) THEN
            CALL open_bcs_2
          ELSE
            CALL open_bcs_3
          END IF

          v0 = SQRT(vx(ix, iy, nz+1)**2 &
              + vy(ix, iy, nz+1)**2 + vz(ix, iy, nz+1)**2)
          v1 = SQRT(vxbc(0)**2 + vybc(0)**2 + vzbc(0)**2)

          IF (ABS(v0 - v1) < ABS(dv_back(ix, iy))) THEN
            rho(ix, iy, nz+1) = rbc(0)
            energy(ix, iy, nz+1) = ebc(0)
            bz(ix, iy, nz+1) = bxbc(0)
            by(ix, iy, nz+1) = bybc(0)
            vx(ix, iy, nz+1) = vzbc(0)
            vy(ix, iy, nz+1) = vybc(0)
            vz(ix, iy, nz+1) = vxbc(0)
            bx(ix, iy, nz+1) = bzbc(0)
          ELSE
            vx(ix, iy, nz+1) = vzbc(1)
            vy(ix, iy, nz+1) = vybc(1)
            vz(ix, iy, nz+1) = vxbc(1)
          END IF
        END DO
      END DO
    END IF

    ! front boundary
    IF (zbc_front == BC_OPEN .AND. front == MPI_PROC_NULL) THEN
      DO iy = -1, ny+1
        DO ix = -1, nx+1
          vxbc(1) = vz(ix, iy, 0)
          vybc(1) = vy(ix, iy, 0)
          vzbc(1) = vx(ix, iy, 0)
          bxbc(1) = bz(ix, iy, 0)
          bybc(1) = by(ix, iy, 1)
          bzbc(1) = bx(ix, iy, 1)
          rbc(1) = rho(ix, iy, 1)
          rbc(0) = rho(ix, iy, 0)
          ebc(1) = energy(ix, iy, 1)

          pbc = (gamma - 1.0_num) * energy(ix, iy, 1) * rho(ix, iy, 1)

          pfar = (gamma - 1.0_num) * energy(ix, iy, -1) * rho(ix, iy, -1)

          rhofar = rho(ix, iy, -1)
          uxfar = vz(ix, iy, -2)
          uyfar = vy(ix, iy, -2)
          uzfar = vx(ix, iy, -2)
          bxfar = bz(ix, iy, -2)
          byfar = by(ix, iy, -1)
          bzfar = bx(ix, iy, -1)

          direction = -1.0_num
          vnorm = vz(ix, iy, 0)
          bperp = SQRT(byfar**2 + bzfar**2)

          IF (ABS(bxfar) <= 0.01_num * bperp) THEN
            CALL open_bcs_1
          ELSE IF (bperp <= 0.01_num * ABS(bxfar)) THEN
            CALL open_bcs_2
          ELSE
            CALL open_bcs_3
          END IF

          v0 = SQRT(vx(ix, iy, -1)**2 + vy(ix, iy, -1)**2 + vz(ix, iy, -1)**2)
          v1 = SQRT(vxbc(0)**2 + vybc(0)**2 + vzbc(0)**2)

          IF (ABS(v0 - v1) < ABS(dv_front(ix, iy))) THEN
            rho(ix, iy, 0) = rbc(0)
            energy(ix, iy, 0) = ebc(0)
            bz(ix, iy, -1) = bxbc(0)
            by(ix, iy, 0) = bybc(0)
            vx(ix, iy, -1) = vzbc(0)
            vy(ix, iy, -1) = vybc(0)
            vz(ix, iy, -1) = vxbc(0)
            bx(ix, iy, 0) = bzbc(0)
          ELSE
            vx(ix, iy, -1) = vzbc(1)
            vy(ix, iy, -1) = vybc(1)
            vz(ix, iy, -1) = vxbc(1)
          END IF
        END DO
      END DO
    END IF

  END SUBROUTINE open_bcs



  SUBROUTINE open_bcs_1

    ! open bc when bx = 0
    REAL(num) :: c0, ct, cf
    REAL(num) :: pg, rhog, cffar, c0far, ctfar
    REAL(num) :: pmagg, uxg
    REAL(num), DIMENSION(3) :: vtest, pstar, vstar, rhostar, pmagstar
    INTEGER :: i

    c0far = SQRT(gamma * pfar / rhofar)
    ctfar = SQRT((byfar**2 + bzfar**2) / rhofar)
    cffar = SQRT(c0far**2 + ctfar**2)

    c0 = SQRT(gamma * (gamma - 1.0_num) * ebc(1))
    ct = SQRT((bybc(1)**2 + bzbc(1)**2) / rbc(1))
    cf = SQRT(c0**2 + ct**2)

    ! Define the speeds of the characteristics to be checked along
    vtest(1) = vnorm + cf
    vtest(2) = vnorm - cf
    vtest(3) = vnorm

    DO i = 1, 3
      IF (direction * vtest(i) >= 0.0_num) THEN
        pstar(i) =  pbc + 0.5_num * (bybc(1)**2 + bzbc(1)**2 - bxbc(1)**2)
        pmagstar(i) =  0.5_num * (bybc(1)**2 + bzbc(1)**2 - bxbc(1)**2)
        vstar(i) = vxbc(1)
        rhostar(i) = rbc(1)
      ELSE
        pstar(i) = pfar + 0.5_num * (byfar**2 + bzfar**2 - bxfar**2)
        pmagstar(i) = 0.5_num * (byfar**2 + bzfar**2 - bxfar**2)
        vstar(i) = uxfar
        rhostar(i) = rhofar
      END IF
    END DO

    bxbc(0) = bxbc(1)
    bybc(0) = bybc(1)
    bzbc(0) = bzbc(1)
    pmagg = 0.5_num * (bybc(0)**2 + bzbc(0)**2 - bxbc(0)**2)
    pg = 0.5_num &
        * (pstar(1) + pstar(2) + rhofar * cffar * (vstar(1) - vstar(2)))
    rhog = ((pg - pmagg) - (pstar(3) - pmagstar(3))) / c0far**2 + rhostar(3)
    rbc(0) = MAX(rhog, none_zero)
    ebc(0) = MAX(pg - pmagg, none_zero) / (gamma - 1.0_num) / rbc(0)
    uxg = 0.5_num &
        * (vstar(1) + vstar(2) + (pstar(1) - pstar(2)) / (rhofar * cffar))
    vxbc(0) = uxg
    vybc(0) = vybc(1)
    vzbc(0) = vzbc(1)

  END SUBROUTINE open_bcs_1



  SUBROUTINE open_bcs_2

    ! open bc when bperp = 0
    REAL(num) :: lambdayfar, lambdazfar
    REAL(num) :: c0, cx
    REAL(num) :: pg, rhog, c0far, cxfar
    REAL(num) :: pmagg, uxg, uyg, uzg, lambdag
    REAL(num), DIMENSION(5) :: vtest, pstar, uxstar, rhostar, pmagstar
    REAL(num), DIMENSION(5) :: uystar, lambdaystar, lambdazstar, uzstar
    INTEGER :: i

    lambdayfar = -bxfar * byfar
    lambdazfar = -bxfar * bzfar
    c0far = SQRT(gamma * pfar / rhofar)
    cxfar = SQRT(bxfar**2 / rhofar)

    c0 = SQRT(gamma * (gamma - 1.0_num) * ebc(1))
    cx = SQRT(bxbc(1)**2 / rbc(1))

    ! Define the speeds of the characteristics to be checked along
    vtest(1) = vnorm + c0
    vtest(2) = vnorm - c0
    vtest(3) = vnorm + cx
    vtest(4) = vnorm - cx
    vtest(5) = vnorm

    DO i = 1, 5
      IF (direction * vtest(i) >= 0.0_num) THEN
        pstar(i) = pbc + 0.5_num * (bybc(1)**2 + bzbc(1)**2 - bxbc(1)**2)
        pmagstar(i) = 0.5_num * (bybc(1)**2 + bzbc(1)**2 - bxbc(1)**2)
        lambdaystar(i) = -bybc(1) * bxbc(1)
        lambdazstar(i) = -bzbc(1) * bxbc(1)
        uxstar(i) = vxbc(1)
        uystar(i) = vybc(1)
        uzstar(i) = vzbc(1)
        rhostar(i) = rbc(1)
      ELSE
        pstar(i) = pfar + 0.5_num * (byfar**2 + bzfar**2 - bxfar**2)
        pmagstar(i) = 0.5_num * (byfar**2 + bzfar**2 - bxfar**2)
        lambdaystar(i) = lambdayfar
        lambdazstar(i) = lambdazfar
        uxstar(i) = uxfar
        uystar(i) = uyfar
        uzstar(i) = uzfar
        rhostar(i) = rhofar
      END IF
    END DO

    bxbc(0) = bxbc(1)
    lambdag = 0.5_num * (lambdaystar(3) + lambdaystar(4) &
        + rhofar * cxfar * (uystar(3) - uystar(4)))

    bybc(0) = -lambdag / bxbc(0)
    lambdag = 0.5_num * (lambdazstar(3) + lambdazstar(4) &
        + rhofar * cxfar * (uzstar(3) - uzstar(4)))

    bzbc(0) = -lambdag  / bxbc(0)
    uyg = 0.5_num * (uystar(3) + uystar(4) &
        + (lambdaystar(3) - lambdaystar(4)) / (rhofar * cxfar))

    vybc(0) = uyg
    uzg = 0.5_num * (uzstar(3) + uzstar(4) &
        + (lambdazstar(3) - lambdazstar(4)) / (rhofar * cxfar))

    vzbc(0) = uzg
    pmagg = 0.5_num * (bybc(0)**2 + bzbc(0)**2 - bxbc(0)**2)

    pg = 0.5_num &
        * (pstar(1) + pstar(2) + rhofar * c0far * (uxstar(1) - uxstar(2)))

    rhog = ((pg - pmagg) - (pstar(5) - pmagstar(5))) / c0far**2 + rhostar(5)
    rbc(0) = MAX(rhog, none_zero)
    ebc(0) = MAX(pg - pmagg, none_zero) / (gamma - 1.0_num) / rbc(0)

    uxg = 0.5_num &
        * (uxstar(1) + uxstar(2) + (pstar(1) - pstar(2)) / (rhofar * c0far))
    vxbc(0) = uxg

  END SUBROUTINE open_bcs_2



  SUBROUTINE open_bcs_3
    ! Solve for when bx and bperp are non zero. Solves in the coordinate system
    ! such that y-axis points along by_farfield
    REAL(num), DIMENSION(7) :: vtest
    INTEGER :: i
    REAL(num) :: a, b, c, d, e, f, g
    REAL(num) :: pmagg, pmagfar, phi, theta
    REAL(num) :: c0, cx, ct, cf, cs
    REAL(num) :: lambdafar, byfar2
    REAL(num) :: c0far, cxfar, ctfar, cffar, csfar
    REAL(num) :: pg, rhog, uxg, uyg, uzg, lambdag, byg, bxg, bzg
    REAL(num), DIMENSION(7) :: pstar, uxstar, uystar, uzstar, rhostar
    REAL(num), DIMENSION(7) :: lambdastar, pmagstar, bzstar

    ! Setup the far field variables
    byfar2 = SQRT(byfar**2 + bzfar**2)
    phi = ATAN2(bzfar, byfar)
    pmagfar = 0.5_num * (byfar2**2 - bxfar**2)
    pfar = pfar + pmagfar
    lambdafar = -bxfar * byfar2
    c0far = SQRT(gamma * (pfar - pmagfar) / rhofar)
    cxfar = SQRT(bxfar**2 / rhofar)
    ctfar = SQRT(byfar2**2 / rhofar)
    cffar = SQRT(0.5_num * ((c0far**2 + cxfar**2 + ctfar**2) &
        + SQRT((c0far**2 + cxfar**2 + ctfar**2)**2 &
        - 4.0_num * c0far**2 * cxfar**2)))
    csfar = SQRT(0.5_num * ((c0far**2 + cxfar**2 + ctfar**2) &
        - SQRT((c0far**2 + cxfar**2 + ctfar**2)**2 &
        - 4.0_num * c0far**2 * cxfar**2)))

    ! Setup the speeds
    c0 = SQRT(gamma * (gamma - 1.0_num) * ebc(1))
    cx = SQRT(bxbc(1)**2 / rbc(1))
    ct = SQRT((bybc(1)**2 + bzbc(1)**2) / rbc(1))
    cf = SQRT(0.5_num * ((c0**2 + cx**2 + ct**2) &
        + SQRT((c0**2 + cx**2 + ct**2)**2 - 4.0_num * c0**2 * cx**2)))
    cs = SQRT(0.5_num * ((c0**2 + cx**2 + ct**2) &
        - SQRT((c0**2 + cx**2 + ct**2)**2 - 4.0_num * c0**2 * cx**2)))

    ! Define the speeds of the characteristics to be checked along
    vtest(1) = vnorm + cf
    vtest(2) = vnorm - cf
    vtest(3) = vnorm - cs
    vtest(4) = vnorm + cs
    vtest(5) = vnorm
    vtest(6) = vnorm + cx
    vtest(7) = vnorm - cx

    ! Now check which characteristics are inflowing, outflowing, non-moving
    DO i = 1, 7
      IF (direction * vtest(i)  >= 0.0_num) THEN
        pstar(i) = pbc + 0.5_num * (bybc(1)**2 + bzbc(1)**2 - bxbc(1)**2)
        pmagstar(i) = 0.5_num * (bybc(1)**2 + bzbc(1)**2 - bxbc(1)**2)
        rhostar(i) = rbc(1)
        lambdastar(i) = -bxbc(1) * (bybc(1) * byfar + bzbc(1) * bzfar) / byfar2
        theta = ATAN2(bzbc(1), bybc(1))
        bzstar(i) = SQRT(bybc(1)**2 + bzbc(1)**2) * SIN(theta - phi)
        uystar(i) = (vybc(1) * byfar + vzbc(1) * bzfar) / byfar2
        uzstar(i) = SQRT(vybc(1)**2 + vzbc(1)**2) * SIN(theta - phi)
        uxstar(i) = vxbc(1)
      ELSE
        pstar(i) = pfar
        pmagstar(i) = pmagfar
        rhostar(i) = rhofar
        lambdastar(i) = lambdafar
        bzstar(i) = 0.0_num
        uystar(i) = uyfar
        uzstar(i) = uzfar
        uxstar(i) = uxfar
      END IF

    END DO

    ! Now setup the constants that are defined in the solution
    a = (cffar**2 - cxfar**2)
    b = (lambdafar / rhofar)
    c = (csfar**2 - cxfar**2)

    d = pstar(1) + pstar(2) + rhofar * cffar * (uxstar(1) - uxstar(2))
    e = lambdastar(1) + lambdastar(2) + rhofar * cffar * (uystar(1) - uystar(2))
    f = pstar(3) + pstar(4) - rhofar * csfar * (uxstar(3) - uxstar(4))
    g = lambdastar(3) + lambdastar(4) - rhofar * csfar * (uystar(3) - uystar(4))
    pg = 0.5_num * (a*d + b*e - c*f - b*g) / (a - c)
    lambdag = 0.5_num * (c * (a*d + b*e) - a * (c*f + b*g)) / (b * (c - a))

    d = (pstar(1) - pstar(2)) / (rhofar * cffar) + (uxstar(1) + uxstar(2))
    e = (lambdastar(1) - lambdastar(2)) &
        / (rhofar * cffar) + (uystar(1) + uystar(2))
    f = (pstar(4) - pstar(3)) / (rhofar * csfar) + (uxstar(3) + uxstar(4))
    g = (lambdastar(4) - lambdastar(3)) &
        / (rhofar * csfar) + (uystar(3) + uystar(4))
    uxg = 0.5_num * (a*d + b*e - c*f - b*g) / (a - c)
    uyg = 0.5_num * (c * (a*d + b*e) - a * (c*f + b*g)) / (b * (c - a))

    a = cxfar * rhofar / bxfar
    bzg = 0.5_num * (bzstar(6) + bzstar(7) + a * (uzstar(7) - uzstar(6)))
    uzg = 0.5_num * (uzstar(6) + uzstar(7) + (bzstar(7) - bzstar(6)) / a)

    bxg = bxbc(1)
    byg = -lambdag / bxg

    pmagg = 0.5_num * (byg**2 + bzg**2 - bxg**2)
    rhog = ((pg - pmagg) - (pstar(5) - pmagstar(5))) / c0**2 + rhostar(5)
    rhog = MAX(rhog, none_zero)
    rbc(0) = rhog
    ebc(0) = MAX(pg - pmagg, none_zero) / ((gamma - 1.0_num) * rhog)

    ! rotate back to grid coordinate system
    bxbc(0) = bxg
    bybc(0) = byg * COS(phi) - bzg * SIN(phi)
    bzbc(0) = byg * SIN(phi) + bzg * COS(phi)
    vxbc(0) = uxg
    vybc(0) = uyg * COS(phi) - uzg * SIN(phi)
    vzbc(0) = uyg * SIN(phi) + uzg * COS(phi)

  END SUBROUTINE open_bcs_3

END MODULE openboundary
