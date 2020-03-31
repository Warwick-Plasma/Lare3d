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

MODULE openboundary

  USE shared_data

  IMPLICIT NONE

  REAL(num) :: pfar, rhofar, efar, uxfar, uyfar, uzfar
  REAL(num) :: bxfar, byfar, bzfar, pbc, vnorm

  REAL(num), DIMENSION(0:1) :: vxbc, vybc, vzbc
  REAL(num), DIMENSION(0:1) :: bxbc, bybc, bzbc
  REAL(num), DIMENSION(0:1) :: rbc, ebc

CONTAINS

  SUBROUTINE open_bcs

    REAL(num) :: bperp

    ! x_min boundary
    IF (xbc_min == BC_OPEN .AND. proc_x_min == MPI_PROC_NULL) THEN
      DO iz = 0, nz+1
        DO iy = 0, ny+1
          ! Variables carried out of domain by Riemann invariants
          vxbc(1) = -vx(0,iy,iz)
          vybc(1) =  vy(0,iy,iz)
          vzbc(1) =  vz(0,iy,iz)
          bxbc(1) = -bx(0,iy,iz)
          bybc(1) =  by(1,iy,iz)
          bzbc(1) =  bz(1,iy,iz)
          rbc (1) = rho(1,iy,iz)
          ebc (1) = energy(1,iy,iz)

          pbc  = (gamma - 1.0_num) * energy( 1,iy,iz) * rho( 1,iy,iz)

          ! Farfield values carried into domain
          pfar = (gamma - 1.0_num) * energy(-1,iy,iz) * rho(-1,iy,iz)

          uxfar  = -vx(-2,iy,iz)
          uyfar  =  vy(-2,iy,iz)
          uzfar  =  vz(-2,iy,iz)
          bxfar  = -bx(-2,iy,iz)
          byfar  =  by(-1,iy,iz)
          bzfar  =  bz(-1,iy,iz)
          rhofar = rho(-1,iy,iz)
          efar = energy(-1,iy,iz)

          vnorm = -vx(0,iy,iz)

          ! Select correct open bc solver
          bperp = SQRT(byfar**2 + bzfar**2)

          IF (ABS(bperp) <= none_zero) THEN
            CALL open_bcs_alfven
          ELSE
            CALL open_bcs_fast
          END IF

          bx (-1,iy,iz) = -bxbc(0)
          by ( 0,iy,iz) =  bybc(0)
          bz ( 0,iy,iz) =  bzbc(0)
          rho( 0,iy,iz) =   rbc(0)
          energy(0,iy,iz) = ebc(0)
        END DO
      END DO
    END IF

    ! x_max boundary
    IF (xbc_max == BC_OPEN .AND. proc_x_max == MPI_PROC_NULL) THEN
      DO iz = 0, nz+1
        DO iy = 0, ny+1
          ! Variables carried out of domain by Riemann invariants
          vxbc(1) =  vx(nx,iy,iz)
          vybc(1) =  vy(nx,iy,iz)
          vzbc(1) =  vz(nx,iy,iz)
          bxbc(1) =  bx(nx,iy,iz)
          bybc(1) =  by(nx,iy,iz)
          bzbc(1) =  bz(nx,iy,iz)
          rbc (1) = rho(nx,iy,iz)
          ebc (1) = energy(nx,iy,iz)

          pbc  = (gamma - 1.0_num) * energy(nx  ,iy,iz) * rho(nx  ,iy,iz)

          ! Farfield values carried into domain
          pfar = (gamma - 1.0_num) * energy(nx+2,iy,iz) * rho(nx+2,iy,iz)

          uxfar  =  vx(nx+2,iy,iz)
          uyfar  =  vy(nx+2,iy,iz)
          uzfar  =  vz(nx+2,iy,iz)
          bxfar  =  bx(nx+2,iy,iz)
          byfar  =  by(nx+2,iy,iz)
          bzfar  =  bz(nx+2,iy,iz)
          rhofar = rho(nx+2,iy,iz)
          efar = energy(nx+2,iy,iz)

          vnorm = vx(nx,iy,iz)

          ! Select correct open bc solver
          bperp = SQRT(byfar**2 + bzfar**2)

          IF (bperp <= none_zero) THEN
            CALL open_bcs_alfven
          ELSE
            CALL open_bcs_fast
          END IF

          bx (nx+1,iy,iz) = bxbc(0)
          by (nx+1,iy,iz) = bybc(0)
          bz (nx+1,iy,iz) = bzbc(0)
          rho(nx+1,iy,iz) =  rbc(0)
          energy(nx+1,iy,iz) = ebc(0)
        END DO
      END DO
    END IF

    ! y_min boundary
    IF (ybc_min == BC_OPEN .AND. proc_y_min == MPI_PROC_NULL) THEN
      DO iz = 0, nz+1
        DO ix = 0, nx+1
          ! Variables carried out of domain by Riemann invariants
          vxbc(1) = -vy(ix,0,iz)
          vybc(1) =  vx(ix,0,iz)
          vzbc(1) =  vz(ix,0,iz)
          bxbc(1) = -by(ix,0,iz)
          bybc(1) =  bx(ix,1,iz)
          bzbc(1) =  bz(ix,1,iz)
          rbc (1) = rho(ix,1,iz)
          ebc (1) = energy(ix,1,iz)

          pbc  = (gamma - 1.0_num) * energy(ix, 1,iz) * rho(ix, 1,iz)

          ! Farfield values carried into domain
          pfar = (gamma - 1.0_num) * energy(ix,-1,iz) * rho(ix,-1,iz)

          uxfar  = -vy(ix,-2,iz)
          uyfar  =  vx(ix,-2,iz)
          uzfar  =  vz(ix,-2,iz)
          bxfar  = -by(ix,-2,iz)
          byfar  =  bx(ix,-1,iz)
          bzfar  =  bz(ix,-1,iz)
          rhofar = rho(ix,-1,iz)
          efar = energy(ix,-1,iz)

          vnorm = -vy(ix,0,iz)

          ! Select correct open bc solver
          bperp = SQRT(byfar**2 + bzfar**2)

          IF (bperp <= none_zero) THEN
            CALL open_bcs_alfven
          ELSE
            CALL open_bcs_fast
          END IF

          bx (ix, 0,iz) =  bybc(0)
          by (ix,-1,iz) = -bxbc(0)
          bz (ix, 0,iz) =  bzbc(0)
          rho(ix, 0,iz) =   rbc(0)
          energy(ix,0,iz) = ebc(0)
        END DO
      END DO
    END IF

    ! y_max boundary
    IF (ybc_max == BC_OPEN .AND. proc_y_max == MPI_PROC_NULL) THEN
      DO iz = 0, nz+1
        DO ix = 0, nx+1
          ! Variables carried out of domain by Riemann invariants
          vxbc(1) =  vy(ix,ny,iz)
          vybc(1) =  vx(ix,ny,iz)
          vzbc(1) =  vz(ix,ny,iz)
          bxbc(1) =  by(ix,ny,iz)
          bybc(1) =  bx(ix,ny,iz)
          bzbc(1) =  bz(ix,ny,iz)
          rbc (1) = rho(ix,ny,iz)
          ebc (1) = energy(ix,ny,iz)

          pbc  = (gamma - 1.0_num) * energy(ix,ny  ,iz) * rho(ix,ny  ,iz)

          ! Farfield values carried into domain
          pfar = (gamma - 1.0_num) * energy(ix,ny+2,iz) * rho(ix,ny+2,iz)

          uxfar  =  vy(ix,ny+2,iz)
          uyfar  =  vx(ix,ny+2,iz)
          uzfar  =  vz(ix,ny+2,iz)
          bxfar  =  by(ix,ny+2,iz)
          byfar  =  bx(ix,ny+2,iz)
          bzfar  =  bz(ix,ny+2,iz)
          rhofar = rho(ix,ny+2,iz)
          efar = energy(ix,ny+2,iz)

          vnorm = vy(ix,ny,iz)

          ! Select correct open bc solver
          bperp = SQRT(byfar**2 + bzfar**2)

          IF (bperp <= none_zero) THEN
            CALL open_bcs_alfven
          ELSE
            CALL open_bcs_fast
          END IF

          bx (ix,ny+1,iz) = bybc(0)
          by (ix,ny+1,iz) = bxbc(0)
          bz (ix,ny+1,iz) = bzbc(0)
          rho(ix,ny+1,iz) =  rbc(0)
          energy(ix,ny+1,iz) = ebc(0)
        END DO
      END DO
    END IF

    ! z_min boundary
    IF (zbc_min == BC_OPEN .AND. proc_z_min == MPI_PROC_NULL) THEN
      DO iy = 0, ny+1
        DO ix = 0, nx+1
          ! Variables carried out of domain by Riemann invariants
          vxbc(1) = -vz(ix,iy,0)
          vybc(1) =  vy(ix,iy,0)
          vzbc(1) =  vx(ix,iy,0)
          bxbc(1) = -bz(ix,iy,0)
          bybc(1) =  by(ix,iy,1)
          bzbc(1) =  bx(ix,iy,1)
          rbc (1) = rho(ix,iy,1)
          rbc (0) = rho(ix,iy,0)
          ebc (1) = energy(ix,iy,1)

          pbc  = (gamma - 1.0_num) * energy(ix,iy, 1) * rho(ix,iy, 1)

          ! Farfield values carried into domain
          pfar = (gamma - 1.0_num) * energy(ix,iy,-1) * rho(ix,iy,-1)

          uxfar  = -vz(ix,iy,-2)
          uyfar  =  vy(ix,iy,-2)
          uzfar  =  vx(ix,iy,-2)
          bxfar  = -bz(ix,iy,-2)
          byfar  =  by(ix,iy,-1)
          bzfar  =  bx(ix,iy,-1)
          rhofar = rho(ix,iy,-1)
          efar = energy(ix,iy,-1)

          vnorm = -vz(ix,iy,0)

          ! Select correct open bc solver
          bperp = SQRT(byfar**2 + bzfar**2)

          IF (bperp <= none_zero) THEN
            CALL open_bcs_alfven
          ELSE
            CALL open_bcs_fast
          END IF

          bx (ix,iy, 0) =  bzbc(0)
          by (ix,iy, 0) =  bybc(0)
          bz (ix,iy,-1) = -bxbc(0)
          rho(ix,iy, 0) =   rbc(0)
          energy(ix,iy,0) = ebc(0)
        END DO
      END DO
    END IF

    ! z_max boundary
    IF (zbc_max == BC_OPEN .AND. proc_z_max == MPI_PROC_NULL) THEN
      DO iy = 0, ny+1
        DO ix = 0, nx+1
          ! Variables carried out of domain by Riemann invariants
          vxbc(1) = vz(ix, iy, nz)
          vybc(1) = vy(ix, iy, nz)
          vzbc(1) = vx(ix, iy, nz)
          bxbc(1) = bz(ix, iy, nz)
          bybc(1) = by(ix, iy, nz)
          bzbc(1) = bx(ix, iy, nz)
          rbc(1) = rho(ix, iy, nz)
          ebc(1) = energy(ix, iy, nz)

          pbc  = (gamma - 1.0_num) * energy(ix,iy,nz  ) * rho(ix,iy,nz  )

          ! Farfield values carried into domain
          pfar = (gamma - 1.0_num) * energy(ix,iy,nz+2) * rho(ix,iy,nz+2)

          uxfar  =  vz(ix,iy,nz+2)
          uyfar  =  vy(ix,iy,nz+2)
          uzfar  =  vx(ix,iy,nz+2)
          bxfar  =  bz(ix,iy,nz+2)
          byfar  =  by(ix,iy,nz+2)
          bzfar  =  bx(ix,iy,nz+2)
          rhofar = rho(ix,iy,nz+2)
          efar = energy(ix,iy,nz+2)

          vnorm = vz(ix,iy,nz)

          ! Select correct open bc solver
          bperp = SQRT(byfar**2 + bzfar**2)

          IF (bperp <= none_zero) THEN
            CALL open_bcs_alfven
          ELSE
            CALL open_bcs_fast
          END IF

          bx (ix,iy,nz+1) = bzbc(0)
          by (ix,iy,nz+1) = bybc(0)
          bz (ix,iy,nz+1) = bxbc(0)
          rho(ix,iy,nz+1) =  rbc(0)
          energy(ix,iy,nz+1) = ebc(0)
        END DO
      END DO
    END IF

  END SUBROUTINE open_bcs



  SUBROUTINE open_bcs_fast

    ! Open bc when bx = 0

    REAL(num) :: c0, ct, cf
    REAL(num) :: pg, rhog, cffar, c0far, ctfar
    REAL(num), DIMENSION(3) :: vtest, pstar, vstar, rhostar, pmagstar
    REAL(num), DIMENSION(3) :: bystar, bzstar
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
      IF (vtest(i) >= 0.0_num) THEN
        pstar(i) =  pbc + 0.5_num * (bybc(1)**2 + bzbc(1)**2)
        pmagstar(i) =  0.5_num * (bybc(1)**2 + bzbc(1)**2)
        vstar(i) = vxbc(1)
        rhostar(i) = rbc(1)
        bystar(i) = bybc(1)
        bzstar(i) = bzbc(1)
      ELSE
        pstar(i) = pfar + 0.5_num * (byfar**2 + bzfar**2)
        pmagstar(i) = 0.5_num * (byfar**2 + bzfar**2)
        vstar(i) = uxfar
        rhostar(i) = rhofar
        bystar(i) = byfar
        bzstar(i) = bzfar
      END IF
    END DO

    bxbc(0) = 0.5_num * (bxbc(1) + bxfar)
    bybc(0) = 0.5_num * (bystar(1) + bystar(2))
    bzbc(0) = 0.5_num * (bzstar(1) + bzstar(2))

    pg = 0.5_num &
        * (pstar(1) + pstar(2) + rhofar * cffar * (vstar(1) - vstar(2))) &
        - 0.5_num * (pmagstar(1) + pmagstar(2))
    pg = MAX(pg, 0.0_num)

    rhog = rhostar(3) + (pg - (pstar(3) - pmagstar(3))) / c0far**2
    rbc(0) = MAX(rhog, none_zero)

    ebc(0) = pg  / (gamma - 1.0_num) / rbc(0)

  END SUBROUTINE open_bcs_fast



  SUBROUTINE open_bcs_alfven

    ! Open bc when bperp = 0

    REAL(num) :: lambdayfar, lambdazfar
    REAL(num) :: c0, cx
    REAL(num) :: pg, rhog, c0far, cxfar
    REAL(num) :: lambdag
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
      IF (vtest(i) >= 0.0_num) THEN
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

    bxbc(0) = 0.5_num * (bxbc(1) + bxfar)
    IF (ABS(bxbc(0)) <= none_zero) bxbc(0) = none_zero

    lambdag = 0.5_num * (lambdaystar(3) + lambdaystar(4) &
        + rhofar * cxfar * (uystar(3) - uystar(4)))
    bybc(0) = -lambdag / bxbc(0)

    lambdag = 0.5_num * (lambdazstar(3) + lambdazstar(4) &
        + rhofar * cxfar * (uzstar(3) - uzstar(4)))
    bzbc(0) = -lambdag  / bxbc(0)

    pg = 0.5_num &
        * (pstar(1) + pstar(2) + rhofar * c0far * (uxstar(1) - uxstar(2))) &
      - 0.5_num * (bybc(0)**2 + bzbc(0)**2 - bxbc(0)**2)
    pg = MAX(pg, none_zero)

    rhog = rhostar(5) + (pg - (pstar(5) - pmagstar(5))) / c0far**2 
    rbc(0) = MAX(rhog, none_zero)
    ebc(0) = pg / (gamma - 1.0_num) / rbc(0)

  END SUBROUTINE open_bcs_alfven



END MODULE openboundary
