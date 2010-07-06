!-------------------------------------------------------------------------
! mass coordinate based Van Leer limited remap.
! See Bram van Leer, JCP, vol 135, p229, (1997)
! Now rewritten to allow compiler vectorizing of loops
! See notes in code
!-------------------------------------------------------------------------
MODULE xremap

  USE shared_data; USE boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: remap_x

  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: rho1, dm, cv2, flux, dxb1
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: rho_v, rho_v1

CONTAINS

  SUBROUTINE remap_x ! remap onto original Eulerian grid

    REAL(num) :: vxb, vxbm, vyb, vybm, vzb, vzbm, dv

    ALLOCATE (rho1(-1:nx+2, -1:ny+2, -1:nz+2), dm(-1:nx+2, -1:ny+2, -1:nz+2), &
        cv2(-1:nx+2, -1:ny+2, -1:nz+2), flux(-2:nx+2, -1:ny+2, -1:nz+2), &
        dxb1(-1:nx+2, -1:ny+2, -1:nz+2), rho_v(-1:nx+2, -1:ny+2, -1:nz+2), &
        rho_v1(-1:nx+2, -1:ny+2, -1:nz+2))

    dm = 0.0_num
    rho1 = rho ! store initial density in rho1

    DO iz = -1, nz+2
      izm = iz - 1
      DO iy = -1, ny+2
        iym = iy - 1
        DO ix = -1, nx+2
          ixm = ix - 1

          ! vx at Sx(i, j, k)
          vxb = (vx1(ix, iy, iz) + vx1(ix, iym, iz) &
              + vx1(ix, iy, izm) + vx1(ix, iym, izm)) / 4.0_num

          ! vx at Sx(i-1, j, k)
          vxbm = (vx1(ixm, iy, iz) + vx1(ixm, iym, iz) &
              + vx1(ixm, iy, izm) + vx1(ixm, iym, izm)) / 4.0_num

          ! vy at Sy(i, j, k)
          vyb = (vy1(ix, iy, iz) + vy1(ixm, iy, iz) &
              + vy1(ix, iy, izm) + vy1(ixm, iy, izm)) / 4.0_num

          ! vy at Sy(i, j-1, k)
          vybm = (vy1(ix, iym, iz) + vy1(ixm, iym, iz) &
              + vy1(ix, iym, izm) + vy1(ixm, iym, izm)) / 4.0_num

          ! vz at Sz(i, j, k)
          vzb = (vz1(ix, iy, iz) + vz1(ixm, iy, iz) &
              + vz1(ix, iym, iz) + vz1(ixm, iym, iz)) / 4.0_num

          ! vz at Sz(i, j, k-1)
          vzbm = (vz1(ix, iy, izm) + vz1(ixm, iy, izm) &
              + vz1(ix, iym, izm) + vz1(ixm, iym, izm)) / 4.0_num

          dv = (REAL(ypass, num) * (vyb - vybm) / dyb(iy) &
              + REAL(zpass, num) * (vzb - vzbm) / dzb(iz) &
              + (vxb - vxbm) / dxb(ix)) * dt

          ! control volume before remap
          cv1(ix, iy, iz) = cv(ix, iy, iz) * (1.0_num + dv)

          dv = (REAL(ypass, num) * (vyb - vybm) / dyb(iy) &
              + REAL(zpass, num) * (vzb - vzbm) / dzb(iz)) * dt

          ! control volume after remap
          cv2(ix, iy, iz) = cv(ix, iy, iz) * (1.0_num + dv)

          ! dxb before remap
          dxb1(ix, iy, iz) = dxb(ix) + (vxb - vxbm) * dt
        END DO
      END DO
    END DO

    ! Evans and Hawley (ApJ, vol 332, p650, (1988))
    ! constrained transport remap of magnetic fluxes
    CALL vx_by_flux
    DO iz = 1, nz
      DO iy = 0, ny
        DO ix = 1, nx
          ixm = ix - 1
          by(ix, iy, iz) = by(ix, iy, iz) - flux(ix, iy, iz) + flux(ixm, iy, iz)
        END DO
      END DO
    END DO

    DO iz = 1, nz
      DO iy = 1, ny
        iym = iy - 1
        DO ix = 0, nx
          bx(ix, iy, iz) = bx(ix, iy, iz) + flux(ix, iy, iz) - flux(ix, iym, iz)
        END DO
      END DO
    END DO

    CALL vx_bz_flux
    DO iz = 0, nz
      DO iy = 1, ny
        DO ix = 1, nx
          ixm = ix - 1
          bz(ix, iy, iz) = bz(ix, iy, iz) - flux(ix, iy, iz) + flux(ixm, iy, iz)
        END DO
      END DO
    END DO

    DO iz = 1, nz
      izm = iz - 1
      DO iy = 1, ny
        DO ix = 0, nx
          bx(ix, iy, iz) = bx(ix, iy, iz) + flux(ix, iy, iz) - flux(ix, iy, izm)
        END DO
      END DO
    END DO

    ! remap of mass + calculation of mass fluxes (dm) needed for later remaps
    CALL x_mass_flux ! calculates dm(0:nx, 0:ny+1)
    CALL dm_x_bcs    ! need dm(0:nx+1, 0:ny+1) for velocity remap
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          ixm = ix - 1
          rho(ix, iy, iz) = (rho1(ix, iy, iz) * cv1(ix, iy, iz) &
              + dm(ixm, iy, iz) - dm(ix, iy, iz)) / cv2(ix, iy, iz)
        END DO
      END DO
    END DO

    ! remap specific energy density using mass coordinates
    CALL x_energy_flux
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          ixm = ix - 1
          energy(ix, iy, iz) = (energy(ix, iy, iz) * cv1(ix, iy, iz) &
              * rho1(ix, iy, iz) + flux(ixm, iy, iz) - flux(ix, iy, iz)) &
              / (cv2(ix, iy, iz) * rho(ix, iy, iz))
        END DO
      END DO
    END DO

    ! redefine dxb1, cv1, cv2, dm and vx1 for velocity (vertex) cells
    ! in some of these calculations the flux variable is used as a
    ! temporary array
    DO iz = 0, nz
      izp = iz + 1
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = -1, nx+1
          ixp = ix + 1

          ! vertex density before remap
          rho_v(ix, iy, iz) = rho1(ix, iy, iz) * cv1(ix, iy, iz) &
              + rho1(ixp, iy , iz ) * cv1(ixp, iy , iz ) &
              + rho1(ix , iyp, iz ) * cv1(ix , iyp, iz ) &
              + rho1(ixp, iyp, iz ) * cv1(ixp, iyp, iz ) &
              + rho1(ix , iy , izp) * cv1(ix , iy , izp) &
              + rho1(ixp, iy , izp) * cv1(ixp, iy , izp) &
              + rho1(ix , iyp, izp) * cv1(ix , iyp, izp) &
              + rho1(ixp, iyp, izp) * cv1(ixp, iyp, izp)

          rho_v(ix, iy, iz) = rho_v(ix, iy, iz) &
              / (cv1(ix, iy , iz ) + cv1(ixp, iy , iz ) &
              + cv1(ix, iyp, iz ) + cv1(ixp, iyp, iz ) &
              + cv1(ix, iy , izp) + cv1(ixp, iy , izp) &
              + cv1(ix, iyp, izp) + cv1(ixp, iyp, izp))
        END DO
      END DO
    END DO

    DO iz = 0, nz
      izp = iz + 1
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = 0, nx
          ixp = ix + 1
          flux(ix, iy, iz) = cv1(ix, iy, iz) + cv1(ixp, iy, iz) &
              + cv1(ix, iyp, iz ) + cv1(ixp, iyp, iz ) &
              + cv1(ix, iy , izp) + cv1(ixp, iy , izp) &
              + cv1(ix, iyp, izp) + cv1(ixp, iyp, izp)
        END DO
      END DO
    END DO
    ! cv1 = vertex CV before remap
    cv1(0:nx, 0:ny, 0:nz) = flux(0:nx, 0:ny, 0:nz) / 8.0_num

    DO iz = 0, nz
      izp = iz + 1
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = 0, nx
          ixp = ix + 1
          flux(ix, iy, iz) = cv2(ix, iy, iz) + cv2(ixp, iy, iz) &
              + cv2(ix, iyp, iz ) + cv2(ixp, iyp, iz ) &
              + cv2(ix, iy , izp) + cv2(ixp, iy , izp) &
              + cv2(ix, iyp, izp) + cv2(ixp, iyp, izp)
        END DO
      END DO
    END DO
    ! cv2 = vertex CV after remap
    cv2(0:nx, 0:ny, 0:nz) = flux(0:nx, 0:ny, 0:nz) / 8.0_num

    DO iz = 0, nz
      DO iy = 0, ny
        DO ix = -2, nx+1
          ixp = ix + 1
          flux(ix, iy, iz) = (vx1(ix, iy, iz) + vx1(ixp, iy, iz)) / 2.0_num
        END DO
      END DO
    END DO
    ! vertex boundary velocity used in remap
    vx1(-2:nx+1, 0:ny, 0:nz) = flux(-2:nx+1, 0:ny, 0:nz)

    DO iz = 0, nz
      DO iy = 0, ny
        DO ix = -1, nx+1
          ixm = ix - 1
          ! dxb1 = width of vertex CV before remap
          dxb1(ix, iy, iz) = dxc(ix) + (vx1(ix, iy, iz) - vx1(ixm, iy, iz)) * dt
        END DO
      END DO
    END DO

    DO iz = 0, nz
      izp = iz + 1
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = -1, nx
          ixp = ix + 1
          flux(ix, iy, iz) = dm(ix, iy, iz) + dm(ixp, iy, iz) &
              + dm(ix, iyp, iz ) + dm(ixp, iyp, iz ) &
              + dm(ix, iy , izp) + dm(ixp, iy , izp) &
              + dm(ix, iyp, izp) + dm(ixp, iyp, izp)
        END DO
      END DO
    END DO
    ! mass flux out of vertex CV
    dm(-1:nx, 0:ny, 0:nz) = flux(-1:nx, 0:ny, 0:nz) / 8.0_num

    DO iz = 0, nz
      DO iy = 0, ny
        DO ix = 0, nx
          ixm = ix - 1
          ! vertex density after remap
          rho_v1(ix, iy, iz) = (rho_v(ix, iy, iz) * cv1(ix, iy, iz) &
              + dm(ixm, iy, iz) - dm(ix, iy, iz)) / cv2(ix, iy, iz)
        END DO
      END DO
    END DO

    CALL x_momy_flux
    DO iz = 0, nz
      DO iy = 0, ny
        DO ix = 0, nx
          ixm = ix - 1
          vy(ix, iy, iz) = (rho_v(ix, iy, iz) * vy(ix, iy, iz) &
              * cv1(ix, iy, iz) + flux(ixm, iy, iz) - flux(ix, iy, iz)) &
              / (cv2(ix, iy, iz) * rho_v1(ix, iy, iz))
        END DO
      END DO
    END DO

    CALL x_momz_flux
    DO iz = 0, nz
      DO iy = 0, ny
        DO ix = 0, nx
          ixm = ix - 1
          vz(ix, iy, iz) = (rho_v(ix, iy, iz) * vz(ix, iy, iz) &
              * cv1(ix, iy, iz) + flux(ixm, iy, iz) - flux(ix, iy, iz)) &
              / (cv2(ix, iy, iz) * rho_v1(ix, iy, iz))
        END DO
      END DO
    END DO

    CALL x_momx_flux
    DO iz = 0, nz
      DO iy = 0, ny
        DO ix = 0, nx
          ixm = ix - 1
          vx(ix, iy, iz) = (rho_v(ix, iy, iz) * vx(ix, iy, iz) &
              * cv1(ix, iy, iz) + flux(ixm, iy, iz) - flux(ix, iy, iz)) &
              / (cv2(ix, iy, iz) * rho_v1(ix, iy, iz))
        END DO
      END DO
    END DO

    CALL boundary_conditions

    DEALLOCATE (rho1, dm, cv2, flux, dxb1, rho_v, rho_v1)
    xpass = 0

  END SUBROUTINE remap_x



  SUBROUTINE vx_by_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: db, dbxp, dbxp2, dbxm
    INTEGER :: ixp2

    DO iz = 0, nz
      izm  = iz - 1
      DO iy = 0, ny
        iyp  = iy + 1
        DO ix = 0, nx
          ixm  = ix - 1
          ixp  = ix + 1
          ixp2 = ix + 2

          v_advect = (vx1(ix, iy, iz) + vx1(ix, iy, izm)) / 2.0_num

          db    = (dxb1(ix  , iy, iz) + dxb1(ix  , iyp, iz)) / 2.0_num
          dbxp  = (dxb1(ixp , iy, iz) + dxb1(ixp , iyp, iz)) / 2.0_num
          dbxp2 = (dxb1(ixp2, iy, iz) + dxb1(ixp2, iyp, iz)) / 2.0_num
          dbxm  = (dxb1(ixm , iy, iz) + dxb1(ixm , iyp, iz)) / 2.0_num

          w4 = by(ix , iy, iz) / db
          w5 = by(ixp, iy, iz) / dbxp

          flux(ix, iy, iz) = (MAX(0.0_num, v_advect) * w4 &
              + MIN(0.0_num, v_advect) * w5) * dt

          w1 = by(ixp , iy, iz) / dbxp  - by(ix , iy, iz) / db
          w2 = by(ix  , iy, iz) / db    - by(ixm, iy, iz) / dbxm
          w3 = by(ixp2, iy, iz) / dbxp2 - by(ixp, iy, iz) / dbxp

          ! vad_p and vad_m are logical switches which determine v_advect>=0
          ! and v_advect<0 respectively. It's written this way to allow vector
          ! optimization

          vad_p = -MIN(SIGN(1.0_num, -v_advect), 0.0_num)
          vad_m =  MAX(SIGN(1.0_num, -v_advect), 0.0_num)

          w5 = ABS(v_advect) * dt / (db * vad_p + dbxp * vad_m)
          w4 = (2.0_num - w5) * ABS(w1) / dxc(ix) &
              + (1.0_num + w5) * ABS(w2 * vad_p + w3 * vad_m) &
              / (dxc(ixm) * vad_p + dxc(ixp) * vad_m)

          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) &
              + SIGN(1.0_num, w2 * vad_p + w3 * vad_m))

          w6 = SIGN(1.0_num, v_advect) * w8 &
              * MIN(ABS(w4) * (dxb(ix) * vad_p + dxb(ixp) * vad_m), &
              ABS(w1), ABS(w2 * vad_p + w3 * vad_m))

          flux(ix, iy, iz) = flux(ix, iy, iz) &
              + v_advect * dt * w6 * (1.0_num - w5)
        END DO
      END DO
    END DO

  END SUBROUTINE vx_by_flux



  SUBROUTINE vx_bz_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: db, dbxp, dbxp2, dbxm
    INTEGER :: ixp2

    DO iz = 0, nz
      izp  = iz + 1
      DO iy = 0, ny
        iym  = iy - 1
        DO ix = 0, nx
          ixm  = ix - 1
          ixp  = ix + 1
          ixp2 = ix + 2

          v_advect = (vx1(ix, iy, iz) + vx1(ix, iym, iz)) / 2.0_num

          db    = (dxb1(ix  , iy, iz) + dxb1(ix  , iy, izp)) / 2.0_num
          dbxp  = (dxb1(ixp , iy, iz) + dxb1(ixp , iy, izp)) / 2.0_num
          dbxp2 = (dxb1(ixp2, iy, iz) + dxb1(ixp2, iy, izp)) / 2.0_num
          dbxm  = (dxb1(ixm , iy, iz) + dxb1(ixm , iy, izp)) / 2.0_num

          w4 = bz(ix , iy, iz) / db
          w5 = bz(ixp, iy, iz) / dbxp

          flux(ix, iy, iz) = (MAX(0.0_num, v_advect) * w4 &
              + MIN(0.0_num, v_advect) * w5) * dt

          w1 = bz(ixp , iy, iz) / dbxp  - bz(ix , iy, iz) / db
          w2 = bz(ix  , iy, iz) / db    - bz(ixm, iy, iz) / dbxm
          w3 = bz(ixp2, iy, iz) / dbxp2 - bz(ixp, iy, iz) / dbxp

          ! vad_p and vad_m are logical switches which determine v_advect>=0
          ! and v_advect<0 respectively. It's written this way to allow vector
          ! optimization. See example in **SECTION 2**

          vad_p = -MIN(SIGN(1.0_num, -v_advect), 0.0_num)
          vad_m =  MAX(SIGN(1.0_num, -v_advect), 0.0_num)

          ! This code is the vectorizable replacement for the code
          ! in **SECTION 2**
          w5 = ABS(v_advect) * dt / (db * vad_p + dbxp * vad_m)
          w4 = (2.0_num - w5) * ABS(w1) / dxc(ix) &
              + (1.0_num + w5) * ABS(w2 * vad_p + w3 * vad_m) &
              / (dxc(ixm) * vad_p + dxc(ixp) * vad_m)

          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) &
              + SIGN(1.0_num, w2 * vad_p + w3 * vad_m))

          w6 = SIGN(1.0_num, v_advect) * w8 &
              * MIN(ABS(w4) * (dxb(ix) * vad_p + dxb(ixp) * vad_m), &
              ABS(w1), ABS(w2 * vad_p + w3 * vad_m))

          flux(ix, iy, iz) = flux(ix, iy, iz) &
              + v_advect * dt * w6 * (1.0_num - w5)

!!$          !**SECTION 2**
!!$          IF (v_advect > 0.0) THEN
!!$            w5 = ABS(v_advect) * dt / db
!!$            w4 = (2.0_num - w5) * ABS(w1) / dxc(ix) &
!!$                + (1.0_num + w5) * ABS(w2) / dxc(ixm)
!!$            w4 = w4 / 6.0_num
!!$            w4 = ABS(w1) / dxc(ix)
!!$            w4 = w4 / 2.0_num
!!$            w8 = 0.5_num * (SIGN(1.0_num, w1) + SIGN(1.0_num, w2))
!!$            w6 = w8 * MIN(ABS(w4)*dxb(ix), ABS(w1), ABS(w2))
!!$            flux(ix, iy, iz) = flux(ix, iy, iz) &
!!$                + v_advect * dt * w6 * (1.0_num - w5)
!!$          ELSE
!!$            w5 = ABS(v_advect) * dt / dbxp
!!$            w4 = (2.0_num - w5) * ABS(w1) / dxc(ix) &
!!$                + (1.0_num + w5) * ABS(w3) / dxc(ixp)
!!$            w4 = w4 / 6.0_num
!!$            w8 = 0.5_num * (SIGN(1.0_num, w1) + SIGN(1.0_num, w3))
!!$            w6 = -w8 * MIN(ABS(w4)*dxb(ixp), ABS(w1), ABS(w3))
!!$            flux(ix, iy, iz) = flux(ix, iy, iz) &
!!$                + v_advect * dt * w6 * (1.0_num - w5)
!!$          END IF
!!$          !**END SECTION 2**

        END DO
      END DO
    END DO

  END SUBROUTINE vx_bz_flux



  SUBROUTINE x_mass_flux

    REAL(num) :: v_advect, flux_rho, vad_p, vad_m
    INTEGER :: ixp2

    DO iz = 0, nz+1
      izm  = iz - 1
      DO iy = 0, ny+1
        iym  = iy - 1
        DO ix = 0, nx
          ixm  = ix - 1
          ixp  = ix + 1
          ixp2 = ix + 2

          v_advect = (vx1(ix, iy, iz) + vx1(ix, iym, iz) &
              + vx1(ix, iy, izm) + vx1(ix, iym, izm)) / 4.0_num

          dm(ix, iy, iz) = (MAX(0.0_num, v_advect) * rho(ix, iy, iz) &
              + MIN(0.0_num, v_advect) * rho(ixp, iy, iz)) * dt

          w1 = rho(ixp , iy, iz) - rho(ix , iy, iz)
          w2 = rho(ix  , iy, iz) - rho(ixm, iy, iz)
          w3 = rho(ixp2, iy, iz) - rho(ixp, iy, iz)

          ! vad_p and vad_m are logical switches which determine v_advect>=0
          ! and v_advect<0 respectively. It's written this way to allow vector
          ! optimization

          vad_p = -MIN(SIGN(1.0_num, -v_advect), 0.0_num)
          vad_m =  MAX(SIGN(1.0_num, -v_advect), 0.0_num)

          w5 = ABS(v_advect) * dt &
              / (dxb1(ix, iy, iz) * vad_p + dxb1(ixp, iy, iz) * vad_m)

          w4 = (2.0_num - w5) * ABS(w1) / dxc(ix) &
              + (1.0_num + w5) * ABS(w2 * vad_p + w3 * vad_m) &
              / (dxc(ixm) * vad_p + dxc(ixp) * vad_m)

          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) &
              + SIGN(1.0_num, w2 * vad_p + w3 * vad_m))

          w6 = SIGN(1.0_num, v_advect) * w8 &
              * MIN(ABS(w4) * (dxb(ix) * vad_p + dxb(ixp) * vad_m), &
              ABS(w1), ABS(w2 * vad_p + w3 * vad_m))

          flux_rho = v_advect * dt * w6 * (1.0_num - w5)

          dm(ix, iy, iz) = (flux_rho + dm(ix, iy, iz)) * dyb(iy) * dzb(iz)
        END DO
      END DO
    END DO

  END SUBROUTINE x_mass_flux



  SUBROUTINE x_energy_flux ! energy remap in mass coordinates

    REAL(num) :: v_advect, vad_p, vad_m
    INTEGER :: ixp2

    DO iz = 0, nz
      izm  = iz - 1
      DO iy = 0, ny
        iym  = iy - 1
        DO ix = 0, nx
          ixm  = ix - 1
          ixp  = ix + 1
          ixp2 = ix + 2

          v_advect = (vx1(ix, iy, iz) + vx1(ix, iym, iz) &
              + vx1(ix, iy, izm) + vx1(ix, iym, izm)) / 4.0_num

          w1 = energy(ixp , iy, iz) - energy(ix , iy, iz)
          w2 = energy(ix  , iy, iz) - energy(ixm, iy, iz)
          w3 = energy(ixp2, iy, iz) - energy(ixp, iy, iz)

          ! vad_p and vad_m are logical switches which determine v_advect>=0
          ! and v_advect<0 respectively. It's written this way to allow vector
          ! optimization

          vad_p = -MIN(SIGN(1.0_num, -v_advect), 0.0_num)
          vad_m =  MAX(SIGN(1.0_num, -v_advect), 0.0_num)

          w5 = ABS(v_advect) * dt &
              / (dxb1(ix, iy, iz) * vad_p + dxb1(ixp, iy, iz) * vad_m)

          w7 = energy(ix, iy, iz) * vad_p + energy(ixp, iy, iz) * vad_m

          w6 = ABS(dm(ix, iy, iz)) / dyb(iy) / dzb(iz) &
              / (rho1(ix , iy, iz) * dxb1(ix , iy, iz) * vad_p &
              +  rho1(ixp, iy, iz) * dxb1(ixp, iy, iz) * vad_m)

          w4 = (2.0_num - w5) * ABS(w1) / dxc(ix) &
              + (1.0_num + w5) * ABS(w2 * vad_p + w3 * vad_m) &
              / (dxc(ixm) * vad_p + dxc(ixp) * vad_m)

          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) &
              + SIGN(1.0_num, w2 * vad_p + w3 * vad_m))

          w5 = SIGN(1.0_num, v_advect) * w8 &
              * MIN(ABS(w4) * (dxb(ix) * vad_p + dxb(ixp) * vad_m), &
              ABS(w1), ABS(w2 * vad_p + w3 * vad_m))

          flux(ix, iy, iz) = dm(ix, iy, iz) * (w7 + w5 * (1.0_num - w6))
        END DO
      END DO
    END DO

  END SUBROUTINE x_energy_flux



  SUBROUTINE x_momy_flux ! energy remap in mass coordinates

    REAL(num) :: v_advect, m, mp, ai, aip, dk
    REAL(num) :: vad_p, vad_m
    INTEGER :: ixp2

    DO iz = 0, nz
      DO iy = 0, ny
        DO ix = -1, nx
          ixm  = ix - 1
          ixp  = ix + 1
          ixp2 = ix + 2

          v_advect = vx1(ix, iy, iz)

          w1 = vy(ixp , iy, iz) - vy(ix , iy, iz)
          w2 = vy(ix  , iy, iz) - vy(ixm, iy, iz)
          w3 = vy(ixp2, iy, iz) - vy(ixp, iy, iz)

          ! vad_p and vad_m are logical switches which determine v_advect>=0
          ! and v_advect<0 respectively. It's written this way to allow vector
          ! optimization

          vad_p = -MIN(SIGN(1.0_num, -v_advect), 0.0_num)
          vad_m =  MAX(SIGN(1.0_num, -v_advect), 0.0_num)

          w5 = ABS(v_advect) * dt &
              / (dxb1(ix, iy, iz) * vad_p + dxb1(ixp, iy, iz) * vad_m)

          w7 = vy(ix, iy, iz) * vad_p + vy(ixp, iy, iz) * vad_m

          w6 = ABS(dm(ix, iy, iz)) / dyc(iy) / dzc(iz) &
              / (rho_v(ix , iy, iz) * dxb1(ix , iy, iz) * vad_p &
              +  rho_v(ixp, iy, iz) * dxb1(ixp, iy, iz) * vad_m)

          w4 = (2.0_num - w5) * ABS(w1) / dxb(ixp) &
              + (1.0_num + w5) * ABS(w2 * vad_p + w3 * vad_m) &
              / (dxb(ix) * vad_p + dxb(ixp2) * vad_m)

          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) &
              + SIGN(1.0_num, w2 * vad_p + w3 * vad_m))

          w5 = SIGN(1.0_num, v_advect) * w8 &
              * MIN(ABS(w4) * (dxc(ix) * vad_p + dxc(ixp) * vad_m), &
              ABS(w1), ABS(w2 * vad_p + w3 * vad_m))

          flux(ix, iy, iz) = w7 + w5 * (1.0_num - w6)
        END DO
      END DO
    END DO

    IF (rke) THEN
      DO iz = 0, nz
        izp = iz + 1
        DO iy = 0, ny
          iyp = iy + 1
          DO ix = 0, nx-1
            ixm = ix - 1
            ixp = ix + 1

            m = rho_v1(ix, iy, iz) * cv2(ix, iy, iz)
            mp = rho_v1(ixp, iy, iz) * cv2(ixp, iy, iz)

            ai = (vy(ix, iy, iz) - flux(ixm, iy, iz)) * dm(ixm, iy, iz) / m &
                + (flux(ix, iy, iz) - vy(ix, iy, iz)) * dm(ix, iy, iz) / m

            aip = (vy(ixp, iy, iz) - flux(ix, iy, iz)) * dm(ix, iy, iz) / mp &
                + (flux(ixp, iy, iz) - vy(ixp, iy, iz)) * dm(ixp, iy, iz) / mp

            dk = (vy(ixp, iy, iz) - vy(ix, iy, iz)) * (flux(ix, iy, iz) &
                - 0.5_num * (vy(ixp, iy, iz) + vy(ix, iy, iz))) &
                + 0.5_num * ai * (flux(ix, iy, iz) - vy(ix, iy, iz)) &
                + 0.5_num * aip * (vy(ixp, iy, iz) - flux(ix, iy, iz))

            dk = dk * dm(ix, iy, iz) / 4.0_num
            delta_ke(ixp, iy , iz ) = delta_ke(ixp, iy , iz ) + dk
            delta_ke(ixp, iyp, iz ) = delta_ke(ixp, iyp, iz ) + dk
            delta_ke(ixp, iy , izp) = delta_ke(ixp, iy , izp) + dk
            delta_ke(ixp, iyp, izp) = delta_ke(ixp, iyp, izp) + dk
          END DO
        END DO
      END DO
    END IF

    flux(-1:nx, 0:ny, 0:nz) = flux(-1:nx, 0:ny, 0:nz) * dm(-1:nx, 0:ny, 0:nz)

  END SUBROUTINE x_momy_flux



  SUBROUTINE x_momz_flux ! energy remap in mass coordinates

    REAL(num) :: v_advect, m, mp, ai, aip, dk
    REAL(num) :: vad_p, vad_m
    INTEGER :: ixp2

    DO iz = 0, nz
      DO iy = 0, ny
        DO ix = -1, nx
          ixm  = ix - 1
          ixp  = ix + 1
          ixp2 = ix + 2

          v_advect = vx1(ix, iy, iz)

          w1 = vz(ixp , iy, iz) - vz(ix , iy, iz)
          w2 = vz(ix  , iy, iz) - vz(ixm, iy, iz)
          w3 = vz(ixp2, iy, iz) - vz(ixp, iy, iz)

          ! vad_p and vad_m are logical switches which determine v_advect>=0
          ! and v_advect<0 respectively. It's written this way to allow vector
          ! optimization

          vad_p = -MIN(SIGN(1.0_num, -v_advect), 0.0_num)
          vad_m =  MAX(SIGN(1.0_num, -v_advect), 0.0_num)

          w5 = ABS(v_advect) * dt &
              / (dxb1(ix, iy, iz) * vad_p + dxb1(ixp, iy, iz) * vad_m)

          w7 = vz(ix, iy, iz) * vad_p + vz(ixp, iy, iz) * vad_m

          w6 = ABS(dm(ix, iy, iz)) / dyc(iy) / dzc(iz) &
              / (rho_v(ix , iy, iz) * dxb1(ix , iy, iz) * vad_p &
              +  rho_v(ixp, iy, iz) * dxb1(ixp, iy, iz) * vad_m)

          w4 = (2.0_num - w5) * ABS(w1) / dxb(ixp) &
              + (1.0_num + w5) * ABS(w2 * vad_p + w3 * vad_m) &
              / (dxb(ix) * vad_p + dxb(ixp2) * vad_m)

          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) &
              + SIGN(1.0_num, w2 * vad_p + w3 * vad_m))

          w5 = SIGN(1.0_num, v_advect) * w8 &
              * MIN(ABS(w4) * (dxc(ix) * vad_p + dxc(ixp) * vad_m), &
              ABS(w1), ABS(w2 * vad_p + w3 * vad_m))

          flux(ix, iy, iz) = w7 + w5 * (1.0_num - w6)
        END DO
      END DO
    END DO

    IF (rke) THEN
      DO iz = 0, nz
        izp = iz + 1
        DO iy = 0, ny
          iyp = iy + 1
          DO ix = 0, nx-1
            ixm = ix - 1
            ixp = ix + 1

            m = rho_v1(ix, iy, iz) * cv2(ix, iy, iz)
            mp = rho_v1(ixp, iy, iz) * cv2(ixp, iy, iz)

            ai = (vz(ix, iy, iz) - flux(ixm, iy, iz)) * dm(ixm, iy, iz) / m &
                + (flux(ix, iy, iz) - vz(ix, iy, iz)) * dm(ix, iy, iz) / m

            aip = (vz(ixp, iy, iz) - flux(ix, iy, iz)) * dm(ix, iy, iz) / mp &
                + (flux(ixp, iy, iz) - vz(ixp, iy, iz)) * dm(ixp, iy, iz) / mp

            dk = (vz(ixp, iy, iz) - vz(ix, iy, iz)) * (flux(ix, iy, iz) &
                - 0.5_num * (vz(ixp, iy, iz) + vz(ix, iy, iz))) &
                + 0.5_num * ai * (flux(ix, iy, iz) - vz(ix, iy, iz)) &
                + 0.5_num * aip * (vz(ixp, iy, iz) - flux(ix, iy, iz))

            dk = dk * dm(ix, iy, iz) / 4.0_num
            delta_ke(ixp, iy , iz ) = delta_ke(ixp, iy , iz ) + dk
            delta_ke(ixp, iyp, iz ) = delta_ke(ixp, iyp, iz ) + dk
            delta_ke(ixp, iy , izp) = delta_ke(ixp, iy , izp) + dk
            delta_ke(ixp, iyp, izp) = delta_ke(ixp, iyp, izp) + dk
          END DO
        END DO
      END DO
    END IF

    flux(-1:nx, 0:ny, 0:nz) = flux(-1:nx, 0:ny, 0:nz) * dm(-1:nx, 0:ny, 0:nz)

  END SUBROUTINE x_momz_flux



  SUBROUTINE x_momx_flux

    REAL(num) :: v_advect, m, mp, ai, aip, dk
    REAL(num) :: vad_p, vad_m
    INTEGER :: ixp2

    DO iz = 0, nz
      DO iy = 0, ny
        DO ix = -1, nx
          ixm  = ix - 1
          ixp  = ix + 1
          ixp2 = ix + 2

          v_advect = vx1(ix, iy, iz)

          w1 = vx(ixp , iy, iz) - vx(ix , iy, iz)
          w2 = vx(ix  , iy, iz) - vx(ixm, iy, iz)
          w3 = vx(ixp2, iy, iz) - vx(ixp, iy, iz)

          ! vad_p and vad_m are logical switches which determine v_advect>=0
          ! and v_advect<0 respectively. It's written this way to allow vector
          ! optimization

          vad_p = -MIN(SIGN(1.0_num, -v_advect), 0.0_num)
          vad_m =  MAX(SIGN(1.0_num, -v_advect), 0.0_num)

          w5 = ABS(v_advect) * dt &
              / (dxb1(ix, iy, iz) * vad_p + dxb1(ixp, iy, iz) * vad_m)

          w7 = vx(ix, iy, iz) * vad_p + vx(ixp, iy, iz) * vad_m

          w6 = ABS(dm(ix, iy, iz)) / dyc(iy) / dzc(iz) &
              / (rho_v(ix , iy, iz) * dxb1(ix , iy, iz) * vad_p &
              +  rho_v(ixp, iy, iz) * dxb1(ixp, iy, iz) * vad_m)

          w4 = (2.0_num - w5) * ABS(w1) / dxb(ixp) &
              + (1.0_num + w5) * ABS(w2 * vad_p + w3 * vad_m) &
              / (dxb(ix) * vad_p + dxb(ixp2) * vad_m)

          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) &
              + SIGN(1.0_num, w2 * vad_p + w3 * vad_m))

          w5 = SIGN(1.0_num, v_advect) * w8 &
              * MIN(ABS(w4) * (dxc(ix) * vad_p + dxc(ixp) * vad_m), &
              ABS(w1), ABS(w2 * vad_p + w3 * vad_m))

          flux(ix, iy, iz) = w7 + w5 * (1.0_num - w6)
        END DO
      END DO
    END DO

    IF (rke) THEN
      DO iz = 0, nz
        izp = iz + 1
        DO iy = 0, ny
          iyp = iy + 1
          DO ix = 0, nx-1
            ixm = ix - 1
            ixp = ix + 1

            m = rho_v1(ix, iy, iz) * cv2(ix, iy, iz)
            mp = rho_v1(ixp, iy, iz) * cv2(ixp, iy, iz)

            ai = (vx(ix, iy, iz) - flux(ixm, iy, iz)) * dm(ixm, iy, iz) / m &
                + (flux(ix, iy, iz) - vx(ix, iy, iz)) * dm(ix, iy, iz) / m

            aip = (vx(ixp, iy, iz) - flux(ix, iy, iz)) * dm(ix, iy, iz) / mp &
                + (flux(ixp, iy, iz) - vx(ixp, iy, iz)) * dm(ixp, iy, iz) / mp

            dk = (vx(ixp, iy, iz) - vx(ix, iy, iz)) * (flux(ix, iy, iz) &
                - 0.5_num * (vx(ixp, iy, iz) + vx(ix, iy, iz))) &
                + 0.5_num * ai * (flux(ix, iy, iz) - vx(ix, iy, iz)) &
                + 0.5_num * aip * (vx(ixp, iy, iz) - flux(ix, iy, iz))

            dk = dk * dm(ix, iy, iz) / 4.0_num
            delta_ke(ixp, iy , iz ) = delta_ke(ixp, iy , iz ) + dk
            delta_ke(ixp, iyp, iz ) = delta_ke(ixp, iyp, iz ) + dk
            delta_ke(ixp, iy , izp) = delta_ke(ixp, iy , izp) + dk
            delta_ke(ixp, iyp, izp) = delta_ke(ixp, iyp, izp) + dk
          END DO
        END DO
      END DO
    END IF

    flux(-1:nx, 0:ny, 0:nz) = flux(-1:nx, 0:ny, 0:nz) * dm(-1:nx, 0:ny, 0:nz)

  END SUBROUTINE x_momx_flux



  SUBROUTINE dm_x_bcs

    CALL MPI_SENDRECV(dm(1, 0:ny+1, 0:nz+1), (ny+2)*(nz+2), mpireal, &
        left, tag, dm(nx+1, 0:ny+1, 0:nz+1), (ny+2)*(nz+2), mpireal, &
        right, tag, comm, status, errcode)

    IF (right == MPI_PROC_NULL) &
        dm(nx+1, 0:ny+1, 0:nz+1) = dm(nx, 0:ny+1, 0:nz+1)

    CALL MPI_SENDRECV(dm(nx-1, 0:ny+1, 0:nz+1), (ny+2)*(nz+2), mpireal, &
        right, tag, dm(-1, 0:ny+1, 0:nz+1), (ny+2)*(nz+2), mpireal, &
        left, tag, comm, status, errcode)

    IF (left == MPI_PROC_NULL) &
        dm(-1, 0:ny+1, 0:nz+1) = dm(0, 0:ny+1, 0:nz+1)

  END SUBROUTINE dm_x_bcs

END MODULE xremap
