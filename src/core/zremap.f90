!-------------------------------------------------------------------------
! mass coordinate based Van Leer limited remap.
! See Bram van Leer, JCP, vol 135, p229, (1997)
! Now rewritten to allow compiler vectorizing of loops
! See notes in code
!-------------------------------------------------------------------------
MODULE zremap

  USE shared_data; USE boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: remap_z

  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: rho1, dm, cv2, flux, dzb1
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: rho_v, rho_v1

CONTAINS

  SUBROUTINE remap_z ! remap onto original Eulerian grid

    REAL(num) :: vxb, vxbm, vyb, vybm, vzb, vzbm, dv

    ALLOCATE (rho1(-1:nx+2, -1:ny+2, -1:nz+2), dm(-1:nx+2, -1:ny+2, -1:nz+2), &
        cv2(-1:nx+2, -1:ny+2, -1:nz+2), flux(-1:nx+2, -1:ny+2, -2:nz+2), &
        dzb1(-1:nx+2, -1:ny+2, -1:nz+2), rho_v(-1:nx+2, -1:ny+2, -1:nz+2), &
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

          dv = (REAL(xpass, num) * (vxb - vxbm) / dxb(ix) &
              + REAL(ypass, num) * (vyb - vybm) / dyb(iy) &
              + (vzb - vzbm) / dzb(iz)) * dt

          ! control volume before remap
          cv1(ix, iy, iz) = cv(ix, iy, iz) * (1.0_num + dv)

          dv = (REAL(xpass, num) * (vxb - vxbm) / dxb(ix) &
              + REAL(ypass, num) * (vyb - vybm) / dyb(iy)) * dt

          ! control volume after remap
          cv2(ix, iy, iz) = cv(ix, iy, iz) * (1.0_num + dv)

          ! dzb before remap
          dzb1(ix, iy, iz) = dzb(iz) + (vzb - vzbm) * dt
        END DO
      END DO
    END DO

    ! Evans and Hawley (ApJ, vol 332, p650, (1988))
    ! constrained transport remap of magnetic fluxes
    CALL vz_bx_flux
    DO iz = 1, nz
      izm = iz - 1
      DO iy = 1, ny
        DO ix = 0, nx
          bx(ix, iy, iz) = bx(ix, iy, iz) - flux(ix, iy, iz) + flux(ix, iy, izm)
        END DO
      END DO
    END DO

    DO iz = 0, nz
      DO iy = 1, ny
        DO ix = 1, nx
          ixm = ix - 1
          bz(ix, iy, iz) = bz(ix, iy, iz) + flux(ix, iy, iz) - flux(ixm, iy, iz)
        END DO
      END DO
    END DO

    CALL vz_by_flux
    DO iz = 1, nz
      izm = iz - 1
      DO iy = 0, ny
        DO ix = 1, nx
          by(ix, iy, iz) = by(ix, iy, iz) - flux(ix, iy, iz) + flux(ix, iy, izm)
        END DO
      END DO
    END DO

    DO iz = 0, nz
      DO iy = 1, ny
        iym = iy - 1
        DO ix = 1, nx
          bz(ix, iy, iz) = bz(ix, iy, iz) + flux(ix, iy, iz) - flux(ix, iym, iz)
        END DO
      END DO
    END DO

    ! remap of mass + calculation of mass fluxes (dm) needed for later remaps
    CALL z_mass_flux ! calculates dm(0:nx, 0:ny+1)
    CALL dm_z_bcs    ! need dm(0:nx+1, 0:ny+1) for velocity remap
    DO iz = 1, nz
      izm = iz - 1
      DO iy = 1, ny
        DO ix = 1, nx
          rho(ix, iy, iz) = (rho1(ix, iy, iz) * cv1(ix, iy, iz) &
              + dm(ix, iy, izm) - dm(ix, iy, iz)) / cv2(ix, iy, iz)
        END DO
      END DO
    END DO

    ! remap specific energy density using mass coordinates
    CALL z_energy_flux
    DO iz = 1, nz
      izm = iz - 1
      DO iy = 1, ny
        DO ix = 1, nx
          energy(ix, iy, iz) = (energy(ix, iy, iz) * cv1(ix, iy, iz) &
              * rho1(ix, iy, iz) + flux(ix, iy, izm) - flux(ix, iy, iz)) &
              / (cv2(ix, iy, iz) * rho(ix, iy, iz))
        END DO
      END DO
    END DO

    ! redefine dzb1, cv1, cv2, dm and vz1 for velocity (vertex) cells
    ! in some of these calculations the flux variable is used as a
    ! temporary array
    DO iz = -1, nz+1
      izp = iz + 1
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = 0, nx
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

    DO iz = -2, nz+1
      izp = iz + 1
      DO iy = 0, ny
        DO ix = 0, nx
          flux(ix, iy, iz) = (vz1(ix, iy, iz) + vz1(ix, iy, izp)) / 2.0_num
        END DO
      END DO
    END DO
    ! vertex boundary velocity used in remap
    vz1(0:nx, 0:ny, -2:nz+1) = flux(0:nx, 0:ny, -2:nz+1)

    DO iz = -1, nz+1
      izm = iz - 1
      DO iy = 0, ny
        DO ix = 0, nx
          ! dzb1 = width of vertex CV before remap
          dzb1(ix, iy, iz) = dzc(iz) + (vz1(ix, iy, iz) - vz1(ix, iy, izm)) * dt
        END DO
      END DO
    END DO

    DO iz = -1, nz
      izp = iz + 1
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = 0, nx
          ixp = ix + 1
          flux(ix, iy, iz) = dm(ix, iy, iz) + dm(ixp, iy, iz) &
              + dm(ix, iyp, iz ) + dm(ixp, iyp, iz ) &
              + dm(ix, iy , izp) + dm(ixp, iy , izp) &
              + dm(ix, iyp, izp) + dm(ixp, iyp, izp)
        END DO
      END DO
    END DO
    ! mass flux out of vertex CV
    dm(0:nx, 0:ny, -1:nz) = flux(0:nx, 0:ny, -1:nz) / 8.0_num

    DO iz = 0, nz
      izm = iz - 1
      DO iy = 0, ny
        DO ix = 0, nx
          ! vertex density after remap
          rho_v1(ix, iy, iz) = (rho_v(ix, iy, iz) * cv1(ix, iy, iz) &
              + dm(ix, iy, izm) - dm(ix, iy, iz)) / cv2(ix, iy, iz)
        END DO
      END DO
    END DO

    CALL z_momy_flux
    DO iz = 0, nz
      izm = iz - 1
      DO iy = 0, ny
        DO ix = 0, nx
          vy(ix, iy, iz) = (rho_v(ix, iy, iz) * vy(ix, iy, iz) &
              * cv1(ix, iy, iz) + flux(ix, iy, izm) - flux(ix, iy, iz)) &
              / (cv2(ix, iy, iz) * rho_v1(ix, iy, iz))
        END DO
      END DO
    END DO

    CALL z_momx_flux
    DO iz = 0, nz
      izm = iz - 1
      DO iy = 0, ny
        DO ix = 0, nx
          vx(ix, iy, iz) = (rho_v(ix, iy, iz) * vx(ix, iy, iz) &
              * cv1(ix, iy, iz) + flux(ix, iy, izm) - flux(ix, iy, iz)) &
              / (cv2(ix, iy, iz) * rho_v1(ix, iy, iz))
        END DO
      END DO
    END DO

    CALL z_momz_flux
    DO iz = 0, nz
      izm = iz - 1
      DO iy = 0, ny
        DO ix = 0, nx
          vz(ix, iy, iz) = (rho_v(ix, iy, iz) * vz(ix, iy, iz) &
              * cv1(ix, iy, iz) + flux(ix, iy, izm) - flux(ix, iy, iz)) &
              / (cv2(ix, iy, iz) * rho_v1(ix, iy, iz))
        END DO
      END DO
    END DO

    CALL boundary_conditions

    DEALLOCATE (rho1, dm, cv2, flux, dzb1, rho_v, rho_v1)
    zpass = 0

  END SUBROUTINE remap_z



  SUBROUTINE vz_bx_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: db, dbzp, dbzp2, dbzm
    INTEGER :: izp2

    DO iz = 0, nz
      izm  = iz - 1
      izp  = iz + 1
      izp2 = iz + 2
      DO iy = 0, ny
        iym  = iy - 1
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = 0, nx
          ixp  = ix + 1

          v_advect = (vz1(ix, iy, iz) + vz1(ix, iym, iz)) / 2.0_num

          db    = (dzb1(ix, iy, iz  ) + dzb1(ixp, iy, iz  )) / 2.0_num
          dbzp  = (dzb1(ix, iy, izp ) + dzb1(ixp, iy, izp )) / 2.0_num
          dbzp2 = (dzb1(ix, iy, izp2) + dzb1(ixp, iy, izp2)) / 2.0_num
          dbzm  = (dzb1(ix, iy, izm ) + dzb1(ixp, iy, izm )) / 2.0_num

          w4 = bx(ix, iy, iz ) / db
          w5 = bx(ix, iy, izp) / dbzp

          flux(ix, iy, iz) = (MAX(0.0_num, v_advect) * w4 &
              + MIN(0.0_num, v_advect) * w5) * dt

          w1 = bx(ix, iy, izp ) / dbzp  - bx(ix, iy, iz ) / db
          w2 = bx(ix, iy, iz  ) / db    - bx(ix, iy, izm) / dbzm
          w3 = bx(ix, iy, izp2) / dbzp2 - bx(ix, iy, izp) / dbzp

          ! vad_p and vad_m are logical switches which determine v_advect>=0
          ! and v_advect<0 respectively. It's written this way to allow vector
          ! optimization

          vad_p = -MIN(SIGN(1.0_num, -v_advect), 0.0_num)
          vad_m =  MAX(SIGN(1.0_num, -v_advect), 0.0_num)

          w5 = ABS(v_advect) * dt / (db * vad_p + dbzp * vad_m)
          w4 = (2.0_num - w5) * ABS(w1) / dzc(iz) &
              + (1.0_num + w5) * ABS(w2 * vad_p + w3 * vad_m) &
              / (dzc(izm) * vad_p + dzc(izp) * vad_m)

          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) &
              + SIGN(1.0_num, w2 * vad_p + w3 * vad_m))

          w6 = SIGN(1.0_num, v_advect) * w8 &
              * MIN(ABS(w4) * (dzb(iz) * vad_p + dzb(izp) * vad_m), &
              ABS(w1), ABS(w2 * vad_p + w3 * vad_m))

          flux(ix, iy, iz) = flux(ix, iy, iz) &
              + v_advect * dt * w6 * (1.0_num - w5)
        END DO
      END DO
    END DO

  END SUBROUTINE vz_bx_flux



  SUBROUTINE vz_by_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: db, dbzp, dbzp2, dbzm
    INTEGER :: izp2

    DO iz = 0, nz
      izm  = iz - 1
      izp  = iz + 1
      izp2 = iz + 2
      DO iy = 0, ny
        iyp  = iy + 1
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = 0, nx
          ixm  = ix - 1

          v_advect = (vz1(ix, iy, iz) + vz1(ixm, iy, iz)) / 2.0_num

          db    = (dzb1(ix, iy, iz  ) + dzb1(ix, iyp, iz  )) / 2.0_num
          dbzp  = (dzb1(ix, iy, izp ) + dzb1(ix, iyp, izp )) / 2.0_num
          dbzp2 = (dzb1(ix, iy, izp2) + dzb1(ix, iyp, izp2)) / 2.0_num
          dbzm  = (dzb1(ix, iy, izm ) + dzb1(ix, iyp, izm )) / 2.0_num

          w4 = by(ix, iy, iz ) / db
          w5 = by(ix, iy, izp) / dbzp

          flux(ix, iy, iz) = (MAX(0.0_num, v_advect) * w4 &
              + MIN(0.0_num, v_advect) * w5) * dt

          w1 = by(ix, iy, izp ) / dbzp  - by(ix, iy, iz ) / db
          w2 = by(ix, iy, iz  ) / db    - by(ix, iy, izm) / dbzm
          w3 = by(ix, iy, izp2) / dbzp2 - by(ix, iy, izp) / dbzp

          ! vad_p and vad_m are logical switches which determine v_advect>=0
          ! and v_advect<0 respectively. It's written this way to allow vector
          ! optimization

          vad_p = -MIN(SIGN(1.0_num, -v_advect), 0.0_num)
          vad_m =  MAX(SIGN(1.0_num, -v_advect), 0.0_num)

          w5 = ABS(v_advect) * dt / (db * vad_p + dbzp + vad_m)
          w4 = (2.0_num - w5) * ABS(w1) / dzc(iz) &
              + (1.0_num + w5) * ABS(w2 * vad_p + w3 * vad_m) &
              / (dzc(izm) * vad_p + dzc(izp) * vad_m)

          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) &
              + SIGN(1.0_num, w2 * vad_p + w3 * vad_m))

          w6 = SIGN(1.0_num, v_advect) * w8 &
              * MIN(ABS(w4) * (dzb(iz) * vad_p + dzb(izp) * vad_m), &
              ABS(w1), ABS(w2 * vad_p + w3 * vad_m))

          flux(ix, iy, iz) = flux(ix, iy, iz) &
              + v_advect * dt * w6 * (1.0_num - w5)
        END DO
      END DO
    END DO

  END SUBROUTINE vz_by_flux



  SUBROUTINE z_mass_flux

    REAL(num) :: v_advect, flux_rho, vad_p, vad_m
    INTEGER :: izp2

    DO iz = 0, nz
      izm  = iz - 1
      izp  = iz + 1
      izp2 = iz + 2
      DO iy = 0, ny+1
        iym  = iy - 1
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = 0, nx+1
          ixm  = ix - 1

          v_advect = (vz1(ix, iy, iz) + vz1(ix, iym, iz) &
              + vz1(ixm, iy, iz) + vz1(ixm, iym, iz)) / 4.0_num

          dm(ix, iy, iz) = (MAX(0.0_num, v_advect) * rho(ix, iy, iz) &
              + MIN(0.0_num, v_advect) * rho(ix, iy, izp)) * dt

          w1 = rho(ix, iy, izp ) - rho(ix, iy, iz )
          w2 = rho(ix, iy, iz  ) - rho(ix, iy, izm)
          w3 = rho(ix, iy, izp2) - rho(ix, iy, izp)

          ! vad_p and vad_m are logical switches which determine v_advect>=0
          ! and v_advect<0 respectively. It's written this way to allow vector
          ! optimization

          vad_p = -MIN(SIGN(1.0_num, -v_advect), 0.0_num)
          vad_m =  MAX(SIGN(1.0_num, -v_advect), 0.0_num)

          w5 = ABS(v_advect) * dt &
              / (dzb1(ix, iy, iz) * vad_p + dzb1(ix, iy, izp) * vad_m)

          w4 = (2.0_num - w5) * ABS(w1) / dzc(iz) &
              + (1.0_num + w5) * ABS(w2 * vad_p + w3 * vad_m) &
              / (dzc(izm) * vad_p + dzc(izp) * vad_m)

          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) &
              + SIGN(1.0_num, w2 * vad_p + w3 * vad_m))

          w6 = SIGN(1.0_num, v_advect) * w8 &
              * MIN(ABS(w4) * (dzb(iz) * vad_p + dzb(izp) * vad_m), &
              ABS(w1), ABS(w2 * vad_p + w3 * vad_m))

          flux_rho = v_advect * dt * w6 * (1.0_num - w5)

          dm(ix, iy, iz) = (flux_rho + dm(ix, iy, iz)) * dxb(ix) * dyb(iy)
        END DO
      END DO
    END DO

  END SUBROUTINE z_mass_flux



  SUBROUTINE z_energy_flux ! energy remap in mass coordinates

    REAL(num) :: v_advect, vad_p, vad_m
    INTEGER :: izp2

    DO iz = 0, nz
      izm  = iz - 1
      izp  = iz + 1
      izp2 = iz + 2
      DO iy = 0, ny
        iym  = iy - 1
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = 0, nx
          ixm  = ix - 1

          v_advect = (vz1(ix, iy, iz) + vz1(ix, iym, iz) &
              + vz1(ixm, iy, iz) + vz1(ixm, iym, iz)) / 4.0_num

          w1 = energy(ix, iy, izp ) - energy(ix, iy, iz )
          w2 = energy(ix, iy, iz  ) - energy(ix, iy, izm)
          w3 = energy(ix, iy, izp2) - energy(ix, iy, izp)

          ! vad_p and vad_m are logical switches which determine v_advect>=0
          ! and v_advect<0 respectively. It's written this way to allow vector
          ! optimization

          vad_p = -MIN(SIGN(1.0_num, -v_advect), 0.0_num)
          vad_m =  MAX(SIGN(1.0_num, -v_advect), 0.0_num)

          w5 = ABS(v_advect) * dt &
              / (dzb1(ix, iy, iz) * vad_p + dzb1(ix, iy, izp) * vad_m)

          w7 = energy(ix, iy, iz) * vad_p + energy(ix, iy, izp) * vad_m

          w6 = ABS(dm(ix, iy, iz)) / dxb(ix) / dyb(iy) &
              / (rho1(ix, iy, iz ) * dzb1(ix, iy, iz ) * vad_p &
              +  rho1(ix, iy, izp) * dzb1(ix, iy, izp) * vad_m)

          w4 = (2.0_num - w5) * ABS(w1) / dzc(iz) &
              + (1.0_num + w5) * ABS(w2 * vad_p + w3 * vad_m) &
              / (dzc(izm) * vad_p + dzc(izp) * vad_m)

          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) &
              + SIGN(1.0_num, w2 * vad_p + w3 * vad_m))

          w5 = SIGN(1.0_num, v_advect) * w8 &
              * MIN(ABS(w4) * (dzb(iz) * vad_p + dzb(izp) * vad_m), &
              ABS(w1), ABS(w2 * vad_p + w3 * vad_m))

          flux(ix, iy, iz) = dm(ix, iy, iz) * (w7 + w5 * (1.0_num - w6))
        END DO
      END DO
    END DO

  END SUBROUTINE z_energy_flux



  SUBROUTINE z_momy_flux ! energy remap in mass coordinates

    REAL(num) :: v_advect, m, mp, ai, aip, dk
    REAL(num) :: vad_p, vad_m
    INTEGER :: izp2

    DO iz = -1, nz
      izm  = iz - 1
      izp  = iz + 1
      izp2 = iz + 2
      DO iy = 0, ny
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = 0, nx

          v_advect = vz1(ix, iy, iz)

          w1 = vy(ix, iy, izp ) - vy(ix, iy, iz )
          w2 = vy(ix, iy, iz  ) - vy(ix, iy, izm)
          w3 = vy(ix, iy, izp2) - vy(ix, iy, izp)

          ! vad_p and vad_m are logical switches which determine v_advect>=0
          ! and v_advect<0 respectively. It's written this way to allow vector
          ! optimization

          vad_p = -MIN(SIGN(1.0_num, -v_advect), 0.0_num)
          vad_m =  MAX(SIGN(1.0_num, -v_advect), 0.0_num)

          w5 = ABS(v_advect) * dt &
              / (dzb1(ix, iy, iz) * vad_p + dzb1(ix, iy, izp) * vad_m)

          w7 = vy(ix, iy, iz) * vad_p + vy(ix, iy, izp) * vad_m

          w6 = ABS(dm(ix, iy, iz)) / dyc(iy) / dxc(ix) &
              / (rho_v(ix, iy, iz ) * dzb1(ix, iy, iz ) * vad_p &
              +  rho_v(ix, iy, izp) * dzb1(ix, iy, izp) * vad_m)

          w4 = (2.0_num - w5) * ABS(w1) / dzb(izp) &
              + (1.0_num + w5) * ABS(w2 * vad_p + w3 * vad_m) &
              / (dzb(iz) * vad_p + dzb(izp2) * vad_m)

          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) &
              + SIGN(1.0_num, w2 * vad_p + w3 * vad_m))

          w5 = SIGN(1.0_num, v_advect) * w8 &
              * MIN(ABS(w4) * (dzc(iz) * vad_p + dzc(izp) * vad_m), &
              ABS(w1), ABS(w2 * vad_p + w3 * vad_m))

          flux(ix, iy, iz) = w7 + w5 * (1.0_num - w6)
        END DO
      END DO
    END DO

    IF (rke) THEN
      DO iz = 0, nz-1
        izm = iz - 1
        izp = iz + 1
        DO iy = 0, ny
          iyp = iy + 1
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
          DO ix = 0, nx
            ixp = ix + 1

            m = rho_v1(ix, iy, iz) * cv2(ix, iy, iz)
            mp = rho_v1(ix, iy, izp) * cv2(ix, iy, izp)

            ai = (vy(ix, iy, iz) - flux(ix, iy, izm)) * dm(ix, iy, izm) / m &
                + (flux(ix, iy, iz) - vy(ix, iy, iz)) * dm(ix, iy, iz) / m

            aip = (vy(ix, iy, izp) - flux(ix, iy, iz)) * dm(ix, iy, iz) / mp &
                + (flux(ix, iy, izp) - vy(ix, iy, izp)) * dm(ix, iy, izp) / mp

            dk = (vy(ix, iy, izp) - vy(ix, iy, iz)) * (flux(ix, iy, iz) &
                - 0.5_num * (vy(ix, iy, izp) + vy(ix, iy, iz))) &
                + 0.5_num * ai * (flux(ix, iy, iz) - vy(ix, iy, iz)) &
                + 0.5_num * aip * (vy(ix, iy, izp) - flux(ix, iy, iz))

            dk = dk * dm(ix, iy, iz) / 4.0_num
            delta_ke(ix , iy , izp) = delta_ke(ix , iy , izp) + dk
            delta_ke(ixp, iy , izp) = delta_ke(ixp, iy , izp) + dk
            delta_ke(ix , iyp, izp) = delta_ke(ix , iyp, izp) + dk
            delta_ke(ixp, iyp, izp) = delta_ke(ixp, iyp, izp) + dk
          END DO
        END DO
      END DO
    END IF

    flux(0:nx, 0:ny, -1:nz) = flux(0:nx, 0:ny, -1:nz) * dm(0:nx, 0:ny, -1:nz)

  END SUBROUTINE z_momy_flux



  SUBROUTINE z_momx_flux ! energy remap in mass coordinates

    REAL(num) :: v_advect, m, mp, ai, aip, dk
    REAL(num) :: vad_p, vad_m
    INTEGER :: izp2

    DO iz = -1, nz
      izm  = iz - 1
      izp  = iz + 1
      izp2 = iz + 2
      DO iy = 0, ny
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = 0, nx

          v_advect = vz1(ix, iy, iz)

          w1 = vx(ix, iy, izp ) - vx(ix, iy, iz )
          w2 = vx(ix, iy, iz  ) - vx(ix, iy, izm)
          w3 = vx(ix, iy, izp2) - vx(ix, iy, izp)

          ! vad_p and vad_m are logical switches which determine v_advect>=0
          ! and v_advect<0 respectively. It's written this way to allow vector
          ! optimization

          vad_p = -MIN(SIGN(1.0_num, -v_advect), 0.0_num)
          vad_m =  MAX(SIGN(1.0_num, -v_advect), 0.0_num)

          w5 = ABS(v_advect) * dt &
              / (dzb1(ix, iy, iz) * vad_p + dzb1(ix, iy, izp) * vad_m)

          w7 = vx(ix, iy, iz) * vad_p + vx(ix, iy, izp) * vad_m

          w6 = ABS(dm(ix, iy, iz)) / dyc(iy) / dxc(ix) &
              / (rho_v(ix, iy, iz ) * dzb1(ix, iy, iz ) * vad_p &
              +  rho_v(ix, iy, izp) * dzb1(ix, iy, izp) * vad_m)

          w4 = (2.0_num - w5) * ABS(w1) / dzb(izp) &
              + (1.0_num + w5) * ABS(w2 * vad_p + w3 * vad_m) &
              / (dzb(iz) * vad_p + dzb(izp2) * vad_m)

          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) &
              + SIGN(1.0_num, w2 * vad_p + w3 * vad_m))

          w5 = SIGN(1.0_num, v_advect) * w8 &
              * MIN(ABS(w4) * (dzc(iz) * vad_p + dzc(izp) * vad_m), &
              ABS(w1), ABS(w2 * vad_p + w3 * vad_m))

          flux(ix, iy, iz) = w7 + w5 * (1.0_num - w6)
        END DO
      END DO
    END DO

    IF (rke) THEN
      DO iz = 0, nz-1
        izm = iz - 1
        izp = iz + 1
        DO iy = 0, ny
          iyp = iy + 1
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
          DO ix = 0, nx
            ixp = ix + 1

            m = rho_v1(ix, iy, iz) * cv2(ix, iy, iz)
            mp = rho_v1(ix, iy, izp) * cv2(ix, iy, izp)

            ai = (vx(ix, iy, iz) - flux(ix, iy, izm)) * dm(ix, iy, izm) / m &
                + (flux(ix, iy, iz) - vx(ix, iy, iz)) * dm(ix, iy, iz) / m

            aip = (vx(ix, iy, izp) - flux(ix, iy, iz)) * dm(ix, iy, iz) / mp &
                + (flux(ix, iy, izp) - vx(ix, iy, izp)) * dm(ix, iy, izp) / mp

            dk = (vx(ix, iy, izp) - vx(ix, iy, iz)) * (flux(ix, iy, iz) &
                - 0.5_num * (vx(ix, iy, izp) + vx(ix, iy, iz))) &
                + 0.5_num * ai * (flux(ix, iy, iz) - vx(ix, iy, iz)) &
                + 0.5_num * aip * (vx(ix, iy, izp) - flux(ix, iy, iz))

            dk = dk * dm(ix, iy, iz) / 4.0_num
            delta_ke(ix , iy , izp) = delta_ke(ix , iy , izp) + dk
            delta_ke(ixp, iy , izp) = delta_ke(ixp, iy , izp) + dk
            delta_ke(ix , iyp, izp) = delta_ke(ix , iyp, izp) + dk
            delta_ke(ixp, iyp, izp) = delta_ke(ixp, iyp, izp) + dk
          END DO
        END DO
      END DO
    END IF

    flux(0:nx, 0:ny, -1:nz) = flux(0:nx, 0:ny, -1:nz) * dm(0:nx, 0:ny, -1:nz)

  END SUBROUTINE z_momx_flux



  SUBROUTINE z_momz_flux

    REAL(num) :: v_advect, m, mp, ai, aip, dk
    REAL(num) :: vad_p, vad_m
    INTEGER :: izp2

    DO iz = -1, nz
      izm  = iz - 1
      izp  = iz + 1
      izp2 = iz + 2
      DO iy = 0, ny
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = 0, nx

          v_advect = vz1(ix, iy, iz)

          w1 = vz(ix, iy, izp ) - vz(ix, iy, iz )
          w2 = vz(ix, iy, iz  ) - vz(ix, iy, izm)
          w3 = vz(ix, iy, izp2) - vz(ix, iy, izp)

          ! vad_p and vad_m are logical switches which determine v_advect>=0
          ! and v_advect<0 respectively. It's written this way to allow vector
          ! optimization

          vad_p = -MIN(SIGN(1.0_num, -v_advect), 0.0_num)
          vad_m =  MAX(SIGN(1.0_num, -v_advect), 0.0_num)

          w5 = ABS(v_advect) * dt &
              / (dzb1(ix, iy, iz) * vad_p + dzb1(ix, iy, izp) * vad_m)

          w7 = vz(ix, iy, iz) * vad_p + vz(ix, iy, izp) * vad_m

          w6 = ABS(dm(ix, iy, iz)) / dyc(iy) / dxc(ix) &
              / (rho_v(ix, iy, iz ) * dzb1(ix, iy, iz ) * vad_p &
              +  rho_v(ix, iy, izp) * dzb1(ix, iy, izp) * vad_m)

          w4 = (2.0_num - w5) * ABS(w1) / dzb(izp) &
              + (1.0_num + w5) * ABS(w2 * vad_p + w3 * vad_m) &
              / (dzb(iz) * vad_p + dzb(izp2) * vad_m)

          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) &
              + SIGN(1.0_num, w2 * vad_p + w3 * vad_m))

          w5 = SIGN(1.0_num, v_advect) * w8 &
              * MIN(ABS(w4) * (dzc(iz) * vad_p + dzc(izp) * vad_m), &
              ABS(w1), ABS(w2 * vad_p + w3 * vad_m))

          flux(ix, iy, iz) = w7 + w5 * (1.0_num - w6)
        END DO
      END DO
    END DO

    IF (rke) THEN
      DO iz = 0, nz-1
        izm = iz - 1
        izp = iz + 1
        DO iy = 0, ny
          iyp = iy + 1
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
          DO ix = 0, nx
            ixp = ix + 1

            m = rho_v1(ix, iy, iz) * cv2(ix, iy, iz)
            mp = rho_v1(ix, iy, izp) * cv2(ix, iy, izp)

            ai = (vz(ix, iy, iz) - flux(ix, iy, izm)) * dm(ix, iy, izm) / m &
                + (flux(ix, iy, iz) - vz(ix, iy, iz)) * dm(ix, iy, iz) / m

            aip = (vz(ix, iy, izp) - flux(ix, iy, iz)) * dm(ix, iy, iz) / mp &
                + (flux(ix, iy, izp) - vz(ix, iy, izp)) * dm(ix, iy, izp) / mp

            dk = (vz(ix, iy, izp) - vz(ix, iy, iz)) * (flux(ix, iy, iz) &
                - 0.5_num * (vz(ix, iy, izp) + vz(ix, iy, iz))) &
                + 0.5_num * ai * (flux(ix, iy, iz) - vz(ix, iy, iz)) &
                + 0.5_num * aip * (vz(ix, iy, izp) - flux(ix, iy, iz))

            dk = dk * dm(ix, iy, iz) / 4.0_num
            delta_ke(ix , iy , izp) = delta_ke(ix , iy , izp) + dk
            delta_ke(ixp, iy , izp) = delta_ke(ixp, iy , izp) + dk
            delta_ke(ix , iyp, izp) = delta_ke(ix , iyp, izp) + dk
            delta_ke(ixp, iyp, izp) = delta_ke(ixp, iyp, izp) + dk
          END DO
        END DO
      END DO
    END IF

    flux(0:nx, 0:ny, -1:nz) = flux(0:nx, 0:ny, -1:nz) * dm(0:nx, 0:ny, -1:nz)

  END SUBROUTINE z_momz_flux



  SUBROUTINE dm_z_bcs

    CALL MPI_SENDRECV(dm(0:nx+1, 0:ny+1, 1), (nx+2)*(ny+2), mpireal, &
        front, tag, dm(0:nx+1, 0:ny+1, nz+1), (nx+2)*(ny+2), mpireal, &
        back, tag, comm, status, errcode)

    IF (back == MPI_PROC_NULL) &
        dm(0:nx+1, 0:ny+1, nz+1) = dm(0:nx+1, 0:ny+1, nz)

    CALL MPI_SENDRECV(dm(0:nx+1, 0:ny+1, nz-1), (nx+2)*(ny+2), mpireal, &
        back, tag, dm(0:nx+1, 0:ny+1, -1), (nx+2)*(ny+2), mpireal, &
        front, tag, comm, status, errcode)

    IF (front == MPI_PROC_NULL) &
        dm(0:nx+1, 0:ny+1, -1) = dm(0:nx+1, 0:ny+1, 0)

  END SUBROUTINE dm_z_bcs

END MODULE zremap
