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
! Mass coordinate based Van Leer limited remap.
! See Bram van Leer, JCP, vol 135, p229, (1997)
!******************************************************************************

MODULE zremap

  USE shared_data; USE boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: remap_z

  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: dzb1, dzc1, rho_v, rho_v1
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: rho1, dm, cv2, flux

CONTAINS

  ! Remap onto original Eulerian grid

  SUBROUTINE remap_z

    REAL(num) :: vxb, vxbm, vyb, vybm, vzb, vzbm, dv

    IF (predictor_step) THEN
      ALLOCATE(cv2   (-1:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(flux  (-1:nx+2, -1:ny+2, -2:nz+2))
      ALLOCATE(dzb1  (-1:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(dzc1  (-1:nx+2, -1:ny+2, -1:nz+2))
    ELSE
      ALLOCATE(rho1  (-1:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(dm    (-1:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(cv2   (-1:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(flux  (-1:nx+2, -1:ny+2, -2:nz+2))
      ALLOCATE(dzb1  (-1:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(dzc1  (-1:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(rho_v (-1:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(rho_v1(-1:nx+2, -1:ny+2, -1:nz+2))
    END IF

    DO iz = -1, nz + 2
      izm = iz - 1
      DO iy = -1, ny + 2
        iym = iy - 1
        DO ix = -1, nx + 2
          ixm = ix - 1

          ! vx at Bx(i,j,k)
          vxb  = (vx1(ix ,iy ,iz ) + vx1(ix ,iym,iz ) &
                + vx1(ix ,iy ,izm) + vx1(ix ,iym,izm)) * 0.25_num

          ! vx at Bx(i-1,j,k)
          vxbm = (vx1(ixm,iy ,iz ) + vx1(ixm,iym,iz ) &
                + vx1(ixm,iy ,izm) + vx1(ixm,iym,izm)) * 0.25_num

          ! vy at By(i,j,k)
          vyb  = (vy1(ix ,iy ,iz ) + vy1(ixm,iy ,iz ) &
                + vy1(ix ,iy ,izm) + vy1(ixm,iy ,izm)) * 0.25_num

          ! vy at By(i,j-1,k)
          vybm = (vy1(ix ,iym,iz ) + vy1(ixm,iym,iz ) &
                + vy1(ix ,iym,izm) + vy1(ixm,iym,izm)) * 0.25_num

          ! vz at Bz(i,j,k)
          vzb  = (vz1(ix ,iy ,iz ) + vz1(ixm,iy ,iz ) &
                + vz1(ix ,iym,iz ) + vz1(ixm,iym,iz )) * 0.25_num

          ! vz at Bz(i,j,k-1)
          vzbm = (vz1(ix ,iy ,izm) + vz1(ixm,iy ,izm) &
                + vz1(ix ,iym,izm) + vz1(ixm,iym,izm)) * 0.25_num

          dv = (REAL(xpass, num) * (vxb - vxbm) / dxb(ix) &
              + REAL(ypass, num) * (vyb - vybm) / dyb(iy) &
              + (vzb - vzbm) / dzb(iz)) * dt

          ! Control volume before remap
          cv1(ix,iy,iz) = cv(ix,iy,iz) * (1.0_num + dv)

          dv = (REAL(xpass, num) * (vxb - vxbm) / dxb(ix) &
              + REAL(ypass, num) * (vyb - vybm) / dyb(iy)) * dt

          ! Control volume after remap
          cv2(ix,iy,iz) = cv(ix,iy,iz) * (1.0_num + dv)

          ! dzb before remap
          dzb1(ix,iy,iz) = dzb(iz) + (vzb - vzbm) * dt
        END DO
      END DO
    END DO

    DO iz = -1, nz + 2
      izp = MIN(iz + 1, nz + 2)
      DO iy = -1, ny + 2
        DO ix = -1, nx + 2
          ! dzc before remap
          dzc1(ix,iy,iz) = 0.5_num * (dzb1(ix,iy,iz) + dzb1(ix,iy,izp))
        END DO
      END DO
    END DO

    ! Evans and Hawley (ApJ, vol 332, p650, (1988))
    ! constrained transport remap of magnetic fluxes
    CALL vz_bx_flux

    DO iz = 0, nz + 1
      izm = iz - 1
      DO iy = 0, ny + 1
        DO ix = 0, nx + 1
          bx(ix,iy,iz) = bx(ix,iy,iz) - flux(ix,iy,iz) + flux(ix,iy,izm)
        END DO
      END DO
    END DO

    DO iz = 0, nz + 1
      DO iy = 0, ny + 1
        DO ix = 0, nx + 1
          ixm = ix - 1
          bz(ix,iy,iz) = bz(ix,iy,iz) + flux(ix,iy,iz) - flux(ixm,iy,iz)
        END DO
      END DO
    END DO

    CALL vz_by_flux

    DO iz = 0, nz + 1
      izm = iz - 1
      DO iy = 0, ny + 1
        DO ix = 0, nx + 1
          by(ix,iy,iz) = by(ix,iy,iz) - flux(ix,iy,iz) + flux(ix,iy,izm)
        END DO
      END DO
    END DO

    DO iz = 0, nz + 1
      DO iy = 0, ny + 1
        iym = iy - 1
        DO ix = 0, nx + 1
          bz(ix,iy,iz) = bz(ix,iy,iz) + flux(ix,iy,iz) - flux(ix,iym,iz)
        END DO
      END DO
    END DO

    IF (predictor_step) THEN
      DEALLOCATE(cv2, flux, dzb1, dzc1)
      RETURN
    END IF

    dm = 0.0_num
    ! Store initial density in rho1
    rho1(:,:,:) = rho(:,:,:)

    ! Remap of mass + calculation of mass fluxes (dm) needed for later remaps
    ! Calculates dm(0:nx+1,0:ny+1,0:nz)
    CALL z_mass_flux
    ! Need dm(0:nx+1,0:ny+1,-1:nz+1) for velocity remap
    CALL dm_z_bcs

    DO iz = 1, nz
      izm = iz - 1
      DO iy = 1, ny
        DO ix = 1, nx
          rho(ix,iy,iz) = (rho1(ix,iy,iz) * cv1(ix,iy,iz) &
              + dm(ix,iy,izm) - dm(ix,iy,iz)) / cv2(ix,iy,iz)
        END DO
      END DO
    END DO

    ! Remap specific energy density using mass coordinates
    CALL z_energy_flux

    DO iz = 1, nz
      izm = iz - 1
      DO iy = 1, ny
        DO ix = 1, nx
          energy(ix,iy,iz) = (energy(ix,iy,iz) * cv1(ix,iy,iz) &
              * rho1(ix,iy,iz) + flux(ix,iy,izm) - flux(ix,iy,iz)) &
              / (cv2(ix,iy,iz) * rho(ix,iy,iz))
        END DO
      END DO
    END DO

    ! Redefine dzb1, cv1, cv2, dm and vz1 for velocity (vertex) cells.
    ! In some of these calculations the flux variable is used as a
    ! temporary array
    DO iz = -1, nz + 1
      izp = iz + 1
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = 0, nx
          ixp = ix + 1

          ! Vertex density before remap
          rho_v(ix,iy,iz) = rho1(ix,iy,iz) * cv1(ix,iy,iz) &
              + rho1(ixp,iy ,iz ) * cv1(ixp,iy ,iz ) &
              + rho1(ix ,iyp,iz ) * cv1(ix ,iyp,iz ) &
              + rho1(ixp,iyp,iz ) * cv1(ixp,iyp,iz ) &
              + rho1(ix ,iy ,izp) * cv1(ix ,iy ,izp) &
              + rho1(ixp,iy ,izp) * cv1(ixp,iy ,izp) &
              + rho1(ix ,iyp,izp) * cv1(ix ,iyp,izp) &
              + rho1(ixp,iyp,izp) * cv1(ixp,iyp,izp)

          rho_v(ix,iy,iz) = rho_v(ix,iy,iz) &
              / (cv1(ix,iy ,iz ) + cv1(ixp,iy ,iz ) &
              + cv1(ix,iyp,iz ) + cv1(ixp,iyp,iz ) &
              + cv1(ix,iy ,izp) + cv1(ixp,iy ,izp) &
              + cv1(ix,iyp,izp) + cv1(ixp,iyp,izp))
        END DO
      END DO
    END DO

    DO iz = 0, nz
      izp = iz + 1
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = 0, nx
          ixp = ix + 1
          flux(ix,iy,iz) = cv1(ix,iy,iz) + cv1(ixp,iy,iz) &
              + cv1(ix,iyp,iz ) + cv1(ixp,iyp,iz ) &
              + cv1(ix,iy ,izp) + cv1(ixp,iy ,izp) &
              + cv1(ix,iyp,izp) + cv1(ixp,iyp,izp)
        END DO
      END DO
    END DO
    cv1(0:nx,0:ny,0:nz) = flux(0:nx,0:ny,0:nz) * 0.125_num

    DO iz = 0, nz
      izp = iz + 1
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = 0, nx
          ixp = ix + 1
          flux(ix,iy,iz) = cv2(ix,iy,iz) + cv2(ixp,iy,iz) &
              + cv2(ix,iyp,iz ) + cv2(ixp,iyp,iz ) &
              + cv2(ix,iy ,izp) + cv2(ixp,iy ,izp) &
              + cv2(ix,iyp,izp) + cv2(ixp,iyp,izp)
        END DO
      END DO
    END DO
    cv2(0:nx,0:ny,0:nz) = flux(0:nx,0:ny,0:nz) * 0.125_num

    DO iz = -2, nz + 1
      izp = iz + 1
      DO iy = 0, ny
        DO ix = 0, nx
          flux(ix,iy,iz) = (vz1(ix,iy,iz) + vz1(ix,iy,izp)) * 0.5_num
        END DO
      END DO
    END DO
    vz1(0:nx,0:ny,-2:nz+1) = flux(0:nx,0:ny,-2:nz+1)

    DO iz = -1, nz
      izp = iz + 1
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = 0, nx
          ixp = ix + 1
          flux(ix,iy,iz) = dm(ix,iy,iz) + dm(ixp,iy,iz) &
              + dm(ix,iyp,iz ) + dm(ixp,iyp,iz ) &
              + dm(ix,iy ,izp) + dm(ixp,iy ,izp) &
              + dm(ix,iyp,izp) + dm(ixp,iyp,izp)
        END DO
      END DO
    END DO
    dm(0:nx,0:ny,-1:nz) = flux(0:nx,0:ny,-1:nz) * 0.125_num

    DO iz = 0, nz
      izm = iz - 1
      DO iy = 0, ny
        DO ix = 0, nx
          rho_v1(ix,iy,iz) = (rho_v(ix,iy,iz) * cv1(ix,iy,iz) &
              + dm(ix,iy,izm) - dm(ix,iy,iz)) / cv2(ix,iy,iz)
        END DO
      END DO
    END DO

    CALL z_momx_flux

    DO iz = 0, nz
      izm = iz - 1
      DO iy = 0, ny
        DO ix = 0, nx
          vx(ix,iy,iz) = (rho_v(ix,iy,iz) * vx(ix,iy,iz) &
              * cv1(ix,iy,iz) + flux(ix,iy,izm) - flux(ix,iy,iz)) &
              / (cv2(ix,iy,iz) * rho_v1(ix,iy,iz))
        END DO
      END DO
    END DO

    CALL z_momy_flux

    DO iz = 0, nz
      izm = iz - 1
      DO iy = 0, ny
        DO ix = 0, nx
          vy(ix,iy,iz) = (rho_v(ix,iy,iz) * vy(ix,iy,iz) &
              * cv1(ix,iy,iz) + flux(ix,iy,izm) - flux(ix,iy,iz)) &
              / (cv2(ix,iy,iz) * rho_v1(ix,iy,iz))
        END DO
      END DO
    END DO

    CALL z_momz_flux

    DO iz = 0, nz
      izm = iz - 1
      DO iy = 0, ny
        DO ix = 0, nx
          vz(ix,iy,iz) = (rho_v(ix,iy,iz) * vz(ix,iy,iz) &
              * cv1(ix,iy,iz) + flux(ix,iy,izm) - flux(ix,iy,iz)) &
              / (cv2(ix,iy,iz) * rho_v1(ix,iy,iz))
        END DO
      END DO
    END DO

    CALL boundary_conditions

    DEALLOCATE (rho1, dm, cv2, flux, dzb1, dzc1, rho_v, rho_v1)
    zpass = 0

  END SUBROUTINE remap_z



  SUBROUTINE vz_bx_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dzci, dzcu, dzbu, phi, ss, Da, Di
    REAL(num) :: dzu, dbz, dbzp, dbzp2, dbzm
    INTEGER :: izp2

    DO iz = -1, nz + 1
      izm  = MAX(iz - 1, -1)
      izp  = iz + 1
      izp2 = MIN(iz + 2, nz + 2)
      DO iy = -1, ny + 1
        iym = iy - 1
        DO ix = -1, nx + 1
          ixp = ix + 1

          v_advect = (vz1(ix,iy,iz) + vz1(ix,iym,iz)) * 0.5_num

          dbz   = (dzb1(ix,iy,iz  ) + dzb1(ixp,iy,iz  )) * 0.5_num
          dbzp  = (dzb1(ix,iy,izp ) + dzb1(ixp,iy,izp )) * 0.5_num
          dbzp2 = (dzb1(ix,iy,izp2) + dzb1(ixp,iy,izp2)) * 0.5_num
          dbzm  = (dzb1(ix,iy,izm ) + dzb1(ixp,iy,izm )) * 0.5_num

          fm  = bx(ix,iy,izm ) / dbzm
          fi  = bx(ix,iy,iz  ) / dbz
          fp  = bx(ix,iy,izp ) / dbzp
          fp2 = bx(ix,iy,izp2) / dbzp2

          dfm = fi - fm
          dfi = fp - fi
          dfp = fp2 - fp

          ! vad_p and vad_m are logical switches which determine v_advect>=0
          ! and v_advect<0 respectively. It's written this way to allow vector
          ! optimization

          sign_v = SIGN(1.0_num, v_advect)
          vad_p = (sign_v + 1.0_num) * 0.5_num
          vad_m = 1.0_num - vad_p

          fu = fi * vad_p + fp * vad_m
          dfu = dfm * vad_p + dfp * vad_m
          dzci = dzc1(ix,iy,iz )
          dzcu = dzc1(ix,iy,izm) * vad_p + dzc1(ix,iy,izp) * vad_m
          dzbu = dzb1(ix,iy,iz ) * vad_p + dzb1(ix,iy,izp) * vad_m

          dzu = dbz * vad_p + dbzp * vad_m
          phi = ABS(v_advect) * dt / dzu

          Da =  (2.0_num - phi) * ABS(dfi) / dzci &
              + (1.0_num + phi) * ABS(dfu) / dzcu
          Da = Da * sixth

          ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

          Di = sign_v * ss * MIN(ABS(Da) * dzbu, ABS(dfi), ABS(dfu))

          flux(ix,iy,iz) = (fu + Di * (1.0_num - phi)) * v_advect * dt
        END DO
      END DO
    END DO

  END SUBROUTINE vz_bx_flux



  SUBROUTINE vz_by_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dzci, dzcu, dzbu, phi, ss, Da, Di
    REAL(num) :: dzu, dbz, dbzp, dbzp2, dbzm
    INTEGER :: izp2

    DO iz = -1, nz + 1
      izm  = MAX(iz - 1, -1)
      izp  = iz + 1
      izp2 = MIN(iz + 2, nz + 2)
      DO iy = -1, ny + 1
        iyp = iy + 1
        DO ix = -1, nx + 1
          ixm = ix - 1

          v_advect = (vz1(ix,iy,iz) + vz1(ixm,iy,iz)) * 0.5_num

          dbz   = (dzb1(ix,iy,iz  ) + dzb1(ix,iyp,iz  )) * 0.5_num
          dbzp  = (dzb1(ix,iy,izp ) + dzb1(ix,iyp,izp )) * 0.5_num
          dbzp2 = (dzb1(ix,iy,izp2) + dzb1(ix,iyp,izp2)) * 0.5_num
          dbzm  = (dzb1(ix,iy,izm ) + dzb1(ix,iyp,izm )) * 0.5_num

          fm  = by(ix,iy,izm ) / dbzm
          fi  = by(ix,iy,iz  ) / dbz
          fp  = by(ix,iy,izp ) / dbzp
          fp2 = by(ix,iy,izp2) / dbzp2

          dfm = fi - fm
          dfi = fp - fi
          dfp = fp2 - fp

          ! vad_p and vad_m are logical switches which determine v_advect>=0
          ! and v_advect<0 respectively. It's written this way to allow vector
          ! optimization

          sign_v = SIGN(1.0_num, v_advect)
          vad_p = (sign_v + 1.0_num) * 0.5_num
          vad_m = 1.0_num - vad_p

          fu = fi * vad_p + fp * vad_m
          dfu = dfm * vad_p + dfp * vad_m
          dzci = dzc1(ix,iy,iz )
          dzcu = dzc1(ix,iy,izm) * vad_p + dzc1(ix,iy,izp) * vad_m
          dzbu = dzb1(ix,iy,iz ) * vad_p + dzb1(ix,iy,izp) * vad_m

          dzu = dbz * vad_p + dbzp * vad_m
          phi = ABS(v_advect) * dt / dzu

          Da =  (2.0_num - phi) * ABS(dfi) / dzci &
              + (1.0_num + phi) * ABS(dfu) / dzcu
          Da = Da * sixth

          ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

          Di = sign_v * ss * MIN(ABS(Da) * dzbu, ABS(dfi), ABS(dfu))

          flux(ix,iy,iz) = (fu + Di * (1.0_num - phi)) * v_advect * dt
        END DO
      END DO
    END DO

  END SUBROUTINE vz_by_flux



  SUBROUTINE z_mass_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dzci, dzcu, dzbu, phi, ss, Da, Di
    REAL(num) :: area
    INTEGER :: izp2

    DO iz = 0, nz
      izm  = iz - 1
      izp  = iz + 1
      izp2 = iz + 2
      DO iy = 0, ny + 1
        iym = iy - 1
        DO ix = 0, nx + 1
          ixm = ix - 1
          area = dxb(ix) * dyb(iy)

          v_advect = (vz1(ix,iy ,iz) + vz1(ixm,iy ,iz) &
                    + vz1(ix,iym,iz) + vz1(ixm,iym,iz)) * 0.25_num

          fm  = rho(ix,iy,izm )
          fi  = rho(ix,iy,iz  )
          fp  = rho(ix,iy,izp )
          fp2 = rho(ix,iy,izp2)

          dfm = fi - fm
          dfi = fp - fi
          dfp = fp2 - fp

          ! vad_p and vad_m are logical switches which determine v_advect>=0
          ! and v_advect<0 respectively. It's written this way to allow vector
          ! optimization

          sign_v = SIGN(1.0_num, v_advect)
          vad_p = (sign_v + 1.0_num) * 0.5_num
          vad_m = 1.0_num - vad_p

          fu = fi * vad_p + fp * vad_m
          dfu = dfm * vad_p + dfp * vad_m
          dzci = dzc1(ix,iy,iz )
          dzcu = dzc1(ix,iy,izm) * vad_p + dzc1(ix,iy,izp) * vad_m
          dzbu = dzb1(ix,iy,iz ) * vad_p + dzb1(ix,iy,izp) * vad_m

          phi = ABS(v_advect) * dt / dzbu

          Da =  (2.0_num - phi) * ABS(dfi) / dzci &
              + (1.0_num + phi) * ABS(dfu) / dzcu
          Da = Da * sixth

          ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

          Di = sign_v * ss * MIN(ABS(Da) * dzbu, ABS(dfi), ABS(dfu))

          dm(ix,iy,iz) = (fu + Di * (1.0_num - phi)) * v_advect * dt * area
        END DO
      END DO
    END DO

  END SUBROUTINE z_mass_flux



  ! Energy remap in mass coordinates

  SUBROUTINE z_energy_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dzci, dzcu, dzbu, phi, ss, Da, Di
    REAL(num) :: area, rhou, dmu
    INTEGER :: izp2

    DO iz = 0, nz
      izm  = iz - 1
      izp  = iz + 1
      izp2 = iz + 2
      DO iy = 0, ny
        iym = iy - 1
        DO ix = 0, nx
          ixm = ix - 1
          area = dxb(ix) * dyb(iy)

          v_advect = (vz1(ix,iy ,iz) + vz1(ixm,iy ,iz) &
                    + vz1(ix,iym,iz) + vz1(ixm,iym,iz)) * 0.25_num

          fm  = energy(ix,iy,izm )
          fi  = energy(ix,iy,iz  )
          fp  = energy(ix,iy,izp )
          fp2 = energy(ix,iy,izp2)

          dfm = fi - fm
          dfi = fp - fi
          dfp = fp2 - fp

          ! vad_p and vad_m are logical switches which determine v_advect>=0
          ! and v_advect<0 respectively. It's written this way to allow vector
          ! optimization

          sign_v = SIGN(1.0_num, v_advect)
          vad_p = (sign_v + 1.0_num) * 0.5_num
          vad_m = 1.0_num - vad_p

          fu = fi * vad_p + fp * vad_m
          dfu = dfm * vad_p + dfp * vad_m
          dzci = dzc1(ix,iy,iz )
          dzcu = dzc1(ix,iy,izm) * vad_p + dzc1(ix,iy,izp) * vad_m
          dzbu = dzb1(ix,iy,iz ) * vad_p + dzb1(ix,iy,izp) * vad_m

          phi = ABS(v_advect) * dt / dzbu

          Da =  (2.0_num - phi) * ABS(dfi) / dzci &
              + (1.0_num + phi) * ABS(dfu) / dzcu
          Da = Da * sixth

          ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

          Di = sign_v * ss * MIN(ABS(Da) * dzbu, ABS(dfi), ABS(dfu))

          rhou = rho1(ix,iy,iz) * vad_p + rho1(ix,iy,izp) * vad_m
          dmu = ABS(dm(ix,iy,iz)) / area / dzbu / rhou

          flux(ix,iy,iz) = (fu + Di * (1.0_num - dmu)) * dm(ix,iy,iz)
        END DO
      END DO
    END DO

  END SUBROUTINE z_energy_flux



  SUBROUTINE z_momx_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dzci, dzcu, dzbu, phi, ss, Da, Di
    REAL(num) :: area, rhou, dmu, m, mp, ai, aip, dk
    INTEGER :: izp2

    DO iz = -1, nz
      izm  = iz - 1
      izp  = iz + 1
      izp2 = iz + 2
      DO iy = 0, ny
        DO ix = 0, nx
          area = dxc(ix) * dyc(iy)

          v_advect = vz1(ix,iy,iz)

          fm  = vx(ix,iy,izm )
          fi  = vx(ix,iy,iz  )
          fp  = vx(ix,iy,izp )
          fp2 = vx(ix,iy,izp2)

          dfm = fi - fm
          dfi = fp - fi
          dfp = fp2 - fp

          ! vad_p and vad_m are logical switches which determine v_advect>=0
          ! and v_advect<0 respectively. It's written this way to allow vector
          ! optimization

          sign_v = SIGN(1.0_num, v_advect)
          vad_p = (sign_v + 1.0_num) * 0.5_num
          vad_m = 1.0_num - vad_p

          fu = fi * vad_p + fp * vad_m
          dfu = dfm * vad_p + dfp * vad_m
          dzci = dzb1(ix,iy,izp)
          dzcu = dzb1(ix,iy,iz ) * vad_p + dzb1(ix,iy,izp2) * vad_m
          dzbu = dzc1(ix,iy,iz ) * vad_p + dzc1(ix,iy,izp ) * vad_m

          phi = ABS(v_advect) * dt / dzbu

          Da =  (2.0_num - phi) * ABS(dfi) / dzci &
              + (1.0_num + phi) * ABS(dfu) / dzcu
          Da = Da * sixth

          ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

          Di = sign_v * ss * MIN(ABS(Da) * dzbu, ABS(dfi), ABS(dfu))

          rhou = rho_v(ix,iy,iz) * vad_p + rho_v(ix,iy,izp) * vad_m
          dmu = ABS(dm(ix,iy,iz)) / area / dzbu / rhou

          flux(ix,iy,iz) = fu + Di * (1.0_num - dmu)
        END DO
      END DO
    END DO

    IF (rke) THEN
      DO iz = 0, nz - 1
        izm = iz - 1
        izp = iz + 1
        DO iy = 0, ny
          iyp = iy + 1
          DO ix = 0, nx
            ixp = ix + 1

            m  = rho_v1(ix,iy,iz ) * cv2(ix,iy,iz )
            mp = rho_v1(ix,iy,izp) * cv2(ix,iy,izp)

            ai =  (vx(ix,iy,iz ) - flux(ix,iy,izm)) * dm(ix,iy,izm) / m &
                - (vx(ix,iy,iz ) - flux(ix,iy,iz )) * dm(ix,iy,iz ) / m

            aip = (vx(ix,iy,izp) - flux(ix,iy,iz )) * dm(ix,iy,iz ) / mp &
                - (vx(ix,iy,izp) - flux(ix,iy,izp)) * dm(ix,iy,izp) / mp

            dk = (vx(ix,iy,izp) - vx(ix,iy,iz)) * (flux(ix,iy,iz) &
                - 0.5_num * (vx(ix,iy,izp) + vx(ix,iy,iz))) &
                - 0.5_num * ai  * (vx(ix,iy,iz ) - flux(ix,iy,iz)) &
                + 0.5_num * aip * (vx(ix,iy,izp) - flux(ix,iy,iz))

            dk = dk * dm(ix,iy,iz) * 0.5_num
            delta_ke(ix ,iy ,izp) = delta_ke(ix ,iy ,izp) + dk
            delta_ke(ixp,iy ,izp) = delta_ke(ixp,iy ,izp) + dk
            delta_ke(ix ,iyp,izp) = delta_ke(ix ,iyp,izp) + dk
            delta_ke(ixp,iyp,izp) = delta_ke(ixp,iyp,izp) + dk
          END DO
        END DO
      END DO
    END IF

    flux(0:nx,0:ny,-1:nz) = flux(0:nx,0:ny,-1:nz) * dm(0:nx,0:ny,-1:nz)

  END SUBROUTINE z_momx_flux



  SUBROUTINE z_momy_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dzci, dzcu, dzbu, phi, ss, Da, Di
    REAL(num) :: area, rhou, dmu, m, mp, ai, aip, dk
    INTEGER :: izp2

    DO iz = -1, nz
      izm  = iz - 1
      izp  = iz + 1
      izp2 = iz + 2
      DO iy = 0, ny
        DO ix = 0, nx
          area = dxc(ix) * dyc(iy)

          v_advect = vz1(ix,iy,iz)

          fm  = vy(ix,iy,izm )
          fi  = vy(ix,iy,iz  )
          fp  = vy(ix,iy,izp )
          fp2 = vy(ix,iy,izp2)

          dfm = fi - fm
          dfi = fp - fi
          dfp = fp2 - fp

          ! vad_p and vad_m are logical switches which determine v_advect>=0
          ! and v_advect<0 respectively. It's written this way to allow vector
          ! optimization

          sign_v = SIGN(1.0_num, v_advect)
          vad_p = (sign_v + 1.0_num) * 0.5_num
          vad_m = 1.0_num - vad_p

          fu = fi * vad_p + fp * vad_m
          dfu = dfm * vad_p + dfp * vad_m
          dzci = dzb1(ix,iy,izp)
          dzcu = dzb1(ix,iy,iz ) * vad_p + dzb1(ix,iy,izp2) * vad_m
          dzbu = dzc1(ix,iy,iz ) * vad_p + dzc1(ix,iy,izp ) * vad_m

          phi = ABS(v_advect) * dt / dzbu

          Da =  (2.0_num - phi) * ABS(dfi) / dzci &
              + (1.0_num + phi) * ABS(dfu) / dzcu
          Da = Da * sixth

          ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

          Di = sign_v * ss * MIN(ABS(Da) * dzbu, ABS(dfi), ABS(dfu))

          rhou = rho_v(ix,iy,iz) * vad_p + rho_v(ix,iy,izp) * vad_m
          dmu = ABS(dm(ix,iy,iz)) / area / dzbu / rhou

          flux(ix,iy,iz) = fu + Di * (1.0_num - dmu)
        END DO
      END DO
    END DO

    IF (rke) THEN
      DO iz = 0, nz - 1
        izm = iz - 1
        izp = iz + 1
        DO iy = 0, ny
          iyp = iy + 1
          DO ix = 0, nx
            ixp = ix + 1

            m  = rho_v1(ix,iy,iz ) * cv2(ix,iy,iz )
            mp = rho_v1(ix,iy,izp) * cv2(ix,iy,izp)

            ai =  (vy(ix,iy,iz ) - flux(ix,iy,izm)) * dm(ix,iy,izm) / m &
                - (vy(ix,iy,iz ) - flux(ix,iy,iz )) * dm(ix,iy,iz ) / m

            aip = (vy(ix,iy,izp) - flux(ix,iy,iz )) * dm(ix,iy,iz ) / mp &
                - (vy(ix,iy,izp) - flux(ix,iy,izp)) * dm(ix,iy,izp) / mp

            dk = (vy(ix,iy,izp) - vy(ix,iy,iz)) * (flux(ix,iy,iz) &
                - 0.5_num * (vy(ix,iy,izp) + vy(ix,iy,iz))) &
                - 0.5_num * ai  * (vy(ix,iy,iz ) - flux(ix,iy,iz)) &
                + 0.5_num * aip * (vy(ix,iy,izp) - flux(ix,iy,iz))

            dk = dk * dm(ix,iy,iz) * 0.5_num
            delta_ke(ix ,iy ,izp) = delta_ke(ix ,iy ,izp) + dk
            delta_ke(ixp,iy ,izp) = delta_ke(ixp,iy ,izp) + dk
            delta_ke(ix ,iyp,izp) = delta_ke(ix ,iyp,izp) + dk
            delta_ke(ixp,iyp,izp) = delta_ke(ixp,iyp,izp) + dk
          END DO
        END DO
      END DO
    END IF

    flux(0:nx,0:ny,-1:nz) = flux(0:nx,0:ny,-1:nz) * dm(0:nx,0:ny,-1:nz)

  END SUBROUTINE z_momy_flux



  SUBROUTINE z_momz_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dzci, dzcu, dzbu, phi, ss, Da, Di
    REAL(num) :: area, rhou, dmu, m, mp, ai, aip, dk
    INTEGER :: izp2

    DO iz = -1, nz
      izm  = iz - 1
      izp  = iz + 1
      izp2 = iz + 2
      DO iy = 0, ny
        DO ix = 0, nx
          area = dxc(ix) * dyc(iy)

          v_advect = vz1(ix,iy,iz)

          fm  = vz(ix,iy,izm )
          fi  = vz(ix,iy,iz  )
          fp  = vz(ix,iy,izp )
          fp2 = vz(ix,iy,izp2)

          dfm = fi - fm
          dfi = fp - fi
          dfp = fp2 - fp

          ! vad_p and vad_m are logical switches which determine v_advect>=0
          ! and v_advect<0 respectively. It's written this way to allow vector
          ! optimization

          sign_v = SIGN(1.0_num, v_advect)
          vad_p = (sign_v + 1.0_num) * 0.5_num
          vad_m = 1.0_num - vad_p

          fu = fi * vad_p + fp * vad_m
          dfu = dfm * vad_p + dfp * vad_m
          dzci = dzb1(ix,iy,izp)
          dzcu = dzb1(ix,iy,iz ) * vad_p + dzb1(ix,iy,izp2) * vad_m
          dzbu = dzc1(ix,iy,iz ) * vad_p + dzc1(ix,iy,izp ) * vad_m

          phi = ABS(v_advect) * dt / dzbu

          Da =  (2.0_num - phi) * ABS(dfi) / dzci &
              + (1.0_num + phi) * ABS(dfu) / dzcu
          Da = Da * sixth

          ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

          Di = sign_v * ss * MIN(ABS(Da) * dzbu, ABS(dfi), ABS(dfu))

          rhou = rho_v(ix,iy,iz) * vad_p + rho_v(ix,iy,izp) * vad_m
          dmu = ABS(dm(ix,iy,iz)) / area / dzbu / rhou

          flux(ix,iy,iz) = fu + Di * (1.0_num - dmu)
        END DO
      END DO
    END DO

    IF (rke) THEN
      DO iz = 0, nz - 1
        izm = iz - 1
        izp = iz + 1
        DO iy = 0, ny
          iyp = iy + 1
          DO ix = 0, nx
            ixp = ix + 1

            m  = rho_v1(ix,iy,iz ) * cv2(ix,iy,iz )
            mp = rho_v1(ix,iy,izp) * cv2(ix,iy,izp)

            ai =  (vz(ix,iy,iz ) - flux(ix,iy,izm)) * dm(ix,iy,izm) / m &
                - (vz(ix,iy,iz ) - flux(ix,iy,iz )) * dm(ix,iy,iz ) / m

            aip = (vz(ix,iy,izp) - flux(ix,iy,iz )) * dm(ix,iy,iz ) / mp &
                - (vz(ix,iy,izp) - flux(ix,iy,izp)) * dm(ix,iy,izp) / mp

            dk = (vz(ix,iy,izp) - vz(ix,iy,iz)) * (flux(ix,iy,iz) &
                - 0.5_num * (vz(ix,iy,izp) + vz(ix,iy,iz))) &
                - 0.5_num * ai  * (vz(ix,iy,iz ) - flux(ix,iy,iz)) &
                + 0.5_num * aip * (vz(ix,iy,izp) - flux(ix,iy,iz))

            dk = dk * dm(ix,iy,iz) * 0.5_num
            delta_ke(ix ,iy ,izp) = delta_ke(ix ,iy ,izp) + dk
            delta_ke(ixp,iy ,izp) = delta_ke(ixp,iy ,izp) + dk
            delta_ke(ix ,iyp,izp) = delta_ke(ix ,iyp,izp) + dk
            delta_ke(ixp,iyp,izp) = delta_ke(ixp,iyp,izp) + dk
          END DO
        END DO
      END DO
    END IF

    flux(0:nx,0:ny,-1:nz) = flux(0:nx,0:ny,-1:nz) * dm(0:nx,0:ny,-1:nz)

  END SUBROUTINE z_momz_flux



  SUBROUTINE dm_z_bcs

    CALL MPI_SENDRECV(&
        dm(-1,-1,1   ), 1, cell_zface, proc_z_min, tag, &
        dm(-1,-1,nz+1), 1, cell_zface, proc_z_max, tag, &
        comm, status, errcode)

    IF (proc_z_max == MPI_PROC_NULL) &
        dm(0:nx+1,0:ny+1,nz+1) = dm(0:nx+1,0:ny+1,nz)

    CALL MPI_SENDRECV(&
        dm(-1,-1,nz-1), 1, cell_zface, proc_z_max, tag, &
        dm(-1,-1,-1  ), 1, cell_zface, proc_z_min, tag, &
        comm, status, errcode)

    IF (proc_z_min == MPI_PROC_NULL) &
        dm(0:nx+1,0:ny+1,-1) = dm(0:nx+1,0:ny+1,0)

  END SUBROUTINE dm_z_bcs

END MODULE zremap
