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

MODULE xremap

  USE shared_data; USE boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: remap_x

  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: dxb1, dxc1, rho_v, rho_v1
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: rho1, dm, cv2, flux

CONTAINS

  ! Remap onto original Eulerian grid

  SUBROUTINE remap_x

    REAL(num) :: vxb, vxbm, vyb, vybm, vzb, vzbm, dv

    IF (predictor_step) THEN
      ALLOCATE(cv2   (-1:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(flux  (-2:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(dxb1  (-1:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(dxc1  (-1:nx+2, -1:ny+2, -1:nz+2))
    ELSE
      ALLOCATE(rho1  (-1:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(dm    (-1:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(cv2   (-1:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(flux  (-2:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(dxb1  (-1:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(dxc1  (-1:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(rho_v (-1:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(rho_v1(-1:nx+2, -1:ny+2, -1:nz+2))
    ENDIF

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

          dv = (REAL(ypass, num) * (vyb - vybm) / dyb(iy) &
              + REAL(zpass, num) * (vzb - vzbm) / dzb(iz) &
              + (vxb - vxbm) / dxb(ix)) * dt

          ! Control volume before remap
          cv1(ix,iy,iz) = cv(ix,iy,iz) * (1.0_num + dv)

          dv = (REAL(ypass, num) * (vyb - vybm) / dyb(iy) &
              + REAL(zpass, num) * (vzb - vzbm) / dzb(iz)) * dt

          ! Control volume after remap
          cv2(ix,iy,iz) = cv(ix,iy,iz) * (1.0_num + dv)

          ! dxb before remap
          dxb1(ix,iy,iz) = dxb(ix) + (vxb - vxbm) * dt
        END DO
      END DO
    END DO

    DO iz = -1, nz + 2
      DO iy = -1, ny + 2
        DO ix = -1, nx + 2
          ixp = MIN(ix + 1, nx + 2)
          ! dxc before remap
          dxc1(ix,iy,iz) = 0.5_num * (dxb1(ix,iy,iz) + dxb1(ixp,iy,iz))
        END DO
      END DO
    END DO

    ! Evans and Hawley (ApJ, vol 332, p650, (1988))
    ! constrained transport remap of magnetic fluxes
    CALL vx_by_flux

    DO iz = 0, nz + 1
      DO iy = 0, ny + 1
        DO ix = 0, nx + 1
          ixm = ix - 1
          by(ix,iy,iz) = by(ix,iy,iz) - flux(ix,iy,iz) + flux(ixm,iy,iz)
        END DO
      END DO
    END DO

    DO iz = 0, nz + 1
      DO iy = 0, ny + 1
        iym = iy - 1
        DO ix = 0, nx + 1
          bx(ix,iy,iz) = bx(ix,iy,iz) + flux(ix,iy,iz) - flux(ix,iym,iz)
        END DO
      END DO
    END DO

    CALL vx_bz_flux

    DO iz = 0, nz
      DO iy = 1, ny
        DO ix = 1, nx
          ixm = ix - 1
          bz(ix,iy,iz) = bz(ix,iy,iz) - flux(ix,iy,iz) + flux(ixm,iy,iz)
        END DO
      END DO
    END DO

    DO iz = 0, nz + 1
      izm = iz - 1
      DO iy = 0, ny + 1
        DO ix = 0, nx
          bx(ix,iy,iz) = bx(ix,iy,iz) + flux(ix,iy,iz) - flux(ix,iy,izm)
        END DO
      END DO
    END DO

    IF (predictor_step) THEN
      DEALLOCATE(cv2, flux, dxb1, dxc1)
      RETURN
    END IF

    dm = 0.0_num
    ! Store initial density in rho1
    rho1(:,:,:) = rho(:,:,:)

    ! Remap of mass + calculation of mass fluxes (dm) needed for later remaps
    ! Calculates dm(0:nx,0:ny+1,0:nz+1)
    CALL x_mass_flux
    ! Need dm(-1:nx+1,0:ny+1,0:nz+1) for velocity remap
    CALL dm_x_bcs

    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          ixm = ix - 1
          rho(ix,iy,iz) = (rho1(ix,iy,iz) * cv1(ix,iy,iz) &
              + dm(ixm,iy,iz) - dm(ix,iy,iz)) / cv2(ix,iy,iz)
        END DO
      END DO
    END DO

    ! Remap specific energy density using mass coordinates
    CALL x_energy_flux

    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          ixm = ix - 1
          energy(ix,iy,iz) = (energy(ix,iy,iz) * cv1(ix,iy,iz) &
              * rho1(ix,iy,iz) + flux(ixm,iy,iz) - flux(ix,iy,iz)) &
              / (cv2(ix,iy,iz) * rho(ix,iy,iz))
        END DO
      END DO
    END DO

    ! Redefine dxb1, cv1, cv2, dm and vx1 for velocity (vertex) cells.
    ! In some of these calculations the flux variable is used as a
    ! temporary array
    DO iz = 0, nz
      izp = iz + 1
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = -1, nx + 1
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

    DO iz = 0, nz
      DO iy = 0, ny
        DO ix = -2, nx + 1
          ixp = ix + 1
          flux(ix,iy,iz) = (vx1(ix,iy,iz) + vx1(ixp,iy,iz)) * 0.5_num
        END DO
      END DO
    END DO
    vx1(-2:nx+1,0:ny,0:nz) = flux(-2:nx+1,0:ny,0:nz)

    DO iz = 0, nz
      izp = iz + 1
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = -1, nx
          ixp = ix + 1
          flux(ix,iy,iz) = dm(ix,iy,iz) + dm(ixp,iy,iz) &
              + dm(ix,iyp,iz ) + dm(ixp,iyp,iz ) &
              + dm(ix,iy ,izp) + dm(ixp,iy ,izp) &
              + dm(ix,iyp,izp) + dm(ixp,iyp,izp)
        END DO
      END DO
    END DO
    dm(-1:nx,0:ny,0:nz) = flux(-1:nx,0:ny,0:nz) * 0.125_num

    DO iz = 0, nz
      DO iy = 0, ny
        DO ix = 0, nx
          ixm = ix - 1
          rho_v1(ix,iy,iz) = (rho_v(ix,iy,iz) * cv1(ix,iy,iz) &
              + dm(ixm,iy,iz) - dm(ix,iy,iz)) / cv2(ix,iy,iz)
        END DO
      END DO
    END DO

    CALL x_momx_flux

    DO iz = 0, nz
      DO iy = 0, ny
        DO ix = 0, nx
          ixm = ix - 1
          vx(ix,iy,iz) = (rho_v(ix,iy,iz) * vx(ix,iy,iz) &
              * cv1(ix,iy,iz) + flux(ixm,iy,iz) - flux(ix,iy,iz)) &
              / (cv2(ix,iy,iz) * rho_v1(ix,iy,iz))
        END DO
      END DO
    END DO

    CALL x_momy_flux

    DO iz = 0, nz
      DO iy = 0, ny
        DO ix = 0, nx
          ixm = ix - 1
          vy(ix,iy,iz) = (rho_v(ix,iy,iz) * vy(ix,iy,iz) &
              * cv1(ix,iy,iz) + flux(ixm,iy,iz) - flux(ix,iy,iz)) &
              / (cv2(ix,iy,iz) * rho_v1(ix,iy,iz))
        END DO
      END DO
    END DO

    CALL x_momz_flux

    DO iz = 0, nz
      DO iy = 0, ny
        DO ix = 0, nx
          ixm = ix - 1
          vz(ix,iy,iz) = (rho_v(ix,iy,iz) * vz(ix,iy,iz) &
              * cv1(ix,iy,iz) + flux(ixm,iy,iz) - flux(ix,iy,iz)) &
              / (cv2(ix,iy,iz) * rho_v1(ix,iy,iz))
        END DO
      END DO
    END DO

    CALL boundary_conditions

    DEALLOCATE (rho1, dm, cv2, flux, dxb1, dxc1, rho_v, rho_v1)
    xpass = 0

  END SUBROUTINE remap_x



  SUBROUTINE vx_by_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dxci, dxcu, dxbu, phi, ss, Da, Di
    REAL(num) :: dxu, dbx, dbxp, dbxp2, dbxm
    INTEGER :: ixp2

    DO iz = -1, nz + 1
      izm = iz - 1
      DO iy = -1, ny + 1
        iyp = iy + 1
        DO ix = -1, nx + 1
          ixm  = MAX(ix - 1, -1)
          ixp  = ix + 1
          ixp2 = MIN(ix + 2, nx + 2)

          v_advect = (vx1(ix,iy,iz) + vx1(ix,iy,izm)) * 0.5_num

          dbx   = (dxb1(ix  ,iy,iz) + dxb1(ix  ,iyp,iz)) * 0.5_num
          dbxp  = (dxb1(ixp ,iy,iz) + dxb1(ixp ,iyp,iz)) * 0.5_num
          dbxp2 = (dxb1(ixp2,iy,iz) + dxb1(ixp2,iyp,iz)) * 0.5_num
          dbxm  = (dxb1(ixm ,iy,iz) + dxb1(ixm ,iyp,iz)) * 0.5_num

          fm  = by(ixm,iy,iz) / dbxm
          fi  = by(ix ,iy,iz) / dbx
          fp  = by(ixp,iy,iz) / dbxp
          fp2 = by(ixp2,iy,iz) / dbxp2

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
          dxci = dxc1(ix ,iy,iz)
          dxcu = dxc1(ixm,iy,iz) * vad_p + dxc1(ixp,iy,iz) * vad_m
          dxbu = dxb1(ix ,iy,iz) * vad_p + dxb1(ixp,iy,iz) * vad_m

          dxu = dbx * vad_p + dbxp * vad_m
          phi = ABS(v_advect) * dt / dxu

          Da =  (2.0_num - phi) * ABS(dfi) / dxci &
              + (1.0_num + phi) * ABS(dfu) / dxcu
          Da = Da * sixth

          ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

          Di = sign_v * ss * MIN(ABS(Da) * dxbu, ABS(dfi), ABS(dfu))

          flux(ix,iy,iz) = (fu + Di * (1.0_num - phi)) * v_advect * dt
        END DO
      END DO
    END DO

  END SUBROUTINE vx_by_flux



  SUBROUTINE vx_bz_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dxci, dxcu, dxbu, phi, ss, Da, Di
    REAL(num) :: dxu, dbx, dbxp, dbxp2, dbxm
    INTEGER :: ixp2

    DO iz = -1, nz + 1
      izp = iz + 1
      DO iy = -1, ny + 1
        iym = iy - 1
        DO ix = -1 , nx + 1
          ixm  = MAX(ix - 1, -1)
          ixp  = ix + 1
          ixp2 = MIN(ix + 2, nx+2)

          v_advect = (vx1(ix,iy,iz) + vx1(ix,iym,iz)) * 0.5_num

          dbx   = (dxb1(ix  ,iy,iz) + dxb1(ix  ,iy,izp)) * 0.5_num
          dbxp  = (dxb1(ixp ,iy,iz) + dxb1(ixp ,iy,izp)) * 0.5_num
          dbxp2 = (dxb1(ixp2,iy,iz) + dxb1(ixp2,iy,izp)) * 0.5_num
          dbxm  = (dxb1(ixm ,iy,iz) + dxb1(ixm ,iy,izp)) * 0.5_num

          fm  = bz(ixm,iy,iz) / dbxm
          fi  = bz(ix ,iy,iz) / dbx
          fp  = bz(ixp,iy,iz) / dbxp
          fp2 = bz(ixp2,iy,iz) / dbxp2

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
          dxci = dxc1(ix ,iy,iz)
          dxcu = dxc1(ixm,iy,iz) * vad_p + dxc1(ixp,iy,iz) * vad_m
          dxbu = dxb1(ix ,iy,iz) * vad_p + dxb1(ixp,iy,iz) * vad_m

          dxu = dbx * vad_p + dbxp * vad_m
          phi = ABS(v_advect) * dt / dxu

          Da =  (2.0_num - phi) * ABS(dfi) / dxci &
              + (1.0_num + phi) * ABS(dfu) / dxcu
          Da = Da * sixth

          ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

          Di = sign_v * ss * MIN(ABS(Da) * dxbu, ABS(dfi), ABS(dfu))

          flux(ix,iy,iz) = (fu + Di * (1.0_num - phi)) * v_advect * dt
        END DO
      END DO
    END DO

  END SUBROUTINE vx_bz_flux



  SUBROUTINE x_mass_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dxci, dxcu, dxbu, phi, ss, Da, Di
    REAL(num) :: area
    INTEGER :: ixp2

    DO iz = 0, nz + 1
      izm = iz - 1
      DO iy = 0, ny + 1
        iym = iy - 1
        area = dyb(iy) * dzb(iz)
        DO ix = 0, nx
          ixm  = ix - 1
          ixp  = ix + 1
          ixp2 = ix + 2

          v_advect = (vx1(ix,iy,iz ) + vx1(ix,iym,iz ) &
                    + vx1(ix,iy,izm) + vx1(ix,iym,izm)) * 0.25_num

          fm  = rho(ixm ,iy,iz)
          fi  = rho(ix  ,iy,iz)
          fp  = rho(ixp ,iy,iz)
          fp2 = rho(ixp2,iy,iz)

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
          dxci = dxc1(ix ,iy,iz)
          dxcu = dxc1(ixm,iy,iz) * vad_p + dxc1(ixp,iy,iz) * vad_m
          dxbu = dxb1(ix ,iy,iz) * vad_p + dxb1(ixp,iy,iz) * vad_m

          phi = ABS(v_advect) * dt / dxbu

          Da =  (2.0_num - phi) * ABS(dfi) / dxci &
              + (1.0_num + phi) * ABS(dfu) / dxcu
          Da = Da * sixth

          ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

          Di = sign_v * ss * MIN(ABS(Da) * dxbu, ABS(dfi), ABS(dfu))

          dm(ix,iy,iz) = (fu + Di * (1.0_num - phi)) * v_advect * dt * area
        END DO
      END DO
    END DO

  END SUBROUTINE x_mass_flux



  ! Energy remap in mass coordinates

  SUBROUTINE x_energy_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dxci, dxcu, dxbu, phi, ss, Da, Di
    REAL(num) :: area, rhou, dmu
    INTEGER :: ixp2

    DO iz = 0, nz
      izm = iz - 1
      DO iy = 0, ny
        iym = iy - 1
        area = dyb(iy) * dzb(iz)
        DO ix = 0, nx
          ixm  = ix - 1
          ixp  = ix + 1
          ixp2 = ix + 2

          v_advect = (vx1(ix,iy,iz ) + vx1(ix,iym,iz ) &
                    + vx1(ix,iy,izm) + vx1(ix,iym,izm)) * 0.25_num

          fm  = energy(ixm ,iy,iz)
          fi  = energy(ix  ,iy,iz)
          fp  = energy(ixp ,iy,iz)
          fp2 = energy(ixp2,iy,iz)

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
          dxci = dxc1(ix ,iy,iz)
          dxcu = dxc1(ixm,iy,iz) * vad_p + dxc1(ixp,iy,iz) * vad_m
          dxbu = dxb1(ix ,iy,iz) * vad_p + dxb1(ixp,iy,iz) * vad_m

          phi = ABS(v_advect) * dt / dxbu

          Da =  (2.0_num - phi) * ABS(dfi) / dxci &
              + (1.0_num + phi) * ABS(dfu) / dxcu
          Da = Da * sixth

          ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

          Di = sign_v * ss * MIN(ABS(Da) * dxbu, ABS(dfi), ABS(dfu))

          rhou = rho1(ix,iy,iz) * vad_p + rho1(ixp,iy,iz) * vad_m
          dmu = ABS(dm(ix,iy,iz)) / area / dxbu / rhou

          flux(ix,iy,iz) = (fu + Di * (1.0_num - dmu)) * dm(ix,iy,iz)
        END DO
      END DO
    END DO

  END SUBROUTINE x_energy_flux



  SUBROUTINE x_momx_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dxci, dxcu, dxbu, phi, ss, Da, Di
    REAL(num) :: area, rhou, dmu, m, mp, ai, aip, dk
    INTEGER :: ixp2

    DO iz = 0, nz
      DO iy = 0, ny
        area = dyc(iy) * dzc(iz)
        DO ix = -1, nx
          ixm  = ix - 1
          ixp  = ix + 1
          ixp2 = ix + 2

          v_advect = vx1(ix,iy,iz)

          fm  = vx(ixm ,iy,iz)
          fi  = vx(ix  ,iy,iz)
          fp  = vx(ixp ,iy,iz)
          fp2 = vx(ixp2,iy,iz)

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
          dxci = dxb1(ixp,iy,iz)
          dxcu = dxb1(ix ,iy,iz) * vad_p + dxb1(ixp2,iy,iz) * vad_m
          dxbu = dxc1(ix ,iy,iz) * vad_p + dxc1(ixp ,iy,iz) * vad_m

          phi = ABS(v_advect) * dt / dxbu

          Da =  (2.0_num - phi) * ABS(dfi) / dxci &
              + (1.0_num + phi) * ABS(dfu) / dxcu
          Da = Da * sixth

          ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

          Di = sign_v * ss * MIN(ABS(Da) * dxbu, ABS(dfi), ABS(dfu))

          rhou = rho_v(ix,iy,iz) * vad_p + rho_v(ixp,iy,iz) * vad_m
          dmu = ABS(dm(ix,iy,iz)) / area / dxbu / rhou

          flux(ix,iy,iz) = fu + Di * (1.0_num - dmu)
        END DO
      END DO
    END DO

    IF (rke) THEN
      DO iz = 0, nz
        izp = iz + 1
        DO iy = 0, ny
          iyp = iy + 1
          DO ix = 0, nx - 1
            ixm = ix - 1
            ixp = ix + 1

            m  = rho_v1(ix ,iy,iz) * cv2(ix ,iy,iz)
            mp = rho_v1(ixp,iy,iz) * cv2(ixp,iy,iz)

            ai =  (vx(ix ,iy,iz) - flux(ixm,iy,iz)) * dm(ixm,iy,iz) / m &
                - (vx(ix ,iy,iz) - flux(ix ,iy,iz)) * dm(ix ,iy,iz) / m

            aip = (vx(ixp,iy,iz) - flux(ix ,iy,iz)) * dm(ix ,iy,iz) / mp &
                - (vx(ixp,iy,iz) - flux(ixp,iy,iz)) * dm(ixp,iy,iz) / mp

            dk = (vx(ixp,iy,iz) - vx(ix,iy,iz)) * (flux(ix,iy,iz) &
                - 0.5_num * (vx(ixp,iy,iz) + vx(ix,iy,iz))) &
                - 0.5_num * ai  * (vx(ix ,iy,iz) - flux(ix,iy,iz)) &
                + 0.5_num * aip * (vx(ixp,iy,iz) - flux(ix,iy,iz))

            dk = dk * dm(ix,iy,iz) * 0.5_num
            delta_ke(ixp,iy ,iz ) = delta_ke(ixp,iy ,iz ) + dk
            delta_ke(ixp,iyp,iz ) = delta_ke(ixp,iyp,iz ) + dk
            delta_ke(ixp,iy ,izp) = delta_ke(ixp,iy ,izp) + dk
            delta_ke(ixp,iyp,izp) = delta_ke(ixp,iyp,izp) + dk
          END DO
        END DO
      END DO
    END IF

    flux(-1:nx,0:ny,0:nz) = flux(-1:nx,0:ny,0:nz) * dm(-1:nx,0:ny,0:nz)

  END SUBROUTINE x_momx_flux



  SUBROUTINE x_momy_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dxci, dxcu, dxbu, phi, ss, Da, Di
    REAL(num) :: area, rhou, dmu, m, mp, ai, aip, dk
    INTEGER :: ixp2

    DO iz = 0, nz
      DO iy = 0, ny
        area = dyc(iy) * dzc(iz)
        DO ix = -1, nx
          ixm  = ix - 1
          ixp  = ix + 1
          ixp2 = ix + 2

          v_advect = vx1(ix,iy,iz)

          fm  = vy(ixm ,iy,iz)
          fi  = vy(ix  ,iy,iz)
          fp  = vy(ixp ,iy,iz)
          fp2 = vy(ixp2,iy,iz)

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
          dxci = dxb1(ixp,iy,iz)
          dxcu = dxb1(ix ,iy,iz) * vad_p + dxb1(ixp2,iy,iz) * vad_m
          dxbu = dxc1(ix ,iy,iz) * vad_p + dxc1(ixp ,iy,iz) * vad_m

          phi = ABS(v_advect) * dt / dxbu

          Da =  (2.0_num - phi) * ABS(dfi) / dxci &
              + (1.0_num + phi) * ABS(dfu) / dxcu
          Da = Da * sixth

          ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

          Di = sign_v * ss * MIN(ABS(Da) * dxbu, ABS(dfi), ABS(dfu))

          rhou = rho_v(ix,iy,iz) * vad_p + rho_v(ixp,iy,iz) * vad_m
          dmu = ABS(dm(ix,iy,iz)) / area / dxbu / rhou

          flux(ix,iy,iz) = fu + Di * (1.0_num - dmu)
        END DO
      END DO
    END DO

    IF (rke) THEN
      DO iz = 0, nz
        izp = iz + 1
        DO iy = 0, ny
          iyp = iy + 1
          DO ix = 0, nx - 1
            ixm = ix - 1
            ixp = ix + 1

            m  = rho_v1(ix ,iy,iz) * cv2(ix ,iy,iz)
            mp = rho_v1(ixp,iy,iz) * cv2(ixp,iy,iz)

            ai =  (vy(ix ,iy,iz) - flux(ixm,iy,iz)) * dm(ixm,iy,iz) / m &
                - (vy(ix ,iy,iz) - flux(ix ,iy,iz)) * dm(ix ,iy,iz) / m

            aip = (vy(ixp,iy,iz) - flux(ix ,iy,iz)) * dm(ix ,iy,iz) / mp &
                - (vy(ixp,iy,iz) - flux(ixp,iy,iz)) * dm(ixp,iy,iz) / mp

            dk = (vy(ixp,iy,iz) - vy(ix,iy,iz)) * (flux(ix,iy,iz) &
                - 0.5_num * (vy(ixp,iy,iz) + vy(ix,iy,iz))) &
                - 0.5_num * ai  * (vy(ix ,iy,iz) - flux(ix,iy,iz)) &
                + 0.5_num * aip * (vy(ixp,iy,iz) - flux(ix,iy,iz))

            dk = dk * dm(ix,iy,iz) * 0.5_num
            delta_ke(ixp,iy ,iz ) = delta_ke(ixp,iy ,iz ) + dk
            delta_ke(ixp,iyp,iz ) = delta_ke(ixp,iyp,iz ) + dk
            delta_ke(ixp,iy ,izp) = delta_ke(ixp,iy ,izp) + dk
            delta_ke(ixp,iyp,izp) = delta_ke(ixp,iyp,izp) + dk
          END DO
        END DO
      END DO
    END IF

    flux(-1:nx,0:ny,0:nz) = flux(-1:nx,0:ny,0:nz) * dm(-1:nx,0:ny,0:nz)

  END SUBROUTINE x_momy_flux



  SUBROUTINE x_momz_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dxci, dxcu, dxbu, phi, ss, Da, Di
    REAL(num) :: area, rhou, dmu, m, mp, ai, aip, dk
    INTEGER :: ixp2

    DO iz = 0, nz
      DO iy = 0, ny
        area = dyc(iy) * dzc(iz)
        DO ix = -1, nx
          ixm  = ix - 1
          ixp  = ix + 1
          ixp2 = ix + 2

          v_advect = vx1(ix,iy,iz)

          fm  = vz(ixm ,iy,iz)
          fi  = vz(ix  ,iy,iz)
          fp  = vz(ixp ,iy,iz)
          fp2 = vz(ixp2,iy,iz)

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
          dxci = dxb1(ixp,iy,iz)
          dxcu = dxb1(ix ,iy,iz) * vad_p + dxb1(ixp2,iy,iz) * vad_m
          dxbu = dxc1(ix ,iy,iz) * vad_p + dxc1(ixp ,iy,iz) * vad_m

          phi = ABS(v_advect) * dt / dxbu

          Da =  (2.0_num - phi) * ABS(dfi) / dxci &
              + (1.0_num + phi) * ABS(dfu) / dxcu
          Da = Da * sixth

          ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

          Di = sign_v * ss * MIN(ABS(Da) * dxbu, ABS(dfi), ABS(dfu))

          rhou = rho_v(ix,iy,iz) * vad_p + rho_v(ixp,iy,iz) * vad_m
          dmu = ABS(dm(ix,iy,iz)) / area / dxbu / rhou

          flux(ix,iy,iz) = fu + Di * (1.0_num - dmu)
        END DO
      END DO
    END DO

    IF (rke) THEN
      DO iz = 0, nz
        izp = iz + 1
        DO iy = 0, ny
          iyp = iy + 1
          DO ix = 0, nx - 1
            ixm = ix - 1
            ixp = ix + 1

            m  = rho_v1(ix ,iy,iz) * cv2(ix ,iy,iz)
            mp = rho_v1(ixp,iy,iz) * cv2(ixp,iy,iz)

            ai =  (vz(ix ,iy,iz) - flux(ixm,iy,iz)) * dm(ixm,iy,iz) / m &
                - (vz(ix ,iy,iz) - flux(ix ,iy,iz)) * dm(ix ,iy,iz) / m

            aip = (vz(ixp,iy,iz) - flux(ix ,iy,iz)) * dm(ix ,iy,iz) / mp &
                - (vz(ixp,iy,iz) - flux(ixp,iy,iz)) * dm(ixp,iy,iz) / mp

            dk = (vz(ixp,iy,iz) - vz(ix,iy,iz)) * (flux(ix,iy,iz) &
                - 0.5_num * (vz(ixp,iy,iz) + vz(ix,iy,iz))) &
                - 0.5_num * ai  * (vz(ix ,iy,iz) - flux(ix,iy,iz)) &
                + 0.5_num * aip * (vz(ixp,iy,iz) - flux(ix,iy,iz))

            dk = dk * dm(ix,iy,iz) * 0.5_num
            delta_ke(ixp,iy ,iz ) = delta_ke(ixp,iy ,iz ) + dk
            delta_ke(ixp,iyp,iz ) = delta_ke(ixp,iyp,iz ) + dk
            delta_ke(ixp,iy ,izp) = delta_ke(ixp,iy ,izp) + dk
            delta_ke(ixp,iyp,izp) = delta_ke(ixp,iyp,izp) + dk
          END DO
        END DO
      END DO
    END IF

    flux(-1:nx,0:ny,0:nz) = flux(-1:nx,0:ny,0:nz) * dm(-1:nx,0:ny,0:nz)

  END SUBROUTINE x_momz_flux



  SUBROUTINE dm_x_bcs

    CALL MPI_SENDRECV(&
        dm(1   ,-1,-1), 1, cell_xface, proc_x_min, tag, &
        dm(nx+1,-1,-1), 1, cell_xface, proc_x_max, tag, &
        comm, status, errcode)

    IF (proc_x_max == MPI_PROC_NULL) &
        dm(nx+1,0:ny+1,0:nz+1) = dm(nx,0:ny+1,0:nz+1)

    CALL MPI_SENDRECV(&
        dm(nx-1,-1,-1), 1, cell_xface, proc_x_max, tag, &
        dm(-1  ,-1,-1), 1, cell_xface, proc_x_min, tag, &
        comm, status, errcode)

    IF (proc_x_min == MPI_PROC_NULL) &
        dm(-1,0:ny+1,0:nz+1) = dm(0,0:ny+1,0:nz+1)

  END SUBROUTINE dm_x_bcs

END MODULE xremap
