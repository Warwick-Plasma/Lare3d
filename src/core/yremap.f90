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

MODULE yremap

  USE shared_data; USE boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: remap_y

  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: dyb1, dyc1, rho_v, rho_v1
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: rho1, dm, cv2, flux

CONTAINS

  ! Remap onto original Eulerian grid

  SUBROUTINE remap_y

    REAL(num) :: vxb, vxbm, vyb, vybm, vzb, vzbm, dv

    IF (predictor_step) THEN
      ALLOCATE(cv2   (-1:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(flux  (-1:nx+2, -2:ny+2, -1:nz+2))
      ALLOCATE(dyb1  (-1:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(dyc1  (-1:nx+2, -1:ny+2, -1:nz+2))
    ELSE
      ALLOCATE(rho1  (-1:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(dm    (-1:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(cv2   (-1:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(flux  (-1:nx+2, -2:ny+2, -1:nz+2))
      ALLOCATE(dyb1  (-1:nx+2, -1:ny+2, -1:nz+2))
      ALLOCATE(dyc1  (-1:nx+2, -1:ny+2, -1:nz+2))
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
              + REAL(zpass, num) * (vzb - vzbm) / dzb(iz) &
              + (vyb - vybm) / dyb(iy)) * dt

          ! Control volume before remap
          cv1(ix,iy,iz) = cv(ix,iy,iz) * (1.0_num + dv)

          dv = (REAL(xpass, num) * (vxb - vxbm) / dxb(ix) &
              + REAL(zpass, num) * (vzb - vzbm) / dzb(iz)) * dt

          ! Control volume after remap
          cv2(ix,iy,iz) = cv(ix,iy,iz) * (1.0_num + dv)

          ! dyb before remap
          dyb1(ix,iy,iz) = dyb(iy) + (vyb - vybm) * dt
        END DO
      END DO
    END DO

    DO iz = -1, nz + 2
      DO iy = -1, ny + 2
        iyp = MIN(iy + 1, ny + 2)
        DO ix = -1, nx + 2
          ! dyc before remap
          dyc1(ix,iy,iz) = 0.5_num * (dyb1(ix,iy,iz) + dyb1(ix,iyp,iz))
        END DO
      END DO
    END DO

    ! Evans and Hawley (ApJ, vol 332, p650, (1988))
    ! constrained transport remap of magnetic fluxes
    CALL vy_bx_flux

    DO iz = 0, nz + 1
      DO iy = 0, ny + 1
        iym = iy - 1
        DO ix = 0, nx + 1
          bx(ix,iy,iz) = bx(ix,iy,iz) - flux(ix,iy,iz) + flux(ix,iym,iz)
        END DO
      END DO
    END DO

    DO iz = 0, nz + 1
      DO iy = 0, ny + 1
        DO ix = 0, nx + 1
          ixm = ix - 1
          by(ix,iy,iz) = by(ix,iy,iz) + flux(ix,iy,iz) - flux(ixm,iy,iz)
        END DO
      END DO
    END DO

    CALL vy_bz_flux

    DO iz = 0, nz + 1
      DO iy = 0, ny + 1
        iym = iy - 1
        DO ix = 0, nx + 1
          bz(ix,iy,iz) = bz(ix,iy,iz) - flux(ix,iy,iz) + flux(ix,iym,iz)
        END DO
      END DO
    END DO

    DO iz = 1, nz
      izm = iz - 1
      DO iy = 0, ny
        DO ix = 1, nx
          by(ix,iy,iz) = by(ix,iy,iz) + flux(ix,iy,iz) - flux(ix,iy,izm)
        END DO
      END DO
    END DO

    IF (predictor_step) THEN
      DEALLOCATE(cv2, flux, dyb1, dyc1)
      RETURN
    END IF

    dm = 0.0_num
    ! Store initial density in rho1
    rho1(:,:,:) = rho(:,:,:)

    ! Remap of mass + calculation of mass fluxes (dm) needed for later remaps
    ! Calculates dm(0:nx+1,0:ny,0:nz+1)
    CALL y_mass_flux
    ! Need dm(0:nx+1,-1:ny+1,0:nz+1) for velocity remap
    CALL dm_y_bcs

    DO iz = 1, nz
      DO iy = 1, ny
        iym = iy - 1
        DO ix = 1, nx
          rho(ix,iy,iz) = (rho1(ix,iy,iz) * cv1(ix,iy,iz) &
              + dm(ix,iym,iz) - dm(ix,iy,iz)) / cv2(ix,iy,iz)
        END DO
      END DO
    END DO

    ! Remap specific energy density using mass coordinates
    CALL y_energy_flux

    DO iz = 1, nz
      DO iy = 1, ny
        iym = iy - 1
        DO ix = 1, nx
          energy(ix,iy,iz) = (energy(ix,iy,iz) * cv1(ix,iy,iz) &
              * rho1(ix,iy,iz) + flux(ix,iym,iz) - flux(ix,iy,iz)) &
              / (cv2(ix,iy,iz) * rho(ix,iy,iz))
        END DO
      END DO
    END DO

    ! Redefine dyb1, cv1, cv2, dm and vy1 for velocity (vertex) cells.
    ! In some of these calculations the flux variable is used as a
    ! temporary array
    DO iz = 0, nz
      izp = iz + 1
      DO iy = -1, ny + 1
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

    DO iz = 0, nz
      DO iy = -2, ny + 1
        iyp = iy + 1
        DO ix = 0, nx
          flux(ix,iy,iz) = (vy1(ix,iy,iz) + vy1(ix,iyp,iz)) * 0.5_num
        END DO
      END DO
    END DO
    vy1(0:nx,-2:ny+1,0:nz) = flux(0:nx,-2:ny+1,0:nz)

    DO iz = 0, nz
      izp = iz + 1
      DO iy = -1, ny
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
    dm(0:nx,-1:ny,0:nz) = flux(0:nx,-1:ny,0:nz) * 0.125_num

    DO iz = 0, nz
      DO iy = 0, ny
        iym = iy - 1
        DO ix = 0, nx
          rho_v1(ix,iy,iz) = (rho_v(ix,iy,iz) * cv1(ix,iy,iz) &
              + dm(ix,iym,iz) - dm(ix,iy,iz)) / cv2(ix,iy,iz)
        END DO
      END DO
    END DO

    CALL y_momx_flux

    DO iz = 0, nz
      DO iy = 0, ny
        iym = iy - 1
        DO ix = 0, nx
          vx(ix,iy,iz) = (rho_v(ix,iy,iz) * vx(ix,iy,iz) &
              * cv1(ix,iy,iz) + flux(ix,iym,iz) - flux(ix,iy,iz)) &
              / (cv2(ix,iy,iz) * rho_v1(ix,iy,iz))
        END DO
      END DO
    END DO

    CALL y_momy_flux

    DO iz = 0, nz
      DO iy = 0, ny
        iym = iy - 1
        DO ix = 0, nx
          vy(ix,iy,iz) = (rho_v(ix,iy,iz) * vy(ix,iy,iz) &
              * cv1(ix,iy,iz) + flux(ix,iym,iz) - flux(ix,iy,iz)) &
              / (cv2(ix,iy,iz) * rho_v1(ix,iy,iz))
        END DO
      END DO
    END DO

    CALL y_momz_flux

    DO iz = 0, nz
      DO iy = 0, ny
        iym = iy - 1
        DO ix = 0, nx
          vz(ix,iy,iz) = (rho_v(ix,iy,iz) * vz(ix,iy,iz) &
              * cv1(ix,iy,iz) + flux(ix,iym,iz) - flux(ix,iy,iz)) &
              / (cv2(ix,iy,iz) * rho_v1(ix,iy,iz))
        END DO
      END DO
    END DO

    CALL boundary_conditions

    DEALLOCATE (rho1, dm, cv2, flux, dyb1, dyc1, rho_v, rho_v1)
    ypass = 0

  END SUBROUTINE remap_y



  SUBROUTINE vy_bx_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dyci, dycu, dybu, phi, ss, Da, Di
    REAL(num) :: dyu, dby, dbyp, dbyp2, dbym
    INTEGER :: iyp2

    DO iz = -1, nz + 1
      izm = iz - 1
      DO iy = -1, ny +1
        iym  = MAX(iy - 1, -1)
        iyp  = iy + 1
        iyp2 = MIN(iy + 2, ny + 2)
        DO ix = -1, nx + 1
          ixp = ix + 1

          v_advect = (vy1(ix,iy,iz) + vy1(ix,iy,izm)) * 0.5_num

          dby   = (dyb1(ix,iy  ,iz) + dyb1(ixp,iy  ,iz)) * 0.5_num
          dbyp  = (dyb1(ix,iyp ,iz) + dyb1(ixp,iyp ,iz)) * 0.5_num
          dbyp2 = (dyb1(ix,iyp2,iz) + dyb1(ixp,iyp2,iz)) * 0.5_num
          dbym  = (dyb1(ix,iym ,iz) + dyb1(ixp,iym ,iz)) * 0.5_num

          fm  = bx(ix,iym ,iz) / dbym
          fi  = bx(ix,iy  ,iz) / dby
          fp  = bx(ix,iyp ,iz) / dbyp
          fp2 = bx(ix,iyp2,iz) / dbyp2

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
          dyci = dyc1(ix,iy ,iz)
          dycu = dyc1(ix,iym,iz) * vad_p + dyc1(ix,iyp,iz) * vad_m
          dybu = dyb1(ix,iy ,iz) * vad_p + dyb1(ix,iyp,iz) * vad_m

          dyu = dby * vad_p + dbyp * vad_m
          phi = ABS(v_advect) * dt / dyu

          Da =  (2.0_num - phi) * ABS(dfi) / dyci &
              + (1.0_num + phi) * ABS(dfu) / dycu
          Da = Da * sixth

          ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

          Di = sign_v * ss * MIN(ABS(Da) * dybu, ABS(dfi), ABS(dfu))

          flux(ix,iy,iz) = (fu + Di * (1.0_num - phi)) * v_advect * dt
        END DO
      END DO
    END DO

  END SUBROUTINE vy_bx_flux



  SUBROUTINE vy_bz_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dyci, dycu, dybu, phi, ss, Da, Di
    REAL(num) :: dyu, dby, dbyp, dbyp2, dbym
    INTEGER :: iyp2

    DO iz = -1, nz + 1
      izp = iz + 1
      DO iy = -1, ny + 1
        iym  = MAX(iy - 1, -1)
        iyp  = iy + 1
        iyp2 = MIN(iy + 2, ny + 2)
        DO ix = -1, nx + 1
          ixm = ix - 1

          v_advect = (vy1(ix,iy,iz) + vy1(ixm,iy,iz)) * 0.5_num

          dby   = (dyb1(ix,iy  ,iz) + dyb1(ix,iy  ,izp)) * 0.5_num
          dbyp  = (dyb1(ix,iyp ,iz) + dyb1(ix,iyp ,izp)) * 0.5_num
          dbyp2 = (dyb1(ix,iyp2,iz) + dyb1(ix,iyp2,izp)) * 0.5_num
          dbym  = (dyb1(ix,iym ,iz) + dyb1(ix,iym ,izp)) * 0.5_num

          fm  = bz(ix,iym ,iz) / dbym
          fi  = bz(ix,iy  ,iz) / dby
          fp  = bz(ix,iyp ,iz) / dbyp
          fp2 = bz(ix,iyp2,iz) / dbyp2

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
          dyci = dyc1(ix,iy ,iz)
          dycu = dyc1(ix,iym,iz) * vad_p + dyc1(ix,iyp,iz) * vad_m
          dybu = dyb1(ix,iy ,iz) * vad_p + dyb1(ix,iyp,iz) * vad_m

          dyu = dby * vad_p + dbyp * vad_m
          phi = ABS(v_advect) * dt / dyu

          Da =  (2.0_num - phi) * ABS(dfi) / dyci &
              + (1.0_num + phi) * ABS(dfu) / dycu
          Da = Da * sixth

          ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

          Di = sign_v * ss * MIN(ABS(Da) * dybu, ABS(dfi), ABS(dfu))

          flux(ix,iy,iz) = (fu + Di * (1.0_num - phi)) * v_advect * dt
        END DO
      END DO
    END DO

  END SUBROUTINE vy_bz_flux



  SUBROUTINE y_mass_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dyci, dycu, dybu, phi, ss, Da, Di
    REAL(num) :: area
    INTEGER :: iyp2

    DO iz = 0, nz + 1
      izm = iz - 1
      DO iy = 0, ny
        iym  = iy - 1
        iyp  = iy + 1
        iyp2 = iy + 2
        DO ix = 0, nx + 1
          ixm = ix - 1
          area = dxb(ix) * dzb(iz)

          v_advect = (vy1(ix,iy,iz ) + vy1(ixm,iy,iz ) &
                    + vy1(ix,iy,izm) + vy1(ixm,iy,izm)) * 0.25_num

          fm  = rho(ix,iym ,iz)
          fi  = rho(ix,iy  ,iz)
          fp  = rho(ix,iyp ,iz)
          fp2 = rho(ix,iyp2,iz)

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
          dyci = dyc1(ix,iy ,iz)
          dycu = dyc1(ix,iym,iz) * vad_p + dyc1(ix,iyp,iz) * vad_m
          dybu = dyb1(ix,iy ,iz) * vad_p + dyb1(ix,iyp,iz) * vad_m

          phi = ABS(v_advect) * dt / dybu

          Da =  (2.0_num - phi) * ABS(dfi) / dyci &
              + (1.0_num + phi) * ABS(dfu) / dycu
          Da = Da * sixth

          ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

          Di = sign_v * ss * MIN(ABS(Da) * dybu, ABS(dfi), ABS(dfu))

          dm(ix,iy,iz) = (fu + Di * (1.0_num - phi)) * v_advect * dt * area
        END DO
      END DO
    END DO

  END SUBROUTINE y_mass_flux



  ! Energy remap in mass coordinates

  SUBROUTINE y_energy_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dyci, dycu, dybu, phi, ss, Da, Di
    REAL(num) :: area, rhou, dmu
    INTEGER :: iyp2

    DO iz = 0, nz
      izm = iz - 1
      DO iy = 0, ny
        iym  = iy - 1
        iyp  = iy + 1
        iyp2 = iy + 2
        DO ix = 0, nx
          ixm = ix - 1
          area = dxb(ix) * dzb(iz)

          v_advect = (vy1(ix,iy,iz ) + vy1(ixm,iy,iz ) &
                    + vy1(ix,iy,izm) + vy1(ixm,iy,izm)) * 0.25_num

          fm  = energy(ix,iym ,iz)
          fi  = energy(ix,iy  ,iz)
          fp  = energy(ix,iyp ,iz)
          fp2 = energy(ix,iyp2,iz)

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
          dyci = dyc1(ix,iy ,iz)
          dycu = dyc1(ix,iym,iz) * vad_p + dyc1(ix,iyp,iz) * vad_m
          dybu = dyb1(ix,iy ,iz) * vad_p + dyb1(ix,iyp,iz) * vad_m

          phi = ABS(v_advect) * dt / dybu

          Da =  (2.0_num - phi) * ABS(dfi) / dyci &
              + (1.0_num + phi) * ABS(dfu) / dycu
          Da = Da * sixth

          ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

          Di = sign_v * ss * MIN(ABS(Da) * dybu, ABS(dfi), ABS(dfu))

          rhou = rho1(ix,iy,iz) * vad_p + rho1(ix,iyp,iz) * vad_m
          dmu = ABS(dm(ix,iy,iz)) / area / dybu / rhou

          flux(ix,iy,iz) = (fu + Di * (1.0_num - dmu)) * dm(ix,iy,iz)
        END DO
      END DO
    END DO

  END SUBROUTINE y_energy_flux



  SUBROUTINE y_momx_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dyci, dycu, dybu, phi, ss, Da, Di
    REAL(num) :: area, rhou, dmu, m, mp, ai, aip, dk
    INTEGER :: iyp2

    DO iz = 0, nz
      DO iy = -1, ny
        iym  = iy - 1
        iyp  = iy + 1
        iyp2 = iy + 2
        DO ix = 0, nx
          area = dxc(ix) * dzc(iz)

          v_advect = vy1(ix,iy,iz)

          fm  = vx(ix,iym ,iz)
          fi  = vx(ix,iy  ,iz)
          fp  = vx(ix,iyp ,iz)
          fp2 = vx(ix,iyp2,iz)

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
          dyci = dyb1(ix,iyp,iz)
          dycu = dyb1(ix,iy ,iz) * vad_p + dyb1(ix,iyp2,iz) * vad_m
          dybu = dyc1(ix,iy ,iz) * vad_p + dyc1(ix,iyp ,iz) * vad_m

          phi = ABS(v_advect) * dt / dybu

          Da =  (2.0_num - phi) * ABS(dfi) / dyci &
              + (1.0_num + phi) * ABS(dfu) / dycu
          Da = Da * sixth

          ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

          Di = sign_v * ss * MIN(ABS(Da) * dybu, ABS(dfi), ABS(dfu))

          rhou = rho_v(ix,iy,iz) * vad_p + rho_v(ix,iyp,iz) * vad_m
          dmu = ABS(dm(ix,iy,iz)) / area / dybu / rhou

          flux(ix,iy,iz) = fu + Di * (1.0_num - dmu)
        END DO
      END DO
    END DO

    IF (rke) THEN
      DO iz = 0, nz
        izp = iz + 1
        DO iy = 0, ny - 1
          iym = iy - 1
          iyp = iy + 1
          DO ix = 0, nx
            ixp = ix + 1

            m  = rho_v1(ix,iy ,iz) * cv2(ix,iy ,iz)
            mp = rho_v1(ix,iyp,iz) * cv2(ix,iyp,iz)

            ai =  (vx(ix,iy ,iz) - flux(ix,iym,iz)) * dm(ix,iym,iz) / m &
                - (vx(ix,iy ,iz) - flux(ix,iy ,iz)) * dm(ix,iy ,iz) / m

            aip = (vx(ix,iyp,iz) - flux(ix,iy ,iz)) * dm(ix,iy ,iz) / mp &
                - (vx(ix,iyp,iz) - flux(ix,iyp,iz)) * dm(ix,iyp,iz) / mp

            dk = (vx(ix,iyp,iz) - vx(ix,iy,iz)) * (flux(ix,iy,iz) &
                - 0.5_num * (vx(ix,iyp,iz) + vx(ix,iy,iz))) &
                - 0.5_num * ai  * (vx(ix,iy ,iz) - flux(ix,iy,iz)) &
                + 0.5_num * aip * (vx(ix,iyp,iz) - flux(ix,iy,iz))

            dk = dk * dm(ix,iy,iz) * 0.5_num
            delta_ke(ix ,iyp,iz ) = delta_ke(ix ,iyp,iz ) + dk
            delta_ke(ixp,iyp,iz ) = delta_ke(ixp,iyp,iz ) + dk
            delta_ke(ix ,iyp,izp) = delta_ke(ix ,iyp,izp) + dk
            delta_ke(ixp,iyp,izp) = delta_ke(ixp,iyp,izp) + dk
          END DO
        END DO
      END DO
    END IF

    flux(0:nx,-1:ny,0:nz) = flux(0:nx,-1:ny,0:nz) * dm(0:nx,-1:ny,0:nz)

  END SUBROUTINE y_momx_flux



  SUBROUTINE y_momy_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dyci, dycu, dybu, phi, ss, Da, Di
    REAL(num) :: area, rhou, dmu, m, mp, ai, aip, dk
    INTEGER :: iyp2

    DO iz = 0, nz
      DO iy = -1, ny
        iym  = iy - 1
        iyp  = iy + 1
        iyp2 = iy + 2
        DO ix = 0, nx
          area = dxc(ix) * dzc(iz)

          v_advect = vy1(ix,iy,iz)

          fm  = vy(ix,iym ,iz)
          fi  = vy(ix,iy  ,iz)
          fp  = vy(ix,iyp ,iz)
          fp2 = vy(ix,iyp2,iz)

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
          dyci = dyb1(ix,iyp,iz)
          dycu = dyb1(ix,iy ,iz) * vad_p + dyb1(ix,iyp2,iz) * vad_m
          dybu = dyc1(ix,iy ,iz) * vad_p + dyc1(ix,iyp ,iz) * vad_m

          phi = ABS(v_advect) * dt / dybu

          Da =  (2.0_num - phi) * ABS(dfi) / dyci &
              + (1.0_num + phi) * ABS(dfu) / dycu
          Da = Da * sixth

          ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

          Di = sign_v * ss * MIN(ABS(Da) * dybu, ABS(dfi), ABS(dfu))

          rhou = rho_v(ix,iy,iz) * vad_p + rho_v(ix,iyp,iz) * vad_m
          dmu = ABS(dm(ix,iy,iz)) / area / dybu / rhou

          flux(ix,iy,iz) = fu + Di * (1.0_num - dmu)
        END DO
      END DO
    END DO

    IF (rke) THEN
      DO iz = 0, nz
        izp = iz + 1
        DO iy = 0, ny - 1
          iym = iy - 1
          iyp = iy + 1
          DO ix = 0, nx
            ixp = ix + 1

            m  = rho_v1(ix,iy ,iz) * cv2(ix,iy ,iz)
            mp = rho_v1(ix,iyp,iz) * cv2(ix,iyp,iz)

            ai =  (vy(ix,iy ,iz) - flux(ix,iym,iz)) * dm(ix,iym,iz) / m &
                - (vy(ix,iy ,iz) - flux(ix,iy ,iz)) * dm(ix,iy ,iz) / m

            aip = (vy(ix,iyp,iz) - flux(ix,iy ,iz)) * dm(ix,iy ,iz) / mp &
                - (vy(ix,iyp,iz) - flux(ix,iyp,iz)) * dm(ix,iyp,iz) / mp

            dk = (vy(ix,iyp,iz) - vy(ix,iy,iz)) * (flux(ix,iy,iz) &
                - 0.5_num * (vy(ix,iyp,iz) + vy(ix,iy,iz))) &
                - 0.5_num * ai  * (vy(ix,iy ,iz) - flux(ix,iy,iz)) &
                + 0.5_num * aip * (vy(ix,iyp,iz) - flux(ix,iy,iz))

            dk = dk * dm(ix,iy,iz) * 0.5_num
            delta_ke(ix ,iyp,iz ) = delta_ke(ix ,iyp,iz ) + dk
            delta_ke(ixp,iyp,iz ) = delta_ke(ixp,iyp,iz ) + dk
            delta_ke(ix ,iyp,izp) = delta_ke(ix ,iyp,izp) + dk
            delta_ke(ixp,iyp,izp) = delta_ke(ixp,iyp,izp) + dk
          END DO
        END DO
      END DO
    END IF

    flux(0:nx,-1:ny,0:nz) = flux(0:nx,-1:ny,0:nz) * dm(0:nx,-1:ny,0:nz)

  END SUBROUTINE y_momy_flux



  SUBROUTINE y_momz_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dyci, dycu, dybu, phi, ss, Da, Di
    REAL(num) :: area, rhou, dmu, m, mp, ai, aip, dk
    INTEGER :: iyp2

    DO iz = 0, nz
      DO iy = -1, ny
        iym  = iy - 1
        iyp  = iy + 1
        iyp2 = iy + 2
        DO ix = 0, nx
          area = dxc(ix) * dzc(iz)

          v_advect = vy1(ix,iy,iz)

          fm  = vz(ix,iym ,iz)
          fi  = vz(ix,iy  ,iz)
          fp  = vz(ix,iyp ,iz)
          fp2 = vz(ix,iyp2,iz)

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
          dyci = dyb1(ix,iyp,iz)
          dycu = dyb1(ix,iy ,iz) * vad_p + dyb1(ix,iyp2,iz) * vad_m
          dybu = dyc1(ix,iy ,iz) * vad_p + dyc1(ix,iyp ,iz) * vad_m

          phi = ABS(v_advect) * dt / dybu

          Da =  (2.0_num - phi) * ABS(dfi) / dyci &
              + (1.0_num + phi) * ABS(dfu) / dycu
          Da = Da * sixth

          ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

          Di = sign_v * ss * MIN(ABS(Da) * dybu, ABS(dfi), ABS(dfu))

          rhou = rho_v(ix,iy,iz) * vad_p + rho_v(ix,iyp,iz) * vad_m
          dmu = ABS(dm(ix,iy,iz)) / area / dybu / rhou

          flux(ix,iy,iz) = fu + Di * (1.0_num - dmu)
        END DO
      END DO
    END DO

    IF (rke) THEN
      DO iz = 0, nz
        izp = iz + 1
        DO iy = 0, ny - 1
          iym = iy - 1
          iyp = iy + 1
          DO ix = 0, nx
            ixp = ix + 1

            m  = rho_v1(ix,iy ,iz) * cv2(ix,iy ,iz)
            mp = rho_v1(ix,iyp,iz) * cv2(ix,iyp,iz)

            ai =  (vz(ix,iy ,iz) - flux(ix,iym,iz)) * dm(ix,iym,iz) / m &
                - (vz(ix,iy ,iz) - flux(ix,iy ,iz)) * dm(ix,iy ,iz) / m

            aip = (vz(ix,iyp,iz) - flux(ix,iy ,iz)) * dm(ix,iy ,iz) / mp &
                - (vz(ix,iyp,iz) - flux(ix,iyp,iz)) * dm(ix,iyp,iz) / mp

            dk = (vz(ix,iyp,iz) - vz(ix,iy,iz)) * (flux(ix,iy,iz) &
                - 0.5_num * (vz(ix,iyp,iz) + vz(ix,iy,iz))) &
                - 0.5_num * ai  * (vz(ix,iy ,iz) - flux(ix,iy,iz)) &
                + 0.5_num * aip * (vz(ix,iyp,iz) - flux(ix,iy,iz))

            dk = dk * dm(ix,iy,iz) * 0.5_num
            delta_ke(ix ,iyp,iz ) = delta_ke(ix ,iyp,iz ) + dk
            delta_ke(ixp,iyp,iz ) = delta_ke(ixp,iyp,iz ) + dk
            delta_ke(ix ,iyp,izp) = delta_ke(ix ,iyp,izp) + dk
            delta_ke(ixp,iyp,izp) = delta_ke(ixp,iyp,izp) + dk
          END DO
        END DO
      END DO
    END IF

    flux(0:nx,-1:ny,0:nz) = flux(0:nx,-1:ny,0:nz) * dm(0:nx,-1:ny,0:nz)

  END SUBROUTINE y_momz_flux



  SUBROUTINE dm_y_bcs

    CALL MPI_SENDRECV(&
        dm(-1,1   ,-1), 1, cell_yface, proc_y_min, tag, &
        dm(-1,ny+1,-1), 1, cell_yface, proc_y_max, tag, &
        comm, status, errcode)

    IF (proc_y_max == MPI_PROC_NULL) &
        dm(0:nx+1,ny+1,0:nz+1) = dm(0:nx+1,ny,0:nz+1)

    CALL MPI_SENDRECV(&
        dm(-1,ny-1,-1), 1, cell_yface, proc_y_max, tag, &
        dm(-1,-1  ,-1), 1, cell_yface, proc_y_min, tag, &
        comm, status, errcode)

    IF (proc_y_min == MPI_PROC_NULL) &
        dm(0:nx+1,-1,0:nz+1) = dm(0:nx+1,0,0:nz+1)

  END SUBROUTINE dm_y_bcs

END MODULE yremap
