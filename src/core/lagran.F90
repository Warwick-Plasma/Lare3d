MODULE lagran

  !-----------------------------------------------------------------
  ! This subroutine performs the Lagrangian step
  !-----------------------------------------------------------------

  USE shared_data
  USE boundary
  USE neutral
  USE diagnostics
  USE conduct

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lagrangian_step, eta_calc

  ! only used inside lagran.f90
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: qxy, qxz, qyz, pressure
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: qxx, qyy, qzz, visc_heat
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: flux_x, flux_y, flux_z, curlb

CONTAINS

  SUBROUTINE lagrangian_step

    INTEGER :: substeps, subcycle
    REAL(num) :: actual_dt, dt_sub

    ALLOCATE (bx1(0:nx+1, 0:ny+1, 0:nz+1), by1(0:nx+1, 0:ny+1, 0:nz+1), &
        bz1(0:nx+1, 0:ny+1, 0:nz+1), qxy(0:nx+1, 0:ny+1, 0:nz+1), &
        qxz(0:nx+1, 0:ny+1, 0:nz+1), qyz(0:nx+1, 0:ny+1, 0:nz+1), &
        visc_heat(0:nx+1, 0:ny+1, 0:nz+1), &
        pressure(-1:nx+2, -1:ny+2, -1:nz+2), qxx(0:nx+1, 0:ny+1, 0:nz+1), &
        qyy(0:nx+1, 0:ny+1, 0:nz+1), qzz(0:nx+1, 0:ny+1, 0:nz+1), &
        flux_x(0:nx, 0:ny, 0:nz), flux_y(0:nx, 0:ny, 0:nz), &
        flux_z(0:nx, 0:ny, 0:nz), curlb(0:nx, 0:ny, 0:nz))

    IF (resistive_mhd) THEN
      ! if subcycling isn't wanted set dt = dtr in set_dt, don't just
      ! set substeps to 1.
      dt_sub = dtr

      substeps = INT(dt / dt_sub) + 1

      IF (substeps > peak_substeps) peak_substeps = substeps
      actual_dt = dt
      dt = dt / REAL(substeps, num)

      DO subcycle = 1, substeps
        CALL eta_calc
        IF (eos_number /= EOS_IDEAL) CALL neutral_fraction
        IF (cowling_resistivity) CALL perpendicular_resistivity
        ! IF (hall_mhd) CALL hall_effects
        IF (resistive_mhd) CALL resistive_effects
      END DO

      dt = actual_dt
    END IF

    IF (conduction) CALL conduct_heat 

    DO iz = 0, nz+1
      izm = iz - 1
      DO iy = 0, ny+1
        iym = iy - 1
        DO ix = 0, nx+1
          ixm = ix - 1
          bx1(ix, iy, iz) = (bx(ix, iy, iz) + bx(ixm, iy, iz)) * 0.5_num
          by1(ix, iy, iz) = (by(ix, iy, iz) + by(ix, iym, iz)) * 0.5_num
          bz1(ix, iy, iz) = (bz(ix, iy, iz) + bz(ix, iy, izm)) * 0.5_num
        END DO
      END DO
    END DO

    CALL predictor_corrector_step

    DEALLOCATE (bx1, by1, bz1, qxy, qxz, qyz, visc_heat, pressure, qxx, qyy, &
        qzz, flux_x, flux_y, flux_z, curlb)

    CALL energy_bcs
    CALL density_bcs
    CALL velocity_bcs

  END SUBROUTINE lagrangian_step



  SUBROUTINE predictor_corrector_step

    REAL(num) :: p, pxp, pyp, pxpyp
    REAL(num) :: pzp, pxpzp, pypzp, pxpypzp
    REAL(num) :: e1, rho_v
    REAL(num) :: fx, fy, fz
    REAL(num) :: vxb, vxbm, vyb, vybm, vzb, vzbm
    REAL(num) :: bxv, byv, bzv, jx, jy, jz
    REAL(num) :: cvx, cvxp, cvy, cvyp, cvz, cvzp
    REAL(num) :: dv

    dt2 = dt * 0.5_num
    CALL viscosity_and_b_update

    bx1 = bx1 * cv1(0:nx+1, 0:ny+1, 0:nz+1)
    by1 = by1 * cv1(0:nx+1, 0:ny+1, 0:nz+1)
    bz1 = bz1 * cv1(0:nx+1, 0:ny+1, 0:nz+1)

    DO iz = 0, nz+1
      DO iy = 0, ny+1
        DO ix = 0, nx+1
          dv = cv1(ix, iy, iz) / cv(ix, iy, iz) - 1.0_num
          ! predictor energy
#ifdef Q_MONO
            ! add shock viscosity
            pressure(ix, iy, iz) = pressure(ix, iy, iz) + p_visc(ix, iy, iz)
#endif
          e1 = energy(ix, iy, iz) - pressure(ix, iy, iz) * dv / rho(ix, iy, iz)
          e1 = e1  + visc_heat(ix, iy, iz) * dt2 / rho(ix, iy, iz) 

          ! now define the predictor step pressures 
          pressure(ix,iy,iz) = (e1 - (1.0_num - xi_n(ix,iy,iz)) * ionise_pot) &
                  * (gamma - 1.0_num) * rho(ix,iy,iz) * cv(ix,iy,iz) / cv1(ix,iy,iz)
#ifdef Q_MONO
          ! add shock viscosity
          pressure(ix, iy, iz) = pressure(ix, iy, iz) + p_visc(ix, iy, iz)
#endif
        END DO
      END DO
    END DO

    DO iz = 0, nz
      izp = iz + 1
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = 0, nx
          ixp = ix + 1

          p       = pressure(ix , iy , iz )
          pxp     = pressure(ixp, iy , iz )
          pyp     = pressure(ix , iyp, iz )
          pxpyp   = pressure(ixp, iyp, iz )
          pzp     = pressure(ix , iy , izp)
          pxpzp   = pressure(ixp, iy , izp)
          pypzp   = pressure(ix , iyp, izp)
          pxpypzp = pressure(ixp, iyp, izp)

          w1 = (p + pyp + pzp + pypzp) * 0.25_num
          w2 = (pxp + pxpyp + pxpzp + pxpypzp) * 0.25_num
          fx = -(w2 - w1) / dxc(ix)
          w1 = (p + pxp + pzp + pxpzp) * 0.25_num
          w2 = (pyp + pxpyp + pypzp + pxpypzp) * 0.25_num
          fy = -(w2 - w1) / dyc(iy)
          w1 = (p + pxp + pyp + pxpyp) * 0.25_num
          w2 = (pzp + pxpzp + pypzp + pxpypzp) * 0.25_num
          fz = -(w2 - w1) / dzc(iz)

          ! add diagonal components
          w1 = (qxx(ix, iy, iz) + qxx(ix, iyp, iz) &
              + qxx(ix, iy, izp) + qxx(ix, iyp, izp)) * 0.25_num
          w2 = (qxx(ixp, iy, iz) + qxx(ixp, iyp, iz) &
              + qxx(ixp, iy, izp) + qxx(ixp, iyp, izp)) * 0.25_num
          fx = fx + (w2 - w1) / dxc(ix)
  
          w1 = (qyy(ix, iy, iz) + qyy(ixp, iy, iz) &
              + qyy(ix, iy, izp) + qyy(ixp, iy, izp)) * 0.25_num
          w2 = (qyy(ix, iyp, iz) + qyy(ixp, iyp, iz) &
              + qyy(ix, iyp, izp) + qyy(ixp, iyp, izp)) * 0.25_num
          fy = fy + (w2 - w1) / dyc(iy)
  
          w1 = (qzz(ix, iy, iz) + qzz(ixp, iy, iz) &
              + qzz(ix, iyp, iz) + qzz(ixp, iyp, iz)) * 0.25_num
          w2 = (qzz(ix, iy, izp) + qzz(ixp, iy, izp) &
              + qzz(ix, iyp, izp) + qzz(ixp, iyp, izp)) * 0.25_num
          fz = fz + (w2 - w1) / dzc(iz)
  
          ! add shear viscosity forces
          ! fx
          w1 = (qxy(ix, iy, iz) + qxy(ixp, iy, iz) &
              + qxy(ix, iy, izp) + qxy(ixp, iy, izp)) * 0.25_num
          w2 = (qxy(ix, iyp, iz) + qxy(ixp, iyp, iz) &
              + qxy(ix, iyp, izp) + qxy(ixp, iyp, izp)) * 0.25_num
          fx = fx + (w2 - w1) / dyc(iy)
  
          w1 = (qxz(ix, iy, iz) + qxz(ixp, iy, iz) &
              + qxz(ix, iyp, iz) + qxz(ixp, iyp, iz)) * 0.25_num
          w2 = (qxz(ix, iy, izp) + qxz(ixp, iy, izp) &
              + qxz(ix, iyp, izp) + qxz(ixp, iyp, izp)) * 0.25_num
          fx = fx + (w2 - w1) / dzc(iz)
  
          ! fy
          w1 = (qxy(ix, iy, iz) + qxy(ix, iyp, iz) &
              + qxy(ix, iy, izp) + qxy(ix, iyp, izp)) * 0.25_num
          w2 = (qxy(ixp, iy, iz) + qxy(ixp, iyp, iz) &
              + qxy(ixp, iy, izp) + qxy(ixp, iyp, izp)) * 0.25_num
          fy = fy + (w2 - w1) / dxc(ix)
  
          w1 = (qyz(ix, iy, iz) + qyz(ixp, iy, iz) &
              + qyz(ix, iyp, iz) + qyz(ixp, iyp, iz)) * 0.25_num
          w2 = (qyz(ix, iy, izp) + qyz(ixp, iy, izp) &
              + qyz(ix, iyp, izp) + qyz(ixp, iyp, izp)) * 0.25_num
          fy = fy + (w2 - w1) / dzc(iz)
 
          ! fz
          w1 = (qxz(ix, iy, iz) + qxz(ix, iyp, iz) &
              + qxz(ix, iy, izp) + qxz(ix, iyp, izp)) * 0.25_num
          w2 = (qxz(ixp, iy, iz) + qxz(ixp, iyp, iz) &
              + qxz(ixp, iy, izp) + qxz(ixp, iyp, izp)) * 0.25_num
          fz = fz + (w2 - w1) / dxc(ix)
 
          w1 = (qyz(ix, iy, iz) + qyz(ixp, iy, iz) &
              + qyz(ix, iy, izp) + qyz(ixp, iy, izp)) * 0.25_num
          w2 = (qyz(ix, iyp, iz) + qyz(ixp, iyp, iz) &
              + qyz(ix, iyp, izp) + qyz(ixp, iyp, izp)) * 0.25_num
          fz = fz + (w2 - w1) / dyc(iy)

          cvx = cv1(ix, iy, iz) + cv1(ix, iyp, iz) &
              + cv1(ix, iy, izp) + cv1(ix, iyp, izp)
          cvxp = cv1(ixp, iy, iz) + cv1(ixp, iyp, iz) &
              + cv1(ixp, iy, izp) + cv1(ixp, iyp, izp)
          cvy = cv1(ix, iy, iz) + cv1(ixp, iy, iz) &
              + cv1(ix, iy, izp) + cv1(ixp, iy, izp)
          cvyp = cv1(ix, iyp, iz) + cv1(ixp, iyp, iz) &
              + cv1(ix, iyp, izp) + cv1(ixp, iyp, izp)
          cvz = cv1(ix, iy, iz) + cv1(ixp, iy, iz) &
              + cv1(ix, iyp, iz) + cv1(ixp, iyp, iz)
          cvzp = cv1(ix, iy, izp) + cv1(ixp, iy, izp) &
              + cv1(ix, iyp, izp) + cv1(ixp, iyp, izp)

          w1 = (bz1(ix, iy, iz) + bz1(ixp, iy, iz) &
              + bz1(ix, iy, izp) + bz1(ixp, iy, izp)) / cvy
          w2 = (bz1(ix, iyp, iz) + bz1(ixp, iyp, iz) &
              + bz1(ix, iyp, izp) + bz1(ixp, iyp, izp)) / cvyp
          jx = (w2 - w1) / dyc(iy)
          w1 = (by1(ix, iy, iz) + by1(ixp, iy, iz) &
              + by1(ix, iyp, iz) + by1(ixp, iyp, iz)) / cvz
          w2 = (by1(ix, iy, izp) + by1(ixp, iy, izp) &
              + by1(ix, iyp, izp) + by1(ixp, iyp, izp)) / cvzp
          jx = jx - (w2 - w1) / dzc(iz)

          w1 = (bz1(ix, iy, iz) + bz1(ix, iyp, iz) &
              + bz1(ix, iy, izp) + bz1(ix, iyp, izp)) / cvx
          w2 = (bz1(ixp, iy, iz) + bz1(ixp, iyp, iz) &
              + bz1(ixp, iy, izp) + bz1(ixp, iyp, izp)) / cvxp
          jy = -(w2 - w1) / dxc(ix)
          w1 = (bx1(ix, iy, iz) + bx1(ixp, iy, iz) &
              + bx1(ix, iyp, iz) + bx1(ixp, iyp, iz)) / cvz
          w2 = (bx1(ix, iy, izp) + bx1(ixp, iy, izp) &
              + bx1(ix, iyp, izp) + bx1(ixp, iyp, izp)) / cvzp
          jy = jy + (w2 - w1) / dzc(iz)

          w1 = (by1(ix, iy, iz) + by1(ix, iyp, iz) &
              + by1(ix, iy, izp) + by1(ix, iyp, izp)) / cvx
          w2 = (by1(ixp, iy, iz) + by1(ixp, iyp, iz) &
              + by1(ixp, iy, izp) + by1(ixp, iyp, izp)) / cvxp
          jz = (w2 - w1) / dxc(ix)
          w1 = (bx1(ix, iy, iz) + bx1(ixp, iy, iz) &
              + bx1(ix, iy, izp) + bx1(ixp, iy, izp)) / cvy
          w2 = (bx1(ix, iyp, iz) + bx1(ixp, iyp, iz) &
              + bx1(ix, iyp, izp) + bx1(ixp, iyp, izp)) / cvyp
          jz = jz - (w2 - w1) / dyc(iy)

          bxv = (bx1(ix , iy , iz ) + bx1(ixp, iy , iz ) &
              + bx1(ix , iy , izp) + bx1(ixp, iy , izp) &
              + bx1(ix , iyp, iz ) + bx1(ixp, iyp, iz ) &
              + bx1(ix , iyp, izp) + bx1(ixp, iyp, izp)) &
              / (cvx + cvxp)

          byv = (by1(ix , iy , iz ) + by1(ixp, iy , iz ) &
              + by1(ix , iy , izp) + by1(ixp, iy , izp) &
              + by1(ix , iyp, iz ) + by1(ixp, iyp, iz ) &
              + by1(ix , iyp, izp) + by1(ixp, iyp, izp)) &
              / (cvx + cvxp)

          bzv = (bz1(ix , iy , iz ) + bz1(ixp, iy , iz ) &
              + bz1(ix , iy , izp) + bz1(ixp, iy , izp) &
              + bz1(ix , iyp, iz ) + bz1(ixp, iyp, iz ) &
              + bz1(ix , iyp, izp) + bz1(ixp, iyp, izp)) &
              / (cvx + cvxp)

          fx = fx + (jy * bzv - jz * byv)
          fy = fy - (jx * bzv - jz * bxv)
          fz = fz + (jx * byv - jy * bxv)

          rho_v = rho(ix , iy , iz ) * cv(ix , iy , iz ) &
              + rho(ixp, iy , iz ) * cv(ixp, iy , iz ) &
              + rho(ix , iyp, iz ) * cv(ix , iyp, iz ) &
              + rho(ixp, iyp, iz ) * cv(ixp, iyp, iz ) &
              + rho(ix , iy , izp) * cv(ix , iy , izp) &
              + rho(ixp, iy , izp) * cv(ixp, iy , izp) &
              + rho(ix , iyp, izp) * cv(ix , iyp, izp) &
              + rho(ixp, iyp, izp) * cv(ixp, iyp, izp)

          rho_v = rho_v / (cv(ix, iy , iz ) + cv(ixp, iy , iz ) &
              + cv(ix, iyp, iz ) + cv(ixp, iyp, iz ) &
              + cv(ix, iy , izp) + cv(ixp, iy , izp) &
              + cv(ix, iyp, izp) + cv(ixp, iyp, izp))

          fz = fz - (rho_v * grav(iz))

          ! find half step velocity needed for remap
          vx1(ix, iy, iz) = vx(ix, iy, iz) + dt2 * fx / rho_v
          vy1(ix, iy, iz) = vy(ix, iy, iz) + dt2 * fy / rho_v
          vz1(ix, iy, iz) = vz(ix, iy, iz) + dt2 * fz / rho_v

          ! velocity at the end of the Lagrangian step
          vx(ix, iy, iz) = vx(ix, iy, iz) + dt * fx / rho_v
          vy(ix, iy, iz) = vy(ix, iy, iz) + dt * fy / rho_v
          vz(ix, iy, iz) = vz(ix, iy, iz) + dt * fz / rho_v

        END DO
      END DO
    END DO

    CALL remap_v_bcs

    CALL visc_heating

    ! finally correct density and energy to final values
    DO iz = 1, nz
      izm = iz - 1
      DO iy = 1, ny
        iym = iy - 1
        DO ix = 1, nx
          ixm = ix - 1

          vxb = (vx1(ix, iy, iz) + vx1(ix, iym, iz) &
              + vx1(ix, iy, izm) + vx1(ix, iym, izm)) * 0.25_num

          vxbm = (vx1(ixm, iy, iz) + vx1(ixm, iym, iz) &
              + vx1(ixm, iy, izm) + vx1(ixm, iym, izm)) * 0.25_num

          vyb = (vy1(ix, iy, iz) + vy1(ixm, iy, iz) &
              + vy1(ix, iy, izm) + vy1(ixm, iy, izm)) * 0.25_num

          vybm = (vy1(ix, iym, iz) + vy1(ixm, iym, iz) &
              + vy1(ix, iym, izm) + vy1(ixm, iym, izm)) * 0.25_num

          vzb = (vz1(ix, iy, iz) + vz1(ixm, iy, iz) &
              + vz1(ix, iym, iz) + vz1(ixm, iym, iz)) * 0.25_num

          vzbm = (vz1(ix, iy, izm) + vz1(ixm, iy, izm) &
              + vz1(ix, iym, izm) + vz1(ixm, iym, izm)) * 0.25_num

          dv = ((vxb - vxbm) / dxb(ix) + (vyb - vybm) / dyb(iy) &
              + (vzb - vzbm) / dzb(iz)) * dt
          ! it is possible that dv has changed sign since the predictor step
          ! in this case p_visc * dv ought to be removed from the heating
          ! if QMONO is set - not done for simplicity since this is a rare
          ! combination

          cv1(ix, iy, iz) = cv(ix, iy, iz) * (1.0_num + dv)

          ! energy at end of Lagrangian step

          energy(ix, iy, iz) = energy(ix, iy, iz) &
              - pressure(ix, iy, iz) * dv / rho(ix, iy, iz) &
              + dt * visc_heat(ix, iy, iz) / rho(ix, iy, iz)

          ! density at end of Lagrangian step
          rho(ix, iy, iz) = rho(ix, iy, iz) / (1.0_num + dv)

          total_visc_heating = total_visc_heating &
                + dt * visc_heat(ix, iy, iz) * cv(ix, iy, iz) 
                
#ifdef Q_MONO
          total_visc_heating = total_visc_heating &
              - p_visc(ix, iy, iz) * dv * cv(ix, iy, iz)
#endif
        END DO
      END DO
    END DO

  END SUBROUTINE predictor_corrector_step



  SUBROUTINE viscosity_and_b_update
    ! wilkins viscosity and B field update
    REAL(num) :: vxb, vxbm, vyb, vybm, vzb, vzbm
    REAL(num) :: p, pxp, pxm, pyp, pym, pzp, pzm, fx, fy, fz, dv
    REAL(num) :: dvxdx, dvydy, dvzdz, dvxy, dvxz, dvyz, s, L, L2, cf
    REAL(num) :: sxx, syy, szz, sxy, sxz, syz
    REAL(num) :: dvxdy, dvxdz, dvydx, dvydz, dvzdx, dvzdy
    REAL(num) :: cs
    REAL(num) :: w2_1,w2_2,w2_3,w2_4
    REAL(num) :: flag1,flag2,flag3,flag4,sg0,dvg0

    DO iz = -1, nz+2
      DO iy = -1, ny+2
        DO ix = -1, nx+2
           pressure(ix,iy,iz) = (energy(ix,iy,iz) - (1.0_num - xi_n(ix,iy,iz)) * ionise_pot) &
                  * (gamma - 1.0_num) * rho(ix,iy,iz)
        END DO
      END DO
    END DO

    DO iz = 0, nz+1
      izm = iz - 1
      izp = iz + 1
      DO iy = 0, ny+1
        iym = iy - 1
        iyp = iy + 1
        DO ix = 0, nx+1
          ixm = ix - 1
          ixp = ix + 1

          ! vx at Sx(i, j, k)
          vxb = (vx(ix, iy, iz) + vx(ix, iym, iz) &
              + vx(ix, iy, izm) + vx(ix, iym, izm)) * 0.25_num

          ! vx at Sx(i-1, j, k)
          vxbm = (vx(ixm, iy, iz) + vx(ixm, iym, iz) &
              + vx(ixm, iy, izm) + vx(ixm, iym, izm)) * 0.25_num

          ! vy at Sy(i, j, k)
          vyb = (vy(ix, iy, iz) + vy(ixm, iy, iz) &
              + vy(ix, iy, izm) + vy(ixm, iy, izm)) * 0.25_num

          ! vy at Sy(i, j-1, k)
          vybm = (vy(ix, iym, iz) + vy(ixm, iym, iz) &
              + vy(ix, iym, izm) + vy(ixm, iym, izm)) * 0.25_num

          ! vz at Sz(i, j, k)
          vzb = (vz(ix, iy, iz) + vz(ixm, iy, iz) &
              + vz(ix, iym, iz) + vz(ixm, iym, iz)) * 0.25_num

          ! vz at Sz(i, j, k-1)
          vzbm = (vz(ix, iy, izm) + vz(ixm, iy, izm) &
              + vz(ix, iym, izm) + vz(ixm, iym, izm)) * 0.25_num

          dvxdx = (vxb - vxbm) / dxb(ix)
          dvydy = (vyb - vybm) / dyb(iy)
          dvzdz = (vzb - vzbm) / dzb(iz)

          dv = (dvxdx + dvydy + dvzdz) * dt2
          cv1(ix, iy, iz) = cv(ix, iy, iz) * (1.0_num + dv)

          ! vx at Sy(i, j, k)
          vxb = (vx(ix, iy, iz) + vx(ixm, iy, iz) &
              + vx(ix, iy, izm) + vx(ixm, iy, izm)) * 0.25_num

          ! vx at Sy(i, j-1, k)
          vxbm = (vx(ix, iym, iz) + vx(ixm, iym, iz) &
              + vx(ix, iym, izm) + vx(ixm, iym, izm)) * 0.25_num

          ! vy at Sx(i, j, k)
          vyb = (vy(ix, iy, iz) + vy(ix, iym, iz) &
              + vy(ix, iy, izm) + vy(ix, iym, izm)) * 0.25_num

          ! vy at Sx(i-1, j, k)
          vybm = (vy(ixm, iy, iz) + vy(ixm, iym, iz) &
              + vy(ixm, iy, izm) + vy(ixm, iym, izm)) * 0.25_num

          dvxdy = (vxb - vxbm) / dyb(iy)
          dvydx = (vyb - vybm) / dxb(ix)
          dvxy = dvxdy + dvydx
          sxy = dvxy * 0.5_num
          sxx = 2.0_num * dvxdx * third - (dvydy + dvzdz) * third
          syy = 2.0_num * dvydy * third - (dvxdx + dvzdz) * third
          szz = 2.0_num * dvzdz * third - (dvxdx + dvydy) * third

          ! vx at Sz(i, j, k)
          vxb = (vx(ix, iy, iz) + vx(ixm, iy, iz) &
              + vx(ix, iym, iz) + vx(ixm, iym, iz)) * 0.25_num

          ! vx at Sz(i, j, k-1)
          vxbm = (vx(ix, iy, izm) + vx(ixm, iy, izm) &
              + vx(ix, iym, izm) + vx(ixm, iym, izm)) * 0.25_num

          ! vz at Sx(i, j, k)
          vzb = (vz(ix, iy, iz) + vz(ix, iym, iz) &
              + vz(ix, iy, izm) + vz(ix, iym, izm)) * 0.25_num

          ! vz at Sx(i-1, j, k)
          vzbm = (vz(ixm, iy, iz) + vz(ixm, iym, iz) &
              + vz(ixm, iy, izm) + vz(ixm, iym, izm)) * 0.25_num

          dvxdz = (vxb - vxbm) / dzb(iz)
          dvzdx = (vzb - vzbm) / dxb(ix)
          dvxz = dvxdz + dvzdx
          sxz = dvxz * 0.5_num

          ! vy at Sz(i, j, k)
          vyb = (vy(ix, iy, iz) + vy(ixm, iy, iz) &
              + vy(ix, iym, iz) + vy(ixm, iym, iz)) * 0.25_num

          ! vy at Sz(i, j, k-1)
          vybm = (vy(ix, iy, izm) + vy(ixm, iy, izm) &
              + vy(ix, iym, izm) + vy(ixm, iym, izm)) * 0.25_num

          ! vz at Sy(i, j, k)
          vzb = (vz(ix, iy, iz) + vz(ixm, iy, iz) &
              + vz(ix, iy, izm) + vz(ixm, iy, izm)) * 0.25_num

          ! vz at Sy(i, j-1, k)
          vzbm = (vz(ix, iym, iz) + vz(ixm, iym, iz) &
              + vz(ix, iym, izm) + vz(ixm, iym, izm)) * 0.25_num

          dvydz = (vyb - vybm) / dzb(iz)
          dvzdy = (vzb - vzbm) / dyb(iy)
          dvyz = dvydz + dvzdy
          syz = dvyz * 0.5_num

          p = energy(ix, iy, iz) * (gamma - 1.0_num) * rho(ix, iy, iz)
          pxp = energy(ixp, iy, iz) * (gamma - 1.0_num) * rho(ixp, iy, iz)
          pxm = energy(ixm, iy, iz) * (gamma - 1.0_num) * rho(ixm, iy, iz)
          pyp = energy(ix, iyp, iz) * (gamma - 1.0_num) * rho(ix, iyp, iz)
          pym = energy(ix, iym, iz) * (gamma - 1.0_num) * rho(ix, iym, iz)
          pzp = energy(ix, iy, izp) * (gamma - 1.0_num) * rho(ix, iy, izp)
          pzm = energy(ix, iy, izm) * (gamma - 1.0_num) * rho(ix, iy, izm)

          ! should be half this but this cancels later
          fx = -(pxp - pxm) / dxb(ix)
          fy = -(pyp - pym) / dyb(iy)
          fz = -(pzp - pzm) / dzb(iz)

          w1 = fx**2 + fy**2 + fz**2
          s = (dvxdx * fx**2 + dvydy * fy**2 + dvzdz * fz**2 + dvxy * fx * fy &
              + dvxz * fx * fz + dvyz * fy * fz) / MAX(w1, none_zero)

!	     These flags are used to replace the rather clearer code in 
!	     **SECTION 1**. They are used instead to facilitate vector
!	     optimization
             flag1=MAX(SIGN(1.0_num,dyb(iy)*ABS(fx)-dxb(ix)*ABS(fy)),0.0_num)
             flag2=MAX(SIGN(1.0_num,dzb(iz)*ABS(fx)-dxb(ix)*ABS(fz)),0.0_num)
             flag3=MAX(SIGN(1.0_num,dzb(iz)*ABS(fy)-dyb(iy)*ABS(fz)),0.0_num)
             flag4=MAX(SIGN(1.0_num,w1-1.e-6_num),0.0_num)

             w2_1=dxb(ix)**2*w1 / MAX(fx**2, 1.e-20_num)
             w2_2=dzb(iz)**2*w1 / MAX(fz**2, 1.e-20_num)
             w2_3=dyb(iy)**2*w1 / MAX(fy**2, 1.e-20_num)
             w2_4=dzb(iz)**2*w1 / MAX(fz**2, 1.e-20_num)

             w2=w2_1*flag1*flag2 + w2_2*flag1*(1.0_num-flag2)&
                  +w2_3*(1.0_num-flag1)*flag3 + w2_4*(1.0_num-flag1)*(1.0_num-flag3)

             w2=w2*flag4 + MIN(dxb(ix), dyb(iy), dzb(iz))**2 * (1.0_num-flag4)
!!$		!BEGIN **SECTION 1**
!!$		!*******************
!          IF (dxb(ix) * ABS(fy) < dyb(iy) * ABS(fx)) THEN
!            IF (dxb(ix) * ABS(fz) < dzb(iz) * ABS(fx)) THEN
!              w2 = dxb(ix)**2 * w1 / (fx**2 + 1.e-20_num)
!            ELSE
!              w2 = dzb(iz)**2 * w1 / (fz**2 + 1.e-20_num)
!            END IF
!          ELSE
!            IF (dyb(iy) * ABS(fz) < dzb(iz) * ABS(fy)) THEN
!              w2 = dyb(iy)**2 * w1 / (fy**2 + 1.e-20_num)
!            ELSE
!              w2 = dzb(iz)**2 * w1 / (fz**2 + 1.e-20_num)
!            END IF
!          END IF
!          IF (w1 < 1.e-6_num) w2 = MIN(dxb(ix), dyb(iy), dzb(iz))**2
!$		!*****************

          L = SQRT(w2)

          L2 = L
			!This code is equivalent to IF (s > 0 .OR. dv > 0) L=0.0
  			  sg0 = MAX(SIGN(1.0_num,s),0.0_num)
   			  dvg0 = MAX(SIGN(1.0_num,dv),0.0_num)
          L = L * (1.0_num - dvg0) * (1.0_num -sg0) 

          w1 = (bx1(ix, iy, iz)**2 + by1(ix, iy, iz)**2 + bz1(ix, iy, iz)**2) &
                / rho(ix, iy, iz)

          cs = SQRT(gamma*(gamma - 1.0_num) * energy(ix,iy,iz))    
          cf = SQRT(cs**2 + w1)

          p_visc(ix, iy, iz) = visc1 * ABS(s) * L*cf * rho(ix, iy, iz) &
              + visc2 * (s * L)**2 * rho(ix, iy, iz)

          qxx(ix, iy, iz) = 0.0_num
          qxy(ix, iy, iz) = 0.0_num
          qxz(ix, iy, iz) = 0.0_num
          qyy(ix, iy, iz) = 0.0_num
          qyz(ix, iy, iz) = 0.0_num
          qzz(ix, iy, iz) = 0.0_num

#ifndef Q_MONO
            qxy(ix,iy,iz) = sxy * (L2 * rho(ix,iy,iz)  &
                 * (visc1 * cf + L2 * visc2 * ABS(sxy)))
            qxz(ix,iy,iz) = sxz * (L2 * rho(ix,iy,iz)  &
                 * (visc1 * cf + L2 * visc2 * ABS(sxz)))
            qyz(ix,iy,iz) = syz * (L2 * rho(ix,iy,iz)  &
                 * (visc1 * cf + L2 * visc2 * ABS(syz)))
            qxx(ix,iy,iz) = sxx * (L2 * rho(ix,iy,iz)  &
                 * (visc1 * cf + L2 * visc2 * ABS(sxx)))
            qyy(ix,iy,iz) = syy * (L2 * rho(ix,iy,iz)  &
                 * (visc1 * cf + L2 * visc2 * ABS(syy)))
            qzz(ix,iy,iz) = szz * (L2 * rho(ix,iy,iz)  &
                 * (visc1 * cf + L2 * visc2 * ABS(szz)))
#endif
          qxy(ix,iy,iz) = qxy(ix,iy,iz) + 2.0_num * sxy  * rho(ix,iy,iz) * visc3 
          qxz(ix,iy,iz) = qxz(ix,iy,iz) + 2.0_num * sxz  * rho(ix,iy,iz) * visc3 
          qyz(ix,iy,iz) = qyz(ix,iy,iz) + 2.0_num * syz  * rho(ix,iy,iz) * visc3 
          qxx(ix,iy,iz) = qxx(ix,iy,iz) + 2.0_num * sxx  * rho(ix,iy,iz) * visc3 
          qyy(ix,iy,iz) = qyy(ix,iy,iz) + 2.0_num * syy  * rho(ix,iy,iz) * visc3
          qzz(ix,iy,iz) = qzz(ix,iy,iz) + 2.0_num * SZZ  * rho(ix,iy,iz) * visc3 

          visc_heat(ix,iy,iz) = qxy(ix,iy,iz)*dvxy + qxz(ix,iy,iz)*dvxz &
               + qyz(ix,iy,iz)*dvyz + qxx(ix,iy,iz)*dvxdx   &
               + qyy(ix,iy,iz)*dvydy + qzz(ix,iy,iz)*dvzdz 

          w4 = bx1(ix, iy, iz) * dvxdx &
              + by1(ix, iy, iz) * dvxdy + bz1(ix, iy, iz) * dvxdz
          bx1(ix, iy, iz) = (bx1(ix, iy, iz) + w4 * dt2) / (1.0_num + dv)

          w4 = bx1(ix, iy, iz) * dvydx &
              + by1(ix, iy, iz) * dvydy + bz1(ix, iy, iz) * dvydz
          by1(ix, iy, iz) = (by1(ix, iy, iz) + w4 * dt2) / (1.0_num + dv)

          w4 = bx1(ix, iy, iz) * dvzdx &
              + by1(ix, iy, iz) * dvzdy + bz1(ix, iy, iz) * dvzdz
          bz1(ix, iy, iz) = (bz1(ix, iy, iz) + w4 * dt2) / (1.0_num + dv)

        END DO
      END DO
    END DO

  END SUBROUTINE viscosity_and_b_update



  SUBROUTINE visc_heating

    REAL(num) :: vxb, vxbm, vyb, vybm, vzb, vzbm
    REAL(num) :: dvxdx, dvydy, dvzdz, dvxy, dvxz, dvyz

    RETURN

    DO iz = 0, nz+1
      izm = iz - 1
      izp = iz + 1
      DO iy = 0, ny+1
        iym = iy - 1
        iyp = iy + 1
        DO ix = 0, nx+1
          ixm = ix - 1
          ixp = ix + 1

          ! vx at Sx(i, j, k)
          vxb = (vx1(ix, iy, iz) + vx1(ix, iym, iz) &
              + vx1(ix, iy, izm) + vx1(ix, iym, izm)) * 0.25_num

          ! vx at Sx(i-1, j, k)
          vxbm = (vx1(ixm, iy, iz) + vx1(ixm, iym, iz) &
              + vx1(ixm, iy, izm) + vx1(ixm, iym, izm)) * 0.25_num

          ! vy at Sy(i, j, k)
          vyb = (vy1(ix, iy, iz) + vy1(ixm, iy, iz) &
              + vy1(ix, iy, izm) + vy1(ixm, iy, izm)) * 0.25_num

          ! vy at Sy(i, j-1, k)
          vybm = (vy1(ix, iym, iz) + vy1(ixm, iym, iz) &
              + vy1(ix, iym, izm) + vy1(ixm, iym, izm)) * 0.25_num

          ! vz at Sz(i, j, k)
          vzb = (vz1(ix, iy, iz) + vz1(ixm, iy, iz) &
              + vz1(ix, iym, iz) + vz1(ixm, iym, iz)) * 0.25_num

          ! vz at Sz(i, j, k-1)
          vzbm = (vz1(ix, iy, izm) + vz1(ixm, iy, izm) &
              + vz1(ix, iym, izm) + vz1(ixm, iym, izm)) * 0.25_num

          dvxdx = (vxb - vxbm) / dxb(ix)
          dvydy = (vyb - vybm) / dyb(iy)
          dvzdz = (vzb - vzbm) / dzb(iz)

          vxb = (vx1(ix, iy, iz) + vx1(ixm, iy, iz) &
              + vx1(ix, iy, izm) + vx1(ixm, iy, izm)) * 0.25_num
          vxbm = (vx1(ix, iym, iz) + vx1(ixm, iym, iz) &
              + vx1(ix, iym, izm) + vx1(ixm, iym, izm)) * 0.25_num
          vyb = (vy1(ix, iy, iz) + vy1(ix, iym, iz) &
              + vy1(ix, iy, izm) + vy1(ix, iym, izm)) * 0.25_num
          vybm = (vy1(ixm, iy, iz) + vy1(ixm, iym, iz) &
              + vy1(ixm, iy, izm) + vy1(ixm, iym, izm)) * 0.25_num
          dvxy = (vxb - vxbm) / dyb(iy) + (vyb - vybm) / dxb(ix)

          vxb = (vx1(ix, iy, iz) + vx1(ixm, iy, iz) &
              + vx1(ix, iym, iz) + vx1(ixm, iym, iz)) * 0.25_num
          vxbm = (vx1(ix, iy, izm) + vx1(ixm, iy, izm) &
              + vx1(ix, iym, izm) + vx1(ixm, iym, izm)) * 0.25_num
          vzb = (vz1(ix, iy, iz) + vz1(ix, iym, iz) &
              + vz1(ix, iy, izm) + vz1(ix, iym, izm)) * 0.25_num
          vzbm = (vz1(ixm, iy, iz) + vz1(ixm, iym, iz) &
              + vz1(ixm, iy, izm) + vz1(ixm, iym, izm)) * 0.25_num
          dvxz = (vxb - vxbm) / dzb(iz) + (vzb - vzbm) / dxb(ix)

          vyb = (vy1(ix, iy, iz) + vy1(ixm, iy, iz) &
              + vy1(ix, iym, iz) + vy1(ixm, iym, iz)) * 0.25_num
          vybm = (vy1(ix, iy, izm) + vy1(ixm, iy, izm) &
              + vy1(ix, iym, izm) + vy1(ixm, iym, izm)) * 0.25_num
          vzb = (vz1(ix, iy, iz) + vz1(ixm, iy, iz) &
              + vz1(ix, iy, izm) + vz1(ixm, iy, izm)) * 0.25_num
          vzbm = (vz1(ix, iym, iz) + vz1(ixm, iym, iz) &
              + vz1(ix, iym, izm) + vz1(ixm, iym, izm)) * 0.25_num
          dvyz = (vyb - vybm) / dzb(iz) + (vzb - vzbm) / dyb(iy)

          visc_heat(ix, iy, iz) = &
              qxy(ix, iy, iz) * dvxy + qxz(ix, iy, iz) * dvxz &
              + qyz(ix, iy, iz) * dvyz + qxx(ix, iy, iz) * dvxdx &
              + qyy(ix, iy, iz) * dvydy + qzz(ix, iy, iz) * dvzdz
        END DO
      END DO
    END DO

  END SUBROUTINE visc_heating



  SUBROUTINE eta_calc

    REAL(num) :: jx, jy, jz, jxxp, jyyp, jzzp
    REAL(num) :: modj


    IF (resistive_mhd) THEN
      DO iz = -1, nz+1
        izp = iz + 1
        DO iy = -1, ny+1
          iyp = iy + 1
          DO ix = -1, nx+1
            ixp = ix + 1
  
            ! jx at E3(i, j)
            jx = (bz(ix, iyp, iz) - bz(ix, iy, iz)) / dyc(iy) &
                - (by(ix, iy, izp) - by(ix, iy, iz)) / dzc(iz)
  
            ! jx at E3(i+1, j)
            jxxp = (bz(ixp, iyp, iz) - bz(ixp, iy, iz)) / dyc(iy) &
                - (by(ixp, iy, izp) - by(ixp, iy, iz)) / dzc(iz)
  
            ! jy at E2(i, j)
            jy = (bx(ix, iy, izp) - bx(ix, iy, iz)) / dzc(iz) &
                - (bz(ixp, iy, iz) - bz(ix, iy, iz)) / dxc(ix)
  
            ! jy at E2(i, j+1)
            jyyp = (bx(ix, iyp, izp) - bx(ix, iyp, iz)) / dzc(iz) &
                - (bz(ixp, iyp, iz) - bz(ix, iyp, iz)) / dxc(ix)
  
            ! jz at E1(i, j)
            jz = (by(ixp, iy, iz) - by(ix, iy, iz)) / dxc(ix) &
                - (bx(ix, iyp, iz) - bx(ix, iy, iz)) / dyc(iy)
  
            ! jz at E1(i, j)
            jzzp = (by(ixp, iy, izp) - by(ix, iy, izp)) / dxc(ix) &
                - (bx(ix, iyp, izp) - bx(ix, iy, izp)) / dyc(iy)
  
            ! current at V
            w4 = (jx + jxxp) * 0.5_num
            w5 = (jy + jyyp) * 0.5_num
            w6 = (jz + jzzp) * 0.5_num
            modj = SQRT(w4**2 + w5**2 + w6**2)
  
            IF (modj >= j_max) THEN
              eta(ix, iy, iz) = eta_background + eta0
            ELSE
              eta(ix, iy, iz) = eta_background 
            END IF
          END DO
        END DO
      END DO
    ELSE
        eta = 0.0_num   
    END IF 

  END SUBROUTINE eta_calc



  ! Calculate the effect of resistivity on the magnetic field and Ohmic heating
  ! Use the subroutine rkstep
   SUBROUTINE resistive_effects

     REAL(num) :: half_dt, dt6, area
     REAL(num) :: jx1, jx2, jy1, jy2, jz1, jz2
#ifdef Q_FOURTHORDER
     REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: k1x, k2x, k3x, k4x
     REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: k1y, k2y, k3y, k4y
     REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: k1z, k2z, k3z, k4z
     REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: c1, c2, c3, c4

     ALLOCATE(k1x(0:nx, 0:ny, 0:nz), k2x(0:nx, 0:ny, 0:nz))
     ALLOCATE(k3x(0:nx, 0:ny, 0:nz), k4x(0:nx, 0:ny, 0:nz))
     ALLOCATE(k1y(0:nx, 0:ny, 0:nz), k2y(0:nx, 0:ny, 0:nz))
     ALLOCATE(k3y(0:nx, 0:ny, 0:nz), k4y(0:nx, 0:ny, 0:nz))
     ALLOCATE(k1z(0:nx, 0:ny, 0:nz), k2z(0:nx, 0:ny, 0:nz))
     ALLOCATE(k3z(0:nx, 0:ny, 0:nz), k4z(0:nx, 0:ny, 0:nz))
     ALLOCATE(c1(0:nx, 0:ny, 0:nz), c2(0:nx, 0:ny, 0:nz))
     ALLOCATE(c3(0:nx, 0:ny, 0:nz), c4(0:nx, 0:ny, 0:nz))
#endif

     bx1 = bx(0:nx+1, 0:ny+1, 0:nz+1)
     by1 = by(0:nx+1, 0:ny+1, 0:nz+1)
     bz1 = bz(0:nx+1, 0:ny+1, 0:nz+1)

     ! step 1
     CALL rkstep  

!default is first order in time   
#ifndef Q_FOURTHORDER  
     DO iz = 1, nz  
       izm = iz - 1
       DO iy = 1, ny
         iym = iy - 1
         DO ix = 0, nx
           ixm = ix - 1
           area = dzb(iz) * dyb(iy)
           bx(ix, iy, iz) = bx1(ix, iy, iz) &
            + (flux_z(ix,iy,iz) - flux_z(ix,iym,iz) + flux_z(ix,iy,izm) - flux_z(ix,iym,izm)) * dt / area &
            - (flux_y(ix,iy,iz) - flux_y(ix,iy,izm) + flux_y(ix,iym,iz) - flux_y(ix,iym,izm)) * dt / area
         END DO
       END DO
     END DO
     
     DO iz = 1, nz 
       izm = iz - 1
       DO iy = 0, ny
         iym = iy - 1
         DO ix = 1, nx
           ixm = ix - 1
           area = dxb(ix) * dzb(iz)
           by(ix, iy, iz) = by1(ix, iy, iz) &
            + (-flux_z(ix,iy,iz) + flux_z(ixm,iy,iz) - flux_z(ix,iy,izm) + flux_z(ixm,iy,izm)) * dt / area &
            + (flux_x(ix,iy,iz) - flux_x(ix,iy,izm) + flux_x(ixm,iy,iz) - flux_x(ixm,iy,izm)) * dt / area
         END DO
       END DO
     END DO
     
     DO iz = 0, nz
       izm = iz - 1
       DO iy = 1, ny
         iym = iy - 1
         DO ix = 1, nx
           ixm = ix - 1  
           area = dxb(ix) * dyb(iy)
           bz(ix, iy, iz) = bz1(ix, iy, iz) &
            + (flux_y(ix,iy,iz) - flux_y(ixm,iy,iz) + flux_y(ix,iym,iz) - flux_y(ixm,iym,iz)) * dt / area &
            - (flux_x(ix,iy,iz) - flux_x(ix,iym,iz) + flux_x(ixm,iy,iz) - flux_x(ixm,iym,iz)) * dt / area
         END DO
       END DO
     END DO
     CALL bfield_bcs

     DO iz = 1, nz
       izm = iz - 1
       DO iy = 1, ny
         iym = iy - 1
         DO ix = 1, nx
           ixm = ix - 1
           energy(ix, iy, iz) = energy(ix, iy, iz) &
               + (curlb(ix, iy, iz) + curlb(ixm, iy, iz) &
               + curlb(ix, iym, iz) + curlb(ix, iy, izm) &
               + curlb(ixm, iym, iz) + curlb(ixm, iy, izm) &
               + curlb(ix, iym, izm) + curlb(ixm, iym, izm)) &
               * dt / (8.0_num * rho(ix, iy, iz))
         END DO
       END DO
     END DO    
     
     CALL energy_bcs

     DO iz = 0, nz
       DO iy = 0, ny
         DO ix = 0, nx
           w1 = dt * dxc(ix) * dyc(iy) * dzc(iz) * curlb(ix, iy, iz)
           IF ((ix == 0) .OR. (ix == nx)) THEN
             w1 = w1 * 0.5_num
           END IF

           IF ((iy == 0) .OR. (iy == ny)) THEN
             w1 = w1 * 0.5_num
           END IF

           IF ((iz == 0) .OR. (iz == nz)) THEN
             w1 = w1 * 0.5_num
           END IF

           total_ohmic_heating = total_ohmic_heating + w1
         END DO
       END DO
     END DO
     
!if complier flag set then use 4th order Runge-Kutta     
#else    
     half_dt = dt * 0.5_num
     dt6 = dt * sixth

     k1x = flux_x
     k1y = flux_y
     k1z = flux_z
     c1 = curlb
     ! step 2
     DO iz = 1, nz  
       izm = iz - 1
       DO iy = 1, ny
         iym = iy - 1
         DO ix = 0, nx
           ixm = ix - 1
           area = dzb(iz) * dyb(iy)
           bx(ix, iy, iz) = bx1(ix, iy, iz) &
               + (k1z(ix,iy,iz) - k1z(ix,iym,iz) + k1z(ix,iy,izm) - k1z(ix,iym,izm)) * half_dt / area &
               - (k1y(ix,iy,iz) - k1y(ix,iy,izm) + k1y(ix,iym,iz) - k1y(ix,iym,izm)) * half_dt / area
         END DO
       END DO
     END DO

     DO iz = 1, nz 
       izm = iz - 1
       DO iy = 0, ny
         iym = iy - 1
         DO ix = 1, nx
           ixm = ix - 1
           area = dxb(ix) * dzb(iz)
           by(ix, iy, iz) = by1(ix, iy, iz) &
               + (-k1z(ix,iy,iz) + k1z(ixm,iy,iz) - k1z(ix,iy,izm) + k1z(ixm,iy,izm)) * half_dt / area &
               + (k1x(ix,iy,iz) - k1x(ix,iy,izm) + k1x(ixm,iy,iz) - k1x(ixm,iy,izm)) * half_dt / area
         END DO
       END DO
     END DO

     DO iz = 0, nz
       izm = iz - 1
       DO iy = 1, ny
         iym = iy - 1
         DO ix = 1, nx
           ixm = ix - 1  
           area = dxb(ix) * dyb(iy)
           bz(ix, iy, iz) = bz1(ix, iy, iz) &
               + (k1y(ix,iy,iz) - k1y(ixm,iy,iz) + k1y(ix,iym,iz) - k1y(ixm,iym,iz)) * half_dt / area &
               - (k1x(ix,iy,iz) - k1x(ix,iym,iz) + k1x(ixm,iy,iz) - k1x(ixm,iym,iz)) * half_dt / area
         END DO
       END DO
     END DO

     CALL bfield_bcs
     CALL rkstep
     k2x = flux_x
     k2y = flux_y
     k2z = flux_z
     c2 = curlb

     ! step 3
     DO iz = 1, nz  
       izm = iz - 1
       DO iy = 1, ny
         iym = iy - 1
         DO ix = 0, nx
           ixm = ix - 1
           area = dzb(iz) * dyb(iy)
           bx(ix, iy, iz) = bx1(ix, iy, iz) &
               + (k2z(ix,iy,iz) - k2z(ix,iym,iz) + k2z(ix,iy,izm) - k2z(ix,iym,izm)) * half_dt / area &
               - (k2y(ix,iy,iz) - k2y(ix,iy,izm) + k2y(ix,iym,iz) - k2y(ix,iym,izm)) * half_dt / area
         END DO
       END DO
     END DO

     DO iz = 1, nz 
       izm = iz - 1
       DO iy = 0, ny
         iym = iy - 1
         DO ix = 1, nx
           ixm = ix - 1
           area = dxb(ix) * dzb(iz)
           by(ix, iy, iz) = by1(ix, iy, iz) &
               + (-k2z(ix,iy,iz) + k2z(ixm,iy,iz) - k2z(ix,iy,izm) + k2z(ixm,iy,izm)) * half_dt / area &
               + (k2x(ix,iy,iz) - k2x(ix,iy,izm) + k2x(ixm,iy,iz) - k2x(ixm,iy,izm)) * half_dt / area
         END DO
       END DO
     END DO

     DO iz = 0, nz
       izm = iz - 1
       DO iy = 1, ny
         iym = iy - 1
         DO ix = 1, nx
           ixm = ix - 1  
           area = dxb(ix) * dyb(iy)
           bz(ix, iy, iz) = bz1(ix, iy, iz) &
               + (k2y(ix,iy,iz) - k2y(ixm,iy,iz) + k2y(ix,iym,iz) - k2y(ixm,iym,iz)) * half_dt / area &
               - (k2x(ix,iy,iz) - k2x(ix,iym,iz) + k2x(ixm,iy,iz) - k2x(ixm,iym,iz)) * half_dt / area
         END DO
       END DO
     END DO


     CALL bfield_bcs
     CALL rkstep
     k3x = flux_x
     k3y = flux_y
     k3z = flux_z
     c3 = curlb

     ! step 4
     DO iz = 1, nz  
       izm = iz - 1
       DO iy = 1, ny
         iym = iy - 1
         DO ix = 0, nx
           ixm = ix - 1
           area = dzb(iz) * dyb(iy)
           bx(ix, iy, iz) = bx1(ix, iy, iz) &
               + (k3z(ix,iy,iz) - k3z(ix,iym,iz) + k3z(ix,iy,izm) - k3z(ix,iym,izm)) * dt / area &
               - (k3y(ix,iy,iz) - k3y(ix,iy,izm) + k3y(ix,iym,iz) - k3y(ix,iym,izm)) * dt / area
         END DO
       END DO
     END DO

     DO iz = 1, nz 
       izm = iz - 1
       DO iy = 0, ny
         iym = iy - 1
         DO ix = 1, nx
           ixm = ix - 1
           area = dxb(ix) * dzb(iz)
           by(ix, iy, iz) = by1(ix, iy, iz) &
               + (-k3z(ix,iy,iz) + k3z(ixm,iy,iz) - k3z(ix,iy,izm) + k3z(ixm,iy,izm)) * dt / area &
               + (k3x(ix,iy,iz) - k3x(ix,iy,izm) + k3x(ixm,iy,iz) - k3x(ixm,iy,izm)) * dt / area
         END DO
       END DO
     END DO

     DO iz = 0, nz
       izm = iz - 1
       DO iy = 1, ny
         iym = iy - 1
         DO ix = 1, nx
           ixm = ix - 1  
           area = dxb(ix) * dyb(iy)
           bz(ix, iy, iz) = bz1(ix, iy, iz) &
               + (k3y(ix,iy,iz) - k3y(ixm,iy,iz) + k3y(ix,iym,iz) - k3y(ixm,iym,iz)) * dt / area &
               - (k3x(ix,iy,iz) - k3x(ix,iym,iz) + k3x(ixm,iy,iz) - k3x(ixm,iym,iz)) * dt / area
         END DO
       END DO
     END DO

     CALL bfield_bcs
     CALL rkstep
     k4x = flux_x
     k4y = flux_y
     k4z = flux_z
     c4 = curlb        

     ! full update
     k3x = k1x + 2.0_num * k2x + 2.0_num * k3x + k4x
     k3y = k1y + 2.0_num * k2y + 2.0_num * k3y + k4y
     k3z = k1z + 2.0_num * k2z + 2.0_num * k3z + k4z
     c1 = c1 + 2.0_num * c2 + 2.0_num * c3 + c4

     DO iz = 1, nz  
       izm = iz - 1
       DO iy = 1, ny
         iym = iy - 1
         DO ix = 0, nx
           ixm = ix - 1
           area = dzb(iz) * dyb(iy)
           bx(ix, iy, iz) = bx1(ix, iy, iz) &
               + (k3z(ix,iy,iz) - k3z(ix,iym,iz) + k3z(ix,iy,izm) - k3z(ix,iym,izm)) * dt6 / area &
               - (k3y(ix,iy,iz) - k3y(ix,iy,izm) + k3y(ix,iym,iz) - k3y(ix,iym,izm)) * dt6 / area
         END DO
       END DO
     END DO

     DO iz = 1, nz 
       izm = iz - 1
       DO iy = 0, ny
         iym = iy - 1
         DO ix = 1, nx
           ixm = ix - 1
           area = dxb(ix) * dzb(iz)
           by(ix, iy, iz) = by1(ix, iy, iz) &
               + (-k3z(ix,iy,iz) + k3z(ixm,iy,iz) - k3z(ix,iy,izm) + k3z(ixm,iy,izm)) * dt6 / area &
               + (k3x(ix,iy,iz) - k3x(ix,iy,izm) + k3x(ixm,iy,iz) - k3x(ixm,iy,izm)) * dt6 / area
         END DO
       END DO
     END DO

     DO iz = 0, nz
       izm = iz - 1
       DO iy = 1, ny
         iym = iy - 1
         DO ix = 1, nx
           ixm = ix - 1  
           area = dxb(ix) * dyb(iy)
           bz(ix, iy, iz) = bz1(ix, iy, iz) &
               + (k3y(ix,iy,iz) - k3y(ixm,iy,iz) + k3y(ix,iym,iz) - k3y(ixm,iym,iz)) * dt6 / area &
               - (k3x(ix,iy,iz) - k3x(ix,iym,iz) + k3x(ixm,iy,iz) - k3x(ixm,iym,iz)) * dt6 / area
         END DO
       END DO
     END DO

     CALL bfield_bcs

     DO iz = 1, nz
       izm = iz - 1
       DO iy = 1, ny
         iym = iy - 1
         DO ix = 1, nx
           ixm = ix - 1
           energy(ix, iy, iz) = energy(ix, iy, iz) &
               + (c1(ix, iy, iz) + c1(ixm, iy, iz) &
               + c1(ix, iym, iz) + c1(ix, iy, izm) &
               + c1(ixm, iym, iz) + c1(ixm, iy, izm) &
               + c1(ix, iym, izm) + c1(ixm, iym, izm)) &
               * dt6 / (8.0_num * rho(ix, iy, iz))
         END DO
       END DO
     END DO

     CALL energy_bcs

     DO iz = 0, nz
       DO iy = 0, ny
         DO ix = 0, nx
           w1 = dt6 * dxc(ix) * dyc(iy) * dzc(iz) * c1(ix, iy, iz)
           IF ((ix == 0) .OR. (ix == nx)) THEN
             w1 = w1 * 0.5_num
           END IF

           IF ((iy == 0) .OR. (iy == ny)) THEN
             w1 = w1 * 0.5_num
           END IF

           IF ((iz == 0) .OR. (iz == nz)) THEN
             w1 = w1 * 0.5_num
           END IF

           total_ohmic_heating = total_ohmic_heating + w1
         END DO
       END DO
     END DO

#endif

     DO iz = 0, nz
       izp = iz + 1
       DO iy = 0, ny
         iyp = iy + 1
         DO ix = 0, nx
           ixp = ix + 1

           jx1 = (bz(ix, iyp, iz) - bz(ix, iy, iz)) / dyc(iy) &
               - (by(ix, iy, izp) - by(ix, iy, iz)) / dzc(iz)
           jx2 = (bz(ixp, iyp, iz) - bz(ixp, iy, iz)) / dyc(iy) &
               - (by(ixp, iy, izp) - by(ixp, iy, iz)) / dzc(iz)
           jy1 = (bx(ix, iy, izp) - bx(ix, iy, iz)) / dzc(iz) &
               - (bz(ixp, iy, iz) - bz(ix, iy, iz)) / dxc(ix)
           jy2 = (bx(ix, iyp, izp) - bx(ix, iyp, iz)) / dzc(iz) &
               - (bz(ixp, iyp, iz) - bz(ix, iyp, iz)) / dxc(ix)
           jz1 = (by(ixp, iy, iz) - by(ix, iy, iz)) / dxc(ix) &
               - (bx(ix, iyp, iz) - bx(ix, iy, iz)) / dyc(iy)
           jz2 = (by(ixp, iy, izp) - by(ix, iy, izp)) / dxc(ix) &
               - (bx(ix, iyp, izp) - bx(ix, iy, izp)) / dyc(iy)

           jx_r(ix, iy, iz) = (jx1 + jx2) * 0.5_num
           jy_r(ix, iy, iz) = (jy1 + jy2) * 0.5_num
           jz_r(ix, iy, iz) = (jz1 + jz2) * 0.5_num
         END DO
       END DO
     END DO


#ifdef Q_FOURTHORDER
     DEALLOCATE(k1x, k2x, k3x, k4x, k1y, k2y, k3y, k4y, k1z, k2z, k3z, k4z)
     DEALLOCATE(c1, c2, c3, c4)
#endif

   END SUBROUTINE resistive_effects



   ! calculates 'k' values from b[xyz]1 values
   SUBROUTINE rkstep

     REAL(num) :: jx, jy, jz
     REAL(num) :: jx1, jy1, jz1, jx2, jy2, jz2
     REAL(num) :: bxv, byv, bzv
     REAL(num) :: magn_b
     REAL(num) :: j_par_x, j_par_y, j_par_z
     REAL(num) :: j_perp_x, j_perp_y, j_perp_z
     REAL(num) :: magn_j_perp, magn_j_par

     IF (.NOT. cowling_resistivity) THEN
       ! Use simple flux calculation
       DO iz = 0, nz
         izp = iz + 1
         DO iy = 0, ny
           iyp = iy + 1
           DO ix = 0, nx
             ixp = ix + 1

             jx1 = (bz(ix, iyp, iz) - bz(ix, iy, iz)) / dyc(iy) &
                 - (by(ix, iy, izp) - by(ix, iy, iz)) / dzc(iz)
             jx2 = (bz(ixp, iyp, iz) - bz(ixp, iy, iz)) / dyc(iy) &
                 - (by(ixp, iy, izp) - by(ixp, iy, iz)) / dzc(iz)
             jy1 = (bx(ix, iy, izp) - bx(ix, iy, iz)) / dzc(iz) &
                 - (bz(ixp, iy, iz) - bz(ix, iy, iz)) / dxc(ix)
             jy2 = (bx(ix, iyp, izp) - bx(ix, iyp, iz)) / dzc(iz) &
                 - (bz(ixp, iyp, iz) - bz(ix, iyp, iz)) / dxc(ix)
             jz1 = (by(ixp, iy, iz) - by(ix, iy, iz)) / dxc(ix) &
                 - (bx(ix, iyp, iz) - bx(ix, iy, iz)) / dyc(iy)
             jz2 = (by(ixp, iy, izp) - by(ix, iy, izp)) / dxc(ix) &
                 - (bx(ix, iyp, izp) - bx(ix, iy, izp)) / dyc(iy)

             jx = (jx1 + jx2) * 0.5_num
             jy = (jy1 + jy2) * 0.5_num
             jz = (jz1 + jz2) * 0.5_num 

             flux_x(ix, iy, iz) = -jx * eta(ix,iy,iz) * dxc(ix) * 0.5_num
             flux_y(ix, iy, iz) = -jy * eta(ix,iy,iz) * dyc(iy) * 0.5_num
             flux_z(ix, iy, iz) = -jz * eta(ix,iy,iz) * dzc(iz) * 0.5_num
             ! This isn't really curlb. It's actually heat flux
             curlb(ix, iy, iz) = eta(ix, iy, iz) * (jx**2 + jy**2 + jz**2)
           END DO
         END DO
       END DO
     ELSE
       ! Use partially ionised flux calculation
       DO iz = 0, nz 
         izp = iz + 1
         DO iy = 0, ny
           iyp = iy + 1
           DO ix = 0, nx
             ixp = ix + 1 

             jx1 = (bz(ix, iyp, iz) - bz(ix, iy, iz)) / dyc(iy) &
                 - (by(ix, iy, izp) - by(ix, iy, iz)) / dzc(iz)
             jx2 = (bz(ixp, iyp, iz) - bz(ixp, iy, iz)) / dyc(iy) &
                 - (by(ixp, iy, izp) - by(ixp, iy, iz)) / dzc(iz)
             jy1 = (bx(ix, iy, izp) - bx(ix, iy, iz)) / dzc(iz) &
                 - (bz(ixp, iy, iz) - bz(ix, iy, iz)) / dxc(ix)
             jy2 = (bx(ix, iyp, izp) - bx(ix, iyp, iz)) / dzc(iz) &
                 - (bz(ixp, iyp, iz) - bz(ix, iyp, iz)) / dxc(ix)
             jz1 = (by(ixp, iy, iz) - by(ix, iy, iz)) / dxc(ix) &
                 - (bx(ix, iyp, iz) - bx(ix, iy, iz)) / dyc(iy)
             jz2 = (by(ixp, iy, izp) - by(ix, iy, izp)) / dxc(ix) &
                 - (bx(ix, iyp, izp) - bx(ix, iy, izp)) / dyc(iy)

             jx = (jx1 + jx2) * 0.5_num
             jy = (jy1 + jy2) * 0.5_num
             jz = (jz1 + jz2) * 0.5_num

             ! B at vertices
             bxv = (bx(ix, iy, iz) + bx(ix, iyp, iz) + bx(ix, iy, izp) &
                 + bx(ix, iyp, izp)) * 0.25_num
             byv = (by(ix, iy, iz) + by(ixp, iy, iz) + by(ix, iy, izp) &
                 + by(ixp, iy, izp)) * 0.25_num
             bzv = (bz(ix, iy, iz) + bz(ixp, iy, iz) + bz(ix, iyp, iz) &
                 + bz(ixp, iyp, iz)) * 0.25_num
             magn_b = bxv**2 + byv**2 + bzv**2

             ! Calculate parallel and perpendicular currents
             j_par_x = (jx * bxv + jy * byv &
                 + jz * bzv) * bxv / MAX(magn_b, none_zero)
             j_par_y = (jx * bxv + jy * byv &
                 + jz * bzv) * byv / MAX(magn_b, none_zero)
             j_par_z = (jx * bxv + jy * byv &
                 + jz * bzv) * bzv / MAX(magn_b, none_zero)

             ! If b = 0 then there is no parallel current
             IF (magn_b .LT. none_zero) THEN
               j_par_x = 0.0_num
               j_par_y = 0.0_num
               j_par_z = 0.0_num
             END IF

             ! Calculate perpendicular current
             j_perp_x = jx - j_par_x
             j_perp_y = jy - j_par_y
             j_perp_z = jz - j_par_z

             magn_j_par = SQRT(j_par_x**2 + j_par_y**2 + j_par_z**2)
             magn_j_perp = SQRT(j_perp_x**2 + j_perp_y**2 + j_perp_z**2)

             parallel_current(ix, iy, iz) = magn_j_par
             perp_current(ix, iy, iz) = magn_j_perp

             ! This isn't really curlb. It's actually heat flux
             curlb(ix, iy, iz) = eta(ix, iy, iz) * magn_j_par**2 &
                 + (eta_perp(ix, iy, iz) + eta(ix, iy, iz)) * magn_j_perp**2

             flux_x(ix, iy, iz) = -((j_par_x * eta(ix, iy, iz) &
                 + j_perp_x * (eta_perp(ix, iy, iz) + eta(ix, iy, iz))) &
                 * dxc(ix) * 0.5_num)
             flux_y(ix, iy, iz) = -((j_par_y * eta(ix, iy, iz) &
                 + j_perp_y * (eta_perp(ix, iy, iz) + eta(ix, iy, iz))) &
                 * dyc(iy) * 0.5_num)
             flux_z(ix, iy, iz) = -((j_par_z * eta(ix, iy, iz) &
                 + j_perp_z * (eta_perp(ix, iy, iz) + eta(ix, iy, iz))) &
                 * dzc(iz) * 0.5_num)
           END DO
         END DO
       END DO
     END IF

   END SUBROUTINE rkstep



END MODULE lagran
