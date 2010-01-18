MODULE lagran

  !-----------------------------------------------------------------
  ! This subroutine performs the Lagrangian step
  ! Notes:
  ! There are !#DEC$ directives in this routine
  ! These override compilers vector analysis tools
  !-----------------------------------------------------------------

  USE shared_data
  USE boundary
  USE neutral
  USE diagnostics
  USE conduct
  USE eos

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

    IF (include_neutrals) CALL neutral_fraction(eos_number)
    IF (resistive_mhd .OR. hall_mhd) THEN
      ! if subcycling isn't wanted set dt = dtr in set_dt, don't just
      ! set substeps to 1.
      IF (resistive_mhd) THEN
        dt_sub = dtr
      ELSE
        dt_sub = dth
      END IF

      IF (resistive_mhd .AND. hall_mhd) dt_sub = MIN(dtr, dth)
      substeps = INT(dt / dt_sub) + 1

      IF (substeps > peak_substeps) peak_substeps = substeps
      actual_dt = dt
      dt = dt / REAL(substeps, num)

      DO subcycle = 1, substeps
        CALL eta_calc
        IF (include_neutrals) CALL neutral_fraction(eos_number)
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
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = 0, nx+1
          ixm = ix - 1
          bx1(ix, iy, iz) = (bx(ix, iy, iz) + bx(ixm, iy, iz)) / 2.0_num
          by1(ix, iy, iz) = (by(ix, iy, iz) + by(ix, iym, iz)) / 2.0_num
          bz1(ix, iy, iz) = (bz(ix, iy, iz) + bz(ix, iy, izm)) / 2.0_num
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

    dt2 = dt / 2.0_num
    CALL viscosity_and_b_update

    bx1 = bx1 * cv1(0:nx+1, 0:ny+1, 0:nz+1)
    by1 = by1 * cv1(0:nx+1, 0:ny+1, 0:nz+1)
    bz1 = bz1 * cv1(0:nx+1, 0:ny+1, 0:nz+1)

    DO iz = 0, nz+1
      DO iy = 0, ny+1
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
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
          CALL get_pressure(rho(ix, iy, iz) * cv(ix, iy, iz) / cv1(ix, iy, iz),&
              e1, eos_number, ix, iy, iz, pressure(ix, iy, iz))
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
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
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

          w1 = (p + pyp + pzp + pypzp) / 4.0_num
          w2 = (pxp + pxpyp + pxpzp + pxpypzp) / 4.0_num
          fx = -(w2 - w1) / dxc(ix)
          w1 = (p + pxp + pzp + pxpzp) / 4.0_num
          w2 = (pyp + pxpyp + pypzp + pxpypzp) / 4.0_num
          fy = -(w2 - w1) / dyc(iy)
          w1 = (p + pxp + pyp + pxpyp) / 4.0_num
          w2 = (pzp + pxpzp + pypzp + pxpypzp) / 4.0_num
          fz = -(w2 - w1) / dzc(iz)

          ! add diagonal components
          w1 = (qxx(ix, iy, iz) + qxx(ix, iyp, iz) &
              + qxx(ix, iy, izp) + qxx(ix, iyp, izp)) / 4.0_num
          w2 = (qxx(ixp, iy, iz) + qxx(ixp, iyp, iz) &
              + qxx(ixp, iy, izp) + qxx(ixp, iyp, izp)) / 4.0_num
          fx = fx + (w2 - w1) / dxc(ix)
  
          w1 = (qyy(ix, iy, iz) + qyy(ixp, iy, iz) &
              + qyy(ix, iy, izp) + qyy(ixp, iy, izp)) / 4.0_num
          w2 = (qyy(ix, iyp, iz) + qyy(ixp, iyp, iz) &
              + qyy(ix, iyp, izp) + qyy(ixp, iyp, izp)) / 4.0_num
          fy = fy + (w2 - w1) / dyc(iy)
  
          w1 = (qzz(ix, iy, iz) + qzz(ixp, iy, iz) &
              + qzz(ix, iyp, iz) + qzz(ixp, iyp, iz)) / 4.0_num
          w2 = (qzz(ix, iy, izp) + qzz(ixp, iy, izp) &
              + qzz(ix, iyp, izp) + qzz(ixp, iyp, izp)) / 4.0_num
          fz = fz + (w2 - w1) / dzc(iz)
  
          ! add shear viscosity forces
          ! fx
          w1 = (qxy(ix, iy, iz) + qxy(ixp, iy, iz) &
              + qxy(ix, iy, izp) + qxy(ixp, iy, izp)) / 4.0_num
          w2 = (qxy(ix, iyp, iz) + qxy(ixp, iyp, iz) &
              + qxy(ix, iyp, izp) + qxy(ixp, iyp, izp)) / 4.0_num
          fx = fx + (w2 - w1) / dyc(iy)
  
          w1 = (qxz(ix, iy, iz) + qxz(ixp, iy, iz) &
              + qxz(ix, iyp, iz) + qxz(ixp, iyp, iz)) / 4.0_num
          w2 = (qxz(ix, iy, izp) + qxz(ixp, iy, izp) &
              + qxz(ix, iyp, izp) + qxz(ixp, iyp, izp)) / 4.0_num
          fx = fx + (w2 - w1) / dzc(iz)
  
          ! fy
          w1 = (qxy(ix, iy, iz) + qxy(ix, iyp, iz) &
              + qxy(ix, iy, izp) + qxy(ix, iyp, izp)) / 4.0_num
          w2 = (qxy(ixp, iy, iz) + qxy(ixp, iyp, iz) &
              + qxy(ixp, iy, izp) + qxy(ixp, iyp, izp)) / 4.0_num
          fy = fy + (w2 - w1) / dxc(ix)
  
          w1 = (qyz(ix, iy, iz) + qyz(ixp, iy, iz) &
              + qyz(ix, iyp, iz) + qyz(ixp, iyp, iz)) / 4.0_num
          w2 = (qyz(ix, iy, izp) + qyz(ixp, iy, izp) &
              + qyz(ix, iyp, izp) + qyz(ixp, iyp, izp)) / 4.0_num
          fy = fy + (w2 - w1) / dzc(iz)
 
          ! fz
          w1 = (qxz(ix, iy, iz) + qxz(ix, iyp, iz) &
              + qxz(ix, iy, izp) + qxz(ix, iyp, izp)) / 4.0_num
          w2 = (qxz(ixp, iy, iz) + qxz(ixp, iyp, iz) &
              + qxz(ixp, iy, izp) + qxz(ixp, iyp, izp)) / 4.0_num
          fz = fz + (w2 - w1) / dxc(ix)
 
          w1 = (qyz(ix, iy, iz) + qyz(ixp, iy, iz) &
              + qyz(ix, iy, izp) + qyz(ixp, iy, izp)) / 4.0_num
          w2 = (qyz(ix, iyp, iz) + qyz(ixp, iyp, iz) &
              + qyz(ix, iyp, izp) + qyz(ixp, iyp, izp)) / 4.0_num
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

    IF (any_open) CALL store_boundary_dv

    CALL remap_v_bcs

    CALL visc_heating

    ! finally correct density and energy to final values
    DO iz = 1, nz
      izm = iz - 1
      DO iy = 1, ny
        iym = iy - 1
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = 1, nx
          ixm = ix - 1

          vxb = (vx1(ix, iy, iz) + vx1(ix, iym, iz) &
              + vx1(ix, iy, izm) + vx1(ix, iym, izm)) / 4.0_num

          vxbm = (vx1(ixm, iy, iz) + vx1(ixm, iym, iz) &
              + vx1(ixm, iy, izm) + vx1(ixm, iym, izm)) / 4.0_num

          vyb = (vy1(ix, iy, iz) + vy1(ixm, iy, iz) &
              + vy1(ix, iy, izm) + vy1(ixm, iy, izm)) / 4.0_num

          vybm = (vy1(ix, iym, iz) + vy1(ixm, iym, iz) &
              + vy1(ix, iym, izm) + vy1(ixm, iym, izm)) / 4.0_num

          vzb = (vz1(ix, iy, iz) + vz1(ixm, iy, iz) &
              + vz1(ix, iym, iz) + vz1(ixm, iym, iz)) / 4.0_num

          vzbm = (vz1(ix, iy, izm) + vz1(ixm, iy, izm) &
              + vz1(ix, iym, izm) + vz1(ixm, iym, izm)) / 4.0_num

          dv = ((vxb - vxbm) / dxb(ix) + (vyb - vybm) / dyb(iy) &
              + (vzb - vzbm) / dzb(iz)) * dt
          ! it is possible that dv has changed sign since the predictor step
          ! in this case p_visc * dv ought to be removed from the heating
          ! if QMONO is set - not done for simplicity since this is a rare
          ! combination

          cv1(ix, iy, iz) = cv(ix, iy, iz) * (1.0_num + dv)

          ! energy at end of Lagrangian step

          energy(ix, iy, iz) = energy(ix, iy, iz) &
              - pressure(ix, iy, iz) * dv / rho(ix, iy, iz)

          energy(ix, iy, iz) = energy(ix, iy, iz) &
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
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = -1, nx+2
          CALL get_pressure(rho(ix, iy, iz), energy(ix, iy, iz), eos_number, &
              ix, iy, iz, pressure(ix, iy, iz))
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
              + vx(ix, iy, izm) + vx(ix, iym, izm)) / 4.0_num

          ! vx at Sx(i-1, j, k)
          vxbm = (vx(ixm, iy, iz) + vx(ixm, iym, iz) &
              + vx(ixm, iy, izm) + vx(ixm, iym, izm)) / 4.0_num

          ! vy at Sy(i, j, k)
          vyb = (vy(ix, iy, iz) + vy(ixm, iy, iz) &
              + vy(ix, iy, izm) + vy(ixm, iy, izm)) / 4.0_num

          ! vy at Sy(i, j-1, k)
          vybm = (vy(ix, iym, iz) + vy(ixm, iym, iz) &
              + vy(ix, iym, izm) + vy(ixm, iym, izm)) / 4.0_num

          ! vz at Sz(i, j, k)
          vzb = (vz(ix, iy, iz) + vz(ixm, iy, iz) &
              + vz(ix, iym, iz) + vz(ixm, iym, iz)) / 4.0_num

          ! vz at Sz(i, j, k-1)
          vzbm = (vz(ix, iy, izm) + vz(ixm, iy, izm) &
              + vz(ix, iym, izm) + vz(ixm, iym, izm)) / 4.0_num

          dvxdx = (vxb - vxbm) / dxb(ix)
          dvydy = (vyb - vybm) / dyb(iy)
          dvzdz = (vzb - vzbm) / dzb(iz)

          dv = (dvxdx + dvydy + dvzdz) * dt2
          cv1(ix, iy, iz) = cv(ix, iy, iz) * (1.0_num + dv)

          ! vx at Sy(i, j, k)
          vxb = (vx(ix, iy, iz) + vx(ixm, iy, iz) &
              + vx(ix, iy, izm) + vx(ixm, iy, izm)) / 4.0_num

          ! vx at Sy(i, j-1, k)
          vxbm = (vx(ix, iym, iz) + vx(ixm, iym, iz) &
              + vx(ix, iym, izm) + vx(ixm, iym, izm)) / 4.0_num

          ! vy at Sx(i, j, k)
          vyb = (vy(ix, iy, iz) + vy(ix, iym, iz) &
              + vy(ix, iy, izm) + vy(ix, iym, izm)) / 4.0_num

          ! vy at Sx(i-1, j, k)
          vybm = (vy(ixm, iy, iz) + vy(ixm, iym, iz) &
              + vy(ixm, iy, izm) + vy(ixm, iym, izm)) / 4.0_num

          dvxdy = (vxb - vxbm) / dyb(iy)
          dvydx = (vyb - vybm) / dxb(ix)
          dvxy = dvxdy + dvydx
          sxy = dvxy / 2.0_num
          sxx = 2.0_num * dvxdx / 3.0_num - (dvydy + dvzdz) / 3.0_num
          syy = 2.0_num * dvydy / 3.0_num - (dvxdx + dvzdz) / 3.0_num
          szz = 2.0_num * dvzdz / 3.0_num - (dvxdx + dvydy) / 3.0_num

          ! vx at Sz(i, j, k)
          vxb = (vx(ix, iy, iz) + vx(ixm, iy, iz) &
              + vx(ix, iym, iz) + vx(ixm, iym, iz)) / 4.0_num

          ! vx at Sz(i, j, k-1)
          vxbm = (vx(ix, iy, izm) + vx(ixm, iy, izm) &
              + vx(ix, iym, izm) + vx(ixm, iym, izm)) / 4.0_num

          ! vz at Sx(i, j, k)
          vzb = (vz(ix, iy, iz) + vz(ix, iym, iz) &
              + vz(ix, iy, izm) + vz(ix, iym, izm)) / 4.0_num

          ! vz at Sx(i-1, j, k)
          vzbm = (vz(ixm, iy, iz) + vz(ixm, iym, iz) &
              + vz(ixm, iy, izm) + vz(ixm, iym, izm)) / 4.0_num

          dvxdz = (vxb - vxbm) / dzb(iz)
          dvzdx = (vzb - vzbm) / dxb(ix)
          dvxz = dvxdz + dvzdx
          sxz = dvxz / 2.0_num

          ! vy at Sz(i, j, k)
          vyb = (vy(ix, iy, iz) + vy(ixm, iy, iz) &
              + vy(ix, iym, iz) + vy(ixm, iym, iz)) / 4.0_num

          ! vy at Sz(i, j, k-1)
          vybm = (vy(ix, iy, izm) + vy(ixm, iy, izm) &
              + vy(ix, iym, izm) + vy(ixm, iym, izm)) / 4.0_num

          ! vz at Sy(i, j, k)
          vzb = (vz(ix, iy, iz) + vz(ixm, iy, iz) &
              + vz(ix, iy, izm) + vz(ixm, iy, izm)) / 4.0_num

          ! vz at Sy(i, j-1, k)
          vzbm = (vz(ix, iym, iz) + vz(ixm, iym, iz) &
              + vz(ix, iym, izm) + vz(ixm, iym, izm)) / 4.0_num

          dvydz = (vyb - vybm) / dzb(iz)
          dvzdy = (vzb - vzbm) / dyb(iy)
          dvyz = dvydz + dvzdy
          syz = dvyz / 2.0_num

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

             w2_1=dxb(ix)**2*w1 / MAX(fx**2, none_zero)
             w2_2=dzb(iz)**2*w1 / MAX(fz**2,none_zero)
             w2_3=dyb(iy)**2*w1 / MAX(fy**2,none_zero)
             w2_4=dzb(iz)**2*w1 / MAX(fz**2,none_zero)

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
          L = L * sg0 * (1.0_num - dvg0) + L * (1.0_num -sg0) * dvg0 +&
  				L * sg0 * dvg0

          w1 = (bx1(ix, iy, iz)**2 + by1(ix, iy, iz)**2 + bz1(ix, iy, iz)**2) &
                / rho(ix, iy, iz)

          CALL get_cs(rho(ix, iy, iz), energy(ix, iy, iz), eos_number, &
              ix, iy, iz, cs)
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
                 * (visc1 * cf + L2 * visc2 * ABS(s)))
            qxz(ix,iy,iz) = sxz * (L2 * rho(ix,iy,iz)  &
                 * (visc1 * cf + L2 * visc2 * ABS(s)))
            qyz(ix,iy,iz) = syz * (L2 * rho(ix,iy,iz)  &
                 * (visc1 * cf + L2 * visc2 * ABS(s)))
            qxx(ix,iy,iz) = sxx * (L2 * rho(ix,iy,iz)  &
                 * (visc1 * cf + L2 * visc2 * ABS(s)))
            qyy(ix,iy,iz) = syy * (L2 * rho(ix,iy,iz)  &
                 * (visc1 * cf + L2 * visc2 * ABS(s)))
            qzz(ix,iy,iz) = szz * (L2 * rho(ix,iy,iz)  &
                 * (visc1 * cf + L2 * visc2 * ABS(s)))
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
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = 0, nx+1
          ixm = ix - 1
          ixp = ix + 1

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

          dvxdx = (vxb - vxbm) / dxb(ix)
          dvydy = (vyb - vybm) / dyb(iy)
          dvzdz = (vzb - vzbm) / dzb(iz)

          vxb = (vx1(ix, iy, iz) + vx1(ixm, iy, iz) &
              + vx1(ix, iy, izm) + vx1(ixm, iy, izm)) / 4.0_num
          vxbm = (vx1(ix, iym, iz) + vx1(ixm, iym, iz) &
              + vx1(ix, iym, izm) + vx1(ixm, iym, izm)) / 4.0_num
          vyb = (vy1(ix, iy, iz) + vy1(ix, iym, iz) &
              + vy1(ix, iy, izm) + vy1(ix, iym, izm)) / 4.0_num
          vybm = (vy1(ixm, iy, iz) + vy1(ixm, iym, iz) &
              + vy1(ixm, iy, izm) + vy1(ixm, iym, izm)) / 4.0_num
          dvxy = (vxb - vxbm) / dyb(iy) + (vyb - vybm) / dxb(ix)

          vxb = (vx1(ix, iy, iz) + vx1(ixm, iy, iz) &
              + vx1(ix, iym, iz) + vx1(ixm, iym, iz)) / 4.0_num
          vxbm = (vx1(ix, iy, izm) + vx1(ixm, iy, izm) &
              + vx1(ix, iym, izm) + vx1(ixm, iym, izm)) / 4.0_num
          vzb = (vz1(ix, iy, iz) + vz1(ix, iym, iz) &
              + vz1(ix, iy, izm) + vz1(ix, iym, izm)) / 4.0_num
          vzbm = (vz1(ixm, iy, iz) + vz1(ixm, iym, iz) &
              + vz1(ixm, iy, izm) + vz1(ixm, iym, izm)) / 4.0_num
          dvxz = (vxb - vxbm) / dzb(iz) + (vzb - vzbm) / dxb(ix)

          vyb = (vy1(ix, iy, iz) + vy1(ixm, iy, iz) &
              + vy1(ix, iym, iz) + vy1(ixm, iym, iz)) / 4.0_num
          vybm = (vy1(ix, iy, izm) + vy1(ixm, iy, izm) &
              + vy1(ix, iym, izm) + vy1(ixm, iym, izm)) / 4.0_num
          vzb = (vz1(ix, iy, iz) + vz1(ixm, iy, iz) &
              + vz1(ix, iy, izm) + vz1(ixm, iy, izm)) / 4.0_num
          vzbm = (vz1(ix, iym, iz) + vz1(ixm, iym, iz) &
              + vz1(ix, iym, izm) + vz1(ixm, iym, izm)) / 4.0_num
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

    eta = 0.0_num

    DO iz = -1, nz+1
      izp = iz + 1
      DO iy = -1, ny+1
        iyp = iy + 1
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
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
          w4 = (jx + jxxp) / 2.0_num
          w5 = (jy + jyyp) / 2.0_num
          w6 = (jz + jzzp) / 2.0_num
          modj = SQRT(w4**2 + w5**2 + w6**2)

          IF (modj >= j_max) THEN
            eta(ix, iy, iz) = eta_background + eta0
          ELSE
            eta(ix, iy, iz) = eta_background 
          END IF
        END DO
      END DO
    END DO

    IF (.NOT. resistive_mhd) eta = 0.0_num

  END SUBROUTINE eta_calc



  ! Calculate the effect of resistivity on the magnetic field and Ohmic heating
  ! Use the subroutine rkstep
  SUBROUTINE resistive_effects

    REAL(num) :: dt6
    REAL(num) :: jx1, jx2, jy1, jy2, jz1, jz2
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

    dt = dt / 2.0_num

    bx1 = bx(0:nx+1, 0:ny+1, 0:nz+1)
    by1 = by(0:nx+1, 0:ny+1, 0:nz+1)
    bz1 = bz(0:nx+1, 0:ny+1, 0:nz+1)

    ! step 1
    CALL rkstep
    k1x = flux_x
    k1y = flux_y
    k1z = flux_z
    c1 = curlb

#ifdef Q_FIRSTORDER
    dt6 = dt 
    k3x = k1x 
    k3y = k1y 
    k3z = k1z 
#else
    ! step 2
    DO iz = 1, nz
      DO iy = 1, ny
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = 0, nx
          bx(ix, iy, iz) = bx1(ix, iy, iz) &
              + (k1z(ix, iy, iz) - k1z(ix, iy-1, iz)) * dt / dyb(iy) &
              - (k1y(ix, iy, iz) - k1y(ix, iy, iz-1)) * dt / dzb(iz)
        END DO
      END DO
    END DO

    DO iz = 1, nz
      DO iy = 0, ny
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = 1, nx
          by(ix, iy, iz) = by1(ix, iy, iz) &
              + (-k1z(ix, iy, iz) + k1z(ix-1, iy, iz)) * dt / dxb(ix) &
              + (k1x(ix, iy, iz) - k1x(ix, iy, iz-1)) * dt / dzb(iz)
        END DO
      END DO
    END DO

    DO iz = 0, nz
      DO iy = 1, ny
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = 1, nx
          bz(ix, iy, iz) = bz1(ix, iy, iz) &
              + (k1y(ix, iy, iz) - k1z(ix-1, iy, iz)) * dt / dxb(ix) &
              - (k1x(ix, iy, iz) - k1x(ix, iy-1, iz)) * dt / dyb(iy)
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
      DO iy = 1, ny
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = 0, nx
          bx(ix, iy, iz) = bx1(ix, iy, iz) &
              + (k2z(ix, iy, iz) - k2z(ix, iy-1, iz)) * dt / dyb(iy) &
              - (k2y(ix, iy, iz) - k2y(ix, iy, iz-1)) * dt / dzb(iz)
        END DO
      END DO
    END DO

    DO iz = 1, nz
      DO iy = 0, ny
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = 1, nx
          by(ix, iy, iz) = by1(ix, iy, iz) &
              + (-k2z(ix, iy, iz) + k2z(ix-1, iy, iz)) * dt / dxb(ix) &
              + (k2x(ix, iy, iz) - k2x(ix, iy, iz-1)) * dt / dzb(iz)
        END DO
      END DO
    END DO

    DO iz = 0, nz
      DO iy = 1, ny
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = 1, nx
          bz(ix, iy, iz) = bz1(ix, iy, iz) &
              + (k2y(ix, iy, iz) - k2z(ix-1, iy, iz)) * dt / dxb(ix) &
              - (k2x(ix, iy, iz) - k2x(ix, iy-1, iz)) * dt / dyb(iy)
        END DO
      END DO
    END DO

    CALL bfield_bcs
    CALL rkstep
    k3x = flux_x
    k3y = flux_y
    k3z = flux_z
    c3 = curlb

    dt = 2.0_num * dt

    ! step 4
    DO iz = 1, nz
      DO iy = 1, ny
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = 0, nx
          bx(ix, iy, iz) = bx1(ix, iy, iz) &
              + (k3z(ix, iy, iz) - k3z(ix, iy-1, iz)) * dt / dyb(iy) &
              - (k3y(ix, iy, iz) - k3y(ix, iy, iz-1)) * dt / dzb(iz)
        END DO
      END DO
    END DO

    DO iz = 1, nz
      DO iy = 0, ny
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = 1, nx
          by(ix, iy, iz) = by1(ix, iy, iz) &
              + (-k3z(ix, iy, iz) + k3z(ix-1, iy, iz)) * dt / dxb(ix) &
              + (k3x(ix, iy, iz) - k3x(ix, iy, iz-1)) * dt / dzb(iz)
        END DO
      END DO
    END DO

    DO iz = 0, nz
      DO iy = 1, ny
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = 1, nx
          bz(ix, iy, iz) = bz1(ix, iy, iz) &
              + (k3y(ix, iy, iz) - k3z(ix-1, iy, iz)) * dt / dxb(ix) &
              - (k3x(ix, iy, iz) - k3x(ix, iy-1, iz)) * dt / dyb(iy)
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
    dt6 = dt / 6.0_num
    k3x = k1x + 2.0_num * k2x + 2.0_num * k3x + k4x
    k3y = k1y + 2.0_num * k2y + 2.0_num * k3y + k4y
    k3z = k1z + 2.0_num * k2z + 2.0_num * k3z + k4z
    c1 = c1 + 2.0_num * c2 + 2.0_num * c3 + c4
#endif

    DO iz = 1, nz
      DO iy = 1, ny
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = 0, nx
          bx(ix, iy, iz) = bx1(ix, iy, iz) &
              + (k3z(ix, iy, iz) - k3z(ix, iy-1, iz)) * dt6 / dyb(iy) &
              - (k3y(ix, iy, iz) - k3y(ix, iy, iz-1)) * dt6 / dzb(iz)
        END DO
      END DO
    END DO

    DO iz = 1, nz
      DO iy = 0, ny
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = 1, nx
          by(ix, iy, iz) = by1(ix, iy, iz) &
              + (-k3z(ix, iy, iz) + k3z(ix-1, iy, iz)) * dt6 / dxb(ix) &
              + (k3x(ix, iy, iz) - k3x(ix, iy, iz-1)) * dt6 / dzb(iz)
        END DO
      END DO
    END DO

    DO iz = 0, nz
      DO iy = 1, ny
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = 1, nx
          bz(ix, iy, iz) = bz1(ix, iy, iz) &
              + (k3y(ix, iy, iz) - k3z(ix-1, iy, iz)) * dt6 / dxb(ix) &
              - (k3x(ix, iy, iz) - k3x(ix, iy-1, iz)) * dt6 / dyb(iy)
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

    DO iz = 0, nz
      izp = iz + 1
      DO iy = 0, ny
        iyp = iy + 1
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = 0, nx
          ixp = ix + 1

          ! jx at E3(i, j)
          jx1 = (bz(ix, iyp, iz) - bz(ix, iy, iz)) / dyc(iy) &
              - (by(ix, iy, izp) - by(ix, iy, iz)) / dzc(iz)

          ! jx at E3(i+1, j)
          jx2 = (bz(ixp, iyp, iz) - bz(ixp, iy, iz)) / dyc(iy) &
              - (by(ixp, iy, izp) - by(ixp, iy, iz)) / dzc(iz)

          ! jy at E2(i, j)
          jy1 = (bx(ix, iy, izp) - bx(ix, iy, iz)) / dzc(iz) &
              - (bz(ixp, iy, iz) - bz(ix, iy, iz)) / dxc(ix)

          ! jy at E2(i, j+1)
          jy2 = (bx(ix, iyp, izp) - bx(ix, iyp, iz)) / dzc(iz) &
              - (bz(ixp, iyp, iz) - bz(ix, iyp, iz)) / dxc(ix)

          ! jz at E1(i, j)
          jz1 = (by(ixp, iy, iz) - by(ix, iy, iz)) / dxc(ix) &
              - (bx(ix, iyp, iz) - bx(ix, iy, iz)) / dyc(iy)

          ! jz at E1(i, j)
          jz2 = (by(ixp, iy, izp) - by(ix, iy, izp)) / dxc(ix) &
              - (bx(ix, iyp, izp) - bx(ix, iy, izp)) / dyc(iy)

          jx_r(ix, iy, iz) = (jx1 + jx2) / 2.0_num
          jy_r(ix, iy, iz) = (jy1 + jy2) / 2.0_num
          jz_r(ix, iy, iz) = (jz1 + jz2) / 2.0_num
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

    ! Once more to get j_perp and j_par correct
    CALL rkstep

    DEALLOCATE(k1x, k2x, k3x, k4x, k1y, k2y, k3y, k4y, k1z, k2z, k3z, k4z)
    DEALLOCATE(c1, c2, c3, c4)

  END SUBROUTINE resistive_effects



  ! calculates 'k' values from b[xyz]1 values
  SUBROUTINE rkstep

    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: jx, jy, jz
    REAL(num) :: jx1, jy1, jz1, jx2, jy2, jz2
    REAL(num) :: bxv, byv, bzv
    REAL(num) :: magn_b
    REAL(num) :: j_par_x, j_par_y, j_par_z
    REAL(num) :: j_perp_x, j_perp_y, j_perp_z
    REAL(num) :: magn_j_perp, magn_j_par

    ALLOCATE(jx(-1:nx+1, -1:ny+1, -1:nz+1), &
        jy(-1:nx+1, -1:ny+1, -1:nz+1), jz(-1:nx+1, -1:ny+1, -1:nz+1))

    DO iz = -1, nz+1
      izp = iz + 1
      DO iy = -1, ny+1
        iyp = iy + 1
        !DEC$ IVDEP
        !DEC$ VECTOR ALWAYS
        DO ix = -1, nx+1
          ixp = ix + 1

          ! jx at E3(i, j)
          jx1 = (bz(ix, iyp, iz) - bz(ix, iy, iz)) / dyc(iy) &
              - (by(ix, iy, izp) - by(ix, iy, iz)) / dzc(iz)

          ! jx at E3(i+1, j)
          jx2 = (bz(ixp, iyp, iz) - bz(ixp, iy, iz)) / dyc(iy) &
              - (by(ixp, iy, izp) - by(ixp, iy, iz)) / dzc(iz)

          ! jy at E2(i, j)
          jy1 = (bx(ix, iy, izp) - bx(ix, iy, iz)) / dzc(iz) &
              - (bz(ixp, iy, iz) - bz(ix, iy, iz)) / dxc(ix)

          ! jy at E2(i, j+1)
          jy2 = (bx(ix, iyp, izp) - bx(ix, iyp, iz)) / dzc(iz) &
              - (bz(ixp, iyp, iz) - bz(ix, iyp, iz)) / dxc(ix)

          ! jz at E1(i, j)
          jz1 = (by(ixp, iy, iz) - by(ix, iy, iz)) / dxc(ix) &
              - (bx(ix, iyp, iz) - bx(ix, iy, iz)) / dyc(iy)

          ! jz at E1(i, j)
          jz2 = (by(ixp, iy, izp) - by(ix, iy, izp)) / dxc(ix) &
              - (bx(ix, iyp, izp) - bx(ix, iy, izp)) / dyc(iy)

          jx(ix, iy, iz) = (jx1 + jx2) / 2.0_num
          jy(ix, iy, iz) = (jy1 + jy2) / 2.0_num
          jz(ix, iy, iz) = (jz1 + jz2) / 2.0_num
        END DO
      END DO
    END DO

    IF (.NOT. cowling_resistivity) THEN
      ! Use simple flux calculation
      DO iz = 0, nz
        DO iy = 0, ny
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
          DO ix = 0, nx
            flux_x(ix, iy, iz) = -jx(ix, iy, iz) &
                * eta(ix, iy, iz) * dxc(ix) / 2.0_num
            flux_y(ix, iy, iz) = -jy(ix, iy, iz) &
                * eta(ix, iy, iz) * dyc(iy) / 2.0_num
            flux_z(ix, iy, iz) = -jz(ix, iy, iz) &
                * eta(ix, iy, iz) * dzc(iz) / 2.0_num
            ! This isn't really curlb. It's actually heat flux
            curlb(ix, iy, iz) = eta(ix, iy, iz) &
                * (jx(ix, iy, iz)**2 + jy(ix, iy, iz)**2 + jz(ix, iy, iz)**2)
          END DO
        END DO
      END DO
    ELSE
      ! Use partially ionised flux calculation
      DO iz = 0, nz
        DO iy = 0, ny
          iyp = iy + 1
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
          DO ix = 0, nx
            ixp = ix + 1
            ! B at vertices
            bxv = (bx(ix, iy, iz) + bx(ix, iyp, iz) + bx(ix, iy, izp) &
                + bx(ix, iyp, izp)) / 4.0_num
            byv = (by(ix, iy, iz) + by(ixp, iy, iz) + by(ix, iy, izp) &
                + by(ixp, iy, izp)) / 4.0_num
            bzv = (bz(ix, iy, iz) + bz(ixp, iy, iz) + bz(ix, iyp, iz) &
                + bz(ixp, iyp, iz)) / 4.0_num
            magn_b = bxv**2 + byv**2 + bzv**2

            ! Calculate parallel and perpendicular currents
            j_par_x = (jx(ix, iy, iz) * bxv + jy(ix, iy, iz) * byv &
                + jz(ix, iy, iz) * bzv) * bxv / MAX(magn_b, none_zero)
            j_par_y = (jx(ix, iy, iz) * bxv + jy(ix, iy, iz) * byv &
                + jz(ix, iy, iz) * bzv) * byv / MAX(magn_b, none_zero)
            j_par_z = (jx(ix, iy, iz) * bxv + jy(ix, iy, iz) * byv &
                + jz(ix, iy, iz) * bzv) * bzv / MAX(magn_b, none_zero)

            ! If b = 0 then there is no parallel current
            IF (magn_b .LT. none_zero) THEN
              j_par_x = 0.0_num
              j_par_y = 0.0_num
              j_par_z = 0.0_num
            END IF

            ! Calculate perpendicular current
            j_perp_x = jx(ix, iy, iz) - j_par_x
            j_perp_y = jy(ix, iy, iz) - j_par_y
            j_perp_z = jz(ix, iy, iz) - j_par_z

            magn_j_par = SQRT(j_par_x**2 + j_par_y**2 + j_par_z**2)
            magn_j_perp = SQRT(j_perp_x**2 + j_perp_y**2 + j_perp_z**2)

            parallel_current(ix, iy, iz) = magn_j_par
            perp_current(ix, iy, iz) = magn_j_perp

            ! This isn't really curlb. It's actually heat flux
            curlb(ix, iy, iz) = eta(ix, iy, iz) * magn_j_par**2 &
                + (eta_perp(ix, iy, iz) + eta(ix, iy, iz)) * magn_j_perp**2

            flux_x(ix, iy, iz) = -((j_par_x * eta(ix, iy, iz) &
                + j_perp_x * (eta_perp(ix, iy, iz) + eta(ix, iy, iz))) &
                * dxc(ix) / 2.0_num)
            flux_y(ix, iy, iz) = -((j_par_y * eta(ix, iy, iz) &
                + j_perp_y * (eta_perp(ix, iy, iz) + eta(ix, iy, iz))) &
                * dyc(iy) / 2.0_num)
            flux_z(ix, iy, iz) = -((j_par_z * eta(ix, iy, iz) &
                + j_perp_z * (eta_perp(ix, iy, iz) + eta(ix, iy, iz))))
          END DO
        END DO
      END DO
    END IF

    DEALLOCATE (jx, jy, jz)

  END SUBROUTINE rkstep



  SUBROUTINE store_boundary_dv

    REAL(num) :: dvx, dvy, dvz

    IF (xbc_right == BC_OPEN .AND. right == MPI_PROC_NULL) THEN
      DO iz = -2, nz+2
        DO iy = -2, ny+2
          dvx = 2.0_num * (vx(nx, iy, iz) - vx1(nx, iy, iz))
          dvy = 2.0_num * (vy(nx, iy, iz) - vy1(nx, iy, iz))
          dvz = 2.0_num * (vz(nx, iy, iz) - vz1(nx, iy, iz))
          dv_right(iy, iz) = SQRT(dvx**2 + dvy**2 + dvz**2)
        END DO
      END DO
    END IF

    IF (xbc_left == BC_OPEN .AND. left == MPI_PROC_NULL) THEN
      DO iz = -2, nz+2
        DO iy = -2, ny+2
          dvx = 2.0_num * (vx(0, iy, iz) - vx1(0, iy, iz))
          dvy = 2.0_num * (vy(0, iy, iz) - vy1(0, iy, iz))
          dvz = 2.0_num * (vz(0, iy, iz) - vz1(0, iy, iz))
          dv_left(iy, iz) = SQRT(dvx**2 + dvy**2 + dvz**2)
        END DO
      END DO
    END IF

    IF (ybc_up == BC_OPEN .AND. up == MPI_PROC_NULL) THEN
      DO iz = -2, nz+2
        DO ix = -2, nx+2
          dvx = 2.0_num * (vx(ix, ny, iz) - vx1(ix, ny, iz))
          dvy = 2.0_num * (vy(ix, ny, iz) - vy1(ix, ny, iz))
          dvz = 2.0_num * (vz(ix, ny, iz) - vz1(ix, ny, iz))
          dv_up(ix, iz) = SQRT(dvx**2 + dvy**2 + dvz**2)
        END DO
      END DO
    END IF

    IF (ybc_down == BC_OPEN .AND. down == MPI_PROC_NULL) THEN
      DO iz = -2, nz+2
        DO ix = -2, nx+2
          dvx = 2.0_num * (vx(ix, 0, iz) - vx1(ix, 0, iz))
          dvy = 2.0_num * (vy(ix, 0, iz) - vy1(ix, 0, iz))
          dvz = 2.0_num * (vz(ix, 0, iz) - vz1(ix, 0, iz))
          dv_down(ix, iz) = SQRT(dvx**2 + dvy**2 + dvz**2)
        END DO
      END DO
    END IF

    IF (zbc_back == BC_OPEN .AND. back == MPI_PROC_NULL) THEN
      DO iy = -2, ny+2
        DO ix = -2, nx+2
          dvx = 2.0_num * (vx(ix, iy, nz) - vx1(ix, iy, nz))
          dvy = 2.0_num * (vy(ix, iy, nz) - vy1(ix, iy, nz))
          dvz = 2.0_num * (vz(ix, iy, nz) - vz1(ix, iy, nz))
          dv_back(ix, iy) = SQRT(dvx**2 + dvy**2 + dvz**2)
        END DO
      END DO
    END IF

    IF (zbc_front == BC_OPEN .AND. front == MPI_PROC_NULL) THEN
      DO iy = -2, ny+2
        DO ix = -2, nx+2
          dvx = 2.0_num * (vx(ix, iy, 0) - vx1(ix, iy, 0))
          dvy = 2.0_num * (vy(ix, iy, 0) - vy1(ix, iy, 0))
          dvz = 2.0_num * (vz(ix, iy, 0) - vz1(ix, iy, 0))
          dv_front(ix, iy) = SQRT(dvx**2 + dvy**2 + dvz**2)
        END DO
      END DO
    END IF

  END SUBROUTINE store_boundary_dv

END MODULE lagran
