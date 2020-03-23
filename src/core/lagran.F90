!******************************************************************************
! Lagrangian step routines
!******************************************************************************

MODULE lagran

  USE shared_data
  USE boundary
  USE neutral
  USE conduct
  USE radiative
  USE openboundary
  USE remap

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lagrangian_step, eta_calc, set_dt

  ! Only used inside lagran.f90
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: alpha1, alpha2, alpha3
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: visc_heat, pressure, rho_v, cv_v
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: flux_x, flux_y, flux_z, curlb
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: fx_visc, fy_visc, fz_visc
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: bx0, by0, bz0
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: qx, qy, qz
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: energy0, delta_energy

CONTAINS

  !****************************************************************************
  ! This subroutine manages the progress of the lagrangian step
  !****************************************************************************

  SUBROUTINE lagrangian_step

    INTEGER :: substeps, subcycle
    REAL(num) :: actual_dt, dt_sub

#ifndef CAUCHY
    ALLOCATE(bx0(-2:nx+2,-1:ny+2,-1:nz+2))
    ALLOCATE(by0(-1:nx+2,-2:ny+2,-1:nz+2))
    ALLOCATE(bz0(-1:nx+2,-1:ny+2,-2:nz+2))
#endif
    ALLOCATE(bx1(-1:nx+2,-1:ny+2,-1:nz+2))
    ALLOCATE(by1(-1:nx+2,-1:ny+2,-1:nz+2))
    ALLOCATE(bz1(-1:nx+2,-1:ny+2,-1:nz+2))
    ALLOCATE(alpha1(0:nx+1,0:ny+2,0:nz+2))
    ALLOCATE(alpha2(-1:nx+1,0:ny+1,0:nz+2))
    ALLOCATE(alpha3(-1:nx+1,-1:ny+1,0:nz+1))
    ALLOCATE(visc_heat(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(pressure(-1:nx+2,-1:ny+2,-1:nz+2))
    ALLOCATE(rho_v(-1:nx+1,-1:ny+1,-1:nz+1))
    ALLOCATE(cv_v(-1:nx+1,-1:ny+1,-1:nz+1))
    ALLOCATE(fx_visc(0:nx,0:ny,0:nz))
    ALLOCATE(fy_visc(0:nx,0:ny,0:nz))
    ALLOCATE(fz_visc(0:nx,0:ny,0:nz))
    ALLOCATE(flux_x(0:nx,0:ny,0:nz))
    ALLOCATE(flux_y(0:nx,0:ny,0:nz))
    ALLOCATE(flux_z(0:nx,0:ny,0:nz))
    ALLOCATE(curlb (0:nx,0:ny,0:nz))
    ALLOCATE(qx(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(qy(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(qz(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(energy0(-1:nx+2,-1:ny+2,-1:nz+2))
    ALLOCATE(delta_energy(-1:nx+2,-1:ny+2,-1:nz+2))

    DO iz = -1, nz + 2
      izm = iz - 1
      DO iy = -1, ny + 2
        iym = iy - 1
        DO ix = -1, nx + 2
          ixm = ix - 1
          bx1(ix,iy,iz) = (bx(ix,iy,iz) + bx(ixm,iy ,iz )) * 0.5_num
          by1(ix,iy,iz) = (by(ix,iy,iz) + by(ix ,iym,iz )) * 0.5_num
          bz1(ix,iy,iz) = (bz(ix,iy,iz) + bz(ix ,iy ,izm)) * 0.5_num

          pressure(ix,iy,iz) = (gamma - 1.0_num) * rho(ix,iy,iz) &
                * (energy(ix,iy,iz) - (1.0_num - xi_n(ix,iy,iz)) * ionise_pot)
        END DO
      END DO
    END DO

    DO iz = -1, nz + 1
      izp = iz + 1
      DO iy = -1, ny + 1
        iyp = iy + 1
        DO ix = -1, nx + 1
          ixp = ix + 1
          rho_v(ix,iy,iz) = rho(ix ,iy ,iz ) * cv(ix ,iy ,iz ) &
              +   rho(ixp,iy ,iz ) * cv(ixp,iy ,iz ) &
              +   rho(ix ,iyp,iz ) * cv(ix ,iyp,iz ) &
              +   rho(ixp,iyp,iz ) * cv(ixp,iyp,iz ) &
              +   rho(ix ,iy ,izp) * cv(ix ,iy ,izp) &
              +   rho(ixp,iy ,izp) * cv(ixp,iy ,izp) &
              +   rho(ix ,iyp,izp) * cv(ix ,iyp,izp) &
              +   rho(ixp,iyp,izp) * cv(ixp,iyp,izp)

         cv_v(ix,iy,iz) = cv(ix,iy ,iz ) + cv(ixp,iy ,iz ) &
              + cv(ix,iyp,iz ) + cv(ixp,iyp,iz ) &
              + cv(ix,iy ,izp) + cv(ixp,iy ,izp) &
              + cv(ix,iyp,izp) + cv(ixp,iyp,izp)

         rho_v(ix,iy,iz) = rho_v(ix ,iy ,iz ) / cv_v(ix ,iy ,iz ) 

         cv_v(ix,iy,iz) = 0.125_num * cv_v(ix,iy,iz) 

        END DO
      END DO
    END DO

    IF (coronal_heating) CALL user_defined_heating

    IF (use_viscous_damping) CALL viscous_damping
    CALL shock_viscosity
    CALL set_dt
    dt2 = 0.5_num * dt

    IF (resistive_mhd) THEN
      ! If subcycling isn't wanted set dt = dtr in set_dt, don't just
      ! set substeps to 1.
      dt_sub = dtr

      substeps = INT(dt / dt_sub) + 1

      IF (substeps > peak_substeps) peak_substeps = substeps
      actual_dt = dt
      dt = dt / REAL(substeps, num)

      DO subcycle = 1, substeps
        IF (resistive_mhd) CALL eta_calc
        IF (eos_number /= EOS_IDEAL) CALL neutral_fraction
        IF (cowling_resistivity) CALL perpendicular_resistivity
        IF (resistive_mhd) CALL resistive_effects
      END DO

      DO iz = -1, nz + 2
        izm = iz - 1
        DO iy = -1, ny + 2
          iym = iy - 1
          DO ix = -1, nx + 2
            ixm = ix - 1
            bx1(ix,iy,iz) = (bx(ix,iy,iz) + bx(ixm,iy ,iz )) * 0.5_num
            by1(ix,iy,iz) = (by(ix,iy,iz) + by(ix ,iym,iz )) * 0.5_num
            bz1(ix,iy,iz) = (bz(ix,iy,iz) + bz(ix ,iy ,izm)) * 0.5_num
          END DO
        END DO
      END DO

      dt = actual_dt
    END IF

    delta_energy = 0.0_num
    IF (conduction) THEN
      energy0 = energy
      CALL conduct_heat
      delta_energy = energy - energy0
    END IF
    IF (radiation) THEN
      energy0 = energy
      CALL rad_losses
      delta_energy = delta_energy + energy - energy0
    END IF
    energy = energy + delta_energy
    energy = MAX(energy, 0.0_num)

    CALL predictor_corrector_step

    DEALLOCATE(bx1, by1, bz1, alpha1, alpha2, alpha3)
    DEALLOCATE(visc_heat, pressure, rho_v, cv_v)
    DEALLOCATE(fx_visc, fy_visc, fz_visc)
    DEALLOCATE(flux_x, flux_y, flux_z)
    DEALLOCATE(curlb)
    DEALLOCATE(qx, qy, qz)
    DEALLOCATE(energy0, delta_energy)
#ifndef CAUCHY
    DEALLOCATE(bx0, by0, bz0)
#endif

    CALL energy_bcs
    CALL density_bcs
    CALL velocity_bcs

  END SUBROUTINE lagrangian_step



  !****************************************************************************
  ! The main predictor / corrector step which advances the momentum equation
  !****************************************************************************

  SUBROUTINE predictor_corrector_step

    REAL(num) :: pp, ppx, ppy, ppxy
    REAL(num) :: ppz, ppxz, ppyz, ppxyz
    REAL(num) :: e1
    REAL(num) :: vxb, vxbm, vyb, vybm, vzb, vzbm
    REAL(num) :: bxv, byv, bzv, jx, jy, jz
    REAL(num) :: jx1, jx2, jy1, jy2, jz1, jz2
#ifdef CAUCHY
    REAL(num) :: cvx, cvxp, cvy, cvyp, cvz, cvzp
#endif
    REAL(num) :: dv
    REAL(num) :: fx, fy, fz

#ifndef CAUCHY
    bx0(:,:,:) = bx(:,:,:)
    by0(:,:,:) = by(:,:,:)
    bz0(:,:,:) = bz(:,:,:)
#endif

    CALL b_field_and_cv1_update

    bx1(:,:,:) = bx1(:,:,:) * cv1(:,:,:)
    by1(:,:,:) = by1(:,:,:) * cv1(:,:,:)
    bz1(:,:,:) = bz1(:,:,:) * cv1(:,:,:)

    DO iz = 0, nz + 1
      DO iy = 0, ny + 1
        DO ix = 0, nx + 1
          dv = cv1(ix,iy,iz) / cv(ix,iy,iz) - 1.0_num
          ! Predictor energy

          e1 = energy(ix,iy,iz) - pressure(ix,iy,iz) * dv / rho(ix,iy,iz)
          e1 = e1 + visc_heat(ix,iy,iz) * dt2 / rho(ix,iy,iz)

          ! Now define the predictor step pressures
          pressure(ix,iy,iz) = (e1 - (1.0_num - xi_n(ix,iy,iz)) * ionise_pot) &
              * (gamma - 1.0_num) * rho(ix,iy,iz) * cv(ix,iy,iz) / cv1(ix,iy,iz)

        END DO
      END DO
    END DO

    DO iz = 0, nz
      izp = iz + 1
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = 0, nx
          ixp = ix + 1

          pp    = pressure(ix ,iy ,iz )
          ppx   = pressure(ixp,iy ,iz )
          ppy   = pressure(ix ,iyp,iz )
          ppz   = pressure(ix ,iy ,izp)
          ppxy  = pressure(ixp,iyp,iz )
          ppxz  = pressure(ixp,iy ,izp)
          ppyz  = pressure(ix ,iyp,izp)
          ppxyz = pressure(ixp,iyp,izp)

          ! P total at Ex(i,j,k)
          w1 = (pp + ppy + ppz + ppyz) * 0.25_num
          ! P total at Ex(i+1,j,k)
          w2 = (ppx + ppxy + ppxz + ppxyz) * 0.25_num
          fx = -(w2 - w1) / dxc(ix)

          ! P total at Ey(i,j,k)
          w1 = (pp + ppx + ppz + ppxz) * 0.25_num
          ! P total at Ey(i,j+1,k)
          w2 = (ppy + ppxy + ppyz + ppxyz) * 0.25_num
          fy = -(w2 - w1) / dyc(iy)

          ! P total at Ez(i,j,k)
          w1 = (pp + ppx + ppy + ppxy) * 0.25_num
          ! P total at Ez(i,j,k+1)
          w2 = (ppz + ppxz + ppyz + ppxyz) * 0.25_num
          fz = -(w2 - w1) / dzc(iz)

#ifdef CAUCHY
          cvx  = cv1(ix ,iy ,iz ) + cv1(ix ,iyp,iz ) &
              +  cv1(ix ,iy ,izp) + cv1(ix ,iyp,izp)
          cvxp = cv1(ixp,iy ,iz ) + cv1(ixp,iyp,iz ) &
              +  cv1(ixp,iy ,izp) + cv1(ixp,iyp,izp)
          cvy  = cv1(ix ,iy ,iz ) + cv1(ixp,iy ,iz ) &
              +  cv1(ix ,iy ,izp) + cv1(ixp,iy ,izp)
          cvyp = cv1(ix ,iyp,iz ) + cv1(ixp,iyp,iz ) &
              +  cv1(ix ,iyp,izp) + cv1(ixp,iyp,izp)
          cvz  = cv1(ix ,iy ,iz ) + cv1(ixp,iy ,iz ) &
              +  cv1(ix ,iyp,iz ) + cv1(ixp,iyp,iz )
          cvzp = cv1(ix ,iy ,izp) + cv1(ixp,iy ,izp) &
              +  cv1(ix ,iyp,izp) + cv1(ixp,iyp,izp)

          w1 = (bz1(ix ,iy ,iz ) + bz1(ixp,iy ,iz ) &
              + bz1(ix ,iy ,izp) + bz1(ixp,iy ,izp)) / cvy
          w2 = (bz1(ix ,iyp,iz ) + bz1(ixp,iyp,iz ) &
              + bz1(ix ,iyp,izp) + bz1(ixp,iyp,izp)) / cvyp
          jx = (w2 - w1) / dyc(iy)

          w1 = (by1(ix ,iy ,iz ) + by1(ixp,iy ,iz ) &
              + by1(ix ,iyp,iz ) + by1(ixp,iyp,iz )) / cvz
          w2 = (by1(ix ,iy ,izp) + by1(ixp,iy ,izp) &
              + by1(ix ,iyp,izp) + by1(ixp,iyp,izp)) / cvzp
          jx = jx - (w2 - w1) / dzc(iz)

          w1 = (bz1(ix ,iy ,iz ) + bz1(ix ,iyp,iz ) &
              + bz1(ix ,iy ,izp) + bz1(ix ,iyp,izp)) / cvx
          w2 = (bz1(ixp,iy ,iz ) + bz1(ixp,iyp,iz ) &
              + bz1(ixp,iy ,izp) + bz1(ixp,iyp,izp)) / cvxp
          jy = -(w2 - w1) / dxc(ix)

          w1 = (bx1(ix ,iy ,iz ) + bx1(ixp,iy ,iz ) &
              + bx1(ix ,iyp,iz ) + bx1(ixp,iyp,iz )) / cvz
          w2 = (bx1(ix ,iy ,izp) + bx1(ixp,iy ,izp) &
              + bx1(ix ,iyp,izp) + bx1(ixp,iyp,izp)) / cvzp
          jy = jy + (w2 - w1) / dzc(iz)

          w1 = (by1(ix ,iy ,iz ) + by1(ix ,iyp,iz ) &
              + by1(ix ,iy ,izp) + by1(ix ,iyp,izp)) / cvx
          w2 = (by1(ixp,iy ,iz ) + by1(ixp,iyp,iz ) &
              + by1(ixp,iy ,izp) + by1(ixp,iyp,izp)) / cvxp
          jz = (w2 - w1) / dxc(ix)

          w1 = (bx1(ix ,iy ,iz ) + bx1(ixp,iy ,iz ) &
              + bx1(ix ,iy ,izp) + bx1(ixp,iy ,izp)) / cvy
          w2 = (bx1(ix ,iyp,iz ) + bx1(ixp,iyp,iz ) &
              + bx1(ix ,iyp,izp) + bx1(ixp,iyp,izp)) / cvyp
          jz = jz - (w2 - w1) / dyc(iy)

          bxv = (bx1(ix,iy ,iz ) + bx1(ixp,iy ,iz )  &
              +  bx1(ix,iyp,iz ) + bx1(ixp,iyp,iz )  &
              +  bx1(ix,iy ,izp) + bx1(ixp,iy ,izp)  &
              +  bx1(ix,iyp,izp) + bx1(ixp,iyp,izp)) &
              / (cvx + cvxp)

          byv = (by1(ix,iy ,iz ) + by1(ixp,iy ,iz )  &
              +  by1(ix,iyp,iz ) + by1(ixp,iyp,iz )  &
              +  by1(ix,iy ,izp) + by1(ixp,iy ,izp)  &
              +  by1(ix,iyp,izp) + by1(ixp,iyp,izp)) &
              / (cvx + cvxp)

          bzv = (bz1(ix,iy ,iz ) + bz1(ixp,iy ,iz )  &
              +  bz1(ix,iyp,iz ) + bz1(ixp,iyp,iz )  &
              +  bz1(ix,iy ,izp) + bz1(ixp,iy ,izp)  &
              +  bz1(ix,iyp,izp) + bz1(ixp,iyp,izp)) &
              / (cvx + cvxp)
#else

          bxv = 0.25_num * (bx(ix,iy,iz) + bx(ix,iyp,iz) + bx(ix,iy,izp) + bx(ix,iyp,izp))   
          byv = 0.25_num * (by(ix,iy,iz) + by(ixp,iy,iz) + by(ix,iy,izp) + by(ixp,iy,izp))   
          bzv = 0.25_num * (bz(ix,iy,iz) + bz(ix,iyp,iz) + bz(ixp,iy,iz) + bz(ixp,iyp,iz))   

          jx1 = (bz(ix ,iyp,iz ) - bz(ix ,iy ,iz )) / dyc(iy) &
              - (by(ix ,iy ,izp) - by(ix ,iy ,iz )) / dzc(iz)
          jx2 = (bz(ixp,iyp,iz ) - bz(ixp,iy ,iz )) / dyc(iy) &
              - (by(ixp,iy ,izp) - by(ixp,iy ,iz )) / dzc(iz)
          jy1 = (bx(ix ,iy ,izp) - bx(ix ,iy ,iz )) / dzc(iz) &
              - (bz(ixp,iy ,iz ) - bz(ix ,iy ,iz )) / dxc(ix)
          jy2 = (bx(ix ,iyp,izp) - bx(ix ,iyp,iz )) / dzc(iz) &
              - (bz(ixp,iyp,iz ) - bz(ix ,iyp,iz )) / dxc(ix)
          jz1 = (by(ixp,iy ,iz ) - by(ix ,iy ,iz )) / dxc(ix) &
              - (bx(ix ,iyp,iz ) - bx(ix ,iy ,iz )) / dyc(iy)
          jz2 = (by(ixp,iy ,izp) - by(ix ,iy ,izp)) / dxc(ix) &
              - (bx(ix ,iyp,izp) - bx(ix ,iy ,izp)) / dyc(iy)

          jx = (jx1 + jx2) * 0.5_num
          jy = (jy1 + jy2) * 0.5_num
          jz = (jz1 + jz2) * 0.5_num
#endif
          fx = fx + gamma_boris(ix,iy,iz) *(jy * bzv - jz * byv)
          fy = fy + gamma_boris(ix,iy,iz) *(jz * bxv - jx * bzv)
          fz = fz + gamma_boris(ix,iy,iz) *(jx * byv - jy * bxv)

          fz = fz - rho_v(ix,iy,iz) * grav(iz)

          ! Find half step velocity needed for remap
          vx1(ix,iy,iz) = vx(ix,iy,iz) + dt2 * (fx_visc(ix,iy,iz) + fx) / rho_v(ix,iy,iz)
          vy1(ix,iy,iz) = vy(ix,iy,iz) + dt2 * (fy_visc(ix,iy,iz) + fy) / rho_v(ix,iy,iz)
          vz1(ix,iy,iz) = vz(ix,iy,iz) + dt2 * (fz_visc(ix,iy,iz) + fz) / rho_v(ix,iy,iz)
        END DO
      END DO
    END DO

#ifndef CAUCHY
    bx(:,:,:) = bx0(:,:,:) 
    by(:,:,:) = by0(:,:,:) 
    bz(:,:,:) = bz0(:,:,:) 
#endif

    CALL remap_v_bcs

    CALL shock_heating

    DO iz = 0, nz
      DO iy = 0, ny
        DO ix = 0, nx
          vx(ix,iy,iz) = 2.0_num * vx1(ix,iy,iz) - vx(ix,iy,iz) 
          vy(ix,iy,iz) = 2.0_num * vy1(ix,iy,iz) - vy(ix,iy,iz) 
          vz(ix,iy,iz) = 2.0_num * vz1(ix,iy,iz) - vz(ix,iy,iz) 
        END DO
      END DO
    END DO

    CALL velocity_bcs
    IF (any_open) THEN
      CALL open_bcs         
    END IF 

    !Finally correct density and energy to final values
    DO iz = 1, nz
      izm = iz - 1
      DO iy = 1, ny
        iym = iy - 1
        DO ix = 1, nx
          ixm = ix - 1

          ! vx1 at Bx(i,j,k)
          vxb  = (vx1(ix ,iy ,iz ) + vx1(ix ,iym,iz ) &
              +   vx1(ix ,iy ,izm) + vx1(ix ,iym,izm)) * 0.25_num
          ! vx1 at Bx(i-1,j,k)
          vxbm = (vx1(ixm,iy ,iz ) + vx1(ixm,iym,iz ) &
              +   vx1(ixm,iy ,izm) + vx1(ixm,iym,izm)) * 0.25_num
          ! vy1 at By(i,j,k)
          vyb  = (vy1(ix ,iy ,iz ) + vy1(ixm,iy ,iz ) &
              +   vy1(ix ,iy ,izm) + vy1(ixm,iy ,izm)) * 0.25_num
          ! vy1 at By(i,j-1,k)
          vybm = (vy1(ix ,iym,iz ) + vy1(ixm,iym,iz ) &
              +   vy1(ix ,iym,izm) + vy1(ixm,iym,izm)) * 0.25_num
          ! vz1 at Bz(i,j,k)
          vzb  = (vz1(ix ,iy ,iz ) + vz1(ixm,iy ,iz ) &
              +   vz1(ix ,iym,iz ) + vz1(ixm,iym,iz )) * 0.25_num
          ! vz1 at Bz(i,j,k-1)
          vzbm = (vz1(ix ,iy ,izm) + vz1(ixm,iy ,izm) &
              +   vz1(ix ,iym,izm) + vz1(ixm,iym,izm)) * 0.25_num

          dv = ((vxb - vxbm) / dxb(ix) + (vyb - vybm) / dyb(iy) &
              + (vzb - vzbm) / dzb(iz)) * dt
          ! It is possible that dv has changed sign since the predictor step.
          ! In this case p_visc * dv ought to be removed from the heating

          cv1(ix,iy,iz) = cv(ix,iy,iz) * (1.0_num + dv)

          ! Energy at end of Lagrangian step
          energy(ix,iy,iz) = energy(ix,iy,iz) &
              + (dt * visc_heat(ix,iy,iz) - dv * pressure(ix,iy,iz)) &
              / rho(ix,iy,iz)

          visc_dep(ix,iy,iz) = visc_dep(ix,iy,iz) + dt * visc_heat(ix,iy,iz)

          rho(ix,iy,iz) = rho(ix,iy,iz) / (1.0_num + dv)

          total_visc_heating = total_visc_heating &
              + dt * visc_heat(ix,iy,iz) * cv(ix,iy,iz)

        END DO
      END DO
    END DO

    IF (cooling_term) THEN
      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            cool_term_v(ix,iy,iz) = alpha_av * dt * visc_heat(ix,iy,iz)/rho(ix,iy,iz) + &
                (1.0_num - alpha_av) * cool_term_v(ix,iy,iz)
  
            energy(ix,iy,iz) = energy(ix,iy,iz) - cool_term_v(ix,iy,iz)
          END DO
        END DO
      END DO
      energy = MAX(energy, 0.0_num)
    END IF

  END SUBROUTINE predictor_corrector_step




  SUBROUTINE viscous_damping

    qx = 0.0_num
    qy = 0.0_num
    qz = 0.0_num

    DO iz = 0, nz
      izp = iz + 1
      izm = iz - 1
      DO iy = 0, ny
        iyp = iy + 1
        iym = iy - 1
        DO ix = 0, nx
          ixp = ix + 1
          ixm = ix - 1
          qx(ix,iy,iz) = visc3(ix,iy,iz)  &
            * (((vx(ixp,iy,iz) - vx(ix,iy,iz))/ dxb(ixp) - (vx(ix,iy,iz) - vx(ixm,iy,iz))/ dxb(ix)) / dxc(ix) &
            +  ((vx(ix,iyp,iz) - vx(ix,iy,iz))/ dyb(iyp) - (vx(ix,iy,iz) - vx(ix,iym,iz))/ dyb(iy)) / dyc(iy) &
            +  ((vx(ix,iy,izp) - vx(ix,iy,iz))/ dzb(izp) - (vx(ix,iy,iz) - vx(ix,iy,izm))/ dzb(iz)) / dzc(iz)) 
          qy(ix,iy,iz) = visc3(ix,iy,iz)  &
            * (((vy(ixp,iy,iz) - vy(ix,iy,iz))/ dxb(ixp) - (vy(ix,iy,iz) - vy(ixm,iy,iz))/ dxb(ix)) / dxc(ix) &
            +  ((vy(ix,iyp,iz) - vy(ix,iy,iz))/ dyb(iyp) - (vy(ix,iy,iz) - vy(ix,iym,iz))/ dyb(iy)) / dyc(iy) &
            +  ((vy(ix,iy,izp) - vy(ix,iy,iz))/ dzb(izp) - (vy(ix,iy,iz) - vy(ix,iy,izm))/ dzb(iz)) / dzc(iz)) 
          qz(ix,iy,iz) = visc3(ix,iy,iz)  &
            * (((vz(ixp,iy,iz) - vz(ix,iy,iz))/ dxb(ixp) - (vz(ix,iy,iz) - vz(ixm,iy,iz))/ dxb(ix)) / dxc(ix) &
            +  ((vz(ix,iyp,iz) - vz(ix,iy,iz))/ dyb(iyp) - (vz(ix,iy,iz) - vz(ix,iym,iz))/ dyb(iy)) / dyc(iy) &
            +  ((vz(ix,iy,izp) - vz(ix,iy,iz))/ dzb(izp) - (vz(ix,iy,iz) - vz(ix,iy,izm))/ dzb(iz)) / dzc(iz)) 
        END DO
      END DO
    END DO

    DO iz = 0, nz
      DO iy = 0, ny
        DO ix = 0, nx
          vx(ix,iy,iz) = vx(ix,iy,iz) + dt * qx(ix,iy,iz)
          vy(ix,iy,iz) = vy(ix,iy,iz) + dt * qy(ix,iy,iz)
          vz(ix,iy,iz) = vz(ix,iy,iz) + dt * qz(ix,iy,iz)
        END DO
      END DO
    END DO

  END SUBROUTINE viscous_damping


  !****************************************************************************
  ! This subroutine calculates the viscous effects and updates the
  ! magnetic field
  !****************************************************************************

  SUBROUTINE shock_viscosity

    REAL(num) :: dvdots, dx, dxm, dxp
    REAL(num) :: b2, rmin
    REAL(num) :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: cs, cs_v
    INTEGER :: i0, i1, i2, i3, j0, j1, j2, j3, k0, k1, k2, k3
    LOGICAL, SAVE :: first_call = .TRUE.

    ALLOCATE(cs(-1:nx+2,-1:ny+2,-1:nz+2), cs_v(-1:nx+1,-1:ny+1,-1:nz+1))

    IF (first_call) THEN
      first_call = .FALSE.
      visc2_norm = 0.25_num * (gamma + 1.0_num) * visc2
    END IF

    p_visc = 0.0_num
    visc_heat = 0.0_num

    DO iz = -1, nz + 2
      DO ix = -1, nx + 2
        DO iy = -1, ny + 2
          rmin = MAX(rho(ix,iy,iz), none_zero)
          b2 = bx1(ix,iy,iz)**2 + by1(ix,iy,iz)**2 + bz1(ix,iy,iz)**2
          cs(ix,iy,iz) = SQRT((gamma * pressure(ix,iy,iz) + gamma_boris(ix,iy,iz) * b2) / rmin)
        END DO
      END DO
    END DO

    DO iz = -1, nz + 1
      izp = iz + 1
      DO iy = -1, ny + 1
        iyp = iy + 1
        DO ix = -1, nx + 1
          ixp = ix + 1

          cs_v(ix,iy,iz) = cs(ix ,iy ,iz ) * cv(ix ,iy ,iz ) &
              +   cs(ixp,iy ,iz ) * cv(ixp,iy ,iz ) &
              +   cs(ix ,iyp,iz ) * cv(ix ,iyp,iz ) &
              +   cs(ixp,iyp,iz ) * cv(ixp,iyp,iz ) &
              +   cs(ix ,iy ,izp) * cv(ix ,iy ,izp) &
              +   cs(ixp,iy ,izp) * cv(ixp,iy ,izp) &
              +   cs(ix ,iyp,izp) * cv(ix ,iyp,izp) &
              +   cs(ixp,iyp,izp) * cv(ixp,iyp,izp)
         
          cs_v(ix,iy,iz) = 0.125_num * cs_v(ix,iy,iz) / cv_v(ix,iy,iz)

        END DO
      END DO
    END DO

    DO iz = 0, nz + 2
      izm = iz - 1
      izp = iz + 1
      DO iy = 0, ny + 2
        iym = iy - 1
        iyp = iy + 1
        DO ix = 0, nx + 1
          ixm = ix - 1
          ixp = ix + 1

          ! Edge viscosities from Caramana
          i1 = ixm
          j1 = iym
          k1 = izm
          i2 = ix
          j2 = iym
          k2 = izm
          i0 = i1 - 1
          j0 = j1
          k0 = k1
          i3 = i2 + 1
          j3 = j2
          k3 = k2
          dx = dxb(ix)
          dxp = dxb(ixp)
          dxm = dxb(ixm)
          ! dv in direction of dS, i.e. dv.dS / abs(dS)
          dvdots = - (vx(i1,j1,k1) - vx(i2,j2,k2))
          ! Force on node is alpha*dv*ds but store only alpha and convert to force
          ! when needed.  
          alpha1(ix,iy,iz) = edge_viscosity()
        END DO
      END DO
    END DO

    DO iz = 0, nz + 2
      izm = iz - 1
      izp = iz + 1
        DO iy = 0, ny + 1
        iym = iy - 1
        iyp = iy + 1
        DO ix = -1, nx + 1
          ixm = ix - 1
          ixp = ix + 1

          i1 = ix
          j1 = iym
          k1 = izm
          i2 = ix
          j2 = iy
          k2 = izm
          i0 = i1
          j0 = j1 - 1
          k0 = k1
          i3 = i2
          j3 = j2 + 1
          k3 = k2
          dx = dyb(iy)
          dxp = dyb(iyp)
          dxm = dyb(iym)
          dvdots = - (vy(i1,j1,k1) - vy(i2,j2,k2))
          alpha2(ix,iy,iz) = edge_viscosity()
        END DO
      END DO
    END DO

    DO iz = 0, nz + 1
      izm = iz - 1
      izp = iz + 1
        DO iy = -1, ny + 1
        iym = iy - 1
        iyp = iy + 1
        DO ix = -1, nx + 1
          ixm = ix - 1
          ixp = ix + 1

          i1 = ix
          j1 = iy
          k1 = izm
          i2 = ix
          j2 = iy
          k2 = iz
          i0 = i1
          j0 = j1 
          k0 = k1 - 1
          i3 = i2
          j3 = j2 
          k3 = k2 + 1
          dx = dzb(iz)
          dxp = dzb(izp)
          dxm = dzb(izm)
          dvdots = - (vz(i1,j1,k1) - vz(i2,j2,k2))
          alpha3(ix,iy,iz) = edge_viscosity()
        END DO
      END DO
    END DO

    DO iz = 0, nz + 1 
      izm = iz - 1
      izp = iz + 1
      DO iy = 0, ny + 1 
        iym = iy - 1
        iyp = iy + 1
        DO ix = 0, nx + 1 
          ixm = ix - 1
          ixp = ix + 1
          ! Estimate p_visc based on alpha * dv, for timestep control
          a1 = ((vx(ixm,iym,izm) - vx(ix ,iym,izm))**2  &
              + (vy(ixm,iym,izm) - vy(ix ,iym,izm))**2 + (vz(ixm,iym,izm) - vz(ix ,iym,izm))**2) 
          a2 = ((vx(ix ,iym,izm) - vx(ix ,iy, izm))**2  &
              + (vy(ix ,iym,izm) - vy(ix ,iy ,izm))**2 + (vz(ix ,iym,izm) - vz(ix ,iy ,izm))**2)
          a3 = ((vx(ix ,iy ,izm) - vx(ixm,iy ,izm))**2  &
              + (vy(ix ,iy ,izm) - vy(ixm,iy ,izm))**2 + (vz(ix ,iy ,izm) - vz(ixm,iy ,izm))**2) 
          a4 = ((vx(ixm,iy ,izm) - vx(ixm,iym,izm))**2  &
              + (vy(ixm,iy ,izm) - vy(ixm,iym,izm))**2 + (vz(ixm,iy ,izm) - vz(ixm,iym,izm))**2)

          a5 = ((vx(ixm,iym,iz ) - vx(ix ,iym,iz ))**2  &
              + (vy(ixm,iym,iz ) - vy(ix ,iym,iz ))**2 + (vz(ixm,iym,iz ) - vz(ix ,iym,iz ))**2) 
          a6 = ((vx(ix ,iym,iz ) - vx(ix ,iy, iz ))**2  &
              + (vy(ix ,iym,iz ) - vy(ix ,iy ,iz ))**2 + (vz(ix ,iym,iz ) - vz(ix ,iy ,iz ))**2)
          a7 = ((vx(ix ,iy ,iz ) - vx(ixm,iy ,iz ))**2  &
              + (vy(ix ,iy ,iz ) - vy(ixm,iy ,iz ))**2 + (vz(ix ,iy ,iz ) - vz(ixm,iy ,iz ))**2) 
          a8 = ((vx(ixm,iy ,iz ) - vx(ixm,iym,iz ))**2  &
              + (vy(ixm,iy ,iz ) - vy(ixm,iym,iz ))**2 + (vz(ixm,iy ,iz ) - vz(ixm,iym,iz ))**2)

          a9 = ((vx(ix ,iy ,izm) - vx(ix ,iy ,iz ))**2  &
              + (vy(ix ,iy ,izm) - vy(ix ,iy ,iz))**2 + (vz(ix ,iy ,izm) - vz(ix ,iy ,iz))**2) 
          a10= ((vx(ixm,iy ,izm) - vx(ixm,iy ,iz ))**2  &
              + (vy(ixm,iy ,izm) - vy(ixm,iy ,iz))**2 + (vz(ixm,iy ,izm) - vz(ixm,iy ,iz))**2) 
          a11= ((vx(ixm,iym,izm) - vx(ixm,iym,iz ))**2  &
              + (vy(ixm,iym,izm) - vy(ixm,iym,iz))**2 + (vz(ixm,iym,izm) - vz(ixm,iym,iz))**2) 
          a12= ((vx(ix ,iym,izm) - vx(ix ,iym,iz ))**2  &
              + (vy(ix ,iym,izm) - vy(ix ,iym,iz))**2 + (vz(ix ,iym,izm) - vz(ix ,iym,iz))**2) 

  
          p_visc(ix,iy,iz) = MAX(p_visc(ix,iy,iz), - alpha1(ix,iy,iz)*SQRT(a1)) 
          p_visc(ix,iy,iz) = MAX(p_visc(ix,iy,iz), - alpha2(ix,iy,iz)*SQRT(a2)) 
          p_visc(ix,iy,iz) = MAX(p_visc(ix,iy,iz), - alpha3(ix,iy,iz)*SQRT(a9)) 

          visc_heat(ix,iy,iz) = &
              - 0.25_num * dyb(iy) * dzb(iz) * alpha1(ix ,iy ,iz ) * a1  &
              - 0.25_num * dxb(ix) * dzb(iz) * alpha2(ix ,iy ,iz ) * a2  &
              - 0.25_num * dyb(iy) * dzb(iz) * alpha1(ix ,iyp,iz ) * a3  &
              - 0.25_num * dxb(ix) * dzb(iz) * alpha2(ixm,iy ,iz ) * a4  &
              - 0.25_num * dyb(iy) * dzb(iz) * alpha1(ix ,iy ,izp) * a5  &
              - 0.25_num * dxb(ix) * dzb(iz) * alpha2(ix ,iy ,izp) * a6  &
              - 0.25_num * dyb(iy) * dzb(iz) * alpha1(ix ,iyp,izp) * a7  &
              - 0.25_num * dxb(ix) * dzb(iz) * alpha2(ixm,iy ,izp) * a8  &
              - 0.25_num * dxb(ix) * dyb(iy) * alpha3(ix ,iy ,iz ) * a9  &
              - 0.25_num * dxb(ix) * dyb(iy) * alpha3(ixm,iy ,iz ) * a10 &
              - 0.25_num * dxb(ix) * dyb(iy) * alpha3(ixm,iym,iz ) * a11 &
              - 0.25_num * dxb(ix) * dyb(iy) * alpha3(ix ,iym,iz ) * a12
              
          visc_heat(ix,iy,iz) = visc_heat(ix,iy,iz) / cv(ix,iy,iz)
        END DO
      END DO
    END DO

    fx_visc = 0.0_num
    fy_visc = 0.0_num
    fz_visc = 0.0_num

    DO iz = 0, nz
      izm = iz - 1
      izp = iz + 1
      DO iy = 0, ny
        iym = iy - 1
        iyp = iy + 1
        DO ix = 0, nx
          ixm = ix - 1
          ixp = ix + 1
  
          a1 = alpha1(ix ,iyp,izp) * dyc(iy) * dzc(iz)
          a2 = alpha1(ixp,iyp,izp) * dyc(iy) * dzc(iz)
          a3 = alpha2(ix ,iy ,izp) * dxc(ix) * dzc(iz)
          a4 = alpha2(ix ,iyp,izp) * dxc(ix) * dzc(iz)         
          a5 = alpha3(ix ,iy , iz) * dxc(ix) * dyc(iy)
          a6 = alpha3(ix ,iy ,izp) * dxc(ix) * dyc(iy)


          fx_visc(ix,iy,iz) = &
                          +(a1 * (vx(ix,iy,iz) - vx(ixm,iy , iz)) &
                          + a2 * (vx(ix,iy,iz) - vx(ixp,iy , iz)) &
                          + a3 * (vx(ix,iy,iz) - vx(ix ,iym, iz)) &
                          + a4 * (vx(ix,iy,iz) - vx(ix ,iyp, iz)) &
                          + a5 * (vx(ix,iy,iz) - vx(ix , iy,izm)) &
                          + a6 * (vx(ix,iy,iz) - vx(ix , iy,izp)) ) / cv_v(ix,iy,iz)
  
          fy_visc(ix,iy,iz) = &
                          +(a1 * (vy(ix,iy,iz) - vy(ixm,iy , iz)) &
                          + a2 * (vy(ix,iy,iz) - vy(ixp,iy , iz)) &
                          + a3 * (vy(ix,iy,iz) - vy(ix ,iym, iz)) &
                          + a4 * (vy(ix,iy,iz) - vy(ix ,iyp, iz)) &
                          + a5 * (vy(ix,iy,iz) - vy(ix , iy,izm)) &
                          + a6 * (vy(ix,iy,iz) - vy(ix , iy,izp)) ) / cv_v(ix,iy,iz)
  
          fz_visc(ix,iy,iz) = &
                          +(a1 * (vz(ix,iy,iz) - vz(ixm,iy , iz)) &
                          + a2 * (vz(ix,iy,iz) - vz(ixp,iy , iz)) &
                          + a3 * (vz(ix,iy,iz) - vz(ix ,iym, iz)) &
                          + a4 * (vz(ix,iy,iz) - vz(ix ,iyp, iz)) &
                          + a5 * (vz(ix,iy,iz) - vz(ix , iy,izm)) &
                          + a6 * (vz(ix,iy,iz) - vz(ix , iy,izp)) ) / cv_v(ix,iy,iz)
  
        END DO
      END DO
    END DO

    DEALLOCATE(cs, cs_v)

  CONTAINS

    DOUBLE PRECISION FUNCTION edge_viscosity()

      ! Actually returns q_k_bar = q_kur*(1-psi) / abs(dv)
      ! Other symbols follow notation in Caramana

      REAL(num) :: dvx, dvy, dvz, dv, dv2
#ifdef SHOCKLIMITER
      REAL(num) :: dvxm, dvxp, dvym, dvyp, dvzm, dvzp
      REAL(num) :: rl, rr
#endif
      REAL(num) :: psi, rho_edge, cs_edge, q_k_bar

#ifdef SHOCKCOMPRESSION 
      ! Turn off shock viscoity if cell edge expanding
      dvdots = MIN(0.0_num, dvdots)
#else
      ! Allow shock viscoity on expanding edge
      dvdots = -ABS(dvdots)
#endif

      rho_edge = 2.0_num * rho_v(i1,j1,k1) * rho_v(i2,j2,k2) &
          / (rho_v(i1,j1,k1) + rho_v(i2,j2,k2))
      cs_edge = MIN(cs_v(i1,j1,k1), cs_v(i2,j2,k2))

      dvx = vx(i1,j1,k1) - vx(i2,j2,k2)
      dvy = vy(i1,j1,k1) - vy(i2,j2,k2)
      dvz = vz(i1,j1,k1) - vz(i2,j2,k2)
      dv2 = dvx**2 + dvy**2 + dvz**2
      dv = SQRT(dv2)
      psi = 0.0_num
      IF (dv * dt / dx < 1.e-14_num) THEN
        dvdots = 0.0_num
      ELSE
        dvdots = dvdots / dv
      END IF

#ifdef SHOCKLIMITER
      dvxm = vx(i0,j0,k0) - vx(i1,j1,k1)
      dvxp = vx(i2,j2,k2) - vx(i3,j3,k3)
      dvym = vy(i0,j0,k0) - vy(i1,j1,k1)
      dvyp = vy(i2,j2,k2) - vy(i3,j3,k3)
      dvzm = vz(i0,j0,k0) - vz(i1,j1,k1)
      dvzp = vz(i2,j2,k2) - vz(i3,j3,k3)
      IF (dv * dt / dx < 1.e-14_num) THEN
        rl = 1.0_num
        rr = 1.0_num
      ELSE
        rl = (dvxp * dvx + dvyp * dvy + dvzp * dvz) * dx / (dxp * dv2)
        rr = (dvxm * dvx + dvym * dvy + dvzm * dvz) * dx / (dxm * dv2)
      END IF
      psi = MIN(0.5_num * (rr + rl), 2.0_num * rl, 2.0_num * rr, 1.0_num)
      psi = MAX(0.0_num, psi)
#endif

      ! Find q_kur / abs(dv)
      q_k_bar = rho_edge &
          * (visc2_norm * dv + SQRT(visc2_norm**2 * dv2 + (visc1 * cs_edge)**2))

      edge_viscosity = q_k_bar * (1.0_num - psi) * dvdots

    END FUNCTION edge_viscosity

  END SUBROUTINE shock_viscosity



  SUBROUTINE shock_heating

    REAL(num) :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12

    visc_heat = 0.0_num

    DO iz = 0, nz + 1 
      izm = iz - 1
      izp = iz + 1
      DO iy = 0, ny + 1 
        iym = iy - 1
        iyp = iy + 1
        DO ix = 0, nx + 1 
          ixm = ix - 1
          ixp = ix + 1
  
          a1 =  (vx(ixm,iym,izm) - vx(ix ,iym,izm))*(vx1(ixm,iym,izm) - vx1(ix ,iym,izm)) &
              + (vy(ixm,iym,izm) - vy(ix ,iym,izm))*(vy1(ixm,iym,izm) - vy1(ix ,iym,izm)) &
              + (vz(ixm,iym,izm) - vz(ix ,iym,izm))*(vz1(ixm,iym,izm) - vz1(ix ,iym,izm)) 
          a2 =  (vx(ix ,iym,izm) - vx(ix ,iy ,izm))*(vx1(ix ,iym,izm) - vx1(ix ,iy ,izm)) &
              + (vy(ix ,iym,izm) - vy(ix ,iy ,izm))*(vy1(ix ,iym,izm) - vy1(ix ,iy ,izm)) &
              + (vz(ix ,iym,izm) - vz(ix ,iy ,izm))*(vz1(ix ,iym,izm) - vz1(ix ,iy ,izm))
          a3 =  (vx(ix ,iy ,izm) - vx(ixm,iy ,izm))*(vx1(ix ,iy ,izm) - vx1(ixm,iy ,izm)) &
              + (vy(ix ,iy ,izm) - vy(ixm,iy ,izm))*(vy1(ix ,iy ,izm) - vy1(ixm,iy ,izm)) &
              + (vz(ix ,iy ,izm) - vz(ixm,iy ,izm))*(vz1(ix ,iy ,izm) - vz1(ixm,iy ,izm))
          a4 =  (vx(ixm,iy ,izm) - vx(ixm,iym,izm))*(vx1(ixm,iy ,izm) - vx1(ixm,iym,izm)) &
              + (vy(ixm,iy ,izm) - vy(ixm,iym,izm))*(vy1(ixm,iy ,izm) - vy1(ixm,iym,izm)) &
              + (vz(ixm,iy ,izm) - vz(ixm,iym,izm))*(vz1(ixm,iy ,izm) - vz1(ixm,iym,izm))
  
          a5 =  (vx(ixm,iym, iz) - vx(ix ,iym, iz))*(vx1(ixm,iym, iz) - vx1(ix ,iym, iz)) &
              + (vy(ixm,iym, iz) - vy(ix ,iym, iz))*(vy1(ixm,iym, iz) - vy1(ix ,iym, iz)) &
              + (vz(ixm,iym, iz) - vz(ix ,iym, iz))*(vz1(ixm,iym, iz) - vz1(ix ,iym, iz)) 
          a6 =  (vx(ix ,iym, iz) - vx(ix ,iy , iz))*(vx1(ix ,iym, iz) - vx1(ix ,iy , iz)) &
              + (vy(ix ,iym, iz) - vy(ix ,iy , iz))*(vy1(ix ,iym, iz) - vy1(ix ,iy , iz)) &
              + (vz(ix ,iym, iz) - vz(ix ,iy , iz))*(vz1(ix ,iym, iz) - vz1(ix ,iy , iz))
          a7 =  (vx(ix ,iy , iz) - vx(ixm,iy , iz))*(vx1(ix ,iy , iz) - vx1(ixm,iy , iz)) &
              + (vy(ix ,iy , iz) - vy(ixm,iy , iz))*(vy1(ix ,iy , iz) - vy1(ixm,iy , iz)) &
              + (vz(ix ,iy , iz) - vz(ixm,iy , iz))*(vz1(ix ,iy , iz) - vz1(ixm,iy , iz))
          a8 =  (vx(ixm,iy , iz) - vx(ixm,iym, iz))*(vx1(ixm,iy , iz) - vx1(ixm,iym, iz)) &
              + (vy(ixm,iy , iz) - vy(ixm,iym, iz))*(vy1(ixm,iy , iz) - vy1(ixm,iym, iz)) &
              + (vz(ixm,iy , iz) - vz(ixm,iym, iz))*(vz1(ixm,iy , iz) - vz1(ixm,iym, iz))
  
          a9 =  (vx(ix ,iy ,izm) - vx(ix ,iy , iz))*(vx1(ix ,iy ,izm) - vx1(ix ,iy , iz)) &
              + (vy(ix ,iy ,izm) - vy(ix ,iy , iz))*(vy1(ix ,iy ,izm) - vy1(ix ,iy , iz)) &
              + (vz(ix ,iy ,izm) - vz(ix ,iy , iz))*(vz1(ix ,iy ,izm) - vz1(ix ,iy , iz)) 
          a10=  (vx(ixm,iy ,izm) - vx(ixm,iy , iz))*(vx1(ixm,iy ,izm) - vx1(ixm,iy , iz)) &
              + (vy(ixm,iy ,izm) - vy(ixm,iy , iz))*(vy1(ixm,iy ,izm) - vy1(ixm,iy , iz)) &
              + (vz(ixm,iy ,izm) - vz(ixm,iy , iz))*(vz1(ixm,iy ,izm) - vz1(ixm,iy , iz)) 
          a11=  (vx(ixm,iym,izm) - vx(ixm,iym, iz))*(vx1(ixm,iym,izm) - vx1(ixm,iym, iz)) &
              + (vy(ixm,iym,izm) - vy(ixm,iym, iz))*(vy1(ixm,iym,izm) - vy1(ixm,iym, iz)) &
              + (vz(ixm,iym,izm) - vz(ixm,iym, iz))*(vz1(ixm,iym,izm) - vz1(ixm,iym, iz))
          a12=  (vx(ix ,iym,izm) - vx(ix ,iym, iz))*(vx1(ix ,iym,izm) - vx1(ix ,iym, iz)) &
              + (vy(ix ,iym,izm) - vy(ix ,iym, iz))*(vy1(ix ,iym,izm) - vy1(ix ,iym, iz)) &
              + (vz(ix ,iym,izm) - vz(ix ,iym, iz))*(vz1(ix ,iym,izm) - vz1(ix ,iym, iz))
  
          visc_heat(ix,iy,iz) = &
              - 0.25_num * dyb(iy) * dzb(iz) * alpha1(ix ,iy ,iz ) * a1  &
              - 0.25_num * dxb(ix) * dzb(iz) * alpha2(ix ,iy ,iz ) * a2  &
              - 0.25_num * dyb(iy) * dzb(iz) * alpha1(ix ,iyp,iz ) * a3  &
              - 0.25_num * dxb(ix) * dzb(iz) * alpha2(ixm,iy ,iz ) * a4  &
              - 0.25_num * dyb(iy) * dzb(iz) * alpha1(ix ,iy ,izp) * a5  &
              - 0.25_num * dxb(ix) * dzb(iz) * alpha2(ix ,iy ,izp) * a6  &
              - 0.25_num * dyb(iy) * dzb(iz) * alpha1(ix ,iyp,izp) * a7  &
              - 0.25_num * dxb(ix) * dzb(iz) * alpha2(ixm,iy ,izp) * a8  &
              - 0.25_num * dxb(ix) * dyb(iy) * alpha3(ix ,iy ,iz ) * a9  &
              - 0.25_num * dxb(ix) * dyb(iy) * alpha3(ixm,iy ,iz ) * a10 &
              - 0.25_num * dxb(ix) * dyb(iy) * alpha3(ixm,iym,iz ) * a11 &
              - 0.25_num * dxb(ix) * dyb(iy) * alpha3(ix ,iym,iz ) * a12 
  
          visc_heat(ix,iy,iz) = visc_heat(ix,iy,iz) / cv(ix,iy,iz)
        END DO
      END DO
    END DO

    visc_heat = MAX(visc_heat, 0.0_num)

  END SUBROUTINE shock_heating


  !****************************************************************************
  ! This subroutine calculates the viscous effects and updates the
  ! magnetic field
  !****************************************************************************

  SUBROUTINE b_field_and_cv1_update

    REAL(num) :: vxb, vxbm, vyb, vybm, vzb, vzbm, dvxdx, dvydy, dvzdz, dv
#ifdef CAUCHY
    REAL(num) :: dvydx, dvzdx, dvxdy, dvzdy, dvxdz, dvydz
#endif

    DO iz = -1, nz + 2
      izm = iz - 1
      DO iy = -1, ny + 2
        iym = iy - 1
        DO ix = -1, nx + 2
          ixm = ix - 1

          ! vx at Bx(i,j,k)
          vxb  = (vx(ix ,iy ,iz ) + vx(ix ,iym,iz ) &
              +   vx(ix ,iy ,izm) + vx(ix ,iym,izm)) * 0.25_num
          ! vx at Bx(i-1,j,k)
          vxbm = (vx(ixm,iy ,iz ) + vx(ixm,iym,iz ) &
              +   vx(ixm,iy ,izm) + vx(ixm,iym,izm)) * 0.25_num
          ! vy at By(i,j,k)
          vyb  = (vy(ix ,iy ,iz ) + vy(ixm,iy ,iz ) &
              +   vy(ix ,iy ,izm) + vy(ixm,iy ,izm)) * 0.25_num
          ! vy at By(i,j-1,k)
          vybm = (vy(ix ,iym,iz ) + vy(ixm,iym,iz ) &
              +   vy(ix ,iym,izm) + vy(ixm,iym,izm)) * 0.25_num
          ! vz at Bz(i,j,k)
          vzb  = (vz(ix ,iy ,iz ) + vz(ixm,iy ,iz ) &
              +   vz(ix ,iym,iz ) + vz(ixm,iym,iz )) * 0.25_num
          ! vz at Bz(i,j,k-1)
          vzbm = (vz(ix ,iy ,izm) + vz(ixm,iy ,izm) &
              +   vz(ix ,iym,izm) + vz(ixm,iym,izm)) * 0.25_num

          dvxdx = (vxb - vxbm) / dxb(ix)
          dvydy = (vyb - vybm) / dyb(iy)
          dvzdz = (vzb - vzbm) / dzb(iz)

          dv = (dvxdx + dvydy + dvzdz) * dt2
          cv1(ix,iy,iz) = cv(ix,iy,iz) * (1.0_num + dv)

#ifdef CAUCHY
          ! vx at By(i,j,k)
          vxb  = (vx(ix ,iy ,iz ) + vx(ixm,iy ,iz ) &
              +   vx(ix ,iy ,izm) + vx(ixm,iy ,izm)) * 0.25_num
          ! vx at By(i,j-1,k)
          vxbm = (vx(ix ,iym,iz ) + vx(ixm,iym,iz ) &
              +   vx(ix ,iym,izm) + vx(ixm,iym,izm)) * 0.25_num
          ! vy at Bx(i,j,k)
          vyb  = (vy(ix ,iy ,iz ) + vy(ix ,iym,iz ) &
              +   vy(ix ,iy ,izm) + vy(ix ,iym,izm)) * 0.25_num
          ! vy at Bx(i-1,j,k)
          vybm = (vy(ixm,iy ,iz ) + vy(ixm,iym,iz ) &
              +   vy(ixm,iy ,izm) + vy(ixm,iym,izm)) * 0.25_num

          dvxdy = (vxb - vxbm) / dyb(iy)
          dvydx = (vyb - vybm) / dxb(ix)


          ! vx at Bz(i,j,k)
          vxb  = (vx(ix ,iy ,iz ) + vx(ixm,iy ,iz ) &
              +   vx(ix ,iym,iz ) + vx(ixm,iym,iz )) * 0.25_num
          ! vx at Bz(i,j,k-1)
          vxbm = (vx(ix ,iy ,izm) + vx(ixm,iy ,izm) &
              +   vx(ix ,iym,izm) + vx(ixm,iym,izm)) * 0.25_num
          ! vz at Bx(i,j,k)
          vzb  = (vz(ix ,iy ,iz ) + vz(ix ,iym,iz ) &
              +   vz(ix ,iy ,izm) + vz(ix ,iym,izm)) * 0.25_num
          ! vz at Bx(i-1,j,k)
          vzbm = (vz(ixm,iy ,iz ) + vz(ixm,iym,iz ) &
              +   vz(ixm,iy ,izm) + vz(ixm,iym,izm)) * 0.25_num

          dvxdz = (vxb - vxbm) / dzb(iz)
          dvzdx = (vzb - vzbm) / dxb(ix)

          ! vy at Bz(i,j,k)
          vyb  = (vy(ix ,iy ,iz ) + vy(ixm,iy ,iz ) &
              +   vy(ix ,iym,iz ) + vy(ixm,iym,iz )) * 0.25_num
          ! vy at Bz(i,j,k-1)
          vybm = (vy(ix ,iy ,izm) + vy(ixm,iy ,izm) &
              +   vy(ix ,iym,izm) + vy(ixm,iym,izm)) * 0.25_num
          ! vz at By(i,j,k)
          vzb  = (vz(ix ,iy ,iz ) + vz(ixm,iy ,iz ) &
              +   vz(ix ,iy ,izm) + vz(ixm,iy ,izm)) * 0.25_num
          ! vz at By(i,j-1,k)
          vzbm = (vz(ix ,iym,iz ) + vz(ixm,iym,iz ) &
              +   vz(ix ,iym,izm) + vz(ixm,iym,izm)) * 0.25_num

          dvydz = (vyb - vybm) / dzb(iz)
          dvzdy = (vzb - vzbm) / dyb(iy)

          w3 =  bx1(ix,iy,iz) * dvxdx + by1(ix,iy,iz) * dvxdy &
              + bz1(ix,iy,iz) * dvxdz
          w4 =  bx1(ix,iy,iz) * dvydx + by1(ix,iy,iz) * dvydy &
              + bz1(ix,iy,iz) * dvydz
          w5 =  bx1(ix,iy,iz) * dvzdx + by1(ix,iy,iz) * dvzdy &
              + bz1(ix,iy,iz) * dvzdz

          bx1(ix,iy,iz) = (bx1(ix,iy,iz) + w3 * dt2) / (1.0_num + dv)
          by1(ix,iy,iz) = (by1(ix,iy,iz) + w4 * dt2) / (1.0_num + dv)
          bz1(ix,iy,iz) = (bz1(ix,iy,iz) + w5 * dt2) / (1.0_num + dv)
#endif
        END DO
      END DO
    END DO

#ifndef CAUCHY
    vx1(:,:,:) = vx(:,:,:)
    vy1(:,:,:) = vy(:,:,:)
    vz1(:,:,:) = vz(:,:,:)
    
    predictor_step = .TRUE.
    dt = 0.5_num * dt
    CALL eulerian_remap(step)
    dt = 2.0_num * dt
    predictor_step = .FALSE. 
#endif

  END SUBROUTINE b_field_and_cv1_update

 

  !****************************************************************************
  ! Sets CFL limited step
  !****************************************************************************

  SUBROUTINE set_dt

    ! Assumes all variables are defined at the same point. Be careful with
    ! setting 'dt_multiplier' if you expect massive changes across cells.

    REAL(num) :: cs2, c_visc2, rho0, length
    REAL(num) :: dxlocal, dt_local, dtr_local, dt1, ss_limit
    REAL(num) :: dt_locals(2), dt_min(2)
    LOGICAL, SAVE :: first = .TRUE.

    IF (first) THEN
      first = .FALSE.
      IF (restart) THEN
        dt = dt_from_restart
        RETURN
      END IF
    END IF

    dt_local = largest_number
    dtr_local = largest_number

    gamma_boris(:,:,:) = 1.0_num

    DO iz = 0, nz
      izm = iz - 1
      DO iy = 0, ny
        iym = iy - 1
        DO ix = 0, nx
          ixm = ix - 1

          ! Fix dt for Lagrangian step
          w1 = bx(ix,iy,iz)**2 + by(ix,iy,iz)**2 + bz(ix,iy,iz)**2
          ! Sound speed squared
          rho0 = MAX(rho(ix, iy, iz), none_zero)
          w2 = w1 / rho0
          IF (boris .AND. (w2 .GE. va_max2)) THEN
            gamma_boris(ix,iy,iz) = 1.0_num / (1.0_num + w2 / va_max2)
          END IF
          cs2 = (gamma * pressure(ix,iy,iz) + gamma_boris(ix,iy,iz) * w1) / rho0

          !effective speed from viscous pressure
          c_visc2 = p_visc(ix,iy,iz) / rho0

          ! length estimates - could be smoother as in DYNA3D
          length = MIN(dxb(ix), dyb(iy), dzb(iz))

          ! Find ideal MHD CFL limit for Lagrangian step
          dt1 = length / (SQRT(c_visc2) + SQRT(cs2 + c_visc2))
          dt_local = MIN(dt_local, dt1)

          ! Note resistive limits assumes uniform resistivity hence cautious
          ! factor 0.2
          dxlocal = 1.0_num / (1.0_num / dxb(ix)**2 &
              + 1.0_num / dyb(iy)**2 + 1.0_num / dzb(iz)**2)

          dt1 = largest_number
          IF (cowling_resistivity) THEN
            dt1 = 0.2_num * dxlocal &
                / MAX(MAX(eta(ix,iy,iz), eta_perp(ix,iy,iz)), none_zero)
          ELSEIF (resistive_mhd) THEN
            dt1 = 0.2_num * dxlocal / MAX(eta(ix,iy,iz), none_zero)
          END IF

          ! Adjust to accomodate resistive effects
          dtr_local = MIN(dtr_local, dt1)
        END DO
      END DO
    END DO

    dt_locals(1) = dt_local
    dt_locals(2) = dtr_local

    CALL MPI_ALLREDUCE(dt_locals, dt_min, 2, mpireal, MPI_MIN, comm, errcode)

    dt  = dt_multiplier * dt_min(1)
    dtr = dt_multiplier * dt_min(2)

    IF (conduction) THEN
      CALL calc_s_stages(.TRUE.)
      ss_limit = 60
      IF (n_s_stages >= ss_limit) THEN
        dt  = dt  * REAL(2 * ss_limit**2   - 9, num) &
                  / REAL(2 * n_s_stages**2 - 9, num)
        dtr = dtr * REAL(2 * ss_limit**2   - 9, num) &
                  / REAL(2 * n_s_stages**2 - 9, num)
      END IF
    END IF

    time = time + dt

  END SUBROUTINE set_dt



  !****************************************************************************
  ! Calculate the spatial profile of the resistivity at the current timestep
  ! Note that this is a core routine so it works in normalised units
  ! This includes lengths etc.
  !****************************************************************************

  SUBROUTINE eta_calc

    REAL(num) :: jx, jy, jz, jxp, jyp, jzp
    INTEGER :: ixp, iyp, izp

    IF (resistive_mhd) THEN
      DO iz = -1, nz + 1
        izp = iz + 1
        DO iy = -1, ny + 1
          iyp = iy + 1
          DO ix = -1, nx + 1
            ixp = ix + 1

            ! jx at Ex(i,j,k)
            jx  = (bz(ix ,iyp,iz ) - bz(ix ,iy ,iz )) / dyc(iy) &
                - (by(ix ,iy ,izp) - by(ix ,iy ,iz )) / dzc(iz)

            ! jx at Ex(i+1,j,k)
            jxp = (bz(ixp,iyp,iz ) - bz(ixp,iy ,iz )) / dyc(iy) &
                - (by(ixp,iy ,izp) - by(ixp,iy ,iz )) / dzc(iz)

            ! jy at Ey(i,j,k)
            jy  = (bx(ix ,iy ,izp) - bx(ix ,iy ,iz )) / dzc(iz) &
                - (bz(ixp,iy ,iz ) - bz(ix ,iy ,iz )) / dxc(ix)

            ! jy at Ey(i,j+1,k)
            jyp = (bx(ix ,iyp,izp) - bx(ix ,iyp,iz )) / dzc(iz) &
                - (bz(ixp,iyp,iz ) - bz(ix ,iyp,iz )) / dxc(ix)

            ! jz at Ez(i,j,k)
            jz  = (by(ixp,iy ,iz ) - by(ix ,iy ,iz )) / dxc(ix) &
                - (bx(ix ,iyp,iz ) - bx(ix ,iy ,iz )) / dyc(iy)

            ! jz at Ez(i,j,k+1)
            jzp = (by(ixp,iy ,izp) - by(ix ,iy ,izp)) / dxc(ix) &
                - (bx(ix ,iyp,izp) - bx(ix ,iy ,izp)) / dyc(iy)

            ! Current at V
            jx = (jx + jxp) * 0.5_num
            jy = (jy + jyp) * 0.5_num
            jz = (jz + jzp) * 0.5_num

            IF (SQRT(jx**2 + jy**2 + jz**2) > j_max) THEN
              eta(ix,iy,iz) = eta_background + eta0
            ELSE
              eta(ix,iy,iz) = eta_background
            END IF
          END DO
        END DO
      END DO
    ELSE
      eta(:,:,:) = 0.0_num
    END IF

  END SUBROUTINE eta_calc



  !****************************************************************************
  ! Calculate the effect of resistivity on the magnetic field and Ohmic heating
  ! Use the subroutine rkstep
  !****************************************************************************

  SUBROUTINE resistive_effects

    REAL(num) :: jx1, jx2, jy1, jy2, jz1, jz2, local_heating
#ifdef FOURTHORDER
    REAL(num) :: dt6, half_dt
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: k1x, k2x, k3x, k4x
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: k1y, k2y, k3y, k4y
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: k1z, k2z, k3z, k4z
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: c1, c2, c3, c4

    ALLOCATE(k1x(0:nx,0:ny,0:nz), k2x(0:nx,0:ny,0:nz))
    ALLOCATE(k3x(0:nx,0:ny,0:nz), k4x(0:nx,0:ny,0:nz))
    ALLOCATE(k1y(0:nx,0:ny,0:nz), k2y(0:nx,0:ny,0:nz))
    ALLOCATE(k3y(0:nx,0:ny,0:nz), k4y(0:nx,0:ny,0:nz))
    ALLOCATE(k1z(0:nx,0:ny,0:nz), k2z(0:nx,0:ny,0:nz))
    ALLOCATE(k3z(0:nx,0:ny,0:nz), k4z(0:nx,0:ny,0:nz))
    ALLOCATE( c1(0:nx,0:ny,0:nz),  c2(0:nx,0:ny,0:nz))
    ALLOCATE( c3(0:nx,0:ny,0:nz),  c4(0:nx,0:ny,0:nz))
#endif

    bx1(:,:,:) = bx(-1:nx+2,-1:ny+2,-1:nz+2)
    by1(:,:,:) = by(-1:nx+2,-1:ny+2,-1:nz+2)
    bz1(:,:,:) = bz(-1:nx+2,-1:ny+2,-1:nz+2)

    ! Step 1
    CALL rkstep

    ! Default is first order in time
#ifndef FOURTHORDER
    CALL bstep(flux_x, flux_y, flux_z, dt)

    DO iz = 1, nz
      izm = iz - 1
      DO iy = 1, ny
        iym = iy - 1
        DO ix = 1, nx
          ixm = ix - 1
          local_heating = energy(ix,iy,iz) &
              + (curlb(ix ,iy ,iz ) + curlb(ixm,iy ,iz )  &
              +  curlb(ix ,iym,iz ) + curlb(ixm,iym,iz )  &
              +  curlb(ix ,iy ,izm) + curlb(ixm,iy ,izm)  &
              +  curlb(ix ,iym,izm) + curlb(ixm,iym,izm)) &
              * dt / (8.0_num * rho(ix,iy,iz))

          ohmic_dep(ix,iy,iz) = ohmic_dep(ix,iy,iz) + local_heating
          
          energy(ix,iy,iz) = energy(ix,iy,iz) + local_heating 
        END DO
      END DO
    END DO

    IF (cooling_term) THEN 
      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            local_heating = (curlb(ix ,iy ,iz) + curlb(ixm,iy ,iz)  &
                + curlb(ix ,iym,iz) + curlb(ixm,iym,iz)) &
                * dt / (4.0_num * rho(ix,iy,iz))
            cool_term_b(ix,iy,iz) = alpha_av * local_heating  &
              + (1.0_num - alpha_av) * cool_term_b(ix,iy,iz)
  
            energy(ix,iy,iz) = energy(ix,iy,iz) - cool_term_b(ix,iy,iz)
          END DO
        END DO
      END DO
    END IF
    energy = MAX(energy, 0.0_num)

    CALL energy_bcs

    DO iz = 0, nz
      DO iy = 0, ny
        DO ix = 0, nx
          w1 = dt * dxc(ix) * dyc(iy) * dzc(iz) * curlb(ix,iy,iz)
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
#else
    ! If complier flag set then use 4th order Runge-Kutta
    half_dt = dt * 0.5_num
    dt6 = dt * sixth

    k1x(:,:,:) = flux_x(:,:,:)
    k1y(:,:,:) = flux_y(:,:,:)
    k1z(:,:,:) = flux_z(:,:,:)
    c1(:,:,:) = curlb(:,:,:)

    ! Step 2
    CALL bstep(k1x, k1y, k1z, half_dt)

    CALL rkstep

    k2x(:,:,:) = flux_x(:,:,:)
    k2y(:,:,:) = flux_y(:,:,:)
    k2z(:,:,:) = flux_z(:,:,:)
    c2(:,:,:) = curlb(:,:,:)

    ! Step 3
    CALL bstep(k2x, k2y, k2z, half_dt)

    CALL rkstep

    k3x(:,:,:) = flux_x(:,:,:)
    k3y(:,:,:) = flux_y(:,:,:)
    k3z(:,:,:) = flux_z(:,:,:)
    c3(:,:,:)= curlb(:,:,:)

    ! Step 4
    CALL bstep(k3x, k3y, k3z, dt)

    CALL rkstep

    k4x(:,:,:) = flux_x(:,:,:)
    k4y(:,:,:) = flux_y(:,:,:)
    k4z(:,:,:) = flux_z(:,:,:)
    c4(:,:,:) = curlb(:,:,:)

    ! Full update
    k1x(:,:,:) = k1x(:,:,:) + 2.0_num * k2x(:,:,:) + 2.0_num * k3x(:,:,:) + k4x(:,:,:)
    k1y(:,:,:) = k1y(:,:,:) + 2.0_num * k2y(:,:,:) + 2.0_num * k3y(:,:,:) + k4y(:,:,:)
    k1z(:,:,:) = k1z(:,:,:) + 2.0_num * k2z(:,:,:) + 2.0_num * k3z(:,:,:) + k4z(:,:,:)
    c1(:,:,:) = c1(:,:,:) + 2.0_num * c2(:,:,:) + 2.0_num * c3(:,:,:) + c4(:,:,:)

    CALL bstep(k1x, k1y, k1z, dt6)

    DO iz = 1, nz
      izm = iz - 1
      DO iy = 1, ny
        iym = iy - 1
        DO ix = 1, nx
          ixm = ix - 1
          energy(ix,iy,iz) = energy(ix,iy,iz) &
              + (c1(ix ,iy ,iz ) + c1(ixm,iy ,iz )  &
              +  c1(ix ,iym,iz ) + c1(ixm,iym,iz )  &
              +  c1(ix ,iy ,izm) + c1(ixm,iy ,izm)  &
              +  c1(ix ,iym,izm) + c1(ixm,iym,izm)) &
              * dt6 / (8.0_num * rho(ix,iy,iz))
        END DO
      END DO
    END DO

    CALL energy_bcs

    DO iz = 0, nz
      DO iy = 0, ny
        DO ix = 0, nx
          w1 = dt6 * dxc(ix) * dyc(iy) * dzc(iz) * c1(ix,iy,iz)
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

          jx1 = (bz(ix ,iyp,iz ) - bz(ix ,iy ,iz )) / dyc(iy) &
              - (by(ix ,iy ,izp) - by(ix ,iy ,iz )) / dzc(iz)
          jx2 = (bz(ixp,iyp,iz ) - bz(ixp,iy ,iz )) / dyc(iy) &
              - (by(ixp,iy ,izp) - by(ixp,iy ,iz )) / dzc(iz)
          jy1 = (bx(ix ,iy ,izp) - bx(ix ,iy ,iz )) / dzc(iz) &
              - (bz(ixp,iy ,iz ) - bz(ix ,iy ,iz )) / dxc(ix)
          jy2 = (bx(ix ,iyp,izp) - bx(ix ,iyp,iz )) / dzc(iz) &
              - (bz(ixp,iyp,iz ) - bz(ix ,iyp,iz )) / dxc(ix)
          jz1 = (by(ixp,iy ,iz ) - by(ix ,iy ,iz )) / dxc(ix) &
              - (bx(ix ,iyp,iz ) - bx(ix ,iy ,iz )) / dyc(iy)
          jz2 = (by(ixp,iy ,izp) - by(ix ,iy ,izp)) / dxc(ix) &
              - (bx(ix ,iyp,izp) - bx(ix ,iy ,izp)) / dyc(iy)

          jx_r(ix,iy,iz) = (jx1 + jx2) * 0.5_num
          jy_r(ix,iy,iz) = (jy1 + jy2) * 0.5_num
          jz_r(ix,iy,iz) = (jz1 + jz2) * 0.5_num
        END DO
      END DO
    END DO

    ! Once more to get j_perp and j_par correct
    CALL rkstep

#ifdef FOURTHORDER
    DEALLOCATE(k1x, k2x, k3x, k4x, k1y, k2y, k3y, k4y, k1z, k2z, k3z, k4z)
    DEALLOCATE(c1, c2, c3, c4)
#endif

  END SUBROUTINE resistive_effects



  !****************************************************************************
  ! Calculates 'k' values from b[xyz]1 values
  !****************************************************************************

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

            jx1 = (bz(ix ,iyp,iz ) - bz(ix ,iy ,iz )) / dyc(iy) &
                - (by(ix ,iy ,izp) - by(ix ,iy ,iz )) / dzc(iz)
            jx2 = (bz(ixp,iyp,iz ) - bz(ixp,iy ,iz )) / dyc(iy) &
                - (by(ixp,iy ,izp) - by(ixp,iy ,iz )) / dzc(iz)
            jy1 = (bx(ix ,iy ,izp) - bx(ix ,iy ,iz )) / dzc(iz) &
                - (bz(ixp,iy ,iz ) - bz(ix ,iy ,iz )) / dxc(ix)
            jy2 = (bx(ix ,iyp,izp) - bx(ix ,iyp,iz )) / dzc(iz) &
                - (bz(ixp,iyp,iz ) - bz(ix ,iyp,iz )) / dxc(ix)
            jz1 = (by(ixp,iy ,iz ) - by(ix ,iy ,iz )) / dxc(ix) &
                - (bx(ix ,iyp,iz ) - bx(ix ,iy ,iz )) / dyc(iy)
            jz2 = (by(ixp,iy ,izp) - by(ix ,iy ,izp)) / dxc(ix) &
                - (bx(ix ,iyp,izp) - bx(ix ,iy ,izp)) / dyc(iy)

            jx = (jx1 + jx2) * 0.5_num
            jy = (jy1 + jy2) * 0.5_num
            jz = (jz1 + jz2) * 0.5_num

            flux_x(ix,iy,iz) = -jx * eta(ix,iy,iz) * dxc(ix) * 0.5_num
            flux_y(ix,iy,iz) = -jy * eta(ix,iy,iz) * dyc(iy) * 0.5_num
            flux_z(ix,iy,iz) = -jz * eta(ix,iy,iz) * dzc(iz) * 0.5_num
            ! This isn't really curlb. It's actually heat flux
            curlb(ix,iy,iz) = eta(ix,iy,iz) * (jx**2 + jy**2 + jz**2)
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

            jx1 = (bz(ix ,iyp,iz ) - bz(ix ,iy ,iz )) / dyc(iy) &
                - (by(ix ,iy ,izp) - by(ix ,iy ,iz )) / dzc(iz)
            jx2 = (bz(ixp,iyp,iz ) - bz(ixp,iy ,iz )) / dyc(iy) &
                - (by(ixp,iy ,izp) - by(ixp,iy ,iz )) / dzc(iz)
            jy1 = (bx(ix ,iy ,izp) - bx(ix ,iy ,iz )) / dzc(iz) &
                - (bz(ixp,iy ,iz ) - bz(ix ,iy ,iz )) / dxc(ix)
            jy2 = (bx(ix ,iyp,izp) - bx(ix ,iyp,iz )) / dzc(iz) &
                - (bz(ixp,iyp,iz ) - bz(ix ,iyp,iz )) / dxc(ix)
            jz1 = (by(ixp,iy ,iz ) - by(ix ,iy ,iz )) / dxc(ix) &
                - (bx(ix ,iyp,iz ) - bx(ix ,iy ,iz )) / dyc(iy)
            jz2 = (by(ixp,iy ,izp) - by(ix ,iy ,izp)) / dxc(ix) &
                - (bx(ix ,iyp,izp) - bx(ix ,iy ,izp)) / dyc(iy)

            jx = (jx1 + jx2) * 0.5_num
            jy = (jy1 + jy2) * 0.5_num
            jz = (jz1 + jz2) * 0.5_num

            ! B at vertices
            bxv = (bx(ix,iy ,iz ) + bx(ix ,iyp,iz ) &
                 + bx(ix,iy ,izp) + bx(ix ,iyp,izp)) * 0.25_num
            byv = (by(ix,iy ,iz ) + by(ixp,iy ,iz ) &
                 + by(ix,iy ,izp) + by(ixp,iy ,izp)) * 0.25_num
            bzv = (bz(ix,iy ,iz ) + bz(ixp,iy ,iz ) &
                 + bz(ix,iyp,iz ) + bz(ixp,iyp,iz )) * 0.25_num

            magn_b = bxv**2 + byv**2 + bzv**2

            ! Calculate parallel and perpendicular currents
            IF (magn_b > none_zero) THEN
              j_par_x = (jx * bxv + jy * byv + jz * bzv) * bxv / magn_b
              j_par_y = (jx * bxv + jy * byv + jz * bzv) * byv / magn_b
              j_par_z = (jx * bxv + jy * byv + jz * bzv) * bzv / magn_b
            ELSE
              ! If b = 0 then there is no parallel current
              j_par_x = 0.0_num
              j_par_y = 0.0_num
              j_par_z = 0.0_num
            END IF

            ! Calculate perpendicular current
            j_perp_x = jx - j_par_x
            j_perp_y = jy - j_par_y
            j_perp_z = jz - j_par_z

            magn_j_par  = SQRT(j_par_x**2 + j_par_y**2 + j_par_z**2)
            magn_j_perp = SQRT(j_perp_x**2 + j_perp_y**2 + j_perp_z**2)

            parallel_current(ix,iy,iz) = magn_j_par
            perp_current(ix,iy,iz) = magn_j_perp

            ! This isn't really curlb. It's actually heat flux
            curlb(ix,iy,iz) = eta(ix,iy,iz) * magn_j_par**2 &
                + (eta_perp(ix,iy,iz) + eta(ix,iy,iz)) * magn_j_perp**2

            flux_x(ix,iy,iz) = -((j_par_x + j_perp_x) * eta(ix,iy,iz) &
                + j_perp_x * eta_perp(ix,iy,iz)) * dxc(ix) * 0.5_num
            flux_y(ix,iy,iz) = -((j_par_y + j_perp_y) * eta(ix,iy,iz) &
                + j_perp_y * eta_perp(ix,iy,iz)) * dyc(iy) * 0.5_num
            flux_z(ix,iy,iz) = -((j_par_z + j_perp_z) * eta(ix,iy,iz) &
                + j_perp_z * eta_perp(ix,iy,iz)) * dzc(iz) * 0.5_num
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE rkstep



  SUBROUTINE bstep(kx, ky, kz, dt)

    REAL(num), DIMENSION(0:,0:,0:), INTENT(IN) :: kx, ky, kz
    REAL(num), INTENT(IN) :: dt
    REAL(num) :: area

    DO iz = 1, nz
      izm = iz - 1
      DO iy = 1, ny
        iym = iy - 1
        area = dyb(iy) * dzb(iz)
        DO ix = 0, nx
          bx(ix,iy,iz) = bx1(ix,iy,iz) &
              + (kz(ix,iy ,iz ) - kz(ix,iym,iz ) &
              +  kz(ix,iy ,izm) - kz(ix,iym,izm) &
              -  ky(ix,iy ,iz ) + ky(ix,iy ,izm) &
              -  ky(ix,iym,iz ) + ky(ix,iym,izm)) * dt / area
        END DO
      END DO
    END DO

    DO iz = 1, nz
      izm = iz - 1
      DO iy = 0, ny
        DO ix = 1, nx
          ixm = ix - 1
          area = dxb(ix) * dzb(iz)
          by(ix,iy,iz) = by1(ix,iy,iz) &
              + (kx(ix ,iy,iz ) - kx(ix ,iy,izm) &
              +  kx(ixm,iy,iz ) - kx(ixm,iy,izm) &
              -  kz(ix ,iy,iz ) + kz(ixm,iy,iz ) &
              -  kz(ix ,iy,izm) + kz(ixm,iy,izm)) * dt / area
        END DO
      END DO
    END DO

    DO iz = 0, nz
      DO iy = 1, ny
        iym = iy - 1
        DO ix = 1, nx
          ixm = ix - 1
          area = dxb(ix) * dyb(iy)
          bz(ix,iy,iz) = bz1(ix,iy,iz) &
              + (ky(ix ,iy ,iz) - ky(ixm,iy ,iz) &
              +  ky(ix ,iym,iz) - ky(ixm,iym,iz) &
              -  kx(ix ,iy ,iz) + kx(ix ,iym,iz) &
              -  kx(ixm,iy ,iz) + kx(ixm,iym,iz)) * dt / area
        END DO
      END DO
    END DO

    CALL bfield_bcs

  END SUBROUTINE bstep

END MODULE lagran
