MODULE lagran

!-----------------------------------------------------------------
!This subroutine performs the Lagrangian step
!Notes:
!There are !#DEC$ directives in this routine
!These override compilers vector analysis tools
!-----------------------------------------------------------------

  USE shared_data
  USE boundary
  USE diagnostics

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lagrangian_step, eta_calc

  ! only used inside lagran.f90
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: bx1, by1, bz1, qxy, qxz, qyz, pressure
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: qxx, qyy, qzz, visc_heat

CONTAINS


  SUBROUTINE lagrangian_step

    INTEGER :: substeps, subcycle
    REAL(num) :: actual_dt, dt_sub, cumulative_dt
    LOGICAL :: keepcycling

    ALLOCATE (bx1(0:nx+1,0:ny+1,0:nz+1),by1(0:nx+1,0:ny+1,0:nz+1),bz1(0:nx+1,0:ny+1,0:nz+1),&
         qxy(0:nx+1,0:ny+1,0:nz+1), qxz(0:nx+1,0:ny+1,0:nz+1),qyz(0:nx+1,0:ny+1,0:nz+1),  &
         visc_heat(0:nx+1,0:ny+1,0:nz+1),pressure(-1:nx+2,-1:ny+2,-1:nz+2), &
         qxx(0:nx+1,0:ny+1,0:nz+1), qyy(0:nx+1,0:ny+1,0:nz+1),qzz(0:nx+1,0:ny+1,0:nz+1))

    IF (resistiveMHD) THEN
       ! if subcycling isn't wanted set dt = dtr in set_dt, don't just
       ! set substeps to 1.
       keepcycling = .TRUE.
       cumulative_dt = 0.0_num
       substeps = 0
       DO WHILE (keepcycling)
          actual_dt = dt
          dt_sub = dt
          dt_sub = MIN(dt_sub,dtr)
          IF (cumulative_dt + dt_sub >= actual_dt) THEN
             dt_sub = MAX(actual_dt - cumulative_dt, 0.0_num)
	     keepcycling = .FALSE.
          ENDIF
          cumulative_dt = cumulative_dt + dt_sub
          dt = dt_sub
          CALL resistive_effects  !also calculates the current density
          CALL eta_calc  ! check the resisitvity hasn't changed
          CALL set_dt
          substeps = substeps + 1
       ENDDO
       dt = actual_dt
       peak_substeps = MAX(peak_substeps,substeps)
    ENDIF

    time = time + dt

    DO iz = 0, nz+1
       DO iy = 0, ny+1
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
          DO ix = 0, nx+1
             ixm = ix - 1
             iym = iy - 1
             izm = iz - 1
             bx1(ix,iy,iz) = (bx(ix,iy,iz) + bx(ixm,iy,iz)) / 2.0_num
             by1(ix,iy,iz) = (by(ix,iy,iz) + by(ix,iym,iz)) / 2.0_num
             bz1(ix,iy,iz) = (bz(ix,iy,iz) + bz(ix,iy,izm)) / 2.0_num
          END DO
       END DO
    END DO

    CALL predictor_corrector_step

    DEALLOCATE (bx1, by1, bz1, qxy, qxz, qyz, visc_heat, pressure, qxx, qyy, qzz)

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
    REAL(num) :: dv,UV3,NUV3

    dt2 = dt / 2.0_num
    CALL viscosity_and_b_update

    bx1 = bx1 * cv1(0:nx+1,0:ny+1,0:nz+1)
    by1 = by1 * cv1(0:nx+1,0:ny+1,0:nz+1)
    bz1 = bz1 * cv1(0:nx+1,0:ny+1,0:nz+1)

    IF (visc3 > 1.e-6_num) THEN
       UV3=1.0_num
       NUV3=0.0_num
    ELSE
       UV3=0.0_num
       NUV3=1.0_num
    ENDIF

    DO iz = 0, nz+1
       DO iy = 0, ny+1
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
          DO ix = 0, nx+1
             dv = cv1(ix,iy,iz) / cv(ix,iy,iz) - 1.0_num
             e1 = energy(ix,iy,iz) - pressure(ix,iy,iz) * dv/rho(ix,iy,iz)   !predictor energy

#ifndef Q_MONO
             e1 = e1  + visc_heat(ix,iy,iz)*dt2/rho(ix,iy,iz)   
#else
             e1 = e1  + visc_heat(ix,iy,iz)*dt2/rho(ix,iy,iz) * UV3
             e1 = e1 - p_visc(ix,iy,iz) * dv/rho(ix,iy,iz) 
#endif

             ! now define the predictor step pressures
             pressure(ix,iy,iz) = e1 &
                  * (gamma - 1.0_num) * rho(ix,iy,iz) * cv(ix,iy,iz) / cv1(ix,iy,iz)
#ifdef Q_MONO
             ! add shock viscosity
             pressure(ix,iy,iz) = pressure(ix,iy,iz) + p_visc(ix,iy,iz) 
#endif
          END DO
       END DO
    END DO

    DO iz = 0, nz
       DO iy = 0, ny
!DEC$ IVDEP
!DEC$ VECTOR ALWAYS
          DO ix = 0, nx
             ixp = ix + 1
             iyp = iy + 1 
             izp = iz + 1

             p = pressure(ix,iy,iz)
             pxp = pressure(ixp,iy,iz) 
             pyp = pressure(ix,iyp,iz) 
             pxpyp = pressure(ixp,iyp,iz) 
             pzp = pressure(ix,iy,izp) 
             pxpzp = pressure(ixp,iy,izp) 
             pypzp = pressure(ix,iyp,izp) 
             pxpypzp = pressure(ixp,iyp,izp)

             w1 = (p + pyp + pzp + pypzp) / 4.0_num    
             w2 = (pxp + pxpyp + pxpzp + pxpypzp) / 4.0_num       
             fx = - (w2 - w1) / dxc(ix)
             w1 = (p + pxp + pzp + pxpzp) / 4.0_num      
             w2 = (pyp + pxpyp + pypzp + pxpypzp) / 4.0_num        
             fy = - (w2 - w1) / dyc(iy)
             w1 = (p + pxp + pyp + pxpyp) / 4.0_num     
             w2 = (pzp + pxpzp + pypzp + pxpypzp) / 4.0_num      
             fz = - (w2 - w1) / dzc(iz)
#ifdef Q_MONO
             IF ( visc3 > 1.e-6_num) THEN
#endif
                !add diagonal components
                w1 = (qxx(ix,iy,iz)+qxx(ix,iyp,iz)+qxx(ix,iy,izp)+qxx(ix,iyp,izp)) / 4.0_num
                w2 = (qxx(ixp,iy,iz)+qxx(ixp,iyp,iz)+qxx(ixp,iy,izp)+qxx(ixp,iyp,izp)) / 4.0_num
                fx = fx + (w2 - w1) / dxc(ix)
                w1 = (qyy(ix,iy,iz)+qyy(ixp,iy,iz)+qyy(ix,iy,izp)+qyy(ixp,iy,izp)) / 4.0_num
                w2 = (qyy(ix,iyp,iz)+qyy(ixp,iyp,iz)+qyy(ix,iyp,izp)+qyy(ixp,iyp,izp)) / 4.0_num
                fy = fy + (w2 - w1) / dyc(iy)
                w1 = (qzz(ix,iy,iz)+qzz(ixp,iy,iz)+qzz(ix,iyp,iz)+qzz(ixp,iyp,iz)) / 4.0_num
                w2 = (qzz(ix,iy,izp)+qzz(ixp,iy,izp)+qzz(ix,iyp,izp)+qzz(ixp,iyp,izp)) / 4.0_num
                fz = fz + (w2 - w1) / dzc(iz)

                ! add shear viscosity forces
                w1 = (qxy(ix,iy,iz)+qxy(ixp,iy,iz)+qxy(ix,iy,izp)+qxy(ixp,iy,izp)) / 4.0_num
                w2 = (qxy(ix,iyp,iz)+qxy(ixp,iyp,iz)+qxy(ix,iyp,izp)+qxy(ixp,iyp,izp)) / 4.0_num
                fx = fx + (w2 - w1) / dyc(iy)
                w1 = (qxz(ix,iy,iz)+qxz(ixp,iy,iz)+qxz(ix,iyp,iz)+qxz(ixp,iyp,iz)) / 4.0_num
                w2 = (qxz(ix,iy,izp)+qxz(ixp,iy,izp)+qxz(ix,iyp,izp)+qxz(ixp,iyp,izp)) / 4.0_num
                fx = fx + (w2 - w1) / dzc(iz)

                w1 = (qxy(ix,iy,iz)+qxy(ix,iyp,iz)+qxy(ix,iy,izp)+qxy(ix,iyp,izp)) / 4.0_num
                w2 = (qxy(ixp,iy,iz)+qxy(ixp,iyp,iz)+qxy(ixp,iy,izp)+qxy(ixp,iyp,izp)) / 4.0_num
                fy = fy + (w2 - w1) / dxc(ix)
                w1 = (qyz(ix,iy,iz)+qyz(ixp,iy,iz)+qyz(ix,iyp,iz)+qyz(ixp,iyp,iz)) / 4.0_num
                w2 = (qyz(ix,iy,izp)+qyz(ixp,iy,izp)+qyz(ix,iyp,izp)+qyz(ixp,iyp,izp)) / 4.0_num
                fy = fy + (w2 - w1) / dzc(iz)

                w1 = (qxz(ix,iy,iz)+qxz(ix,iyp,iz)+qxz(ix,iy,izp)+qxz(ix,iyp,izp)) / 4.0_num
                w2 = (qxz(ixp,iy,iz)+qxz(ixp,iyp,iz)+qxz(ixp,iy,izp)+qxz(ixp,iyp,izp)) / 4.0_num
                fz = fz + (w2 - w1) / dxc(ix)
                w1 = (qyz(ix,iy,iz)+qyz(ixp,iy,iz)+qyz(ix,iy,izp)+qyz(ixp,iy,izp)) / 4.0_num
                w2 = (qyz(ix,iyp,iz)+qyz(ixp,iyp,iz)+qyz(ix,iyp,izp)+qyz(ixp,iyp,izp)) / 4.0_num
                fz = fz + (w2 - w1) / dyc(iy)
#ifdef Q_MONO
             END IF
#endif

             cvx = cv1(ix,iy,iz)+cv1(ix,iyp,iz)+cv1(ix,iy,izp)+cv1(ix,iyp,izp)
             cvxp = cv1(ixp,iy,iz)+cv1(ixp,iyp,iz)+cv1(ixp,iy,izp)+cv1(ixp,iyp,izp)
             cvy = cv1(ix,iy,iz)+cv1(ixp,iy,iz)+cv1(ix,iy,izp)+cv1(ixp,iy,izp) 
             cvyp = cv1(ix,iyp,iz)+cv1(ixp,iyp,iz)+cv1(ix,iyp,izp)+cv1(ixp,iyp,izp)
             cvz = cv1(ix,iy,iz)+cv1(ixp,iy,iz)+cv1(ix,iyp,iz)+cv1(ixp,iyp,iz)
             cvzp = cv1(ix,iy,izp)+cv1(ixp,iy,izp)+cv1(ix,iyp,izp)+cv1(ixp,iyp,izp)

             w1 = (bz1(ix,iy,iz)+bz1(ixp,iy,iz)+bz1(ix,iy,izp)+bz1(ixp,iy,izp)) / cvy 
             w2 = (bz1(ix,iyp,iz)+bz1(ixp,iyp,iz)+bz1(ix,iyp,izp)+bz1(ixp,iyp,izp)) / cvyp
             jx = (w2 - w1) / dyc(iy)
             w1 = (by1(ix,iy,iz)+by1(ixp,iy,iz)+by1(ix,iyp,iz)+by1(ixp,iyp,iz)) / cvz 
             w2 = (by1(ix,iy,izp)+by1(ixp,iy,izp)+by1(ix,iyp,izp)+by1(ixp,iyp,izp)) / cvzp 
             jx = jx - (w2 - w1) / dzc(iz)

             w1 = (bz1(ix,iy,iz)+bz1(ix,iyp,iz)+bz1(ix,iy,izp)+bz1(ix,iyp,izp)) / cvx 
             w2 = (bz1(ixp,iy,iz)+bz1(ixp,iyp,iz)+bz1(ixp,iy,izp)+bz1(ixp,iyp,izp)) / cvxp 
             jy = - (w2 - w1) / dxc(ix)
             w1 = (bx1(ix,iy,iz)+bx1(ixp,iy,iz)+bx1(ix,iyp,iz)+bx1(ixp,iyp,iz)) / cvz
             w2 = (bx1(ix,iy,izp)+bx1(ixp,iy,izp)+bx1(ix,iyp,izp)+bx1(ixp,iyp,izp)) / cvzp 
             jy = jy + (w2 - w1) / dzc(iz)

             w1 = (by1(ix,iy,iz)+by1(ix,iyp,iz)+by1(ix,iy,izp)+by1(ix,iyp,izp)) / cvx 
             w2 = (by1(ixp,iy,iz)+by1(ixp,iyp,iz)+by1(ixp,iy,izp)+by1(ixp,iyp,izp)) / cvxp 
             jz = (w2 - w1) / dxc(ix)
             w1 = (bx1(ix,iy,iz)+bx1(ixp,iy,iz)+bx1(ix,iy,izp)+bx1(ixp,iy,izp)) / cvy 
             w2 = (bx1(ix,iyp,iz)+bx1(ixp,iyp,iz)+bx1(ix,iyp,izp)+bx1(ixp,iyp,izp)) / cvyp
             jz = jz - (w2 - w1) / dyc(iy)

             bxv = (bx1(ix,iy,iz)+bx1(ixp,iy,iz)+bx1(ix,iy,izp)+bx1(ixp,iy,izp) &
                  + bx1(ix,iyp,iz)+bx1(ixp,iyp,iz)+bx1(ix,iyp,izp)+bx1(ixp,iyp,izp)) &
                  / (cvx + cvxp)
             byv = (by1(ix,iy,iz)+by1(ixp,iy,iz)+by1(ix,iy,izp)+by1(ixp,iy,izp) &
                  + by1(ix,iyp,iz)+by1(ixp,iyp,iz)+by1(ix,iyp,izp)+by1(ixp,iyp,izp)) &
                  / (cvx + cvxp)
             bzv = (bz1(ix,iy,iz)+bz1(ixp,iy,iz)+bz1(ix,iy,izp)+bz1(ixp,iy,izp) &
                  + bz1(ix,iyp,iz)+bz1(ixp,iyp,iz)+bz1(ix,iyp,izp)+bz1(ixp,iyp,izp)) &
                  / (cvx + cvxp)

             fx = fx + (jy*bzv - jz*byv)
             fy = fy - (jx*bzv - jz*bxv)
             fz = fz + (jx*byv - jy*bxv)

             rho_v = rho(ix,iy,iz)*cv(ix,iy,iz) + rho(ixp,iy,iz)*cv(ixp,iy,iz) &
                  + rho(ix,iyp,iz)*cv(ix,iyp,iz) +rho(ixp,iyp,iz)*cv(ixp,iyp,iz) &
                  + rho(ix,iy,izp)*cv(ix,iy,izp) + rho(ixp,iy,izp)*cv(ixp,iy,izp) &
                  + rho(ix,iyp,izp)*cv(ix,iyp,izp) +rho(ixp,iyp,izp)*cv(ixp,iyp,izp)
             rho_v = rho_v / (cv(ix,iy,iz) + cv(ixp,iy,iz) + cv(ix,iyp,iz) + cv(ixp,iyp,iz)  &
                  + cv(ix,iy,izp) + cv(ixp,iy,izp) + cv(ix,iyp,izp) + cv(ixp,iyp,izp)) 

             fz = fz - (rho_v*grav(iz))

             ! find half step velocity needed for remap
             vx1(ix,iy,iz) = vx(ix,iy,iz) + dt2 * fx / rho_v
             vy1(ix,iy,iz) = vy(ix,iy,iz) + dt2 * fy / rho_v 
             vz1(ix,iy,iz) = vz(ix,iy,iz) + dt2 * fz / rho_v

             ! velocity at the end of the Lagrangian step
             vx(ix,iy,iz) = vx(ix,iy,iz) + dt * fx / rho_v
             vy(ix,iy,iz) = vy(ix,iy,iz) + dt * fy / rho_v 
             vz(ix,iy,iz) = vz(ix,iy,iz) + dt * fz / rho_v

          END DO
       END DO
    END DO

    IF (any_open) CALL store_boundary_dv 
    CALL remap_v_bcs
#ifdef Q_MONO
    IF (visc3 > 1.e-6_num) THEN
#endif
       CALL visc_heating
#ifdef Q_MONO
    ENDIF
#endif
    ! finally correct density and energy to final values
    DO iz = 1, nz
       DO iy = 1, ny
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
          DO ix = 1, nx
             ixm = ix - 1
             iym = iy - 1
             izm = iz - 1

             vxb = (vx1(ix,iy,iz)+vx1(ix,iym,iz)+vx1(ix,iy,izm)+vx1(ix,iym,izm)) / 4.0_num     
             vxbm = (vx1(ixm,iy,iz)+vx1(ixm,iym,iz)+vx1(ixm,iy,izm)+vx1(ixm,iym,izm)) / 4.0_num  
             vyb = (vy1(ix,iy,iz)+vy1(ixm,iy,iz)+vy1(ix,iy,izm)+vy1(ixm,iy,izm)) / 4.0_num    
             vybm = (vy1(ix,iym,iz)+vy1(ixm,iym,iz)+vy1(ix,iym,izm)+vy1(ixm,iym,izm)) / 4.0_num  
             vzb = (vz1(ix,iy,iz)+vz1(ixm,iy,iz)+vz1(ix,iym,iz)+vz1(ixm,iym,iz)) / 4.0_num    
             vzbm = (vz1(ix,iy,izm)+vz1(ixm,iy,izm)+vz1(ix,iym,izm)+vz1(ixm,iym,izm)) / 4.0_num  
             dv = ((vxb - vxbm)/dxb(ix) + (vyb - vybm)/dyb(iy) &
                  + (vzb - vzbm)/dzb(iz)) * dt 
             cv1(ix,iy,iz) = cv(ix,iy,iz) * (1.0_num + dv)

             !This line is equivalent to IF (dv > 0.0_num) p_visc(ix,iy,iz) = 0.0_num
             p_visc(ix,iy,iz)=-p_visc(ix,iy,iz)*MIN(SIGN(1.0_num,dv),0.0_num)
             !energy at end of Lagrangian step

             energy(ix,iy,iz) = energy(ix,iy,iz) - pressure(ix,iy,iz) * dv / rho(ix,iy,iz)
#ifdef Q_MONO
             IF (visc3 > 1.e-6_num) THEN
#endif
                energy(ix,iy,iz) = energy(ix,iy,iz) &
                     + dt * visc_heat(ix,iy,iz) / rho(ix,iy,iz)  
#ifdef Q_MONO
             END IF
             energy(ix,iy,iz) = energy(ix,iy,iz) &
                  - p_visc(ix,iy,iz)*dv / rho(ix,iy,iz)
#endif

             rho(ix,iy,iz) = rho(ix,iy,iz) / (1.0_num + dv)!density at end of Lagrangian step
#ifdef Q_MONO
             IF (visc3 > 1.e-6_num) THEN 
#endif
                total_visc_heating = total_visc_heating &
                     + dt * visc_heat(ix,iy,iz) * cv(ix,iy,iz)
#ifdef Q_MONO
             END IF
             total_visc_heating = total_visc_heating &
                  - p_visc(ix,iy,iz) * dv * cv(ix,iy,iz) 
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
    REAL(num) :: sxx, syy, szz, sxy, sxz, syx, syz
    REAL(num) :: dvxdy, dvxdz, dvydx, dvydz, dvzdx, dvzdy
    REAL(num) :: w2_1,w2_2,w2_3,w2_4
    REAL(num) :: flag1,flag2,flag3,flag4

    p_visc = 0.0_num
    pressure = energy(-1:nx+2,-1:ny+2,-1:nz+2) &
         * (gamma - 1.0_num) * rho(-1:nx+2,-1:ny+2,-1:nz+2)

    DO iz = 0, nz+1
       DO iy = 0, ny+1
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
          DO ix = 0, nx+1
             ixp = ix + 1
             iyp = iy + 1
             izp = iz + 1
             ixm = ix - 1
             iym = iy - 1
             izm = iz - 1
             vxb = (vx(ix,iy,iz) + vx(ix,iym,iz) &
                  + vx(ix,iy,izm) + vx(ix,iym,izm)) / 4.0_num     !vx at Sx(i,j,k)
             vxbm = (vx(ixm,iy,iz) + vx(ixm,iym,iz) &
                  + vx(ixm,iy,izm) + vx(ixm,iym,izm)) / 4.0_num   !vx at Sx(i-1,j,k) 
             vyb = (vy(ix,iy,iz) + vy(ixm,iy,iz)  &
                  + vy(ix,iy,izm) + vy(ixm,iy,izm)) / 4.0_num    !vy at Sy(i,j,k) 
             vybm = (vy(ix,iym,iz) + vy(ixm,iym,iz)  &
                  + vy(ix,iym,izm) + vy(ixm,iym,izm)) / 4.0_num  !vy at Sy(i,j-1,k)   
             vzb = (vz(ix,iy,iz) + vz(ixm,iy,iz)  &
                  + vz(ix,iym,iz) + vz(ixm,iym,iz)) / 4.0_num    !vz at Sz(i,j,k) 
             vzbm = (vz(ix,iy,izm) + vz(ixm,iy,izm)  &
                  + vz(ix,iym,izm) + vz(ixm,iym,izm)) / 4.0_num  !vz at Sz(i,j,k-1)   

             dvxdx = (vxb - vxbm) / dxb(ix)
             dvydy = (vyb - vybm) / dyb(iy)
             dvzdz = (vzb - vzbm) / dzb(iz)

             dv = (dvxdx + dvydy + dvzdz) * dt2
             cv1(ix,iy,iz) = cv(ix,iy,iz) * (1.0_num + dv)

             vxb = (vx(ix,iy,iz) + vx(ixm,iy,iz)  &
                  + vx(ix,iy,izm) + vx(ixm,iy,izm)) / 4.0_num    !vx at Sy(i,j,k) 
             vxbm = (vx(ix,iym,iz) + vx(ixm,iym,iz)  &
                  + vx(ix,iym,izm) + vx(ixm,iym,izm)) / 4.0_num  !vx at Sy(i,j-1,k)  
             vyb = (vy(ix,iy,iz) + vy(ix,iym,iz) &
                  + vy(ix,iy,izm) + vy(ix,iym,izm)) / 4.0_num     !vy at Sx(i,j,k)
             vybm = (vy(ixm,iy,iz) + vy(ixm,iym,iz) &
                  + vy(ixm,iy,izm) + vy(ixm,iym,izm)) / 4.0_num   !vy at Sx(i-1,j,k) 
             dvxdy = (vxb - vxbm)/dyb(iy)
             dvydx = (vyb - vybm)/dxb(ix)
             dvxy = dvxdy + dvydx
             sxy = dvxy / 2.0_num
             sxx = 2.0_num * dvxdx / 3.0_num - (dvydy + dvzdz) / 3.0_num 
             syy = 2.0_num * dvydy / 3.0_num - (dvxdx + dvzdz) / 3.0_num 
             szz = 2.0_num * dvzdz / 3.0_num - (dvxdx + dvydy) / 3.0_num 

             vxb = (vx(ix,iy,iz) + vx(ixm,iy,iz)  &
                  + vx(ix,iym,iz) + vx(ixm,iym,iz)) / 4.0_num    !vx at Sz(i,j,k) 
             vxbm = (vx(ix,iy,izm) + vx(ixm,iy,izm)  &
                  + vx(ix,iym,izm) + vx(ixm,iym,izm)) / 4.0_num  !vx at Sz(i,j,k-1)   
             vzb = (vz(ix,iy,iz) + vz(ix,iym,iz) &
                  + vz(ix,iy,izm) + vz(ix,iym,izm)) / 4.0_num     !vz at Sx(i,j,k)
             vzbm = (vz(ixm,iy,iz) + vz(ixm,iym,iz) &
                  + vz(ixm,iy,izm) + vz(ixm,iym,izm)) / 4.0_num   !vz at Sx(i-1,j,k) 
             dvxdz = (vxb - vxbm)/dzb(iz) 
             dvzdx = (vzb - vzbm)/dxb(ix)
             dvxz = dvxdz + dvzdx
             sxz = dvxz / 2.0_num

             vyb = (vy(ix,iy,iz) + vy(ixm,iy,iz)  &
                  + vy(ix,iym,iz) + vy(ixm,iym,iz)) / 4.0_num    !vy at Sz(i,j,k) 
             vybm = (vy(ix,iy,izm) + vy(ixm,iy,izm)  &
                  + vy(ix,iym,izm) + vy(ixm,iym,izm)) / 4.0_num  !vy at Sz(i,j,k-1)   
             vzb = (vz(ix,iy,iz) + vz(ixm,iy,iz)  &
                  + vz(ix,iy,izm) + vz(ixm,iy,izm)) / 4.0_num    !vz at Sy(i,j,k) 
             vzbm = (vz(ix,iym,iz) + vz(ixm,iym,iz)  &
                  + vz(ix,iym,izm) + vz(ixm,iym,izm)) / 4.0_num  !vz at Sy(i,j-1,k)   
             dvydz = (vyb - vybm)/dzb(iz) 
             dvzdy = (vzb - vzbm)/dyb(iy)
             dvyz = dvydz + dvzdy
             syz = dvyz / 2.0_num

             p = pressure(ix,iy,iz) 
             pxp = pressure(ixp,iy,iz) 
             pxm = pressure(ixm,iy,iz) 
             pyp = pressure(ix,iyp,iz) 
             pym = pressure(ix,iym,iz) 
             pzp = pressure(ix,iy,izp) 
             pzm = pressure(ix,iy,izm) 

             fx = - (pxp - pxm) / dxb(ix)  !should be half this but this cancels later
             fy = - (pyp - pym) / dyb(iy)  !should be half this but this cancels later
             fz = - (pzp - pzm) / dzb(iz)  !should be half this but this cancels later

             w1 = fx**2 + fy**2 + fz**2
             s = (dvxdx*fx**2 + dvydy*fy**2 + dvzdz*fz**2 + dvxy*fx*fy &
                  + dvxz*fx*fz + dvyz*fy*fz) / (w1 + 1.e-20_num)

!	     These flags are used to replace the rather clearer code in 
!	     **SECTION 1**. They are used instead to facilitate vector
!	     optimization
             flag1=MAX(SIGN(1.0_num,dyb(iy)*ABS(fx)-dxb(ix)*ABS(fy)),0.0_num)
             flag2=MAX(SIGN(1.0_num,dzb(iz)*ABS(fx)-dxb(ix)*ABS(fz)),0.0_num)
             flag3=MAX(SIGN(1.0_num,dzb(iz)*ABS(fy)-dyb(iy)*ABS(fz)),0.0_num)
             flag4=MAX(SIGN(1.0_num,w1-1.e-6_num),0.0_num)

             w2_1=dxb(ix)**2*w1 / (fx**2 + 1.e-20_num)
             w2_2=dzb(iz)**2*w1 / (fz**2 + 1.e-20_num)
             w2_3=dyb(iy)**2*w1 / (fy**2 + 1.e-20_num)
             w2_4=dzb(iz)**2*w1 / (fz**2 + 1.e-20_num)

             w2=w2_1*flag1*flag2 + w2_2*flag1*(1.0_num-flag2)&
                  +w2_3*(1.0_num-flag1)*flag3 + w2_4*(1.0_num-flag1)*(1.0_num-flag3)

             w2=w2*flag4 + min_grid_spacing**2 * (1.0_num*flag4)
!!$		!BEGIN **SECTION 1**
!!$		!*******************
!!$             IF (dxb(ix)*ABS(fy) < dyb(iy)*ABS(fx)) THEN
!!$                IF (dxb(ix)*ABS(fz) < dzb(iz)*ABS(fx)) THEN
!!$                   w2 = dxb(ix)**2*w1 / (fx**2 + 1.e-20_num)
!!$                ELSE
!!$                   w2 = dzb(iz)**2*w1 / (fz**2 + 1.e-20_num)
!!$                END IF
!!$             ELSE
!!$                IF (dyb(iy)*ABS(fz) < dzb(iz)*ABS(fy)) THEN
!!$                   w2 = dyb(iy)**2*w1 / (fy**2 + 1.e-20_num)
!!$                ELSE
!!$                   w2 = dzb(iz)**2*w1 / (fz**2 + 1.e-20_num)
!!$                END IF
!!$             END IF
!!$             IF (w1 < 1.e-6_num) w2 = min_grid_spacing**2
!!$		!END **SECTION 1**
!!$		!*****************
             L = SQRT(w2)

             L2 = L  
!            These lines are equivalent to : IF (s > 0.0_num .OR. dv > 0.0_num) L = 0.0_num
             L=L*MIN(SIGN(1.0_num,s),0.0_num)
             L=L*MIN(SIGN(1.0_num,dv),0.0_num)

             w1 = (bx1(ix,iy,iz)**2 + by1(ix,iy,iz)**2 + bz1(ix,iy,iz)**2)
             cf = SQRT((w1 + gamma*p)/ rho(ix,iy,iz))

             p_visc(ix,iy,iz) = visc1*ABS(s)*L*cf*rho(ix,iy,iz) &
                  + visc2*(s*L)**2*rho(ix,iy,iz)

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
#ifdef Q_MONO
             IF (visc3 > 1.e-6_num)  THEN
             L = min_grid_spacing
             qxy(ix,iy,iz) = sxy + sxy * L * rho(ix,iy,iz) * visc3
             qxz(ix,iy,iz) = sxz + sxz * L * rho(ix,iy,iz) * visc3 
             qyz(ix,iy,iz) = syz + syz * L * rho(ix,iy,iz) * visc3 
             qxx(ix,iy,iz) = sxx + sxx * L * rho(ix,iy,iz) * visc3 
             qyy(ix,iy,iz) = syy + syy * L * rho(ix,iy,iz) * visc3 
             qzz(ix,iy,iz) = szz + szz * L * rho(ix,iy,iz) * visc3 
          END IF
#endif

#ifdef Q_MONO
          IF (visc3 > 1.e-6_num) THEN
#endif
             visc_heat(ix,iy,iz) = qxy(ix,iy,iz)*dvxy + qxz(ix,iy,iz)*dvxz &
                  + qyz(ix,iy,iz)*dvyz + qxx(ix,iy,iz)*dvxdx   &
                  + qyy(ix,iy,iz)*dvydy + qzz(ix,iy,iz)*dvzdz 
#ifdef Q_MONO
          END IF
#endif

          w4 = bx1(ix,iy,iz)*dvxdx + by1(ix,iy,iz)*dvxdy + bz1(ix,iy,iz)*dvxdz 
          bx1(ix,iy,iz) = (bx1(ix,iy,iz) + w4 * dt2) / (1.0_num + dv)
          w4 = bx1(ix,iy,iz)*dvydx + by1(ix,iy,iz)*dvydy + bz1(ix,iy,iz)*dvydz
          by1(ix,iy,iz) = (by1(ix,iy,iz) + w4 * dt2) / (1.0_num + dv)
          w4 = bx1(ix,iy,iz)*dvzdx + by1(ix,iy,iz)*dvzdy + bz1(ix,iy,iz)*dvzdz 
          bz1(ix,iy,iz) = (bz1(ix,iy,iz) + w4 * dt2) / (1.0_num + dv)

       END DO
    END DO
 END DO

END SUBROUTINE viscosity_and_b_update



  SUBROUTINE visc_heating

    REAL(num) :: vxb, vxbm, vyb, vybm, vzb, vzbm
    REAL(num) :: dvxdx, dvydy, dvzdz, dvxy, dvxz, dvyz

    DO iz = 0, nz+1
       DO iy = 0, ny+1
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
          DO ix = 0, nx+1
             ixp = ix + 1
             iyp = iy + 1
             izp = iz + 1
             ixm = ix - 1
             iym = iy - 1
             izm = iz - 1
             vxb = (vx1(ix,iy,iz) + vx1(ix,iym,iz) &
                  + vx1(ix,iy,izm) + vx1(ix,iym,izm)) / 4.0_num     !vx at Sx(i,j,k)
             vxbm = (vx1(ixm,iy,iz) + vx1(ixm,iym,iz) &
                  + vx1(ixm,iy,izm) + vx1(ixm,iym,izm)) / 4.0_num   !vx at Sx(i-1,j,k) 
             vyb = (vy1(ix,iy,iz) + vy1(ixm,iy,iz)  &
                  + vy1(ix,iy,izm) + vy1(ixm,iy,izm)) / 4.0_num    !vy at Sy(i,j,k) 
             vybm = (vy1(ix,iym,iz) + vy1(ixm,iym,iz)  &
                  + vy1(ix,iym,izm) + vy1(ixm,iym,izm)) / 4.0_num  !vy at Sy(i,j-1,k)   
             vzb = (vz1(ix,iy,iz) + vz1(ixm,iy,iz)  &
                  + vz1(ix,iym,iz) + vz1(ixm,iym,iz)) / 4.0_num    !vz at Sz(i,j,k) 
             vzbm = (vz1(ix,iy,izm) + vz1(ixm,iy,izm)  &
                  + vz1(ix,iym,izm) + vz1(ixm,iym,izm)) / 4.0_num  !vz at Sz(i,j,k-1)   

             dvxdx = (vxb - vxbm) / dxb(ix)
             dvydy = (vyb - vybm) / dyb(iy)
             dvzdz = (vzb - vzbm) / dzb(iz)

             vxb = (vx1(ix,iy,iz) + vx1(ixm,iy,iz)  &
                  + vx1(ix,iy,izm) + vx1(ixm,iy,izm)) / 4.0_num 
             vxbm = (vx1(ix,iym,iz) + vx1(ixm,iym,iz)  &
                  + vx1(ix,iym,izm) + vx1(ixm,iym,izm)) / 4.0_num  
             vyb = (vy1(ix,iy,iz) + vy1(ix,iym,iz) &
                  + vy1(ix,iy,izm) + vy1(ix,iym,izm)) / 4.0_num    
             vybm = (vy1(ixm,iy,iz) + vy1(ixm,iym,iz) &
                  + vy1(ixm,iy,izm) + vy1(ixm,iym,izm)) / 4.0_num  
             dvxy = (vxb - vxbm)/dyb(iy) + (vyb - vybm)/dxb(ix)

             vxb = (vx1(ix,iy,iz) + vx1(ixm,iy,iz)  &
                  + vx1(ix,iym,iz) + vx1(ixm,iym,iz)) / 4.0_num   
             vxbm = (vx1(ix,iy,izm) + vx1(ixm,iy,izm)  &
                  + vx1(ix,iym,izm) + vx1(ixm,iym,izm)) / 4.0_num    
             vzb = (vz1(ix,iy,iz) + vz1(ix,iym,iz) &
                  + vz1(ix,iy,izm) + vz1(ix,iym,izm)) / 4.0_num   
             vzbm = (vz1(ixm,iy,iz) + vz1(ixm,iym,iz) &
                  + vz1(ixm,iy,izm) + vz1(ixm,iym,izm)) / 4.0_num  
             dvxz = (vxb - vxbm)/dzb(iz) + (vzb - vzbm)/dxb(ix) 

             vyb = (vy1(ix,iy,iz) + vy1(ixm,iy,iz)  &
                  + vy1(ix,iym,iz) + vy1(ixm,iym,iz)) / 4.0_num    
             vybm = (vy1(ix,iy,izm) + vy1(ixm,iy,izm)  &
                  + vy1(ix,iym,izm) + vy1(ixm,iym,izm)) / 4.0_num   
             vzb = (vz1(ix,iy,iz) + vz1(ixm,iy,iz)  &
                  + vz1(ix,iy,izm) + vz1(ixm,iy,izm)) / 4.0_num    
             vzbm = (vz1(ix,iym,iz) + vz1(ixm,iym,iz)  &
                  + vz1(ix,iym,izm) + vz1(ixm,iym,izm)) / 4.0_num    
             dvyz = (vyb - vybm)/dzb(iz) + (vzb - vzbm)/dyb(iy)

             visc_heat(ix,iy,iz) = qxy(ix,iy,iz)*dvxy + qxz(ix,iy,iz)*dvxz &
                  + qyz(ix,iy,iz)*dvyz + qxx(ix,iy,iz)*dvxdx   &
                  + qyy(ix,iy,iz)*dvydy + qzz(ix,iy,iz)*dvzdz 
          END DO
       END DO
    END DO

  END SUBROUTINE visc_heating




  SUBROUTINE eta_calc

    REAL(num) :: jx, jy, jz, jxxp, jyyp, jzzp, flux

    eta = 0.0_num
    curlb = 0.0_num

    DO iz = -1, nz+1
       DO iy = -1, ny+1
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
          DO ix = -1, nx+1
             ixp = ix + 1
             iyp = iy + 1
             izp = iz + 1

             jx = (bz(ix,iyp,iz) - bz(ix,iy,iz)) / dyc(iy)  &
                  - (by(ix,iy,izp) - by(ix,iy,iz)) / dzc(iz)       ! jx at E3(i,j)
             jxxp = (bz(ixp,iyp,iz) - bz(ixp,iy,iz)) / dyc(iy)  &
                  - (by(ixp,iy,izp) - by(ixp,iy,iz)) / dzc(iz)    ! jx at E3(i+1,j)
             jy = (bx(ix,iy,izp) - bx(ix,iy,iz)) / dzc(iz)  &
                  - (bz(ixp,iy,iz) - bz(ix,iy,iz)) / dxc(ix)       ! jy at E2(i,j)
             jyyp = (bx(ix,iyp,izp) - bx(ix,iyp,iz)) / dzc(iz)  &
                  - (bz(ixp,iyp,iz) - bz(ix,iyp,iz)) / dxc(ix)    ! jy at E2(i,j+1)
             jz = (by(ixp,iy,iz) - by(ix,iy,iz)) / dxc(ix)  &     
                  - (bx(ix,iyp,iz) - bx(ix,iy,iz)) / dyc(iy)       ! jz at E1(i,j)
             jzzp = (by(ixp,iy,izp) - by(ix,iy,izp)) / dxc(ix)  &     
                  - (bx(ix,iyp,izp) - bx(ix,iy,izp)) / dyc(iy)    ! jz at E1(i,j)

             !current at V
             w4 = (jx+jxxp) / 2.0_num
             w5 = (jy+jyyp) / 2.0_num
             w6 = (jz+jzzp) / 2.0_num
             curlb(ix,iy,iz) = SQRT(w4**2 + w5**2 + w6**2)  
             IF (curlb(ix,iy,iz) >= j_max) THEN
                eta(ix,iy,iz) = eta0 * min_grid_spacing
             ELSE
                eta(ix,iy,iz) = eta_background 
             ENDIF
          END DO
       END DO
    END DO

    IF (.NOT. resistiveMHD) eta = 0.0_num

  END SUBROUTINE eta_calc



  SUBROUTINE resistive_effects

    REAL(num) :: jx, jy, jz, jxxp, jyyp, jzzp, flux
    REAL(num) :: bxv, byv, bzv, j_par_x, j_par_y, j_par_z
    REAL(num) :: j_perp_x, j_perp_y, j_perp_z, magn_b
    REAL(num) :: magn_j_par, magn_j_perp,flip,magn_j
    REAL(num), DIMENSION(:),ALLOCATABLE:: maskx,masky,maskz

    ALLOCATE(maskx(0:nx),masky(0:ny),maskz(0:nz))

    !These mask variables are used later to multiply the
    !edge cells in each direction by 0.5
    maskx=1.0_num
    maskx(0)=0.5_num
    maskx(nx)=0.5_num

    masky=1.0_num
    masky(0)=0.5_num
    masky(ny)=0.5_num

    maskz=1.0_num
    maskz(0)=0.5_num
    maskz(nz)=0.5_num

    DO iz = 0, nz+1
       DO iy = 0, ny+1
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
          DO ix = 0, nx+1
             bx1(ix,iy,iz) = bx(ix,iy,iz) * dyb(iy) * dzb(iz)  
             by1(ix,iy,iz) = by(ix,iy,iz) * dxb(ix) * dzb(iz)
             bz1(ix,iy,iz) = bz(ix,iy,iz) * dxb(ix) * dyb(iy)
          END DO
       END DO
    END DO

    DO iz = 0, nz
       DO iy = 0, ny
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
          DO ix = 0, nx
             ixp = ix + 1
             iyp = iy + 1
             izp = iz + 1

             jx = (bz(ix,iyp,iz) - bz(ix,iy,iz)) / dyc(iy)  &
                  - (by(ix,iy,izp) - by(ix,iy,iz)) / dzc(iz)       ! jx at E3(i,j)
             jxxp = (bz(ixp,iyp,iz) - bz(ixp,iy,iz)) / dyc(iy)  &
                  - (by(ixp,iy,izp) - by(ixp,iy,iz)) / dzc(iz)    ! jx at E3(i+1,j)
             jy = (bx(ix,iy,izp) - bx(ix,iy,iz)) / dzc(iz)  &
                  - (bz(ixp,iy,iz) - bz(ix,iy,iz)) / dxc(ix)       ! jy at E2(i,j)
             jyyp = (bx(ix,iyp,izp) - bx(ix,iyp,iz)) / dzc(iz)  &
                  - (bz(ixp,iyp,iz) - bz(ix,iyp,iz)) / dxc(ix)    ! jy at E2(i,j+1)
             jz = (by(ixp,iy,iz) - by(ix,iy,iz)) / dxc(ix)  &     
                  - (bx(ix,iyp,iz) - bx(ix,iy,iz)) / dyc(iy)       ! jz at E1(i,j)
             jzzp = (by(ixp,iy,izp) - by(ix,iy,izp)) / dxc(ix)  &     
                  - (bx(ix,iyp,izp) - bx(ix,iy,izp)) / dyc(iy)    ! jz at E1(i,j)

             !current at V
             jx = (jx+jxxp) / 2.0_num
             jy = (jy+jyyp) / 2.0_num
             jz = (jz+jzzp) / 2.0_num

             curlb(ix,iy,iz) = SQRT(jx**2 + jy**2 + jz**2)

             ! B at vertices
             bxv = (bx(ix,iy,iz) + bx(ix,iyp,iz) + bx(ix,iy,izp) &
                  + bx(ix,iyp,izp)) / 4.0_num
             byv = (by(ix,iy,iz) + by(ixp,iy,iz) + by(ix,iy,izp) &
                  + by(ixp,iy,izp)) / 4.0_num
             bzv = (bz(ix,iy,iz) + bz(ixp,iy,iz) + bz(ix,iyp,iz) &
                  + bz(ixp,iyp,iz)) / 4.0_num
             magn_b = SQRT(bxv**2 + byv**2 + bzv**2)

             ! it's heat not flux but:
             flux = (dt * ((eta(ix,iy,iz) * magn_j**2))  &
                  / 8.0_num)*PI
             ! it's heat not flux but:
             flux = flux + (dt * eta(ix,iy,iz) * curlb(ix,iy,iz)**2 / 8.0_num) ! eta j^2 at vertex 

             ! Ohmic heating split between 8 cells adjacent to point V
             energy(ix,iy,iz) = energy(ix,iy,iz) +  flux / rho(ix,iy,iz)
             energy(ixp,iy,iz) = energy(ixp,iy,iz) + flux / rho(ixp,iy,iz)
             energy(ix,iyp,iz) = energy(ix,iyp,iz) + flux / rho(ix,iyp,iz)
             energy(ixp,iyp,iz) = energy(ixp,iyp,iz) + flux / rho(ixp,iyp,iz)
             energy(ix,iy,izp) = energy(ix,iy,izp) + flux / rho(ix,iy,izp)
             energy(ixp,iy,izp) = energy(ixp,iy,izp) + flux / rho(ixp,iy,izp)
             energy(ix,iyp,izp) = energy(ix,iyp,izp) + flux / rho(ix,iyp,izp)
             energy(ixp,iyp,izp) = energy(ixp,iyp,izp) + flux / rho(ixp,iyp,izp)

             w2 = flux * cv(ix,iy,iz) * 8.0_num

	     w2 = w2 * maskx(ix)
	     w2 = w2 * masky(iy)
             w2 = w2 * maskz(iz)         

             total_ohmic_heating = total_ohmic_heating + w2

             flux = flux-(jx * eta(ix,iy,iz) * dxc(ix) * dt / 2.0_num)

             bz1(ix,iy,iz) = bz1(ix,iy,iz) - flux
             bz1(ix,iyp,iz) = bz1(ix,iyp,iz) + flux
             bz1(ixp,iy,iz) = bz1(ixp,iy,iz) - flux
             bz1(ixp,iyp,iz) = bz1(ixp,iyp,iz) + flux
             by1(ix,iy,iz) = by1(ix,iy,iz) + flux
             by1(ix,iy,izp) = by1(ix,iy,izp) - flux
             by1(ixp,iy,iz) = by1(ixp,iy,iz) + flux
             by1(ixp,iy,izp) = by1(ixp,iy,izp) - flux

             flux = flux-(jy * eta(ix,iy,iz) * dyc(iy) * dt / 2.0_num )
             bx1(ix,iy,iz) = bx1(ix,iy,iz) - flux
             bx1(ix,iy,izp) = bx1(ix,iy,izp) + flux
             bx1(ix,iyp,iz) = bx1(ix,iyp,iz) - flux
             bx1(ix,iyp,izp) = bx1(ix,iyp,izp) + flux
             bz1(ix,iy,iz) = bz1(ix,iy,iz) + flux
             bz1(ixp,iy,iz) = bz1(ixp,iy,iz) - flux
             bz1(ix,iyp,iz) = bz1(ix,iyp,iz) + flux
             bz1(ixp,iyp,iz) = bz1(ixp,iyp,iz) - flux

             flux = flux-(jz * eta(ix,iy,iz) * dzc(iz) * dt / 2.0_num)
             by1(ix,iy,iz) = by1(ix,iy,iz) - flux
             by1(ixp,iy,iz) = by1(ixp,iy,iz) + flux
             by1(ix,iy,izp) = by1(ix,iy,izp) - flux
             by1(ixp,iy,izp) = by1(ixp,iy,izp) + flux
             bx1(ix,iy,iz) = bx1(ix,iy,iz) + flux
             bx1(ix,iyp,iz) = bx1(ix,iyp,iz) - flux
             bx1(ix,iy,izp) = bx1(ix,iy,izp) + flux
             bx1(ix,iyp,izp) = bx1(ix,iyp,izp) - flux

          END DO
       END DO
    END DO

    DO iz = 0, nz+1
       DO iy = 0, ny+1
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
          DO ix = 0, nx+1
             bx(ix,iy,iz) = bx1(ix,iy,iz) / (dyb(iy) * dzb(iz))  
             by(ix,iy,iz) = by1(ix,iy,iz) / (dxb(ix) * dzb(iz))
             bz(ix,iy,iz) = bz1(ix,iy,iz) / (dxb(ix) * dyb(iy))
          END DO
       END DO
    END DO

    CALL bfield_bcs
    CALL energy_bcs


  END SUBROUTINE resistive_effects



  SUBROUTINE store_boundary_dv

    REAL(num) :: dvx, dvy, dvz

    IF (xbc_right == open .AND. right == MPI_PROC_NULL) THEN
       DO iz = -2,nz+2
          DO iy = -2,ny+2
             dvx=2.0_num*(vx(nx,iy,iz)-vx1(nx,iy,iz))
             dvy=2.0_num*(vy(nx,iy,iz)-vy1(nx,iy,iz))
             dvz=2.0_num*(vz(nx,iy,iz)-vz1(nx,iy,iz))
             dv_right(iy,iz)=SQRT(dvx**2+dvy**2+dvz**2)
          END DO
       END DO
    END IF
    IF (xbc_left == open .AND. left == MPI_PROC_NULL) THEN
       DO iz = -2,nz+2
          DO iy = -2,ny+2
             dvx=2.0_num*(vx(0,iy,iz)-vx1(0,iy,iz))
             dvy=2.0_num*(vy(0,iy,iz)-vy1(0,iy,iz))
             dvz=2.0_num*(vz(0,iy,iz)-vz1(0,iy,iz))
             dv_left(iy,iz)=SQRT(dvx**2+dvy**2+dvz**2)
          END DO
       END DO
    END IF
    IF (ybc_up == open .AND. up == MPI_PROC_NULL) THEN
       DO iz = -2,nz+2
          DO ix = -2,nx+2
             dvx=2.0_num*(vx(ix,ny,iz)-vx1(ix,ny,iz))
             dvy=2.0_num*(vy(ix,ny,iz)-vy1(ix,ny,iz))
             dvz=2.0_num*(vz(ix,ny,iz)-vz1(ix,ny,iz))
             dv_up(ix,iz)=SQRT(dvx**2+dvy**2+dvz**2)
          END DO
       END DO
    END IF
    IF (ybc_down == open .AND. down == MPI_PROC_NULL) THEN
       DO iz = -2,nz+2
          DO ix = -2,nx+2
             dvx=2.0_num*(vx(ix,0,iz)-vx1(ix,0,iz))
             dvy=2.0_num*(vy(ix,0,iz)-vy1(ix,0,iz))
             dvz=2.0_num*(vz(ix,0,iz)-vz1(ix,0,iz))
             dv_down(ix,iz)=SQRT(dvx**2+dvy**2+dvz**2)
          END DO
       END DO
    END IF
    IF (zbc_back == open .AND. back == MPI_PROC_NULL) THEN
       DO iy = -2,ny+2
          DO ix = -2,nx+2
             dvx=2.0_num*(vx(ix,iy,nz)-vx1(ix,iy,nz))
             dvy=2.0_num*(vy(ix,iy,nz)-vy1(ix,iy,nz))
             dvz=2.0_num*(vz(ix,iy,nz)-vz1(ix,iy,nz))
             dv_back(iy,iz)=SQRT(dvx**2+dvy**2+dvz**2)
          END DO
       END DO
    END IF
    IF (zbc_front == open .AND. front == MPI_PROC_NULL) THEN
       DO iy = -2,ny+2
          DO ix = -2,nx+2
             dvx=2.0_num*(vx(ix,iy,0)-vx1(ix,iy,0))
             dvy=2.0_num*(vy(ix,iy,0)-vy1(ix,iy,0))
             dvz=2.0_num*(vz(ix,iy,0)-vz1(ix,iy,0))
             dv_front(iy,iz)=SQRT(dvx**2+dvy**2+dvz**2)
          END DO
       END DO
    END IF

  END SUBROUTINE store_boundary_dv

END MODULE lagran














