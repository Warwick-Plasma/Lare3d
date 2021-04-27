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

MODULE conduct

  USE shared_data
  USE boundary
  USE neutral

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: conduct_heat, calc_s_stages

  REAL(num),PARAMETER :: pow=5.0_num/2.0_num
  REAL(num),PARAMETER :: min_b=1.0e-10_num


CONTAINS


  !****************************************************************************
  ! Subroutine implementing Braginskii parallel thermal conduction.
  ! Notation and algorithm in Appendix of Manual
  !****************************************************************************

  SUBROUTINE conduct_heat

    ALLOCATE(larsen_factor(0:nx,0:ny,0:nz))

    CALL calc_s_stages(.FALSE.)
    CALL heat_conduct_sts2

    DEALLOCATE(larsen_factor)

  END SUBROUTINE conduct_heat



  !****************************************************************************
  ! Function calculating the number of stages needed for the RKL2
  !****************************************************************************

  SUBROUTINE calc_s_stages(lagrangian_call)

    LOGICAL, INTENT(IN) :: lagrangian_call

    REAL(num)  ::  stages, dt_parab, dt1, dt2, dt3, temp
    REAL(num) :: kappa1
    REAL(num) :: q_fs, q_fs2, q_spx, q_spy, q_spz, q_sp2
    INTEGER :: n_s_stages_local, nstages

    ! Make sure arrays are allocated if calling this routine just to determine
    ! dt from the lagrangian setp
    IF (lagrangian_call) ALLOCATE(larsen_factor(0:nx,0:ny,0:nz))

    DO iz=-1,nz+2
      DO iy=-1,ny+2
        DO ix=-1,nx+2
          temperature(ix,iy,iz) = (gamma - 1.0_num) / (2.0_num - xi_n(ix,iy,iz)) &
              * (energy(ix,iy,iz) - (1.0_num - xi_n(ix,iy,iz)) * ionise_pot)
        ENDDO
      ENDDO
    ENDDO

    ! Include flux limiting through a larson factor correction to the
    ! conductivity
    DO iz = 0, nz
      izm = iz - 1
      izp = iz + 1
      DO iy = 0, ny
        iym = iy - 1
        iyp = iy + 1
        DO ix = 0, nx
          ixm = ix - 1
          ixp = ix + 1
          temp = temperature(ix,iy,iz)**pow
          q_fs = flux_limiter * 42.85_num * rho(ix,iy,iz) &  ! 42.85 = SQRT(m_i/m_e)
              * temperature(ix,iy,iz)**1.5_num
          q_fs2 = q_fs**2
          q_spx = - kappa_0 * temp &
              * (temperature(ixp,iy,iz) - temperature(ixm,iy,iz)) &
              * 0.5_num / dxb(ix)
          q_spy = - kappa_0 * temp &
              * (temperature(ix,iyp,iz) - temperature(ix,iym,iz)) &
              * 0.5_num / dyb(iy)
          q_spz = - kappa_0 * temp &
              * (temperature(ix,iy,izp) - temperature(ix,iy,izm)) &
              * 0.5_num / dzb(iz)
          q_sp2 = q_spx**2 + q_spy**2 + q_spz**2
          larsen_factor(ix,iy,iz) = q_fs / SQRT(q_fs2 + q_sp2)
        END DO
      END DO
    END DO

    IF (.NOT. heat_flux_limiter) larsen_factor = 1.0_num

    dt_parab = 1.e10_num
    DO iz = 1 , nz
      DO iy = 1 , ny
        DO ix = 1 , nx
          ! explicit TC time-step
          kappa1 = kappa_0 * larsen_factor(ix,iy,iz)
          temp = (gamma - 1.0_num) / (2.0_num - xi_n(ix,iy,iz)) &
               * (energy(ix,iy,iz) - (1.0_num - xi_n(ix,iy,iz)) * ionise_pot)
          dt1 = rho(ix,iy,iz)*dxb(ix)**2 / (2.0_num * kappa1 * &
              temp**pow)
          dt2 = rho(ix,iy,iz)*dyb(iy)**2 / (2.0_num * kappa1 * &
              temp**pow)
          dt3 = rho(ix,iy,iz)*dzb(iz)**2 / (2.0_num * kappa1 * &
              temp**pow)
          dt_parab = MIN(dt_parab, dt1, dt2, dt3)
        END DO
      END DO
    END DO

    dt_parab = dt_multiplier * dt_parab / SQRT(3.0_num)

    stages = 0.5_num*(SQRT(9.0_num + 16.0_num * (dt/dt_parab))-1.0_num)

    n_s_stages_local = CEILING(stages)
    IF (MODULO(n_s_stages_local,2) .EQ. 0) THEN
      n_s_stages_local = n_s_stages_local + 1
    ENDIF
    CALL MPI_ALLREDUCE(n_s_stages_local, nstages, 1, &
        MPI_INTEGER, MPI_MAX, comm, errcode)
    n_s_stages = nstages

    IF (lagrangian_call) DEALLOCATE(larsen_factor)

  END SUBROUTINE calc_s_stages




  !****************************************************************************
  ! Subroutine to calculate the heat flux
  !****************************************************************************

  SUBROUTINE heat_flux(temperature, flux)

    REAL(num),INTENT(IN),DIMENSION(-1:,-1:,-1:) :: temperature
    REAL(num),INTENT(OUT),DIMENSION(-1:,-1:,-1:) :: flux
    INTEGER :: ix, ixp, ixm
    INTEGER :: iy, iyp, iym
    INTEGER :: iz, izp, izm
    REAL(num) :: tb, tg, fc_sp
    REAL(num) :: tg_a1, tg_a2, tb_p, tb_m
    REAL(num) :: modb
    REAL(num) :: byf, bxf, bzf

    flux=0.0_num
    DO iz = 0,nz
      DO iy = 0,ny
        DO ix = 0,nx
          ixp = ix + 1
          ixm = ix - 1
          iyp = iy + 1
          iym = iy - 1
          izp = iz + 1
          izm = iz - 1

          !X flux
          byf=0.25_num*(by(ix,iy,iz)+by(ixp,iy,iz)+&
              by(ix,iym,iz)+by(ixp,iym,iz))
          bzf=0.25_num*(bz(ix,iy,iz)+bz(ixp,iy,iz)+&
              bz(ix,iy,izm)+bz(ixp,iy,izm))

          modb = bx(ix,iy,iz)**2 + byf**2 + bzf**2 + min_b

          !Braginskii Conductive Flux in X
          !Temperature at the x boundaries in the current cell
          tb = (temperature(ixp,iy,iz) + temperature(ix,iy,iz))/2.0_num
          !Temperature at the x boundaries in the cell above
          tb_p = (temperature(ixp,iyp,iz)+temperature(ix,iyp,iz))/2.0_num
          !Temperature at the x boundaries in the cell below
          tb_m = (temperature(ixp,iym,iz)+temperature(ix,iym,iz))/2.0_num
          !Y temperature gradient at the x boundaries of the current cell
          !Uses centred difference on averaged values, so likely very smoothed
          tg_a1 = (tb_p-tb_m)/(dyc(iy)+dyc(iym))
          !Temperature at the x boundaries in the cell front
          tb_p = (temperature(ixp,iy,izp)+temperature(ix,iy,izp))/2.0_num
          !Temperature at the x boundaries in the cell back
          tb_m = (temperature(ixp,iy,izm)+temperature(ix,iy,izm))/2.0_num
          !Z temperature gradient at the x boundaries of the current cell
          !Uses centred difference on averaged values, so likely very smoothed
          tg_a2 = (tb_p - tb_m)/(dzc(iz)+dzc(izm))

          !X temperature gradient at the x boundaries of the current cell
          tg = (temperature(ixp,iy,iz) - temperature(ix,iy,iz))/dxc(ix)

          tb = MAX(tb, 0.0_num)

          fc_sp = - larsen_factor(ix,iy,iz) * kappa_0 * tb**pow * (bx(ix,iy,iz) * (tg * bx(ix,iy,iz) + &
              tg_a1 * byf + tg_a2 * bzf)+tg*min_b) / modb

          flux(ix,iy,iz) = flux(ix,iy,iz) - fc_sp/dxb(ix)
          flux(ixp,iy,iz) = flux(ixp,iy,iz) + fc_sp/dxb(ix)

          !Y flux
          bxf=0.25_num*(bx(ix,iy,iz)+bx(ix,iyp,iz)+&
              bx(ixm,iy,iz)+bx(ixm,iyp,iz))
          bzf=0.25_num*(bz(ix,iy,iz)+bz(ix,iyp,iz)+&
              bz(ix,iy,izm)+bz(ix,iyp,izm))

          modb = bxf**2 + by(ix,iy,iz)**2 + bzf**2 + min_b

          !Braginskii Conductive Flux in Y
          !Temperature at the y boundaries in the current cell
          tb = (temperature(ix,iyp,iz) + temperature(ix,iy,iz))/2.0_num
          !Temperature at the y boundaries in the cell right
          tb_p = (temperature(ixp,iy,iz)+temperature(ixp,iyp,iz))/2.0_num
          !Temperature at the y boundaries in the cell left
          tb_m = (temperature(ixm,iyp,iz)+temperature(ixm,iy,iz))/2.0_num
          !X temperature gradient at the y boundaries of the current cell
          !Uses centred difference on averaged values, so likely very smoothed
          tg_a1 = (tb_p-tb_m)/(dxc(ix)+dxc(ixm))

          !Temperature at the y boundaries in the cell front
          tb_p = (temperature(ix,iyp,izp)+temperature(ix,iy,izp))/2.0_num
          !Temperature at the x boundaries in the cell back
          tb_m = (temperature(ix,iyp,izm)+temperature(ix,iy,izm))/2.0_num
          !Z temperature gradient at the y boundaries of the current cell
          !Uses centred difference on averaged values, so likely very smoothed
          tg_a2 = (tb_p-tb_m)/(dzc(iz)+dzc(izm))

          !X temperature gradient at the y boundaries of the current cell
          tg = (temperature(ix,iyp,iz) - temperature(ix,iy,iz))/dyc(iy)

          tb = MAX(tb, 0.0_num)

          fc_sp = - larsen_factor(ix,iy,iz) * kappa_0 * tb**pow * (by(ix,iy,iz) * (tg * by(ix,iy,iz) + &
              tg_a1 * bxf + tg_a2 * bzf)+tg*min_b) / modb

          flux(ix,iy,iz) = flux(ix,iy,iz) - fc_sp/dyb(iy)
          flux(ix,iyp,iz) = flux(ix,iyp,iz) + fc_sp/dyb(iy)

          !Z flux
          bxf=0.25_num*(bx(ix,iy,iz)+bx(ix,iy,izp)+&
              bx(ixm,iy,iz)+bx(ixm,iy,izp))
          byf=0.25_num*(by(ix,iy,iz)+by(ix,iy,izp)+&
              by(ix,iym,iz)+by(ix,iym,izp))

          modb = bxf**2 + byf**2 +bz(ix,iy,iz)**2 + min_b

          !Braginskii Conductive Flux in Z
          !Temperature at the z boundaries in the current cell
          tb = (temperature(ix,iy,izp) + temperature(ix,iy,iz))/2.0_num
          !Temperature at the z boundaries in the cell right
          tb_p = (temperature(ixp,iy,izp)+temperature(ixp,iy,iz))/2.0_num
          !Temperature at the z boundaries in the cell left
          tb_m = (temperature(ixm,iy,izp)+temperature(ixm,iy,iz))/2.0_num
          !X temperature gradient at the z boundaries of the current cell
          !Uses centred difference on averaged values, so likely very smoothed
          tg_a1 = (tb_p-tb_m)/(dxc(ix)+dxc(ixm))

          !Temperature at the z boundaries in the cell up
          tb_p = (temperature(ix,iyp,izp)+temperature(ix,iyp,iz))/2.0_num
          !Temperature at the z boundaries in the cell down
          tb_m = (temperature(ix,iym,izp)+temperature(ix,iym,iz))/2.0_num
          !Y temperature gradient at the z boundaries of the current cell
          !Uses centred difference on averaged values, so likely very smoothed
          tg_a2 = (tb_p-tb_m)/(dyc(iy)+dyc(iym))

          !X temperature gradient at the x boundaries of the current cell
          tg = (temperature(ix,iy,izp) - temperature(ix,iy,iz))/dzc(iz)

          tb = MAX(tb, 0.0_num)

          fc_sp = - larsen_factor(ix,iy,iz) * kappa_0 * tb**pow * (bz(ix,iy,iz) * (tg * bz(ix,iy,iz) + &
              tg_a1 * bxf + tg_a2 * byf)+tg*min_b) / modb

          flux(ix,iy,iz) = flux(ix,iy,iz) - fc_sp/dzb(iz)
          flux(ix,iy,izp) = flux(ix,iy,izp) + fc_sp/dzb(iz)
        END DO
      END DO
    END DO

  END SUBROUTINE heat_flux



  !****************************************************************************
  ! Implementation of the RKL2 scheme
  !****************************************************************************
  
  SUBROUTINE heat_conduct_sts2

    !Superstepping based conduction code
    !2nd order Runge-Kutta-Lagrange (RKL2) scheme
    !Based on Meyer at al. 2012 variables named as in that paper
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE  :: flux
    REAL(num)  ::  omega_1
    REAL(num), DIMENSION(0:n_s_stages)  :: a, b
    REAL(num), DIMENSION(1:n_s_stages)  :: mu_tilde
    REAL(num), DIMENSION(2:n_s_stages)  :: mu, nu, gamma_tilde
    ! intermediate solutions Y=[Y_0, Y_j-2, Y_j-1 , Y_j]
    REAL(num), DIMENSION(:,:,:,:), ALLOCATABLE  :: Y             
    REAL(num)  ::  c0, c1
    REAL(num)  ::  Lc_Yj_1                  ! L^c(Y_j-1)
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE  :: Lc_Y0    ! L^c(Y_0)
    INTEGER :: j

    ALLOCATE(flux(-1:nx+2,-1:ny+2,-1:nz+2))
    ALLOCATE(Y(0:3,-1:nx+2,-1:ny+2,-1:nz+2))
    ALLOCATE(Lc_Y0(1:nx,1:ny,1:nz))
    flux=0.0_num
    Y=0.0_num
    
    omega_1 = 4.0_num/(DBLE(n_s_stages)**2.0_num + DBLE(n_s_stages) - 2.0_num)
    b(0:2) = 1.0_num/3.0_num
    DO j = 3, n_s_stages
       b(j) = (DBLE(j)**2.0_num + DBLE(j)- 2.0_num) / &
           (2.0_num * DBLE(j) *(DBLE(j) + 1.0_num))
    ENDDO 
    a = 1.0_num - b
    mu_tilde(1) = omega_1/3.0_num

    DO j=2,n_s_stages
      mu_tilde(j) = ((2.0_num *DBLE(j) - 1.0_num) /DBLE(j))* &
          omega_1 *(b(j)/b(j-1))
      mu(j) = ((2.0_num *DBLE(j) - 1.0_num) /DBLE(j)) *(b(j)/b(j-1))
      nu(j) = -1.0_num * ((DBLE(j) - 1.0_num) /DBLE(j))*(b(j)/b(j-2))
      gamma_tilde(j) = -1.0_num*a(j-1)*mu_tilde(j)
    ENDDO

    !!!!! RKL2 s stage scheme !!!!!
    !! boundary conditions 
    DO j=0,3
      Y(j,:,:,:) = energy
    END DO

    !! First STS stage
    DO iz=-1,nz+2
      DO iy=-1,ny+2
        DO ix=-1,nx+2
          temperature(ix,iy,iz) = (gamma - 1.0_num) / (2.0_num - xi_n(ix,iy,iz)) &
              * (Y(0,ix,iy,iz) - (1.0_num - xi_n(ix,iy,iz)) * ionise_pot)
        ENDDO
      ENDDO
    ENDDO
    CALL heat_flux(temperature,flux)
    DO iz=1,nz
      DO iy=1,ny
        DO ix=1,nx
          ! Store L^c(Y_0) to use in stages 2-s
          Lc_Y0(ix,iy,iz) = (1.0_num / rho(ix,iy,iz)) * flux(ix,iy,iz)
          c0 =  mu_tilde(1) * dt * Lc_Y0(ix,iy,iz)
          Y(1,ix,iy,iz) = Y(0,ix,iy,iz) + c0
        ENDDO
      ENDDO
    ENDDO

    Y(2,1:nx,1:ny,1:nz) = Y(1,1:nx,1:ny,1:nz)
    Y(1,1:nx,1:ny,1:nz) = Y(0,1:nx,1:ny,1:nz)

    DO j=2,n_s_stages
      DO iz=-1,nz+2
        DO iy=-1,ny+2
          DO ix=-1,nx+2
            temperature(ix,iy,iz) = (gamma - 1.0_num) / &
                (2.0_num - xi_n(ix,iy,iz)) &
                * (Y(2,ix,iy,iz) - (1.0_num - &
                xi_n(ix,iy,iz)) * ionise_pot)
          ENDDO
        ENDDO
      ENDDO

      CALL temperature_bcs
      CALL heat_flux(temperature,flux)

      DO iz=1,nz
        DO iy=1,ny
          DO ix=1,nx
            c0 = gamma_tilde(j) * dt * Lc_Y0(ix,iy,iz)
            Lc_Yj_1 = (1.0_num / rho(ix,iy,iz)) * flux(ix,iy,iz)
            c1 = mu_tilde(j) * dt * Lc_Yj_1
          !Y_j
            Y(3,ix,iy,iz) = mu(j)*Y(2,ix,iy,iz) + nu(j)*Y(1,ix,iy,iz) &
                + (1.0_num -mu(j)-nu(j))* Y(0,ix,iy,iz)
            Y(3,ix,iy,iz) = Y(3,ix,iy,iz) + c1 + c0
          ENDDO
        ENDDO
      ENDDO

      IF ( j .LT. n_s_stages ) THEN
        !for jth stage  Y=[Y_0, Y_j-2, Y_j-1 , Y_j]
        Y(1,:,:,:) = Y(2,:,:,:)
        !This is not ideal, but it allows you to not have special boundary conditions
        energy(:,:,:) = Y(3,:,:,:)
        CALL energy_bcs
        IF (eos_number /= EOS_IDEAL) CALL neutral_fraction
        Y(2,:,:,:) = energy
        Y(3,1:nx,1:ny,1:nz) = 0.0_num
      END IF
    ENDDO
    energy(1:nx,1:ny,1:nz) = Y(3,1:nx,1:ny,1:nz)
    CALL energy_bcs

    DEALLOCATE(Y)
    DEALLOCATE(flux)
    DEALLOCATE(Lc_Y0)

  END SUBROUTINE heat_conduct_sts2




END MODULE conduct
