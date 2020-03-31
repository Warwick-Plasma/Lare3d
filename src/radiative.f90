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

MODULE radiative

  USE shared_data
  USE boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rad_losses, user_defined_heating

  INTEGER, PARAMETER :: n = 7   !set the number of temperature boundaries used in Q(T)
  INTEGER, PARAMETER  :: kmax = n - 1
  REAL(num), DIMENSION(:), ALLOCATABLE :: t_boundary, pow, psi
  INTEGER, DIMENSION(:), ALLOCATABLE :: alpha
  REAL(num), DIMENSION(:), ALLOCATABLE :: yk, qk, ratios
  REAL(num) :: cool, ratio

CONTAINS



  SUBROUTINE setup_loss_function
    ! In this subroutine specify the radiative loss Q(T) = psi_k * T^alpha_k
    ! This must be piecewise polynomial with kmax=n-1 regions
    ! bounded by n temperatures. 
    ! Q(T) must be in S.I.
    ! In the energy equation this would appear as a pressure cooling through
    ! dp/dt = -(gamma-1) n_e n_H Q(T)
    ! Only in this form for T~>10^k so fully ionised and n_H=n_e

    REAL(num) :: frac

    t_boundary(:) = (/0.02_num, 0.0398_num, 0.0794_num, 0.251_num, 0.562_num, 1.995_num, 10.0_num/) * 1e6_num

    !Define power for polynomial fit.
    !Alpha is defined as integer but RTV has one none
    !integer power so have integer alpha array and real pow array defined below.
    !The interger array is needed as alpha=1 is a special case.
    frac = -2.0_num / 3.0_num
    alpha(:) = (/0, 2, 0, -2, 0, 0/)
    pow(:) = REAL(alpha, num) + (/0.0_num, 0.0_num, 0.0_num, 0.0_num, 0.0_num, frac/)

    !Usually specify RTV in cgs
    psi(:) = (/-21.85_num, -31.0_num, -21.2_num, -10.4_num, -21.94_num, -17.73_num/)
    psi = 10**psi
    !Convert to SI
    psi = 1.e-13_num * psi

  END SUBROUTINE setup_loss_function




  SUBROUTINE user_defined_heating
    ! Use this to define any heating needed on each timestep
    ! This is a dumb example so change to the heating needed
    ! This example specifies the heating rate (heat_in) in S.I. then
    ! converts to Lare units

    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: heat_in

    ALLOCATE (heat_in(-1:nx+2, -1:ny+2, -1:nz+2))

    ! Specify heating in S.I. units W/m^3
    DO iz = 1, nz
    DO iy = 1, ny
      DO ix = 1, nx
        heat_in(ix,iy,iz) = 0.1_num * EXP(-xc(ix)**2 -yc(iy)**2 - zc(iz)**2)
      END DO
    END DO
    END DO
    ! Convert to internal Lare units
    heat_in(:,:,:) = heat_in(:,:,:) * time_norm**3 / (rho_norm * L_norm)

    energy(:,:,:) = energy(:,:,:) + heat_in(:,:,:) * dt / rho(:,:,:)
    CALL energy_bcs

  END SUBROUTINE user_defined_heating




  SUBROUTINE rad_losses

    LOGICAL :: first_call = .TRUE.

    IF (first_call) THEN
      first_call = .FALSE.
      ALLOCATE (t_boundary(1:n), alpha(1:kmax), pow(1:kmax), psi(1:kmax))
      ALLOCATE (yk(1:n), qk(1:n), ratios(1:n))
      CALL setup_loss_function
      CALL set_exact_integration_arrays
    END IF

    CALL exact_intergation_method
    CALL energy_bcs

  END SUBROUTINE rad_losses



  SUBROUTINE exact_intergation_method

    REAL(num) :: temp_si, inverse_t_cool, yt, fac
    INTEGER :: i, k

    DO ix = 1, nx
      DO iy = 1, ny
        DO iz = 1, nz

          temp_si = 0.5_num * temp_norm * energy(ix,iy,iz) * (gamma - 1.0_num) 
                
          k = -1
          DO i = 1, kmax
            IF (temp_si > t_boundary(i) .AND. temp_si <= t_boundary(i+1)) THEN
              k = i
              EXIT
            END IF
          END DO
          IF (k .LT. 0) CYCLE
  
          IF (alpha(k) .NE. 1) THEN
            fac = 1.0_num / (1.0 - pow(k))
            yt = yk(k) + fac * ratios(k) * (1.0_num - (t_boundary(k)/temp_si)**(pow(k)-1.0_num))
          ELSE
            yt = yk(k) + ratios(k) * LOG((t_boundary(k)/temp_si))
          END IF
  
          inverse_t_cool = cool * rho(ix,iy,iz)
          yt = yt +  dt * inverse_t_cool *  time_norm
  
          k = -1
          DO i = 1, kmax
            IF ((yt > yk(i) .AND. yt <= yk(i+1)) .OR. (yt < yk(i) .AND. yt >= yk(i+1))) THEN
              k = i
              EXIT
            END IF
          END DO
          IF (k .LT. 0) THEN
            temp_si = t_boundary(1)
            energy(ix,iy,iz) = temp_si * 2.0_num / (gamma - 1.0_num) / temp_norm
            CYCLE
          ENDIF
  
          IF (alpha(k) .NE. 1) THEN
            fac = 1.0_num / (1.0 - pow(k))
            temp_si = t_boundary(k) * (1.0_num - (1.0_num - pow(k)) / ratios(k) * (yt - yk(k)))**fac
          ELSE
            temp_si = t_boundary(k) * EXP((yt - yk(k)) / ratios(k))
          END IF
  
          energy(ix,iy,iz) = temp_si * 2.0_num / (gamma - 1.0_num) / temp_norm

        END DO
      END DO
    END DO 


  END SUBROUTINE exact_intergation_method




  SUBROUTINE set_exact_integration_arrays
    ! Define arrays used in Townsend exact integration method
    ! Here qk = Lambda_k from Townsend

    REAL(num) :: fac
    INTEGER :: k

    DO k = 1, kmax
      qk(k) = psi(k) * t_boundary(k)**pow(k)
    END DO
    qk(n) = qk(n-1) * (t_boundary(n)/t_boundary(n-1))**pow(n-1)

    ratio = qk(n) / t_boundary(n)
    cool = 0.5_num * ratio * (gamma - 1.0_num) * rho_norm / (kb_si * mf * mh_si)
    DO k = 1, n
      ratios(k) = ratio * (t_boundary(k) / qk(k))
    END DO

    yk(n) = 0.0_num
    DO k = kmax, 1, -1
      IF (alpha(k) .NE. 1) THEN
        fac = 1.0_num / (1.0 - pow(k))
        yk(k) = yk(k+1) - fac * ratios(k) * (1.0_num - (t_boundary(k)/t_boundary(k+1))**(pow(k)-1.0_num))
      ELSE
        yk(k) = yk(k+1) - ratios(k) * LOG((t_boundary(k)/t_boundary(k+1)))
      END IF
    END DO
   

  END SUBROUTINE set_exact_integration_arrays




END MODULE radiative
