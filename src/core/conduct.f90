MODULE conduct

  USE shared_data
  USE boundary
	USE eos

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: conduct_heat

CONTAINS

  SUBROUTINE conduct_heat

    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: kx, ky, kz, ux, uy, uz
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: energy0, e2temp, kp
		REAL(num) :: e, T
    REAL(num) :: b, bxc, byc, bzc
    REAL(num) :: pow = 5.0_num / 2.0_num
    REAL(num) :: qpx, qmx, q0x
    REAL(num) :: qpy, qmy, q0y
    REAL(num) :: qpz, qmz, q0z
    REAL(num) :: mpx, mmx, m0x
    REAL(num) :: mpy, mmy, m0y
    REAL(num) :: mpz, mmz, m0z
    REAL(num) :: rx, ry, rz
    REAL(num) :: rxx, ryy, rzz
    REAL(num) :: rxy, rxz, ryz
    REAL(num) :: kxx, kyx, kzx, kpx
    REAL(num) :: kxy, kyy, kzy, kpy
    REAL(num) :: kxz, kyz, kzz, kpz
    REAL(num) :: uxx, uyy, uzz
    REAL(num) :: a1, a2, error, q, errmax, errmax_prev = 0.0_num 
    REAL(num) :: w 

    INTEGER :: loop

    LOGICAL :: converged
    REAL(num), PARAMETER :: fractional_error = 1.0e-5_num
    REAL(num), PARAMETER :: b_min = 1.0e-4_num

    ALLOCATE(kx(0:nx+1, 0:ny+1, 0:nz+1))
    ALLOCATE(ky(0:nx+1, 0:ny+1, 0:nz+1))
    ALLOCATE(kz(0:nx+1, 0:ny+1, 0:nz+1))
    ALLOCATE(ux(0:nx+1, 0:ny+1, 0:nz+1))
    ALLOCATE(uy(0:nx+1, 0:ny+1, 0:nz+1))
    ALLOCATE(uz(0:nx+1, 0:ny+1, 0:nz+1))
    ALLOCATE(energy0(-1:nx+2, -1:ny+2, -1:nz+2))
    ALLOCATE(e2temp(0:nx+1, 0:ny+1, 0:nz+1))
    ALLOCATE(kp(0:nx+1, 0:ny+1, 0:nz+1))

		DO iz = 0, nz+1
			DO iy = 0, ny+1
				DO ix = 0, nx+1
					e = energy(ix, iy, iz)
					CALL get_temp(rho(ix, iy, iz), e, eos_number, ix, iy, iz, T)
					e2temp(ix, iy, iz) = T / e
				END DO
			END DO
		END DO

    DO iz = 0, nz+1
      DO iy = 0, ny+1
        DO ix = 0, nx+1
          bxc = (bx(ix, iy, iz) + bx(ix-1, iy, iz))
          byc = (by(ix, iy, iz) + by(ix, iy-1, iz))
          bzc = (bz(ix, iy, iz) + bz(ix, iy, iz-1))
          b = SQRT(bxc**2 + byc**2 + bzc**2)

          ux(ix, iy, iz) = bxc / SQRT(b**2 + b_min**2)
          uy(ix, iy, iz) = byc / SQRT(b**2 + b_min**2)
          uz(ix, iy, iz) = bzc / SQRT(b**2 + b_min**2)
                                
					T = (e2temp(ix, iy, iz) * energy(ix, iy, iz))**pow 
          kx(ix, iy, iz) = ux(ix, iy, iz) * kappa_0 * T
          ky(ix, iy, iz) = uy(ix, iy, iz) * kappa_0 * T
          kz(ix, iy, iz) = uz(ix, iy, iz) * kappa_0 * T 

					kp(ix, iy, iz) = kappa_0 * T * b_min**2 / (b**2 + b_min**2)
        END DO
      END DO
    END DO

    converged = .FALSE.
		w = 1.9_num
    energy0 = energy
    DO loop = 0, 100
      errmax = 0.0_num
			error = 0.0_num
      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            qpx = dxc(ix-1) / (dxc(ix) * (dxc(ix) + dxc(ix-1)))
            qmx = dxc(ix) / (dxc(ix-1) * (dxc(ix) + dxc(ix-1)))
            q0x = (dxc(ix)**2 - dxc(ix-1)**2) &
                  / (dxc(ix) * dxc(ix-1) * (dxc(ix) + dxc(ix-1)))

            qpy = dyc(iy-1) / (dyc(iy) * (dyc(iy) + dyc(iy-1)))
            qmy = dyc(iy) / (dyc(iy-1) * (dyc(iy) + dyc(iy-1)))
            q0y = (dyc(iy)**2 - dyc(iy-1)**2) &
                  / (dyc(iy) * dyc(iy-1) * (dyc(iy) + dyc(iy-1)))

            qpz = dzc(iz-1) / (dzc(iz) * (dzc(iz) + dzc(iz-1)))
            qmz = dzc(iz) / (dzc(iz-1) * (dzc(iz) + dzc(iz-1)))
            q0z = (dzc(iz)**2 - dzc(iz-1)**2) &
                  / (dzc(iz) * dzc(iz-1) * (dzc(iz) + dzc(iz-1)))

            mpx = 1.0_num / (dxc(ix) * dxb(ix))
	          mmx = 1.0_num / (dxc(ix-1) * dxb(ix))
            m0x = (dxc(ix) + dxc(ix-1)) / (dxc(ix) * dxc(ix-1) * dxb(ix))

            mpy = 1.0_num / (dyc(iy) * dyb(iy))
            mmy = 1.0_num / (dyc(iy-1) * dyb(iy))
            m0y = (dyc(iy) + dyc(iy-1)) / (dyc(iy) * dyc(iy-1) * dyb(iy))

            mpz = 1.0_num / (dzc(iz) * dzb(iz))
            mmz = 1.0_num / (dzc(iz-1) * dzb(iz))
            m0z = (dzc(iz) + dzc(iz-1)) / (dzc(iz) * dzc(iz-1) * dzb(iz))

            rx = qpx * e2temp(ix+1, iy, iz) * energy(ix+1, iy, iz) &
                  - qmx * e2temp(ix-1, iy, iz) * energy(ix-1, iy, iz)
            ry = qpy * e2temp(ix, iy+1, iz) * energy(ix, iy+1, iz) &
                  - qmy * e2temp(ix, iy-1, iz) * energy(ix, iy-1, iz)
            rz = qpz * e2temp(ix, iy, iz+1) * energy(ix, iy, iz+1) &
                  - qmz * e2temp(ix, iy, iz-1) * energy(ix, iy, iz-1)

            rxx = mpx * e2temp(ix+1, iy, iz) * energy(ix+1, iy, iz) &
                  + mmx * e2temp(ix-1, iy, iz) * energy(ix-1, iy, iz)
            ryy = mpy * e2temp(ix, iy+1, iz) * energy(ix, iy+1, iz) &
                  + mmy * e2temp(ix, iy-1, iz) * energy(ix, iy-1, iz)
            rzz = mpz * e2temp(ix, iy, iz+1) * energy(ix, iy, iz+1) &
                  + mmz * e2temp(ix, iy, iz-1) * energy(ix, iy, iz-1)

            rxy = qpy  &
                  * (qpx * e2temp(ix+1, iy+1, iz) * energy(ix+1, iy+1, iz) &
                  - qmx * e2temp(ix-1, iy+1, iz) * energy(ix-1, iy+1, iz) &
                  + q0x * e2temp(ix, iy+1, iz) * energy(ix, iy+1, iz)) &
                - qmy  &
                  * (qpx * e2temp(ix+1, iy-1, iz) * energy(ix+1, iy-1, iz) &
                  - qmx * e2temp(ix-1, iy-1, iz) * energy(ix-1, iy-1, iz) &
                  + q0x * e2temp(ix, iy-1, iz) * energy(ix, iy-1, iz)) &
                + q0y  &
                  * (qpx * e2temp(ix+1, iy, iz) * energy(ix+1, iy, iz) &
                  - qmx * e2temp(ix-1, iy, iz) * energy(ix-1, iy, iz))

            rxz = qpz  &
                 * (qpx * e2temp(ix+1, iy, iz+1) * energy(ix+1, iy, iz+1) &
                 - qmx * e2temp(ix-1, iy, iz+1) * energy(ix-1, iy, iz+1) &
                 + q0x * e2temp(ix, iy, iz+1) * energy(ix, iy, iz+1)) &
               - qmz  &
                 * (qpx * e2temp(ix+1, iy, iz-1) * energy(ix+1, iy, iz-1) &
                 - qmx * e2temp(ix-1, iy, iz-1) * energy(ix-1, iy, iz-1) &
                 + q0x * e2temp(ix, iy, iz-1) * energy(ix, iy, iz-1)) &
               + q0z  &
                 * (qpx * e2temp(ix+1, iy, iz) * energy(ix+1, iy, iz) &
                 - qmx * e2temp(ix-1, iy, iz) * energy(ix-1, iy, iz))

            ryz = qpz  &
                 * (qpy * e2temp(ix, iy+1, iz+1) * energy(ix, iy+1, iz+1) &
                 - qmy * e2temp(ix, iy-1, iz+1) * energy(ix, iy-1, iz+1) &
                 + q0y * e2temp(ix, iy, iz+1) * energy(ix, iy, iz+1)) &
                - qmz &
                 * (qpy * e2temp(ix, iy+1, iz-1) * energy(ix, iy+1, iz-1) &
                 - qmy * e2temp(ix, iy-1, iz-1) * energy(ix, iy-1, iz-1) &
                 + q0y * e2temp(ix, iy, iz-1) * energy(ix, iy, iz-1)) &
                + q0z  &
                 * (qpy * e2temp(ix, iy+1, iz) * energy(ix, iy+1, iz) &
                 - qmy * e2temp(ix, iy-1, iz) * energy(ix, iy-1, iz))

            kxx = qpx * kx(ix+1, iy, iz) &
                - qmx * kx(ix-1, iy, iz) + q0x * kx(ix, iy, iz)
            kyx = qpx * ky(ix+1, iy, iz) &
                - qmx * ky(ix-1, iy, iz) + q0x * ky(ix, iy, iz)
            kzx = qpx * kz(ix+1, iy, iz) &
                - qmx * kz(ix-1, iy, iz) + q0x * kz(ix, iy, iz)
            kpx = qpx * kp(ix+1, iy, iz) &
                - qmx * kp(ix-1, iy, iz) + q0x * kp(ix, iy, iz)

            kxy = qpy * kx(ix, iy+1, iz) &
                - qmy * kx(ix, iy-1, iz) + q0y * kx(ix, iy, iz)
            kyy = qpy * ky(ix, iy+1, iz) &
                - qmy * ky(ix, iy-1, iz) + q0y * ky(ix, iy, iz)
            kzy = qpy * kz(ix, iy+1, iz) &
                - qmy * kz(ix, iy-1, iz) + q0y * kz(ix, iy, iz)
            kpy = qpy * kp(ix, iy+1, iz) &
                - qmy * kp(ix, iy-1, iz) + q0y * kp(ix, iy, iz)

            kxz = qpz * kx(ix, iy, iz+1) &
                - qmz * kx(ix, iy, iz-1) + q0z * kx(ix, iy, iz)
            kyz = qpz * ky(ix, iy, iz+1) &
                - qmz * ky(ix, iy, iz-1) + q0z * ky(ix, iy, iz)
            kzz = qpz * kz(ix, iy, iz+1) &
                - qmz * kz(ix, iy, iz-1) + q0z * kz(ix, iy, iz)
            kpz = qpz * kp(ix, iy, iz+1) &
                - qmz * kp(ix, iy, iz-1) + q0z * kp(ix, iy, iz)

            uxx = qpx * ux(ix+1, iy, iz) &
                - qmx * ux(ix-1, iy, iz) + q0x * ux(ix, iy, iz)
            uyy = qpy * uy(ix, iy+1, iz) &
                - qmy * uy(ix, iy-1, iz) + q0y * uy(ix, iy, iz)
            uzz = qpz * uz(ix, iy, iz+1) &
                - qmz * uz(ix, iy, iz-1) + q0z * uz(ix, iy, iz)

            bxc = (bx(ix, iy, iz) + bx(ix-1, iy, iz)) / 2.0_num
            byc = (by(ix, iy, iz) + by(ix, iy-1, iz)) / 2.0_num
            bzc = (bz(ix, iy, iz) + bz(ix, iy, iz-1)) / 2.0_num

            b = SQRT(bxc**2 + byc**2 + bzc**2)
            ! Second differentials in T
            a1 = m0x * ux(ix, iy, iz) * kx(ix, iy, iz) &
                + m0y * uy(ix, iy, iz) * ky(ix, iy, iz) &
                + m0z * uz(ix, iy, iz) * kz(ix, iy, iz) &
                + q0x * q0y * (ux(ix, iy, iz) * ky(ix, iy, iz) &
                + uy(ix, iy, iz) * kx(ix, iy, iz)) &
                + q0x * q0z * (ux(ix, iy, iz) * kz(ix, iy, iz) &
                + uz(ix, iy, iz) * kx(ix, iy, iz)) &
                + q0y * q0z * (uy(ix, iy, iz) * kz(ix, iy, iz) &
                + uz(ix, iy, iz) * ky(ix, iy, iz))

            ! Differentials in kx, ky, kz
            a1 = a1 &
                - ux(ix, iy, iz) * (kxx * q0x + kyx * q0y + kzx * q0z) &
                - uy(ix, iy, iz) * (kxy * q0x + kyy * q0y + kzy * q0z) &
                - uz(ix, iy, iz) * (kxz * q0x + kyz * q0y + kzz * q0z)

            ! Differentials in ux, uy, uz
            a1 = a1 - q0x * kx(ix, iy, iz) * (uxx + uyy + uzz) &
                - q0y * ky(ix, iy, iz) *(uxx + uyy + uzz) &
                - q0z * kz(ix, iy, iz) *(uxx + uyy + uzz)

            ! Second differentials in T
            a2 = rxx * ux(ix, iy, iz) * kx(ix, iy, iz) &
                + ryy * uy(ix, iy, iz) * ky(ix, iy, iz) &
                + rzz * uz(ix, iy, iz) * kz(ix, iy, iz) &
                + rxy * (ux(ix, iy, iz) * ky(ix, iy, iz) &
                + uy(ix, iy, iz) * kx(ix, iy, iz)) &
                + rxz * (ux(ix, iy, iz) * kz(ix, iy, iz) &
                + uz(ix, iy, iz) * kx(ix, iy, iz)) &
                + ryz * (uy(ix, iy, iz) * kz(ix, iy, iz) &
                + uz(ix, iy, iz) * ky(ix, iy, iz))

            ! Differentials in kx, ky, kz
            a2 = a2 + ux(ix, iy, iz) * (kxx * rx + kyx * ry + kzx * rz) &
                + uy(ix, iy, iz) * (kxy * rx + kyy * ry + kzy * rz) &
                + uz(ix, iy, iz) * (kxz * rx + kyz * ry + kzz * rz)

            ! Differentials in ux, uy, uz
            a2 = a2 + uxx * (kx(ix, iy, iz) * rx &
                + ky(ix, iy, iz) * ry + kz(ix, iy, iz) * rz) &
                + uyy * (kx(ix, iy, iz) * rx &
                + ky(ix, iy, iz) * ry + kz(ix, iy, iz) * rz) &
                + uzz * (kx(ix, iy, iz) * rx &
                + ky(ix, iy, iz) * ry + kz(ix, iy, iz) * rz) 

						! add isoptropic elements
	          a1 = a1 + kp(ix, iy, iz) * (m0x + m0y + m0z) &
	                    - kpx * q0x - kpy * q0y - kpz * q0z 
	          a2 = a2 + kp(ix, iy, iz) * (rxx + ryy + rzz) &
	                    + kpx * rx + kpy * ry	+ kpz * rz						

						a1 = a1 * dt * e2temp(ix, iy, iz) / rho(ix, iy, iz) 
						a2 = a2 * dt / rho(ix, iy, iz) 

 						Q = energy(ix, iy, iz)
 						energy(ix, iy, iz) = (1.0_num-w) * energy(ix, iy, iz) &
 					    + w * (energy0(ix, iy, iz)  + a2) / (1.0_num + a1)
 						Q = (Q - energy(ix, iy, iz)) / Q

	          errmax = MAX(errmax, ABS(Q))  
          END DO
        END DO
      END DO
      CALL energy_bcs

      CALL MPI_ALLREDUCE(errmax, error, 1, mpireal, MPI_MAX, comm, errcode)
      errmax = error

			IF (errmax .GT. errmax_prev) w = (1.0_num + w) / 2.0_num
			errmax_prev = errmax

      IF (errmax .LT. fractional_error) THEN
        converged = .TRUE.
        EXIT
      END IF
    END DO

    IF (.NOT. converged .AND. rank == 0) &
        PRINT *, "***WARNING*** Solution failed to converge during ", &
                 "heat conduction"

    DEALLOCATE(kx, ky, kz)
    DEALLOCATE(ux, uy, uz)
    DEALLOCATE(kp, energy0, e2temp)

  END SUBROUTINE conduct_heat

END MODULE conduct
