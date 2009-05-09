MODULE conduct

  USE shared_data
  USE boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: conduct_heat

CONTAINS

  SUBROUTINE conduct_heat

    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: kx, ky, kz, ux, uy, uz
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: energy0, energytot
    REAL(num) :: B, bxc, byc, bzc
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

    REAL(num) :: kxx, kyx, kzx
    REAL(num) :: kxy, kyy, kzy
    REAL(num) :: kxz, kyz, kzz

    REAL(num) :: uxx, uyy, uzz

    REAL(num) :: A1, A2, mx, Q, errtot, mx1
    REAL(num) :: w = 1.5_num

    INTEGER :: CYCLE, sweep, start_index

    LOGICAL :: converged

    mx = 0.0_num

    ALLOCATE(kx(0:nx+1, 0:ny+1, 0:nz+1))
    ALLOCATE(ky(0:nx+1, 0:ny+1, 0:nz+1))
    ALLOCATE(kz(0:nx+1, 0:ny+1, 0:nz+1))
    ALLOCATE(ux(0:nx+1, 0:ny+1, 0:nz+1))
    ALLOCATE(uy(0:nx+1, 0:ny+1, 0:nz+1))
    ALLOCATE(uz(0:nx+1, 0:ny+1, 0:nz+1))
    ALLOCATE(energy0(-1:nx+1, -1:ny+1, -1:nz+1))
    ALLOCATE(energytot(-1:nx+1, -1:ny+1, -1:nz+1))

    energy0 = energy
    energytot = (gamma - 1.0_num)
    converged = .FALSE.
    DO iz = 0, nz+1
      DO iy = 0, ny+1
        DO ix = 0, nx+1
          bxc = (bx(ix, iy, iz) + bx(ix-1, iy, iz))
          byc = (by(ix, iy, iz) + by(ix, iy-1, iz))
          bzc = (bz(ix, iy, iz) + bz(ix, iy, iz-1))
          B = SQRT(bxc**2 + byc**2 + bzc**2)

          IF (B .GT. 1.0e-6_num) THEN
            ux(ix, iy, iz) = bxc / B
            uy(ix, iy, iz) = byc / B
            uz(ix, iy, iz) = bzc / B
          ELSE
            ux(ix, iy, iz) = 1.0_num
            uy(ix, iy, iz) = 1.0_num
            uz(ix, iy, iz) = 1.0_num
          END IF

          kx(ix, iy, iz) = ux(ix, iy, iz) * kappa &
              * (energytot(ix, iy, iz) * energy0(ix, iy, iz))**pow
          ky(ix, iy, iz) = uy(ix, iy, iz) * kappa &
              * (energytot(ix, iy, iz) * energy0(ix, iy, iz))**pow
          kz(ix, iy, iz) = uz(ix, iy, iz) * kappa &
              * (energytot(ix, iy, iz) * energy0(ix, iy, iz))**pow
        END DO
      END DO
    END DO

    DO CYCLE = 0, 200
      errtot = 0.0_num
      DO sweep = 0, 1
        mx = 0.0_num
        mx1 = 0.0_num
        DO iz = 1, nz
          DO iy = 1, ny
            start_index = MOD(iz + iy, 2) + 2 - sweep
            !DEC$ IVDEP
            !DEC$ VECTOR ALWAYS
            DO ix = start_index, nx, 2

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

              rx = qpx * energytot(ix+1, iy, iz) * energy(ix+1, iy, iz) &
                  - qmx * energytot(ix-1, iy, iz) * energy(ix-1, iy, iz)
              ry = qpy * energytot(ix, iy+1, iz) * energy(ix, iy+1, iz) &
                  - qmy * energytot(ix, iy-1, iz) * energy(ix, iy-1, iz)
              rz = qpz * energytot(ix, iy, iz+1) * energy(ix, iy, iz+1) &
                  - qmz * energytot(ix, iy, iz-1) * energy(ix, iy, iz-1)

              rxx = mpx * energytot(ix+1, iy, iz) * energy(ix+1, iy, iz) &
                  + mmx * energytot(ix-1, iy, iz) * energy(ix-1, iy, iz)
              ryy = mpy * energytot(ix, iy+1, iz) * energy(ix, iy+1, iz) &
                  + mmy * energytot(ix, iy-1, iz) * energy(ix, iy-1, iz)
              rzz = mpz * energytot(ix, iy, iz+1) * energy(ix, iy, iz+1) &
                  + mmz * energytot(ix, iy, iz-1) * energy(ix, iy, iz-1)

              rxy = qpy / 16.0_num &
                  * (qpx * energytot(ix+1, iy+1, iz) * energy(ix+1, iy+1, iz) &
                  - qmx * energytot(ix-1, iy+1, iz) * energy(ix-1, iy+1, iz) &
                  + q0x * energytot(ix, iy+1, iz) * energy(ix, iy+1, iz)) &
                  - qmy / 16.0_num &
                  * (qpx * energytot(ix+1, iy-1, iz) * energy(ix+1, iy-1, iz) &
                  - qmx * energytot(ix-1, iy-1, iz) * energy(ix-1, iy-1, iz) &
                  + q0x * energytot(ix, iy-1, iz) * energy(ix, iy-1, iz)) &
                  + q0y / 16.0_num &
                  * (qpx * energytot(ix+1, iy, iz) * energy(ix+1, iy, iz) &
                  - qmx * energytot(ix-1, iy, iz) * energy(ix-1, iy, iz))

              rxz = qpz / 16.0_num &
                  * (qpx * energytot(ix+1, iy, iz+1) * energy(ix+1, iy, iz+1) &
                  - qmx * energytot(ix-1, iy, iz+1) * energy(ix-1, iy, iz+1) &
                  + q0x * energytot(ix, iy, iz+1) * energy(ix, iy, iz+1)) &
                  - qmz / 16.0_num &
                  * (qpx * energytot(ix+1, iy, iz-1) * energy(ix+1, iy, iz-1) &
                  - qmx * energytot(ix-1, iy, iz-1) * energy(ix-1, iy, iz-1) &
                  + q0x * energytot(ix, iy, iz-1) * energy(ix, iy, iz-1)) &
                  + q0z / 16.0_num &
                  * (qpx * energytot(ix+1, iy, iz) * energy(ix+1, iy, iz) &
                  - qmx * energytot(ix-1, iy, iz) * energy(ix-1, iy, iz))

              ryz = qpz / 16.0_num &
                  * (qpy * energytot(ix, iy+1, iz+1) * energy(ix, iy+1, iz+1) &
                  - qmy * energytot(ix, iy-1, iz+1) * energy(ix, iy-1, iz+1) &
                  + q0y * energytot(ix, iy, iz+1) * energy(ix, iy, iz+1)) &
                  - qmz / 16.0_num &
                  * (qpy * energytot(ix, iy+1, iz-1) * energy(ix, iy+1, iz-1) &
                  - qmy * energytot(ix, iy-1, iz-1) * energy(ix, iy-1, iz-1) &
                  + q0y * energytot(ix, iy, iz-1) * energy(ix, iy, iz-1)) &
                  + q0z / 16.0_num &
                  * (qpy * energytot(ix, iy+1, iz) * energy(ix, iy+1, iz) &
                  - qmy * energytot(ix, iy-1, iz) * energy(ix, iy-1, iz))

              kxx = qpx * kx(ix+1, iy, iz) &
                  - qmx * kx(ix-1, iy, iz) + q0x * kx(ix, iy, iz)
              kyx = qpx * ky(ix+1, iy, iz) &
                  - qmx * ky(ix-1, iy, iz) + q0x * ky(ix, iy, iz)
              kzx = qpx * kz(ix+1, iy, iz) &
                  - qmx * kz(ix-1, iy, iz) + q0x * kz(ix, iy, iz)

              kxy = qpy * kx(ix, iy+1, iz) &
                  - qmy * kx(ix, iy-1, iz) + q0y * kx(ix, iy, iz)
              kyy = qpy * ky(ix, iy+1, iz) &
                  - qmy * ky(ix, iy-1, iz) + q0y * ky(ix, iy, iz)
              kzy = qpy * kz(ix, iy+1, iz) &
                  - qmy * kz(ix, iy-1, iz) + q0y * kz(ix, iy, iz)

              kxz = qpz * kx(ix, iy, iz+1) &
                  - qmz * kx(ix, iy, iz-1) + q0z * kx(ix, iy, iz)
              kyz = qpz * ky(ix, iy, iz+1) &
                  - qmz * ky(ix, iy, iz-1) + q0z * ky(ix, iy, iz)
              kzz = qpz * kz(ix, iy, iz+1) &
                  - qmz * kz(ix, iy, iz-1) + q0z * kz(ix, iy, iz)

              uxx = qpx * ux(ix+1, iy, iz) &
                  - qmx * ux(ix-1, iy, iz) + q0x * ux(ix, iy, iz)
              uyy = qpy * uy(ix, iy+1, iz) &
                  - qmy * uy(ix, iy-1, iz) + q0y * uy(ix, iy, iz)
              uzz = qpz * uz(ix, iy, iz+1) &
                  - qmz * uz(ix, iy, iz-1) + q0z * uz(ix, iy, iz)

              bxc = (bx(ix, iy, iz) + bx(ix-1, iy, iz))
              byc = (by(ix, iy, iz) + by(ix, iy-1, iz))
              bzc = (bz(ix, iy, iz) + bz(ix, iy, iz-1))

              B = SQRT(bxc**2 + byc**2 + bzc**2)
              IF (B .GT. 1.0e-6_num) THEN

                ! Second differentials in T
                A1 = m0x * ux(ix, iy, iz) * kx(ix, iy, iz) &
                    + m0y * uy(ix, iy, iz) * ky(ix, iy, iz) &
                    + m0z * uz(ix, iy, iz) * kz(ix, iy, iz) &
                    + q0x * q0y / 16.0_num * (ux(ix, iy, iz) * ky(ix, iy, iz) &
                    + uy(ix, iy, iz) * kx(ix, iy, iz)) &
                    + q0x * q0z / 16.0_num * (ux(ix, iy, iz) * kz(ix, iy, iz) &
                    + uz(ix, iy, iz) * kx(ix, iy, iz)) &
                    + q0y * q0z / 16.0_num * (uy(ix, iy, iz) * kz(ix, iy, iz) &
                    + uz(ix, iy, iz) * ky(ix, iy, iz))

                ! Differentials in kx, ky, kz
                A1 = A1 &
                    + ux(ix, iy, iz) * (kxx * q0x + kyx * q0y + kzx * q0z) &
                    + uy(ix, iy, iz) * (kxy * q0x + kyy * q0y + kzy * q0z) &
                    + uz(ix, iy, iz) * (kxz * q0x + kyz * q0y + kzz * q0z)

                ! Differentials in ux, uy, uz
                A1 = A1 + q0x * kx(ix, iy, iz) * (uxx + uyy + uzz) &
                    + q0y * ky(ix, iy, iz) *(uxx + uyy + uzz) &
                    + q0z * kz(ix, iy, iz) *(uxx + uyy + uzz)

                ! Second differentials in T
                A2 = rxx * ux(ix, iy, iz) * kx(ix, iy, iz) &
                    + ryy * uy(ix, iy, iz) * ky(ix, iy, iz) &
                    + rzz * uz(ix, iy, iz) * kz(ix, iy, iz) &
                    + rxy * (ux(ix, iy, iz) * ky(ix, iy, iz) &
                    + uy(ix, iy, iz) * kx(ix, iy, iz)) &
                    + rxz * (ux(ix, iy, iz) * kz(ix, iy, iz) &
                    + uz(ix, iy, iz) * kx(ix, iy, iz)) &
                    + ryz * (uy(ix, iy, iz) * kz(ix, iy, iz) &
                    + uz(ix, iy, iz) * ky(ix, iy, iz))

                ! Differentials in kx, ky, kz
                A2 = A2 + ux(ix, iy, iz) * (kxx * rx + kyx * ry + kzx * rz) &
                    + uy(ix, iy, iz) * (kxy * rx + kyy * ry + kzy * rz) &
                    + uz(ix, iy, iz) * (kxz * rx + kyz * ry + kzz * rz)

                ! Differentials in ux, uy, uz
                A2 = A2 + uxx * (kx(ix, iy, iz) * rx &
                    + ky(ix, iy, iz) * ry + kz(ix, iy, iz) * rz) &
                    + uyy * (kx(ix, iy, iz) * rx &
                    + ky(ix, iy, iz) * ry + kz(ix, iy, iz) * rz) &
                    + uzz * (kx(ix, iy, iz) * rx &
                    + ky(ix, iy, iz) * ry + kz(ix, iy, iz) * rz)
              ELSE
                ! Isotropic heat conduction with Braginskii conduction
                ! coefficient
                A1 = kx(ix, iy, iz) * (m0x + m0y + m0z) &
                    + kxx * q0x + kyy * q0y + kzz * q0z
                A2 = kx(ix, iy, iz) * rxx + ky(ix, iy, iz) * ryy &
                    + kz(ix, iy, iz) * rzz &
                    + kxx * rx + kyy * ry + kzz * rz
              END IF

              Q = energy(ix, iy, iz)
              energy(ix, iy, iz) = (1.0_num - w) * energy(ix, iy, iz) &
                  + w / (1.0_num + A1 * dt / rho(ix, iy, iz)) &
                  * (energy0(ix, iy, iz) + A2 * dt / rho(ix, iy, iz))
              Q = Q - energy(ix, iy, iz)

              errtot = errtot + ABS(Q)
            END DO
          END DO
        END DO
        CALL energy_bcs
      END DO
      CALL MPI_ALLREDUCE(errtot, mx, 1, mpireal, MPI_SUM, comm, errcode)
      errtot = mx
      IF (errtot .LT. 1e-6_num) THEN
        converged = .TRUE.
        EXIT
      END IF
    END DO

    IF (.NOT. converged .AND. rank == 0) &
        PRINT *, "***WARNING*** Solution failed to converge during ", &
                 "heat conduction"

    DEALLOCATE(kx, ky, kz)
    DEALLOCATE(ux, uy, uz)
    DEALLOCATE(energy0, energytot)

  END SUBROUTINE conduct_heat

END MODULE conduct
