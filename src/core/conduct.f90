MODULE conduct

  USE shared_data
  USE boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: Conduct_Heat
CONTAINS

  SUBROUTINE Conduct_Heat
    REAL(num),DIMENSION(:,:,:),ALLOCATABLE :: kx,ky,kz, ux,uy,uz,energy0,EnergyToT
    REAL(num) :: B,bxc,byc,bzc
    REAL(num) :: pow=5.0_num/2.0_num
    REAL(num) :: QpX,QmX,Q0X
    REAL(num) :: QpY,QmY,Q0Y
    REAL(num) :: QpZ,QmZ,Q0Z

    REAL(num) :: MpX,MmX,M0X
    REAL(num) :: MpY,MmY,M0Y
    REAL(num) :: MpZ,MmZ,M0Z

    REAL(num) :: Rx,Ry,Rz
    REAL(num) :: Rxx,Ryy,Rzz
    REAL(num) :: Rxy,Rxz,Ryz

    REAL(num) :: KxX,KyX,KzX
    REAL(num) :: KxY,KyY,KzY
    REAL(num) :: KxZ,KyZ,KzZ

    REAL(num) :: uxX,uyY,uzZ

    REAL(num) :: A1,A2,mx,Q,errtot,mx1
    REAL(num) :: w=1.5_num

    REAL(dbl) :: runtime = 0.0_dbl, starttime,endtime

    INTEGER :: cycle,sweep,mx_x,mx_y,mx_z, start_index

    LOGICAL :: converged

    mx=0.0_num

    ALLOCATE(kx(0:nx+1,0:ny+1,0:nz+1),ky(0:nx+1,0:ny+1,0:nz+1),kz(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(ux(0:nx+1,0:ny+1,0:nz+1),uy(0:nx+1,0:ny+1,0:nz+1),uz(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(energy0(-1:nx+1,-1:ny+1,-1:nz+1),EnergyToT(-1:nx+1,-1:ny+1,-1:nz+1))

    energy0=energy
    EnergyToT=(gamma-1.0_num)
    converged=.FALSE.
    DO iz=0,nz+1
       DO iy=0,ny+1
          DO ix=0,nx+1
             bxc=(bx(ix,iy,iz)+bx(ix-1,iy,iz))
             byc=(by(ix,iy,iz)+by(ix,iy-1,iz))
             bzc=(bz(ix,iy,iz)+bz(ix,iy,iz-1))
             B=SQRT(bxc**2+byc**2+bzc**2)
             IF (B .GT. 1.0e-6_num) THEN
                ux(ix,iy,iz)=bxc/B
                uy(ix,iy,iz)=byc/B
                uz(ix,iy,iz)=bzc/B
             ELSE
                ux(ix,iy,iz)=1.0_num
                uy(ix,iy,iz)=1.0_num
                uz(ix,iy,iz)=1.0_num
             ENDIF

             kx(ix,iy,iz)=kappa * (energyToT(ix,iy,iz) * energy0(ix,iy,iz))**pow * ux(ix,iy,iz)
             ky(ix,iy,iz)=kappa * (energyToT(ix,iy,iz) * energy0(ix,iy,iz))**pow * uy(ix,iy,iz)
             kz(ix,iy,iz)=kappa * (energyToT(ix,iy,iz) * energy0(ix,iy,iz))**pow * uz(ix,iy,iz)
          ENDDO
       ENDDO
    ENDDO


    DO cycle=0,200
       errtot=0.0_num
       DO sweep=0,1
          mx=0.0_num
          mx1=0.0_num
          DO iz=1,nz
             DO iy=1,ny
                start_index = MOD(iz + iy,2) + 2 - sweep
                !DEC$ IVDEP
                !DEC$ VECTOR ALWAYS
                DO ix=start_index, nx, 2

                   QpX=dxc(ix-1)/(dxc(ix)*(dxc(ix)+dxc(ix-1)))
                   QmX=dxc(ix)/(dxc(ix-1)*(dxc(ix)+dxc(ix-1)))
                   Q0X=(dxc(ix)**2-dxc(ix-1)**2)/(dxc(ix)*dxc(ix-1)*(dxc(ix)+dxc(ix-1)))

                   QpY=dyc(iy-1)/(dyc(iy)*(dyc(iy)+dyc(iy-1)))
                   QmY=dyc(iy)/(dyc(iy-1)*(dyc(iy)+dyc(iy-1)))
                   Q0Y=(dyc(iy)**2-dyc(iy-1)**2)/(dyc(iy)*dyc(iy-1)*(dyc(iy)+dyc(iy-1)))

                   QpZ=dzc(iz-1)/(dzc(iz)*(dzc(iz)+dzc(iz-1)))
                   QmZ=dzc(iz)/(dzc(iz-1)*(dzc(iz)+dzc(iz-1)))
                   Q0Z=(dzc(iz)**2-dzc(iz-1)**2)/(dzc(iz)*dzc(iz-1)*(dzc(iz)+dzc(iz-1)))

                   MpX=1.0_num/(dxc(ix)*dxb(ix))
                   MmX=1.0_num/(dxc(ix-1)*dxb(ix))
                   M0X=(dxc(ix)+dxc(ix-1))/(dxc(ix)*dxc(ix-1)*dxb(ix))

                   MpY=1.0_num/(dyc(iy)*dyb(iy))
                   MmY=1.0_num/(dyc(iy-1)*dyb(iy))
                   M0Y=(dyc(iy)+dyc(iy-1))/(dyc(iy)*dyc(iy-1)*dyb(iy))

                   MpZ=1.0_num/(dzc(iz)*dzb(iz))
                   MmZ=1.0_num/(dzc(iz-1)*dzb(iz))
                   M0Z=(dzc(iz)+dzc(iz-1))/(dzc(iz)*dzc(iz-1)*dzb(iz))

                   Rx=QpX*EnergyToT(ix+1,iy,iz)*energy(ix+1,iy,iz) - QmX*EnergyToT(ix-1,iy,iz)*energy(ix-1,iy,iz)             
                   Ry=QpY*EnergyToT(ix,iy+1,iz)*energy(ix,iy+1,iz) - QmY*EnergyToT(ix,iy-1,iz)*energy(ix,iy-1,iz)
                   Rz=QpZ*EnergyToT(ix,iy,iz+1)*energy(ix,iy,iz+1) - QmZ*EnergyToT(ix,iy,iz-1)*energy(ix,iy,iz-1)

                   Rxx=MpX*EnergyToT(ix+1,iy,iz)*energy(ix+1,iy,iz) + MmX*EnergyToT(ix-1,iy,iz)*energy(ix-1,iy,iz)
                   Ryy=MpY*EnergyToT(ix,iy+1,iz)*energy(ix,iy+1,iz) + MmY*EnergyToT(ix,iy-1,iz)*energy(ix,iy-1,iz)
                   Rzz=MpZ*EnergyToT(ix,iy,iz+1)*energy(ix,iy,iz+1) + MmZ*EnergyToT(ix,iy,iz-1)*energy(ix,iy,iz-1)


                   Rxy=QpY/16.0_num*(QpX*EnergyToT(ix+1,iy+1,iz)*energy(ix+1,iy+1,iz) - QmX*EnergyToT(ix-1,iy+1,iz)*energy(ix-1,iy+1,iz) &
                        + Q0X*EnergyToT(ix,iy+1,iz)*energy(ix,iy+1,iz)) &
                        - QmY/16.0_num*(QpX*EnergyToT(ix+1,iy-1,iz)*energy(ix+1,iy-1,iz) - QmX*EnergyToT(ix-1,iy-1,iz)*energy(ix-1,iy-1,iz) &
                        + Q0X*EnergyToT(ix,iy-1,iz)*energy(ix,iy-1,iz)) + Q0Y/16.0_num*(QpX*EnergyToT(ix+1,iy,iz)*energy(ix+1,iy,iz)-QmX&
                        *EnergyToT(ix-1,iy,iz)*energy(ix-1,iy,iz))

                   Rxz=QpZ/16.0_num*(QpX*EnergyToT(ix+1,iy,iz+1)*energy(ix+1,iy,iz+1) - QmX*EnergyToT(ix-1,iy,iz+1)*energy(ix-1,iy,iz+1) &
                        + Q0X*EnergyToT(ix,iy,iz+1)*energy(ix,iy,iz+1)) - &
                        QmZ/16.0_num*(QpX*EnergyToT(ix+1,iy,iz-1)*energy(ix+1,iy,iz-1) - QmX*EnergyToT(ix-1,iy,iz-1)*energy(ix-1,iy,iz-1) &
                        + Q0X*EnergyToT(ix,iy,iz-1)*energy(ix,iy,iz-1)) + &
                        Q0Z/16.0_num*(QpX*EnergyToT(ix+1,iy,iz)*energy(ix+1,iy,iz)-QmX*EnergyToT(ix-1,iy,iz)*energy(ix-1,iy,iz))

                   Ryz=QpZ/16.0_num*(QpY*EnergyToT(ix,iy+1,iz+1)*energy(ix,iy+1,iz+1) - QmY*EnergyToT(ix,iy-1,iz+1)*energy(ix,iy-1,iz+1) &
                        + Q0Y*EnergyToT(ix,iy,iz+1)*energy(ix,iy,iz+1)) - &
                        QmZ/16.0_num*(QpY*EnergyToT(ix,iy+1,iz-1)*energy(ix,iy+1,iz-1) - QmY*EnergyToT(ix,iy-1,iz-1)*energy(ix,iy-1,iz-1) &
                        + Q0Y*EnergyToT(ix,iy,iz-1)*energy(ix,iy,iz-1)) + &
                        Q0Z/16.0_num*(QpY*EnergyToT(ix,iy+1,iz)*energy(ix,iy+1,iz)-QmY*EnergyToT(ix,iy-1,iz)*energy(ix,iy-1,iz))


                   KxX=QpX * kx(ix+1,iy,iz) - QmX * kx(ix-1,iy,iz) + Q0X * kx(ix,iy,iz)
                   KyX=QpX * ky(ix+1,iy,iz) - QmX * ky(ix-1,iy,iz) + Q0X * ky(ix,iy,iz)                
                   KzX=QpX * kz(ix+1,iy,iz) - QmX * kz(ix-1,iy,iz) + Q0X * kz(ix,iy,iz)

                   KxY=QpY * kx(ix,iy+1,iz) - QmY * kx(ix,iy-1,iz) + Q0Y * kx(ix,iy,iz)
                   KyY=QpY * ky(ix,iy+1,iz) - QmY * ky(ix,iy-1,iz) + Q0Y * ky(ix,iy,iz)                
                   KzY=QpY * kz(ix,iy+1,iz) - QmY * kz(ix,iy-1,iz) + Q0Y * kz(ix,iy,iz)

                   KxZ=QpZ * kx(ix,iy,iz+1) - QmZ * kx(ix,iy,iz-1) + Q0Z * kx(ix,iy,iz)
                   KyZ=QpZ * ky(ix,iy,iz+1) - QmZ * ky(ix,iy,iz-1) + Q0Z * ky(ix,iy,iz)                
                   KzZ=QpZ * kz(ix,iy,iz+1) - QmZ * kz(ix,iy,iz-1) + Q0Z * kz(ix,iy,iz)

                   uxX=QpX * ux(ix+1,iy,iz) - QmX * ux(ix-1,iy,iz) + Q0X * ux(ix,iy,iz)
                   uyY=QpY * uy(ix,iy+1,iz) - QmY * uy(ix,iy-1,iz) + Q0Y * uy(ix,iy,iz)
                   uzZ=QpZ * uz(ix,iy,iz+1) - QmZ * uz(ix,iy,iz-1) + Q0Z * uz(ix,iy,iz)

                   bxc=(bx(ix,iy,iz)+bx(ix-1,iy,iz))
                   byc=(by(ix,iy,iz)+by(ix,iy-1,iz))
                   bzc=(bz(ix,iy,iz)+bz(ix,iy,iz-1))

                   B=SQRT(bxc**2+byc**2+bzc**2)
                   IF (B .GT. 1.0e-6_num) THEN

                      !Second differentials in T
                      A1=M0X*ux(ix,iy,iz)*kx(ix,iy,iz) + M0Y*uy(ix,iy,iz)*ky(ix,iy,iz) + M0Z*uz(ix,iy,iz)*kz(ix,iy,iz) +&
                           Q0X*Q0Y/16.0_num*(ux(ix,iy,iz)*ky(ix,iy,iz)+uy(ix,iy,iz)*kx(ix,iy,iz)) + &                     
                           Q0X*Q0Z/16.0_num*(ux(ix,iy,iz)*kz(ix,iy,iz)+uz(ix,iy,iz)*kx(ix,iy,iz)) + &
                           Q0Y*Q0Z/16.0_num*(uy(ix,iy,iz)*kz(ix,iy,iz)+uz(ix,iy,iz)*ky(ix,iy,iz))
                      !Differentials in kx,ky,kz
                      A1=A1 + ux(ix,iy,iz) * (KxX*Q0X + KyX*Q0Y + KzX*Q0Z) + &
                           uy(ix,iy,iz) * (KxY*Q0X + KyY*Q0Y + KzY*Q0Z) + &
                           uz(ix,iy,iz) * (KxZ*Q0X + KyZ*Q0Y + KzZ*Q0Z)
                      !Differentials in ux,uy,uz
                      A1=A1 + Q0X*kx(ix,iy,iz)*(uxX+uyY+uzZ)&
                           + Q0Y*ky(ix,iy,iz) *(uxX+uyY+uzZ)&
                           + Q0Z*kz(ix,iy,iz) *(uxX+uyY+uzZ)

                      !Second differentials in T
                      A2=Rxx*ux(ix,iy,iz)*kx(ix,iy,iz) + Ryy*uy(ix,iy,iz)*ky(ix,iy,iz) + &
                           Rzz*uz(ix,iy,iz)*kz(ix,iy,iz) +&
                           Rxy*(ux(ix,iy,iz)*ky(ix,iy,iz)+uy(ix,iy,iz)*kx(ix,iy,iz))+&
                           Rxz*(ux(ix,iy,iz)*kz(ix,iy,iz)+uz(ix,iy,iz)*kx(ix,iy,iz))+&
                           Ryz*(uy(ix,iy,iz)*kz(ix,iy,iz)+uz(ix,iy,iz)*ky(ix,iy,iz))
                      !Differentials in kx,ky,kz
                      A2=A2 + ux(ix,iy,iz) * (KxX*Rx + KyX*Ry + KzX*Rz) + &
                           uy(ix,iy,iz) * (KxY*Rx + KyY*Ry + KzY*Rz) + &
                           uz(ix,iy,iz) * (KxZ*Rx + KyZ*Ry + KzZ*Rz)
                      !Differentials in ux,uy,uz
                      A2=A2 + uxX*(kx(ix,iy,iz)*Rx + ky(ix,iy,iz)*Ry + kz(ix,iy,iz)*Rz) + &
                           uyY*(kx(ix,iy,iz)*Rx + ky(ix,iy,iz)*Ry + kz(ix,iy,iz)*Rz) + &
                           uzZ*(kx(ix,iy,iz)*Rx + ky(ix,iy,iz)*Ry + kz(ix,iy,iz)*Rz)
                   ELSE
                      !Isotropic heat conduction with Braginskii conduction coefficient
                      A1=kx(ix,iy,iz)*(M0X+M0Y+M0Z)+KxX*Q0X+KyY*Q0Y+KzZ*Q0Z
                      A2=kx(ix,iy,iz)*Rxx + ky(ix,iy,iz)*Ryy + kz(ix,iy,iz)*Rzz + &
                           KxX*Rx + KyY*Ry + KzZ*Rz
                   ENDIF


                   Q=energy(ix,iy,iz)
                   energy(ix,iy,iz)=(1.0_num-w)*energy(ix,iy,iz) + w/(1.0_num+A1*dt/rho(ix,iy,iz))*(energy0(ix,iy,iz)+A2*dt/rho(ix,iy,iz))
                   Q=Q-energy(ix,iy,iz)


                   errtot=errtot+ABS(Q)
                ENDDO
             ENDDO
          ENDDO
          CALL energy_bcs
       ENDDO
       CALL MPI_ALLREDUCE(errtot, mx, 1, mpireal, MPI_SUM, comm, errcode)
       errtot=mx
       IF (errtot .LT. 1e-6_num) THEN
          converged=.TRUE.
          EXIT
       ENDIF
    ENDDO

    IF (.NOT. converged .AND. rank == 0) PRINT *,"***WARNING*** Solution failed to converge during heat conduction"

    DEALLOCATE(kx,ky,kz)
    DEALLOCATE(ux,uy,uz)
    DEALLOCATE(energy0,EnergyToT)

  END SUBROUTINE Conduct_Heat

END MODULE conduct
