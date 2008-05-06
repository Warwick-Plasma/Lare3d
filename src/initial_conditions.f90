MODULE initial_conditions

  USE shared_data
  USE strings

  IMPLICIT NONE

  PRIVATE

  PUBLIC:: HandleCustomBlock,Equilibrium


CONTAINS

  !-----------------------------------------------------------------------------
  !This function contains the equilibrium
  !-----------------------------------------------------------------------------

  SUBROUTINE Equilibrium
    INTEGER :: ix, iy, iz
    REAL(num) :: p_dropoff=5.0_num

    !Gaussian pressure excess
    grav=0.0_num
    rho=1.0_num
    energy=1.0_num
    bx=0.0_num
    by=1.0_num
    bz=0.0_num
    vx=0.0_num
    vy=0.0_num
    vz=0.0_num


    DO iz=1,nz
       DO iy=1,ny
          DO ix=1,nx
             energy(ix,iy,iz)=energy(ix,iy,iz)*(1.0_num + 10.0_num*EXP(-(xc(ix)**2+yc(iy)**2+zc(iz)**2)/p_dropoff**2))
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE Equilibrium

  !-----------------------------------------------------------------------------
  !These functions contain the user input deck elements
  !-----------------------------------------------------------------------------

  FUNCTION HandleCustomBlock(blockname,Element,Value)

    CHARACTER(len=30),INTENT(IN)::blockname,Element,Value
    INTEGER :: HandleCustomBlock
    LOGICAL :: Result

    !The following line must always be present
    HandleCustomBlock=ERR_UNKNOWN_BLOCK
  
  END FUNCTION HandleCustomBlock


END MODULE initial_conditions
