!**************************************************************
! This module contains the conversion factors between normalised
! and internal code units
!**************************************************************

MODULE normalise
  USE constants
  USE shared_data
  IMPLICIT NONE
  
CONTAINS

  
  SUBROUTINE normalise_transport
  
    REAL(num) :: eta_bar_0, pressure0, v0, energy0, temp0, mbar, kappa0
  
    pressure0 = B0**2 / mu0_si  
    energy0 = pressure0 / rho0  
    v0 = SQRT(energy0) 
    mbar = mf * mh_si
    temp0 = (mbar / kb_si) * energy0  
    
    ! Normalise tbar, r_bar and eta_bar for including Cowling resistivity and neutrals
    t_bar = t_bar / temp0
    r_bar = r_bar * rho0 / temp0**(3.0_num / 2.0_num) 
    eta_bar_0 = (mu0_si * L0 * v0) * rho0**2 * SQRT(temp0) / B0**2 
    eta_bar = eta_bar / eta_bar_0
  
    ! Normalise ionise_pot 
    ionise_pot = ionise_pot_si / (energy0 * mbar)   
    IF (eos_number /= EOS_ION) ionise_pot = 0.0_num
         
    ! normalise tr required for get_neutral etc.
    tr = tr / temp0                             
    
    ! normalise paralle thermal conductivity
    kappa0 = energy0**(3.0_num / 2.0_num) * rho0 * L0 &
        / (mbar / kb_si * energy0)**(7.0_num / 2.0_num)	
  	kappa_0 = 1.e-11_num / kappa0   
  	
  	!find the normalised temperature corresponding to 100MK
  	temperature_100mk = 1.e8_num / temp0      
  	
		! factor required to convert between normalised energy and normalised temperature 
		! only distinguishes between fully ioned and neutral. Ok as only needed in conduction. 
		e2t = (gamma - 1.0_num) / 2.0_num
		IF (eos_number == EOS_IDEAL .AND. neutral_gas) e2t = e2t * 2.0_num

   	!convertion factor to get temperature in MK from normalised energy
    e2tmk = temp0 * e2t / 1.e6_num     
     
    h_star = L0 / (rho0 * v0**3)   
     
    lr_star = 1.148e-35_num * (rho0 / mbar)**2
  
  END SUBROUTINE normalise_transport
  


END MODULE normalise
