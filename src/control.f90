MODULE control

  USE shared_data
  USE normalise

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: user_normalisation, control_variables, set_output_dumps

CONTAINS

  SUBROUTINE user_normalisation

    ! Gamma is the ratio of specific heat capacities
    gamma = 5.0_num / 3.0_num

    ! Average mass of an ion in proton masses
    ! The code assumes a single ion species with this mass
    mf = 1.0_num

    ! The equations describing the normalisation in LARE
    ! have three free parameters which must be specified by
    ! the end user. 

    ! Magnetic field normalisation in Tesla
    B0 = 0.12_num
    ! Length normalisation in m
    L0 = 150.e3_num
    ! Density normalisation in kg / m^3
    RHO0 = 2.7e-4_num

  END SUBROUTINE user_normalisation



  SUBROUTINE control_variables

    ! Set the number of gridpoints in x and y directions
    nx_global = 64
    ny_global = 64
    nz_global = 64

    ! Set the maximum number of iterations of the core solver before the code
    ! terminates. If nsteps < 0 then the code will run until t = t_end
    nsteps = -1

    ! The maximum runtime of the code
    ! If SI_Input is true then this is in seconds
    t_end = 10.0_num

    ! Shock viscosities as detailed in manual - they are dimensionless
    visc1 = 0.1_num
    visc2 = 0.5_num
    ! Real viscosity expressed as the inverse Reynolds number, i.e. the
    ! same for normalised and SI input
    visc3 = 0.0_num

    ! Set these constants to manually
    ! override the domain decomposition.
    ! If either constant is set to zero
    ! then the code will try to automatically
    ! decompose in this direction
    nprocx = 0
    nprocy = 0
    nprocz = 0

    ! The length of the domain in the x direction
    ! If SI_Input is true then this is in metres
    x_start = -75.0_num
    x_end = 75.0_num
    ! Should the x grid be stretched or uniform
    x_stretch = .FALSE.

    ! The length of the domain in the y direction
    ! If SI_Input is true then this is in metres
    y_start = -75.0_num
    y_end = 75.0_num
    ! Should the y grid be stretched of uniform
    y_stretch = .FALSE.

    z_start = -20.0_num
    z_end = 130.0_num
    z_stretch = .FALSE.

    ! Turn on or off the resistive parts of the MHD equations
    resistive_mhd = .TRUE.

    ! The background resistivity expressed as the inverse Lundquist
    ! number, i.e. the same for normalised and SI input
    eta_background = 0.0_num

    ! The critical current for triggering anomalous resistivity
    ! and the resistivity when above the critical current
    ! The resistivity is expressed as the inverse Lundquist number, i.e. the
    ! same for normalised and SI input, bit the j_max must be in SI
    ! if using SI units
    j_max = 0.0_num
    eta0 = 0.0_num

    ! Turn on or off the hall_mhd term in the MHD equations
    ! Well actually this does nothing as it isn't fully
    ! included yet! Sorry!
    hall_mhd = .FALSE.

    ! Set the ion skin depth. If SI_Input is true then
    ! This is in metres. Note that this should be fixed to
    ! (the speed of light) / (ion plasma frequency) with the
    ! plasma frequncy fixed by the chosen normalisation of
    ! density (through the normalising mass density). Here
    ! it is treated as a free pararmeter so be careful!
    lambda0 = 0.0_num

    ! Turn on or off the Braginskii thermal conduction term in
    ! the MHD equations   
    ! WARNING: this is not robust. It is known to have problems 
    ! with steep temperature gradients and very hot regions with
    ! large thermal conductivity. For many problems it is however
    ! fine.     
    conduction = .FALSE. 
    ! Apply a flux limiter to stop heat flows exceeding free streaming limit
    ! this is an experimental feature
    heat_flux_limiter = .FALSE.
    ! Fraction of free streaming heat flux used in limiter
    flux_limiter = 0.01_num
    
    ! Remap kinetic energy correction. LARE does not
    ! perfectly conserve kinetic energy during the remap step
    ! This missing energy can be added back into the simulation
    ! as a uniform heating. Turning rke to true turns on this
    ! addition
    rke = .FALSE.

    ! The code to choose the initial conditions. The valid choices are
    ! IC_NEW - Use set_initial_conditions in "initial_conditions.f90"
    !         to setup new initial conditions
    ! IC_RESTART - Load the output file with index restart_snapshot and
    ! use it as the initial conditions
    initial = IC_NEW
    restart_snapshot = 1

    ! Turn on or off the physics package dealing with partial ionisation
    ! If include_neutrals is true then the code will calculate the
    ! ionisation fraction of the plasma at each point in the domain
    include_neutrals = .FALSE.

    ! If cowling_resistivity is true then the code calculates and
    ! applies the Cowling Resistivity to the MHD equations
    cowling_resistivity = .FALSE.

    ! Set the boundary conditions on the four edges of the simulation domain
    ! Valid constants are
    ! BC_PERIODIC - Periodic boundary conditions
    ! BC_OPEN - Reimann characteristic boundary conditions
    ! BC_OTHER - User boundary conditions specified in "boundary.f90"
    xbc_left = BC_PERIODIC
    xbc_right = BC_PERIODIC
    ybc_up = BC_PERIODIC
    ybc_down = BC_PERIODIC
    zbc_front = BC_OTHER
    zbc_back = BC_OTHER

    ! set to true to turn on routine for damped boundaries
    ! these routines are in boundary.f90 and you should check these
    ! actually do what you want.
    damping = .FALSE.

    ! Set the equation of state. Valid choices are
    ! EOS_IDEAL - Simple ideal gas for perfectly ionised plasma
    ! EOS_PI - Simple ideal gas for partially ionised plasma
    ! EOS_ION - EOS_PI plus the ionisation potential 
    ! N.B. read the manual for notes on these choices
    eos_number = EOS_IDEAL 
    ! EOS_IDEAL also requires that you specific a reduced mass factor
    ! This is the reduced mass in units of the average ion mass
    ! This has no meaning for other EOS choices
    ! For fully ionised simple cases set to 0.5
    ! For neutral hydrogen set to 1.0
    reduced_mass = 0.5 

  END SUBROUTINE control_variables



  SUBROUTINE set_output_dumps

    ! The output directory for the code
    data_dir = "Data"

    ! The interval between output snapshots. If SI_Input is true
    ! Then this is in seconds
    dt_snapshots = 10.0_num

    ! dump_mask is an array which specifies which quantities the
    ! code should output to disk in a data dump.
    ! The codes are
    ! 1  - rho
    ! 2  - energy
    ! 3  - vx
    ! 4  - vy
    ! 5  - vz
    ! 6  - bx
    ! 7  - by
    ! 8  - bz
    ! 9  - temperature
    ! 10 - pressure
    ! 11 - cs (sound speed)
    ! 12 - parallel_current
    ! 13 - perp_current
    ! 14 - neutral_faction
    ! 15 - eta_perp
    ! 16 - eta
    ! 17 - jx
    ! 18 - jy
    ! 19 - jz
    ! If a given element of dump_mask is true then that field is dumped
    ! If the element is false then the field isn't dumped
    ! N.B. if dump_mask(1:8) not true then the restart will not work
    dump_mask = .FALSE.
    dump_mask(1:10) = .TRUE. 
    IF (include_neutrals)  dump_mask(14) = .TRUE. 
    IF (resistive_mhd)  dump_mask(16) = .TRUE. 
    IF (cowling_resistivity)  dump_mask(15) = .TRUE. 

  END SUBROUTINE set_output_dumps

END MODULE control
