MODULE control

  USE shared_data
  USE normalise

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: User_Normalisation, Control_Variables, Set_Output_Dumps

CONTAINS

  SUBROUTINE User_Normalisation
    !Set the normalising constants for LARE
    !This is needed to allow the use of some physics
    !Modules which are coded in SI units

    !Should the code use SI units
    SI = .TRUE.

    !Gamma is the ratio of specific heat capacities
    gamma = 5.0_num/3.0_num

    !Average mass of an ion in proton masses
    !The code assumes a single ion species with this mass
    mf = 1.4_num

    !The equations describing the normalisation in LARE
    !Have three free parameters which must be specified by
    !The end user. These must be the normailisation used for
    !your initial conditions if not running in SI, otherwise
    !they are arbitrary. 

    !Magnetic field normalisation in Tesla
    B0 = 0.12_num
    !Length normalisation in m
    L0 = 150.0_num
    !Density normalisation in kg/m^3
    RHO0 = 2.7e-6_num

  END SUBROUTINE User_Normalisation



  SUBROUTINE Control_Variables

    !Set the number of gridpoints in x and y directions
    nx_global = 40
    ny_global = 40
    nz_global = 40

    !Set the maximum number of iterations of the core solver before the code
    !terminates. If nsteps < 0 then the code will run until t=t_end
    nsteps = 100

    !The maximum runtime of the code
    !If SI_Input is true then this is in seconds
    t_end = 20.0_num

    !Shock viscosities as detailed in manual - they are dimensionless
    visc1 = 0.1_num
    visc2 = 0.5_num
    !Real viscosity expressed as the inverse Reynolds number, i.e. the
    !same for normalised and SI input
    visc3 = 0.0_num

    !Set these constants to manually
    !override the domain decomposition.
    !If either constant is set to zero 
    !then the code will try to automatically
    !decompose in this direction
    nprocx = 0
    nprocy = 0
    nprocz = 0

    !The length of the domain in the x direction
    !If SI_Input is true then this is in metres
    x_start = -2.0e5_num
    x_end = 2.0e5_num
    !Should the x grid be stretched or uniform
    x_stretch = .FALSE.

    !The length of the domain in the y direction
    !If SI_Input is true then this is in metres
    y_start = -2.0e5_num
    y_end = 2.0e5_num
    !Should the y grid be stretched of uniform
    y_stretch = .FALSE.

    z_start = -2.0e5_num!/50.0_num
    z_end = 2.0e5_num!/50.0_num
    z_stretch = .FALSE.

    !Turn on or off the resistive parts of the MHD equations
    ResistiveMHD = .FALSE.
    !The background resistivity expressed as the inverse Lundquist number, i.e. the
    !same for normalised and SI input
    eta_background = 0.0_num

    !The critical current for triggering anomalous resistivity
    !and the resistivity when above the critical current
    !The resistivity is expressed as the inverse Lundquist number, i.e. the
    !same for normalised and SI input, bit the j_max must be in SI if using SI units
    j_max = 0.0_num
    eta0 = 0.0_num

    !Turn on or off the HallMHD term in the MHD equations
    !Well actually this does nothing as it isn't fully
    !included yet! Sorry!
    HallMHD = .FALSE.
    !Set the ion skin depth. If SI_Input is true then
    !This is in metres. Note that this should be fixed to
    !(the speed of light)/(ion plasma frequency) with the
    !plasma frequncy fixed by the chosen normalisation of
    !density (through the normalising mass density). Here
    !it is treated as a free pararmeter so be careful!
    lambda0 = 0.0_num

    !Turn on or off the Braginskii thermal conduction term in
    !the MHD equations
    Conduction = .FALSE.

    !Remap kinetic energy correction. LARE does not
    !perfectly conserve kinetic energy during the remap step
    !This missing energy can be added back into the simulation
    !as a uniform heating. Turning rke to true turns on this
    !addition
    rke = .FALSE.

    !The code to choose the initial conditions. The valid choices are
    !IC_NEW - Use Set_Initial_Conditions in "initial_conditions.f90" 
    !         to setup new initial conditions
    !IC_RESTART - Load the output file with index Restart_Snapshot and
    !use it as the initial conditions
    initial = IC_NEW
    Restart_Snapshot = 1

    !Turn on or off the physics package dealing with partial ionisation
    !If include_neutrals is true then the code will calculate the 
    !ionisation fraction of the plasma at each point in the domain
    include_neutrals = .FALSE.
    !If cowling_resistivity is true then the code calculates and 
    !applies the Cowling Resistivity to the MHD equations
    cowling_resistivity = .FALSE.

    !Set the boundary conditions on the four edges of the simulation domain
    !Valid constants are 
    !BC_PERIODIC - Periodic boundary conditions
    !BC_OPEN - Reimann characteristic boundary conditions
    !BC_OTHER - Other boundary conditions specified in "boundary.f90"
    xbc_left = BC_PERIODIC
    xbc_right = BC_PERIODIC
    ybc_up = BC_PERIODIC
    ybc_down = BC_PERIODIC
    zbc_front = BC_PERIODIC
    zbc_back = BC_PERIODIC

    !set to true to turn on routine for damped boundaries
    damping = .FALSE.

    !Set the equation of state. Valid choices are
    !EOS_IDEAL - Simple ideal gas for perfectly ionised plasma
    !EOS_PI - Simple ideal gas for partially ionised plasma
    !EOS_ION - EOS_PI plus the ionisation potential
    eos_number = EOS_IDEAL

  END SUBROUTINE Control_Variables



  SUBROUTINE Set_Output_Dumps

    !The output directory for the code
    Data_Dir = "Data"

    !The interval between output snapshots. If SI_Input is true
    !Then this is in seconds
    dt_snapshots = 0.!2_num

    !dumpmask is an array which specifies which quantities the 
    !code should output to disk in a data dump.
    !The codes are
    !1  - rho
    !2  - energy
    !3  - vx
    !4  - vy
    !5  - vz
    !6  - bx
    !7  - by
    !8  - bz
    !9  - temperature
    !10 - pressure
    !11 - cs (sound speed)
    !12 - parallel_current
    !13 - perp_current
    !14 - neutral_faction
    !15 - eta_perp
    !16 - eta
    !17 - jx
    !18 - jy
    !19 - jz
    !If a given element of dumpmask is true then that field is dumped
    !If the element is false then the field isn't dumped
    !N.B. if dumpmask(1:8) not true then the restart will not work
    dumpmask = .FALSE.
    dumpmask(1:10) = .TRUE.
    !dumpmask(14) = .TRUE.

  END SUBROUTINE Set_Output_Dumps



END MODULE control
