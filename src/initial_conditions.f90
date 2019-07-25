MODULE initial_conditions

  USE shared_data
  USE neutral
  USE boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_initial_conditions

  REAL(num), DIMENSION(:), ALLOCATABLE :: axis, temperature_1d
  INTEGER :: table_count

CONTAINS

  !****************************************************************************
  ! This function sets up the initial condition for the code
  ! The variables which must be set are:
  !   rho - density
  !   v{x,y,z} - Velocities in x, y, z
  !   b{x,y,z} - Magnetic fields in x, y, z
  !   energy - Specific internal energy
  !
  ! You may also need the neutral fraction. This can be calculated by a
  ! function call to get_neutral(temperature, rho, z). This routine is in
  ! core/neutral.f90 and requires the local temperature and mass density.
  ! For example to set xi_n to the neutral fraction use:
  !   xi_n = get_neutral(temperature, rho, z)
  !****************************************************************************


    SUBROUTINE set_initial_conditions

      !This is actually based on the default initial conditions shipped with LARE.
      !With a few tweaks to constants, it's a very good fit to the C7 model of 
      !Avrett and Loeser  
      INTEGER :: loop
      INTEGER :: ix, iz
      REAL(num) :: a1, a2, dg
      REAL(num) :: a=2.0_num, Tph=11.8_num
      REAL(num) :: r1, maxerr, xi_v
      REAL(num) :: correct1, correct2, t_bottom
      REAL(num) :: pressure0, energy0, v0, mbar, temp0, time0
      REAL(num), DIMENSION(:), ALLOCATABLE :: zc_global, dzb_global, dzc_global
      REAL(num), DIMENSION(:), ALLOCATABLE :: grav_ref, temp_ref, rho_ref
      REAL(num), DIMENSION(:), ALLOCATABLE :: mu_m

      !This is very inelegant, should do better
      pressure0 = B_norm**2 / mu0_si
      energy0 = pressure0 / rho_norm
      v0 = SQRT(energy0)
      mbar = mf * mh_si
      temp0 = (mbar / kb_si) * energy0
      time0 = L_norm/v0

      vx = 0.0_num
      vy = 0.0_num
      vz = 0.0_num
      bx = 0.0_num
      by = 0.0_num
      bz = 0.0_num
      rho = 1.0_num
      energy = 1.0_num
      grav = 0.0_num

      ALLOCATE(zc_global(-1:nz_global+1))
      ALLOCATE(dzb_global(-1:nz_global+1), dzc_global(-1:nz_global))
      ALLOCATE(grav_ref(-1:nz_global+2), temp_ref(-1:nz_global+2))
      ALLOCATE(rho_ref(-1:nx_global+2))
      ALLOCATE(mu_m(-1:nx_global+2))
      CALL read_al_temperatures
      CALL lookup_temperature(0.0_num,t_bottom)
      Tph = t_bottom
    
      !fill in zc_global with the positions central to the zb_global points
      DO iz = -1,nz_global+1
         zc_global(iz) = 0.5_num * (zb_global(iz-1) + zb_global(iz))
      END DO
    
      !fill in dzb_global and dzc_global
      DO iz = -1,nz_global
         dzb_global(iz) = zb_global(iz) - zb_global(iz-1)
         dzc_global(iz) = zc_global(iz+1) - zc_global(iz)
      END DO
      dzb_global(nx_global+1) = dzb_global(nx_global) 
    
      !fill in the reference gravity array - lowering grav to zero at the top 
      !of the corona smoothly from a1 to grav=0 at a2 and above
      grav_ref = 11.78_num
      a1 = 25.0_num!zb_global(nx_global) - 20.0_num
      a2 = 40.0_num
      DO iz = 0,nx_global+2
         IF (zb_global(iz) > a1) THEN
            grav_ref(iz) = 11.78_num * (1.0_num + COS(pi * (zb_global(iz) - a1) &
                 / (a2-a1))) / 2.0_num
         END IF
         IF (zb_global(iz) > a2) THEN
            grav_ref(iz) = 0.0_num
         END IF
      END DO    
      grav_ref(-1) = grav_ref(0)
      grav_ref(nx_global+1:nx_global+2) = grav_ref(nx_global)

      !calculate the density profile, starting from the refence density at the
      !photosphere and calculating up and down from there including beta
      rho_ref = 1.0_num
      mu_m = 1.0_num
      temp_ref = 1.0_num
      IF (eos_number == EOS_IDEAL .AND. (.NOT. neutral_gas)) mu_m = 0.5_num

      DO loop = 1,1000
         maxerr = 0.0_num
         !Go from photosphere down
         DO iz = -1,nx_global+1
            IF (zc_global(iz) < 0.0_num) THEN
               temp_ref(iz) = Tph - a * (gamma - 1.0_num) &
                    * zc_global(iz) * grav_ref(iz) * mu_m(iz) / gamma 
            END IF
            IF (zc_global(iz) >= 0.0_num) THEN
               temp_ref(iz) = 6600.0/temp0
               correct1=-1980.0_num/temp0*EXP(-((zc_global(iz)-6.0e5/L_norm)/(3.2e5/L_norm))**2)
               correct2=-1050.0_num/temp0*EXP(-((zc_global(iz)-2.5e5/L_norm)/(2.2e5/L_norm))**2)
               temp_ref(iz) = (temp_ref(iz) + correct1 + correct2)
            END IF
            IF(zc_global(iz) .GT. 0.0_num) THEN
              CALL lookup_temperature(zc_global(iz),temp_ref(iz))
            ENDIF
          END DO
          temp_ref(nx_global+1:nx_global+2) = temp_ref(nx_global)
           
         DO iz = nx_global,0,-1
            IF (zc_global(iz) < 0.0_num) THEN  
               dg = 1.0_num / (dzb_global(iz) + dzb_global(iz-1))
               rho_ref(iz-1) = rho_ref(iz) * (temp_ref(iz) &
                    /dzc_global(iz-1)/mu_m(iz)+grav_ref(iz-1)*dzb_global(iz)*dg)
               rho_ref(iz-1) = rho_ref(iz-1) / (temp_ref(iz-1) &
                    /dzc_global(iz-1)/mu_m(iz-1)-grav_ref(iz-1)*dzb_global(iz-1)*dg)
            END IF
         END DO
         !Now move from the photosphere up to the corona
         DO iz = 0,nx_global
            IF (zc_global(iz) >= 0.0_num) THEN
               dg = 1.0_num / (dzb_global(iz)+dzb_global(iz-1))
               rho_ref(iz) = rho_ref(iz-1) * (temp_ref(iz-1) &
                    /dzc_global(iz-1)/mu_m(iz-1)-grav_ref(iz-1)*dzb_global(iz-1)*dg)
               rho_ref(iz) = rho_ref(iz) / (temp_ref(iz) &
                    /dzc_global(iz-1)/mu_m(iz)+grav_ref(iz-1)*dzb_global(iz)*dg)
            END IF
         END DO
         IF (eos_number /= EOS_IDEAL) THEN
            DO iz=0,nx_global,1
               xi_v = get_neutral(temp_ref(iz),rho_ref(iz))
               r1 = mu_m(iz)
               mu_m(iz) = 1.0_num / (2.0_num-xi_v)
               maxerr = MAX(maxerr, ABS(mu_m(iz) - r1))
            END DO
         END IF
         IF (maxerr < 1.e-16_num) EXIT
      END DO
    
      rho_ref(nz_global+1:nz_global+2) = rho_ref(nz_global)   

      !fill in all the final arrays from the ref arrays
      grav(:) = grav_ref(coordinates(1)*ny-1:coordinates(1)*ny+ny+2)
      DO iy = -1,ny+2
       DO ix = -1,nx+2
         rho(ix,iy,:) = rho_ref(coordinates(1)*nz-1:coordinates(1)*nz+nz+2)
         energy(ix,iy,:) = temp_ref(coordinates(1)*nz-1:coordinates(1)*nz+nz+2)
       END DO
      END DO

      DO iz = -1,nz+2,1
       DO iy = -1,ny+2,1
          DO ix = -1,nx+2,1
            IF (eos_number /= EOS_IDEAL) THEN         
              xi_v = get_neutral(energy(ix,iy,iz), rho(ix,iy,iz))
            ELSE  
              IF (neutral_gas) THEN
                xi_v = 1.0_num
              ELSE
                xi_v = 0.0_num
              END IF
            END IF
            energy(ix,iy,iz) = (energy(ix,iy,iz) * (2.0_num - xi_v) &
               + (1.0_num - xi_v) * ionise_pot * (gamma - 1.0_num)) &
               / (gamma - 1.0_num)
          END DO
       END DO
      END DO
      DO iy= -1,ny+2
       DO ix= -1,nx+2
          energy(ix,iy,ny+2) = energy(ix,iy,ny+1)
        END DO
      END DO

      DEALLOCATE(zc_global, dzb_global, dzc_global, mu_m)
      DEALLOCATE(grav_ref, temp_ref, rho_ref)

      IF (IAND(initial, IC_NEW) /= 0) CALL potential_field

    END SUBROUTINE set_initial_conditions
    




    SUBROUTINE potential_field()

        REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: phi
        REAL(num) :: w, errmax, error, residual, fractional_error
        REAL(num) :: by_min, by_min_local
        REAL(num) :: by_max, by_max_local
        REAL(num) :: dx1, dy1
        INTEGER :: loop, x1, y1, redblack
        LOGICAL :: converged

        ALLOCATE(phi(-1:nx+2,-1:ny+2, -1:nz+2))
        phi = 0.0_num

        ! Only works for uniform x,y grids - no strecthed grids
        dx1 = (REAL(nx_global, num) / length_x)**2
        dy1 = (REAL(nx_global, num) / length_y)**2

        converged = .FALSE.
        w = 2.0_num / (1.0_num + SIN(pi / REAL(nx_global,num)))
        fractional_error = 1.e-8_num


        DO iz = 1, nz
         DO iy = 1, ny
           DO ix = 0, nx
             bx(ix,iy,iz) = -(phi(ix+1,iy,iz)-phi(ix,iy,iz))/dxc(ix)
           END DO
         END DO
        END DO

        DO iz = 1, nz
         DO iy = 0, ny
           DO ix = 1, nx
             by(ix,iy,iz) = -(phi(ix,iy+1,iz)-phi(ix,iy,iz))/dyc(iy)
           END DO
         END DO
        END DO

        DO iz = 1, nz
         DO iy = 1, ny
           DO ix = 0, nx
             bz(ix,iy,iz) = -(phi(ix,iy,iz+1)-phi(ix,iy,iz))/dzc(iz)
           END DO
         END DO
        END DO

        CALL bfield_bcs

        DEALLOCATE(phi)



    END SUBROUTINE potential_field




    SUBROUTINE read_al_temperatures()

      REAL(num) :: pressure0, energy0, v0, mbar, temp0
      INTEGER :: iel, reverse, start, end, delta

      !This is very inelegant, should do better
      pressure0 = B_norm**2 / mu0_si
      energy0 = pressure0 / rho_norm
      v0 = SQRT(energy0)
      mbar = mf * mh_si
      temp0 = (mbar / kb_si) * energy0  

      IF (rank .EQ. 0) THEN
        OPEN(UNIT=55,FILE='ExternalFiles/tables/c7.table',FORM='FORMATTED')
        READ(55,*) table_count, reverse
        ALLOCATE(axis(1:table_count), temperature_1d(1:table_count))
        IF (reverse .EQ. 0) THEN
          start=1
          end=table_count
          delta=1
        ELSE
          start=table_count
          end=1
          delta=-1
        ENDIF
        DO iel=start,end,delta
          READ(55,*) axis(iel),temperature_1d(iel)
        ENDDO
        CLOSE(55)

        axis=axis/L_norm
        temperature_1d=temperature_1d/temp0

      ENDIF

      CALL MPI_BCAST(table_count,1,MPI_INTEGER,0,comm,errcode)
      IF (rank .NE. 0) ALLOCATE(axis(1:table_count), &
          temperature_1d(1:table_count))
      CALL MPI_BCAST(axis,table_count,mpireal,0,comm,errcode)
      CALL MPI_BCAST(temperature_1d,table_count,mpireal,0,comm,errcode)

    END SUBROUTINE read_al_temperatures



    SUBROUTINE lookup_temperature(y,temp)
      REAL(num), INTENT(IN) :: y
      REAL(num), INTENT(OUT) :: temp
      INTEGER :: iz
      REAL(num) :: frac
      IF (y .GE. MAXVAL(axis)) THEN
        temp=temperature_1d(table_count)
        RETURN
      ENDIF
      DO iz=1,table_count-1
        IF (axis(iz) .LE. y .AND. axis(iz+1) .GT. y) THEN
          frac=(y-axis(iz))/(axis(iz+1)-axis(iz))
          temp=temperature_1d(iz)+(temperature_1d(iz+1)-temperature_1d(iz))*frac
        ENDIF
      ENDDO
    END SUBROUTINE lookup_temperature


END MODULE initial_conditions
