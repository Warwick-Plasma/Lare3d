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

!******************************************************************************
! Interprocessor boundary calls
!******************************************************************************

MODULE mpiboundary

  USE shared_data

  IMPLICIT NONE

CONTAINS

  !****************************************************************************
  ! Boundary exchange for magnetic field
  !****************************************************************************

  SUBROUTINE bfield_mpi

    CALL MPI_SENDRECV( &
        bx(   1,-1,-1), 1, bx_xface,  proc_x_min, tag, &
        bx(nx+1,-1,-1), 1, bx_xface,  proc_x_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        bx(nx-2,-1,-1), 1, bx_xface1, proc_x_max, tag, &
        bx(  -2,-1,-1), 1, bx_xface1, proc_x_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        by(   1,-2,-1), 1, by_xface,  proc_x_min, tag, &
        by(nx+1,-2,-1), 1, by_xface,  proc_x_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        by(nx-1,-2,-1), 1, by_xface,  proc_x_max, tag, &
        by(  -1,-2,-1), 1, by_xface,  proc_x_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        bz(   1,-1,-2), 1, bz_xface,  proc_x_min, tag, &
        bz(nx+1,-1,-2), 1, bz_xface,  proc_x_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        bz(nx-1,-1,-2), 1, bz_xface,  proc_x_max, tag, &
        bz(  -1,-1,-2), 1, bz_xface,  proc_x_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        bx(-2,   1,-1), 1, bx_yface,  proc_y_min, tag, &
        bx(-2,ny+1,-1), 1, bx_yface,  proc_y_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        bx(-2,ny-1,-1), 1, bx_yface,  proc_y_max, tag, &
        bx(-2,  -1,-1), 1, bx_yface,  proc_y_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        by(-1,   1,-1), 1, by_yface,  proc_y_min, tag, &
        by(-1,ny+1,-1), 1, by_yface,  proc_y_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        by(-1,ny-2,-1), 1, by_yface1, proc_y_max, tag, &
        by(-1,  -2,-1), 1, by_yface1, proc_y_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        bz(-1,   1,-2), 1, bz_yface,  proc_y_min, tag, &
        bz(-1,ny+1,-2), 1, bz_yface,  proc_y_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        bz(-1,ny-1,-2), 1, bz_yface,  proc_y_max, tag, &
        bz(-1,  -1,-2), 1, bz_yface,  proc_y_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        bx(-2,-1,   1), 1, bx_zface,  proc_z_min, tag, &
        bx(-2,-1,nz+1), 1, bx_zface,  proc_z_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        bx(-2,-1,nz-1), 1, bx_zface,  proc_z_max, tag, &
        bx(-2,-1,  -1), 1, bx_zface,  proc_z_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        by(-1,-2,   1), 1, by_zface,  proc_z_min, tag, &
        by(-1,-2,nz+1), 1, by_zface,  proc_z_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        by(-1,-2,nz-1), 1, by_zface,  proc_z_max, tag, &
        by(-1,-2,  -1), 1, by_zface,  proc_z_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        bz(-1,-1,   1), 1, bz_zface,  proc_z_min, tag, &
        bz(-1,-1,nz+1), 1, bz_zface,  proc_z_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        bz(-1,-1,nz-2), 1, bz_zface1, proc_z_max, tag, &
        bz(-1,-1,  -2), 1, bz_zface1, proc_z_min, tag, &
        comm, status, errcode)

  END SUBROUTINE bfield_mpi



  !****************************************************************************
  ! Boundary exchange for specific internal energy
  !****************************************************************************

  SUBROUTINE energy_mpi

    CALL MPI_SENDRECV( &
        energy(   1,-1,-1), 1, cell_xface, proc_x_min, tag, &
        energy(nx+1,-1,-1), 1, cell_xface, proc_x_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        energy(nx-1,-1,-1), 1, cell_xface, proc_x_max, tag, &
        energy(  -1,-1,-1), 1, cell_xface, proc_x_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        energy(-1,   1,-1), 1, cell_yface, proc_y_min, tag, &
        energy(-1,ny+1,-1), 1, cell_yface, proc_y_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        energy(-1,ny-1,-1), 1, cell_yface, proc_y_max, tag, &
        energy(-1,  -1,-1), 1, cell_yface, proc_y_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        energy(-1,-1,   1), 1, cell_zface, proc_z_min, tag, &
        energy(-1,-1,nz+1), 1, cell_zface, proc_z_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        energy(-1,-1,nz-1), 1, cell_zface, proc_z_max, tag, &
        energy(-1,-1,  -1), 1, cell_zface, proc_z_min, tag, &
        comm, status, errcode)

  END SUBROUTINE energy_mpi


  !****************************************************************************
  ! Boundary exchange for density
  !****************************************************************************

  SUBROUTINE density_mpi

    CALL MPI_SENDRECV( &
        rho(   1,-1,-1), 1, cell_xface, proc_x_min, tag, &
        rho(nx+1,-1,-1), 1, cell_xface, proc_x_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        rho(nx-1,-1,-1), 1, cell_xface, proc_x_max, tag, &
        rho(  -1,-1,-1), 1, cell_xface, proc_x_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        rho(-1,   1,-1), 1, cell_yface, proc_y_min, tag, &
        rho(-1,ny+1,-1), 1, cell_yface, proc_y_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        rho(-1,ny-1,-1), 1, cell_yface, proc_y_max, tag, &
        rho(-1,  -1,-1), 1, cell_yface, proc_y_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        rho(-1,-1,   1), 1, cell_zface, proc_z_min, tag, &
        rho(-1,-1,nz+1), 1, cell_zface, proc_z_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        rho(-1,-1,nz-1), 1, cell_zface, proc_z_max, tag, &
        rho(-1,-1,  -1), 1, cell_zface, proc_z_min, tag, &
        comm, status, errcode)

  END SUBROUTINE density_mpi



  !****************************************************************************
  ! Full timestep velocity boundary exchange
  !****************************************************************************

  SUBROUTINE velocity_mpi

    CALL MPI_SENDRECV( &
        vx(   0,-2,-2), 1, node_xface1, proc_x_min, tag, &
        vx(  nx,-2,-2), 1, node_xface1, proc_x_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        vx(nx-2,-2,-2), 1, node_xface,  proc_x_max, tag, &
        vx(  -2,-2,-2), 1, node_xface,  proc_x_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        vy(   0,-2,-2), 1, node_xface1, proc_x_min, tag, &
        vy(  nx,-2,-2), 1, node_xface1, proc_x_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        vy(nx-2,-2,-2), 1, node_xface,  proc_x_max, tag, &
        vy(  -2,-2,-2), 1, node_xface,  proc_x_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        vz(   0,-2,-2), 1, node_xface1, proc_x_min, tag, &
        vz(  nx,-2,-2), 1, node_xface1, proc_x_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        vz(nx-2,-2,-2), 1, node_xface,  proc_x_max, tag, &
        vz(  -2,-2,-2), 1, node_xface,  proc_x_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        vx(-2,   0,-2), 1, node_yface1, proc_y_min, tag, &
        vx(-2,  ny,-2), 1, node_yface1, proc_y_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        vx(-2,ny-2,-2), 1, node_yface,  proc_y_max, tag, &
        vx(-2,  -2,-2), 1, node_yface,  proc_y_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        vy(-2,   0,-2), 1, node_yface1, proc_y_min, tag, &
        vy(-2,  ny,-2), 1, node_yface1, proc_y_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        vy(-2,ny-2,-2), 1, node_yface,  proc_y_max, tag, &
        vy(-2,  -2,-2), 1, node_yface,  proc_y_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        vz(-2,   0,-2), 1, node_yface1, proc_y_min, tag, &
        vz(-2,  ny,-2), 1, node_yface1, proc_y_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        vz(-2,ny-2,-2), 1, node_yface,  proc_y_max, tag, &
        vz(-2,  -2,-2), 1, node_yface,  proc_y_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        vx(-2,-2,   0), 1, node_zface1, proc_z_min, tag, &
        vx(-2,-2,  nz), 1, node_zface1, proc_z_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        vx(-2,-2,nz-2), 1, node_zface,  proc_z_max, tag, &
        vx(-2,-2,  -2), 1, node_zface,  proc_z_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        vy(-2,-2,   0), 1, node_zface1, proc_z_min, tag, &
        vy(-2,-2,  nz), 1, node_zface1, proc_z_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        vy(-2,-2,nz-2), 1, node_zface,  proc_z_max, tag, &
        vy(-2,-2,  -2), 1, node_zface,  proc_z_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        vz(-2,-2,   0), 1, node_zface1, proc_z_min, tag, &
        vz(-2,-2,  nz), 1, node_zface1, proc_z_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        vz(-2,-2,nz-2), 1, node_zface,  proc_z_max, tag, &
        vz(-2,-2,  -2), 1, node_zface,  proc_z_min, tag, &
        comm, status, errcode)

  END SUBROUTINE velocity_mpi


  !****************************************************************************
  ! Half timestep velocity boundary exchange
  !****************************************************************************

  SUBROUTINE remap_v_mpi

    CALL MPI_SENDRECV( &
        vx1(   0,-2,-2), 1, node_xface1, proc_x_min, tag, &
        vx1(  nx,-2,-2), 1, node_xface1, proc_x_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        vx1(nx-2,-2,-2), 1, node_xface,  proc_x_max, tag, &
        vx1(  -2,-2,-2), 1, node_xface,  proc_x_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        vy1(   0,-2,-2), 1, node_xface1, proc_x_min, tag, &
        vy1(  nx,-2,-2), 1, node_xface1, proc_x_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        vy1(nx-2,-2,-2), 1, node_xface,  proc_x_max, tag, &
        vy1(  -2,-2,-2), 1, node_xface,  proc_x_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        vz1(   0,-2,-2), 1, node_xface1, proc_x_min, tag, &
        vz1(  nx,-2,-2), 1, node_xface1, proc_x_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        vz1(nx-2,-2,-2), 1, node_xface,  proc_x_max, tag, &
        vz1(  -2,-2,-2), 1, node_xface,  proc_x_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        vx1(-2,   0,-2), 1, node_yface1, proc_y_min, tag, &
        vx1(-2,  ny,-2), 1, node_yface1, proc_y_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        vx1(-2,ny-2,-2), 1, node_yface,  proc_y_max, tag, &
        vx1(-2,  -2,-2), 1, node_yface,  proc_y_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        vy1(-2,   0,-2), 1, node_yface1, proc_y_min, tag, &
        vy1(-2,  ny,-2), 1, node_yface1, proc_y_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        vy1(-2,ny-2,-2), 1, node_yface,  proc_y_max, tag, &
        vy1(-2,  -2,-2), 1, node_yface,  proc_y_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        vz1(-2,   0,-2), 1, node_yface1, proc_y_min, tag, &
        vz1(-2,  ny,-2), 1, node_yface1, proc_y_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        vz1(-2,ny-2,-2), 1, node_yface,  proc_y_max, tag, &
        vz1(-2,  -2,-2), 1, node_yface,  proc_y_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        vx1(-2,-2,   0), 1, node_zface1, proc_z_min, tag, &
        vx1(-2,-2,  nz), 1, node_zface1, proc_z_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        vx1(-2,-2,nz-2), 1, node_zface,  proc_z_max, tag, &
        vx1(-2,-2,  -2), 1, node_zface,  proc_z_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        vy1(-2,-2,   0), 1, node_zface1, proc_z_min, tag, &
        vy1(-2,-2,  nz), 1, node_zface1, proc_z_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        vy1(-2,-2,nz-2), 1, node_zface,  proc_z_max, tag, &
        vy1(-2,-2,  -2), 1, node_zface,  proc_z_min, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV( &
        vz1(-2,-2,   0), 1, node_zface1, proc_z_min, tag, &
        vz1(-2,-2,  nz), 1, node_zface1, proc_z_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV( &
        vz1(-2,-2,nz-2), 1, node_zface,  proc_z_max, tag, &
        vz1(-2,-2,  -2), 1, node_zface,  proc_z_min, tag, &
        comm, status, errcode)

  END SUBROUTINE remap_v_mpi

END MODULE mpiboundary
