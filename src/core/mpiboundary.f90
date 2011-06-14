MODULE mpiboundary

  USE shared_data
  IMPLICIT NONE

CONTAINS

  SUBROUTINE bfield_MPI

    CALL MPI_SENDRECV(bx(1:nx, 1:ny, 1:2), 2*nx*ny, mpireal, proc_z_min, tag, &
        bx(1:nx, 1:ny, nz+1:nz+2), 2*nx*ny, mpireal, proc_z_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(bx(1:nx, 1:ny, nz-1:nz), 2*nx*ny, mpireal, proc_z_max, tag, &
        bx(1:nx, 1:ny, -1:0), 2*nx*ny, mpireal, proc_z_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(by(1:nx, 1:ny, 1:2), 2*nx*ny, mpireal, proc_z_min, tag, &
        by(1:nx, 1:ny, nz+1:nz+2), 2*nx*ny, mpireal, proc_z_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(by(1:nx, 1:ny, nz-1:nz), 2*nx*ny, mpireal, proc_z_max, tag, &
        by(1:nx, 1:ny, -1:0), 2*nx*ny, mpireal, proc_z_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(bz(1:nx, 1:ny, 1:2), 2*nx*ny, mpireal, proc_z_min, tag, &
        bz(1:nx, 1:ny, nz+1:nz+2), 2*nx*ny, mpireal, proc_z_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(bz(1:nx, 1:ny, nz-2:nz), 3*nx*ny, mpireal, proc_z_max, tag, &
        bz(1:nx, 1:ny, -2:0), 3*nx*ny, mpireal, proc_z_min, tag, comm, &
        status, errcode)

    CALL MPI_SENDRECV(bx(1:2, 1:ny, :), 2*ny*(nz+4), mpireal, proc_x_min, tag, &
        bx(nx+1:nx+2, 1:ny, :), 2*ny*(nz+4), mpireal, proc_x_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(bx(nx-2:nx, 1:ny, :), 3*ny*(nz+4), mpireal, proc_x_max, tag, &
        bx(-2:0, 1:ny, :), 3*ny*(nz+4), mpireal, proc_x_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(by(1:2, 1:ny, :), 2*ny*(nz+4), mpireal, proc_x_min, tag, &
        by(nx+1:nx+2, 1:ny, :), 2*ny*(nz+4), mpireal, proc_x_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(by(nx-1:nx, 1:ny, :), 2*ny*(nz+4), mpireal, proc_x_max, tag, &
        by(-1:0, 1:ny, :), 2*ny*(nz+4), mpireal, proc_x_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(bz(1:2, 1:ny, :), 2*ny*(nz+5), mpireal, proc_x_min, tag, &
        bz(nx+1:nx+2, 1:ny, :), 2*ny*(nz+5), mpireal, proc_x_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(bz(nx-1:nx, 1:ny, :), 2*ny*(nz+5), mpireal, proc_x_max, tag, &
        bz(-1:0, 1:ny, :), 2*ny*(nz+5), mpireal, proc_x_min, tag, comm, &
        status, errcode)

    CALL MPI_SENDRECV(bx(:, ny-1:ny, :), 2*(nx+5)*(nz+4), mpireal, proc_y_max, tag, &
        bx(:, -1:0, :), 2*(nx+5)*(nz+4), mpireal, proc_y_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(bx(:, 1:2, :), 2*(nx+5)*(nz+4), mpireal, proc_y_min, tag, &
        bx(:, ny+1:ny+2, :), 2*(nx+5)*(nz+4), mpireal, proc_y_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(by(:, ny-2:ny, :), 3*(nx+4)*(nz+4), mpireal, proc_y_max, tag, &
        by(:, -2:0, :), 3*(nx+4)*(nz+4), mpireal, proc_y_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(by(:, 1:2, :), 2*(nx+4)*(nz+4), mpireal, proc_y_min, tag, &
        by(:, ny+1:ny+2, :), 2*(nx+4)*(nz+4), mpireal, proc_y_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(bz(:, ny-1:ny, :), 2*(nx+4)*(nz+5), mpireal, proc_y_max, tag, &
        bz(:, -1:0, :), 2*(nx+4)*(nz+5), mpireal, proc_y_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(bz(:, 1:2, :), 2*(nx+4)*(nz+5), mpireal, proc_y_min, tag, &
        bz(:, ny+1:ny+2, :), 2*(nx+4)*(nz+5), mpireal, proc_y_max, tag, comm, &
        status, errcode)

  END SUBROUTINE bfield_MPI



  SUBROUTINE energy_MPI

    CALL MPI_SENDRECV(energy(1:nx, 1:ny, 1:2), 2*nx*ny, mpireal, proc_z_min, tag, &
        energy(1:nx, 1:ny, nz+1:nz+2), 2*nx*ny, mpireal, proc_z_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(energy(1:nx, 1:ny, nz-1:nz), 2*nx*ny, mpireal, proc_z_max, &
        tag, energy(1:nx, 1:ny, -1:0), 2*nx*ny, mpireal, proc_z_min, tag, comm, &
        status, errcode)

    CALL MPI_SENDRECV(energy(1:2, 1:ny, :), 2*ny*(nz+4), mpireal, proc_x_min, tag, &
        energy(nx+1:nx+2, 1:ny, :), 2*ny*(nz+4), mpireal, proc_x_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(energy(nx-1:nx, 1:ny, :), 2*ny*(nz+4), mpireal, proc_x_max, &
        tag, energy(-1:0, 1:ny, :), 2*ny*(nz+4), mpireal, proc_x_min, tag, comm, &
        status, errcode)

    CALL MPI_SENDRECV(energy(:, ny-1:ny, :), 2*(nx+4)*(nz+4), mpireal, proc_y_max, &
        tag, energy(:, -1:0, :), 2*(nx+4)*(nz+4), mpireal, proc_y_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(energy(:, 1:2, :), 2*(nx+4)*(nz+4), mpireal, proc_y_min, tag, &
        energy(:, ny+1:ny+2, :), 2*(nx+4)*(nz+4), mpireal, proc_y_max, tag, comm, &
        status, errcode)

  END SUBROUTINE energy_MPI



  SUBROUTINE velocity_MPI

    INTEGER :: x_extent, y_extent, z_extent
    INTEGER :: xy_buf, xz_buf, yz_buf

    x_extent = nx + 5
    y_extent = ny + 5
    z_extent = nz + 5
    xy_buf = x_extent * y_extent
    xz_buf = x_extent * z_extent
    yz_buf = y_extent * z_extent

    CALL MPI_SENDRECV(vx(:, :, 1:2), 2*xy_buf, mpireal, proc_z_min, tag, &
        vx(:, :, nz+1:nz+2), 2*xy_buf, mpireal, proc_z_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vx(:, :, nz-2:nz), 3*xy_buf, mpireal, proc_z_max, tag, &
        vx(:, :, -2:0), 3*xy_buf, mpireal, proc_z_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vy(:, :, 1:2), 2*xy_buf, mpireal, proc_z_min, tag, &
        vy(:, :, nz+1:nz+2), 2*xy_buf, mpireal, proc_z_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vy(:, :, nz-2:nz), 3*xy_buf, mpireal, proc_z_max, tag, &
        vy(:, :, -2:0), 3*xy_buf, mpireal, proc_z_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vz(:, :, 1:2), 2*xy_buf, mpireal, proc_z_min, tag, &
        vz(:, :, nz+1:nz+2), 2*xy_buf, mpireal, proc_z_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vz(:, :, nz-2:nz), 3*xy_buf, mpireal, proc_z_max, tag, &
        vz(:, :, -2:0), 3*xy_buf, mpireal, proc_z_min, tag, comm, &
        status, errcode)

    CALL MPI_SENDRECV(vx(1:2, :, :), 2*yz_buf, mpireal, proc_x_min, tag, &
        vx(nx+1:nx+2, :, :), 2*yz_buf, mpireal, proc_x_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vx(nx-2:nx, :, :), 3*yz_buf, mpireal, proc_x_max, tag, &
        vx(-2:0, :, :), 3*yz_buf, mpireal, proc_x_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vy(1:2, :, :), 2*yz_buf, mpireal, proc_x_min, tag, &
        vy(nx+1:nx+2, :, :), 2*yz_buf, mpireal, proc_x_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vy(nx-2:nx, :, :), 3*yz_buf, mpireal, proc_x_max, tag, &
        vy(-2:0, :, :), 3*yz_buf, mpireal, proc_x_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vz(1:2, :, :), 2*yz_buf, mpireal, proc_x_min, tag, &
        vz(nx+1:nx+2, :, :), 2*yz_buf, mpireal, proc_x_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vz(nx-2:nx, :, :), 3*yz_buf, mpireal, proc_x_max, tag, &
        vz(-2:0, :, :), 3*yz_buf, mpireal, proc_x_min, tag, comm, &
        status, errcode)

    CALL MPI_SENDRECV(vx(:, ny-2:ny, :), 3*xz_buf, mpireal, proc_y_max, tag, &
        vx(:, -2:0, :), 3*xz_buf, mpireal, proc_y_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vx(:, 1:2, :), 2*xz_buf, mpireal, proc_y_min, tag, &
        vx(:, ny+1:ny+2, :), 2*xz_buf, mpireal, proc_y_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vy(:, ny-2:ny, :), 3*xz_buf, mpireal, proc_y_max, tag, &
        vy(:, -2:0, :), 3*xz_buf, mpireal, proc_y_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vy(:, 1:2, :), 2*xz_buf, mpireal, proc_y_min, tag, &
        vy(:, ny+1:ny+2, :), 2*xz_buf, mpireal, proc_y_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vz(:, ny-2:ny, :), 3*xz_buf, mpireal, proc_y_max, tag, &
        vz(:, -2:0, :), 3*xz_buf, mpireal, proc_y_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vz(:, 1:2, :), 2*xz_buf, mpireal, proc_y_min, tag, &
        vz(:, ny+1:ny+2, :), 2*xz_buf, mpireal, proc_y_max, tag, comm, &
        status, errcode)

  END SUBROUTINE velocity_MPI



  SUBROUTINE remap_v_MPI

    INTEGER :: x_extent, y_extent, z_extent
    INTEGER :: xy_buf, xz_buf, yz_buf

    x_extent = nx + 5
    y_extent = ny + 5
    z_extent = nz + 5
    xy_buf = x_extent * y_extent
    xz_buf = x_extent * z_extent
    yz_buf = y_extent * z_extent

    CALL MPI_SENDRECV(vx1(:, :, 1:2), 2*xy_buf, mpireal, proc_z_min, tag, &
        vx1(:, :, nz+1:nz+2), 2*xy_buf, mpireal, proc_z_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vx1(:, :, nz-2:nz), 3*xy_buf, mpireal, proc_z_max, tag, &
        vx1(:, :, -2:0), 3*xy_buf, mpireal, proc_z_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vy1(:, :, 1:2), 2*xy_buf, mpireal, proc_z_min, tag, &
        vy1(:, :, nz+1:nz+2), 2*xy_buf, mpireal, proc_z_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vy1(:, :, nz-2:nz), 3*xy_buf, mpireal, proc_z_max, tag, &
        vy1(:, :, -2:0), 3*xy_buf, mpireal, proc_z_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vz1(:, :, 1:2), 2*xy_buf, mpireal, proc_z_min, tag, &
        vz1(:, :, nz+1:nz+2), 2*xy_buf, mpireal, proc_z_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vz1(:, :, nz-2:nz), 3*xy_buf, mpireal, proc_z_max, tag, &
        vz1(:, :, -2:0), 3*xy_buf, mpireal, proc_z_min, tag, comm, &
        status, errcode)

    CALL MPI_SENDRECV(vx1(1:2, :, :), 2*yz_buf, mpireal, proc_x_min, tag, &
        vx1(nx+1:nx+2, :, :), 2*yz_buf, mpireal, proc_x_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vx1(nx-2:nx, :, :), 3*yz_buf, mpireal, proc_x_max, tag, &
        vx1(-2:0, :, :), 3*yz_buf, mpireal, proc_x_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vy1(1:2, :, :), 2*yz_buf, mpireal, proc_x_min, tag, &
        vy1(nx+1:nx+2, :, :), 2*yz_buf, mpireal, proc_x_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vy1(nx-2:nx, :, :), 3*yz_buf, mpireal, proc_x_max, tag, &
        vy1(-2:0, :, :), 3*yz_buf, mpireal, proc_x_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vz1(1:2, :, :), 2*yz_buf, mpireal, proc_x_min, tag, &
        vz1(nx+1:nx+2, :, :), 2*yz_buf, mpireal, proc_x_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vz1(nx-2:nx, :, :), 3*yz_buf, mpireal, proc_x_max, tag, &
        vz1(-2:0, :, :), 3*yz_buf, mpireal, proc_x_min, tag, comm, &
        status, errcode)

    CALL MPI_SENDRECV(vx1(:, ny-2:ny, :), 3*xz_buf, mpireal, proc_y_max, tag, &
        vx1(:, -2:0, :), 3*xz_buf, mpireal, proc_y_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vx1(:, 1:2, :), 2*xz_buf, mpireal, proc_y_min, tag, &
        vx1(:, ny+1:ny+2, :), 2*xz_buf, mpireal, proc_y_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vy1(:, ny-2:ny, :), 3*xz_buf, mpireal, proc_y_max, tag, &
        vy1(:, -2:0, :), 3*xz_buf, mpireal, proc_y_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vy1(:, 1:2, :), 2*xz_buf, mpireal, proc_y_min, tag, &
        vy1(:, ny+1:ny+2, :), 2*xz_buf, mpireal, proc_y_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vz1(:, ny-2:ny, :), 3*xz_buf, mpireal, proc_y_max, tag, &
        vz1(:, -2:0, :), 3*xz_buf, mpireal, proc_y_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vz1(:, 1:2, :), 2*xz_buf, mpireal, proc_y_min, tag, &
        vz1(:, ny+1:ny+2, :), 2*xz_buf, mpireal, proc_y_max, tag, comm, &
        status, errcode)

  END SUBROUTINE remap_v_MPI



  SUBROUTINE density_MPI

    CALL MPI_SENDRECV(rho(1:nx, 1:ny, 1:2), 2*nx*ny, mpireal, proc_z_min, tag, &
        rho(1:nx, 1:ny, nz+1:nz+2), 2*nx*ny, mpireal, proc_z_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(rho(1:nx, 1:ny, nz-1:nz), 2*nx*ny, mpireal, proc_z_max, tag, &
        rho(1:nx, 1:ny, -1:0), 2*nx*ny, mpireal, proc_z_min, tag, comm, &
        status, errcode)

    CALL MPI_SENDRECV(rho(1:2, 1:ny, :), 2*ny*(nz+4), mpireal, proc_x_min, tag, &
        rho(nx+1:nx+2, 1:ny, :), 2*ny*(nz+4), mpireal, proc_x_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(rho(nx-1:nx, 1:ny, :), 2*ny*(nz+4), mpireal, proc_x_max, tag, &
        rho(-1:0, 1:ny, :), 2*ny*(nz+4), mpireal, proc_x_min, tag, comm, &
        status, errcode)

    CALL MPI_SENDRECV(rho(:, ny-1:ny, :), 2*(nx+4)*(nz+4), mpireal, proc_y_max, tag, &
        rho(:, -1:0, :), 2*(nx+4)*(nz+4), mpireal, proc_y_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(rho(:, 1:2, :), 2*(nx+4)*(nz+4), mpireal, proc_y_min, tag, &
        rho(:, ny+1:ny+2, :), 2*(nx+4)*(nz+4), mpireal, proc_y_max, tag, comm, &
        status, errcode)

  END SUBROUTINE density_MPI

END MODULE mpiboundary
