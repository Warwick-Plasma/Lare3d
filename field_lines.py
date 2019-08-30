import matplotlib.pyplot as plt
import numpy as np

def load_data(i, data_dir='Data'):
    global data, dx, dy, dz
    import sdf_helper as sh
    data = sh.getdata(i, wkd=data_dir)
    x = data.Grid_Grid.data[0]
    y = data.Grid_Grid.data[1]
    z = data.Grid_Grid.data[2]
    dx = (x[-1] - x[0]) / xc.shape[0]
    dy = (y[-1] - y[0]) / yc.shape[0]
    dz = (z[-1] - z[0]) / zc.shape[0]

def get_loc(x,y,z):
    gridx = data.Grid_Grid.data[0]
    gridy = data.Grid_Grid.data[1]
    gridz = data.Grid_Grid.data[2]
    ix = np.digitize(x, gridx)
    nx = data.Grid_Grid.data[0].shape[0] - 1
    ny = data.Grid_Grid.data[1].shape[0] - 1
    nz = data.Grid_Grid.data[2].shape[0] - 1
    if (ix < 1 or ix > nx):
        return None
    iy = np.digitize(y, gridy)
    if (iy < 1 or iy > ny):
        return None
    iz = np.digitize(z, gridz)
    if (iz < 1 or iz > nz):
        return None
    return ix, iy, iz

def get_dir(x, y, z):
    a = get_loc(x, y, z)
    if (a is None):
        return None
    ix = a[0]
    iy = a[1]
    iz = a[2]
    bxm = bx[ix-1,iy-1,iz-1]
    bxp = bx[ix,iy-1,iz-1]
    bym = by[ix-1,iy-1,iz-1]
    byp = by[ix-1,iy,iz-1]
    bzm = bz[ix-1,iy-1,iz-1]
    bzp = bz[ix-1,iy-1,iz]
    # dx, dy, dz not correct for stretched grid
    deltax = (x - data.Grid_Grid.data[0][ix-1]) / dx
    deltay = (y - data.Grid_Grid.data[1][iy-1]) / dy
    deltaz = (z - data.Grid_Grid.data[2][iz-1]) / dz
    b1 = (1.0 - deltax) * bxm + deltax * bxp
    b2 = (1.0 - deltay) * bym + deltay * byp
    b3 = (1.0 - deltaz) * bzm + deltaz * bzp
    bmag = np.sqrt(b1**2 + b2**2 + b3**2)
    return (b1, b2, b3) / bmag

def calculate_streamline(x0, y0, z0, dt):
    locs = []
    pos = np.array([x0, y0, z0])
    locs.append(pos)
    nsteps = 0
    dir = get_dir(pos[0], pos[1], pos[2])
    while True:
        if (dir is None):
            break
        pos = pos + dt * dir
        locs.append(pos)
        dir = get_dir(pos[0], pos[1], pos[2])
        nsteps = nsteps + 1
        if (nsteps > 10000):
            break
    return locs


def plot_field_lines(x0, y0, r, n, i, dt, clear=True, data_dir='Data'):
    from mpl_toolkits.mplot3d import axes3d
    import random
    load_data(i, data_dir=data_dir)
    if (clear):
        plt.clf()
        fig = plt.figure(1)
        ax = fig.gca(projection='3d')
        b_min = 1e-6
        bslice = np.transpose(bz[:,:,0])
        s_min = np.min(bslice)
        if (s_min < 0):
          s_min = s_min * 1.2
        else:
          s_min = s_min * 0.8
        s_max = np.max(bslice) * 1.2
        import matplotlib.colors as color_mod
        X, Y = np.meshgrid(xc,yc)
        lt = 1e-3
        # Plot a contour of the vertical component at the bottom of the domain
        plt.contourf(X, Y, bslice, norm = color_mod.SymLogNorm(linthresh = lt, linscale=1.0, vmin=s_min, vmax=s_max), levels=np.linspace(s_min,s_max,100), cmap = plt.cm.seismic_r, alpha=0.5, offset=0.0)
        # Could at this stage plot a second contour of something
        bslice = np.transpose(bz[:,:,32])
        plt.contourf(X, Y, bslice, norm = color_mod.SymLogNorm(linthresh = lt, linscale=1.0, vmin=s_min, vmax=s_max), levels=np.linspace(s_min,s_max,100), cmap = plt.cm.seismic_r, alpha=0.5, offset=2.5)
    else:
        ax = plt.gca()
    for i in np.arange(n):
        alpha = 2.0 * np.pi * random.random()
        rad = r * random.random()
        x = rad * np.cos(alpha) + x0
        y = rad * np.sin(alpha) + y0
        a = calculate_streamline(x0=x, y0=y, z0=0.0, dt=dt)
        x, y, z = zip(*a)
        ax.plot3D(x, y, z, color='black')
    #todo make these adjust automatically to the domain size
    ax.set_xlim3d([data.Grid_Grid.data[0][0],data.Grid_Grid.data[0][-1]])
    ax.set_ylim3d([data.Grid_Grid.data[1][0],data.Grid_Grid.data[1][-1]])
    ax.set_zlim3d([-0.01*data.Grid_Grid.data[2][-1],data.Grid_Grid.data[2][-1]])
    plt.show(block=False)

    # Some old plotting things
    #ax.plot_surface(X, Y, np.transpose(bz[:,:,0]))
    #bslice = np.log(np.transpose(np.maximum(bz[:,:,0], b_min)))
    #s_max = np.max(bslice)
    #s_min = np.min(bslice)
    #cset = ax.contourf(X, Y, bslice, cmap = plt.cm.seismic, offset=0, alpha=0.5, levels = np.linspace(0.0, s_max, 50))
    #bslice = np.log(-np.transpose(np.minimum(bz[:,:,0], -b_min)))
    #s_max = np.max(bslice)
    #s_min = np.min(bslice)
    #cset = ax.contourf(X, Y, bslice, cmap = plt.cm.seismic_r, offset=0, alpha=0.5, levels = np.linspace(s_min, 0.0, 50))
