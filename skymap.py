import numpy as np
from matplotlib import pyplot as plt
import yt
import healpy as hp

def appendSpherical_np(xyz):
    ptsnew = np.hstack((xyz, np.zeros(xyz.shape)))
    xy = xyz[:,0]**2 + xyz[:,1]**2
    ptsnew[:,3] = np.sqrt(xy + xyz[:,2]**2)
    ptsnew[:,4] = np.arctan2(np.sqrt(xy), xyz[:,2]) # for elevation angle defined from Z-axis down
    ptsnew[:,5] = np.arctan2(xyz[:,1], xyz[:,0])
    return ptsnew

def sph2cart(az, el, r):
    rcos_theta = r * np.cos(el)
    x = rcos_theta * np.cos(az)
    y = rcos_theta * np.sin(az)
    z = r * np.sin(el)
    return x, y, z

def cart2sph(x, y, z):
    hxy = np.hypot(x, y)
    r = np.hypot(hxy, z)
    el = np.arctan2(z, hxy)
    az = np.arctan2(y, x)
    return az, el, r

def cart2sphA(a):
    return cart2sph(a[0], a[1], a[2])

def DM_ray():
    pass

if __name__ == '__main__':
    test = False

    #pts = np.array([[1,0,0], [0,1,0], [0,0,1]])
    #print(appendSpherical_np(pts))
    #for pt in pts:
    #    print(cart2sphA(pt))

    # Create a grid of positions on a unit sphere
    #thetas, phis = np.mgrid[0:180:20j, 0:360:20j] # outputs are reversed rel to args?

    n_degrees = 10
    n_theta = int(180/n_degrees)
    n_phi = int(360/n_degrees)
    theta = np.linspace(0, 180, n_theta, endpoint=False)
    phi = np.linspace(0, 360, n_phi, endpoint=False)
    thetas = np.empty([n_theta, n_phi])
    phis = np.empty([n_theta, n_phi])
    for i,t in enumerate(theta):
        for j,p in enumerate(phi):
            thetas[i,j] = t
            phis[i,j] = p
    thetas = thetas.ravel()
    phis = phis.ravel()
    rs = np.full_like(thetas, 1)
    
    # convert to cartesian coords
    phis = np.radians(phis)
    thetas = np.radians(thetas)
    xs, ys, zs = sph2cart(phis, thetas, rs)

    # plot contours of unit sphere in cartesian space, just to confirm working.
    if test:
        zs = np.abs(zs)
        h = plt.contourf(xs, ys, zs)
        plt.axis('scaled')
        plt.colorbar()
        plt.savefig('test.png')

    # make mollweide projection
    NSIDE = 2
    NPIX = hp.nside2npix(NSIDE)
    print(NPIX)
    m = np.arange(NPIX)
    hp.mollview(m, nest=True, title="Mollview image NESTED")
    plt.savefig('moll.png')

    #data = np.arange(10)
    data = np.arange(len(phis))
    #theta = np.radians(np.arange(10, 110, 10))
    #phi = np.radians(np.linspace(0, 100, 10))
    theta = thetas
    phi = phis
    #print(data)
    #print(theta)
    print(phi)
    nside = 2
    print(np.degrees(hp.nside2resol(nside)))
    pixel_indices = hp.ang2pix(nside, theta, phi)
    print(pixel_indices)
    m = np.zeros(hp.nside2npix(nside))
    m[pixel_indices] = data
    hp.mollview(m)
    plt.savefig('map.png')
    #import pdb; pdb.set_trace()
