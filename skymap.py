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
    #pts = np.array([[1,0,0], [0,1,0], [0,0,1]])
    #print(appendSpherical_np(pts))
    #for pt in pts:
    #    print(cart2sphA(pt))

    # Create a grid of positions on a unit circle
    phis, thetas = np.mgrid[0:180:3j, 0:360:6j] # outputs are reversed rel to args
    rs = np.full_like(thetas, 1)

    # convert to cartesian coords
    phis = np.radians(phis)
    thetas = np.radians(thetas)
    xs, ys, zs = sph2cart(phis, thetas, rs)

    rands = np.random.randn(*rs.shape)
    zz = np.sqrt(phis**2 + thetas**2)

    #h = plt.contourf(phis, thetas, zz)
    h = plt.contourf(xs, ys, zs)
    plt.axis('scaled')
    plt.colorbar()
    plt.savefig('test.png')

    # make aitoff projection
