import numpy as np
from matplotlib import pyplot as plt
import yt
import healpy as hp

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

def get_cart_coords(n_side, radius):
    """
    Returns the cartesian coordinates for a mollweide projection of n_side n and radius
    r.  Returns arrays in x,y,z 
    """
    n_pix = hp.nside2npix(n_side)
    pix_indices = np.arange(n_pix)
    thetas, phis = hp.pix2ang(n_side, pix_indices)
    rs = np.full_like(thetas, radius)
    return sph2cart(phis, thetas, rs)

def calc_DM(x,y,z):
    """
    Calculate dispersion measure for a ray starting at origin and traversing to x,y,z
    """
    return np.random.rand()
    
if __name__ == '__main__':

    n_side = 1
    xs, ys, zs = get_cart_coords(n_side, 1)
    n_pix = len(xs)
    DMs = np.empty(n_pix)
    for i in range(n_pix):
        DMs[i] = calc_DM(xs[i], ys[i], zs[i])
    res = np.degrees(hp.nside2resol(n_side))
    hp.mollview(DMs, title="Angular Size: %.1f" % res)
    plt.savefig('map.png')
