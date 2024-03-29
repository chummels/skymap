import numpy as np
from matplotlib import pyplot as plt
import yt
import healpy as hp
from unyt import kpc, cm
from yt.utilities.math_utils import ortho_find
from tqdm import tqdm

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

def get_cart_coords(n_side, radius, origin):
    """
    Returns the cartesian coordinates for a mollweide projection of n_side n and radius
    r centered on origin.  Returns arrays in x,y,z 
    """
    n_pix = hp.nside2npix(n_side)
    pix_indices = np.arange(n_pix)
    thetas, phis = hp.pix2ang(n_side, pix_indices)
    rs = np.full_like(thetas, radius)
    xs, ys, zs = sph2cart(phis, thetas, rs)
    xs *= kpc
    ys *= kpc
    zs *= kpc
    xs += origin[0]
    ys += origin[1]
    zs += origin[2]
    return xs, ys, zs

def calc_pixel(ds, field, start, end, length, data_source):
    """
    Calculate field value for a ray in ds with start:end 
    Use data_source to speed things up by a small factor

    Field value is just column density = sum(number density * path length)
    """
    #ray = ds.r[start:end]
    ray = ds.ray(start, end, data_source=data_source)
    return np.sum(ray[field] * ray['dts'].d * length)
    #return np.sum(ds.r[start:end][field] * ds.r[start:end]['dts'].d * length)

if __name__ == '__main__':

    fn = '/Users/chummels/src/yt-data/FIRE_M12i_ref11/snapshot_600.hdf5'
    #fn = '/Users/chummels/scratch/FIRE/m12i_res57000/output/snapshot_570.hdf5'
    ds = yt.load(fn)

    # define center and angular momentum vector
    _, center = ds.find_max(('gas', 'density'))
    sp = ds.sphere(center, (10, 'kpc'))
    L = sp.quantities.angular_momentum_vector()
    L, E1, E2 = ortho_find(L)

    # origin = solar location ~ 10 kpc out from center in disk
    #offset = E1 * 10 * kpc
    origin = center
    origin.convert_to_units('kpc')
    #origin = center + offset

    # define basics of projection
    n_side = 5
    radius = 10 * kpc
    field = ('gas', 'H_p0_number_density')
    #field = ('gas', 'El_number_density')

    xs, ys, zs = get_cart_coords(n_side, radius, origin)
    # rotate the healpy coords to align with the edge-on disk
    #CoM = np.array([0, 0, 0])
    #L = np.array([0, 0, 1])
    #P = np.array([[1, 0, 0], [1, 1, 1], [0, 0, 1], [1, 1, 0]])
    #V = np.array([[0, 1, 10], [-1, -1, -1], [1, 1, 1], [1, -1, -1]])
    #cyl_r = compute_cylindrical_radius(CoM, L, P, V)
    #cyl_r

    n_pix = len(xs)
    print("%i pixels" % n_pix)
    DMs = np.empty(n_pix)
    for i in tqdm(range(n_pix), "Generating Pixels"):
        end = ds.arr([xs[i], ys[i], zs[i]], 'kpc')
        DMs[i] = calc_pixel(ds, field, origin, end, radius, sp)
    #import cProfile
    #cProfile.runctx("calc_pixel(ds, field, origin, end, radius)", globals(), locals())

    res = np.degrees(hp.nside2resol(n_side))
    #DMs /= cm**2
    #DMs.convert_to_units('pc/cm**3')
    #hp.mollview(DMs, title="Angular Size: %.1f" % res, unit='DM [pc cm$^{-3}$]', norm='log')
    hp.mollview(DMs, title="Angular Size: %.1f" % res, unit='HI Column Density [cm$^{-2}$]', norm='log')
    plt.savefig('map_%02d.png' % n_side)

    #sp2 = ds.sphere(origin, radius)
    #p = yt.OffAxisProjectionPlot(ds, E1, field, north_vector=L, center=origin, 
    #                             width=2.5*radius, data_source=sp2)

    #p.set_zlim(field, 1e15, 1e23)

    #p.set_unit(field, 'pc/cm**3')
    #p.set_zlim(field, 1e0, 1e3)
    #p.save('projection.png')

    # rotate so plane aligns with plane of mollweide
    # vectorize ray calculations
    # use existing cython routines for ray calculations?
    # offset to solar position
    # run at higher res
    # investigate zack's code for speed clues

    # confirm behavior with cool gas content to ensure disk is aligned with mollweide projection
    # rotate so plane aligns with plane of mollweide
