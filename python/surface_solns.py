import numpy as np
from dataclasses import dataclass
from scipy.optimize import newton

@dataclass
class SliceParams:
    """Class for storing parameters of a single slice."""
    alpha: float
    beta: float
    gamma: float
    c: float
    k: float
    zp: float
    syx: float
    syz: float
    sxy: float
    sxz: float
    theta: float
    szx: float
    szy: float
    u: float
    
@dataclass
class RayIn:
    """Class for storing input ray parameters."""
    xt: float   # x-coordinate
    yt: float   # y-coordinate
    zt: float   # z-coordinate, usually = 0
    l: float    # Direction cosine x
    m: float    # Direction cosine y
    n: float    # Direction cosine z

@dataclass
class RayOut:
    """Class for storing output ray parameters."""
    xs: float   # Transferred x at the surface
    ys: float   # Transferred y at the surface
    zs: float   # Transferred z at the surface, which is the sag
    t: float    # Transfer distance
    ln: float   # Surface normal x
    mn: float   # Surface normal y
    nn: float   # Surface normal z

def convert_angle(t):
    """Converts an angle to be between -180 and 180 degrees.

    Parameters
    ----------
    t: float
        The angle in degrees.
    """
    t = t % 360
    if t > 180:
        return t - 360
    return t

def no_surface_normal(x, y, pslice, normalize):
    """Placeholder surface normal function when evaluating the sag.
    """
    return np.nan, np.nan, np.nan

def convert_ray_in_to_local(ray_in, pslice, transform_func):
    """Converts an input ray from global to local coordinates.
    """
    coords_global = np.array([ray_in.xt, ray_in.yt, ray_in.zt])
    cosines_global = np.array([ray_in.l, ray_in.m, ray_in.n])
    coords_local = transform_func(coords_global, pslice, -1, translate=True)
    cosines_local = transform_func(cosines_global, pslice, -1, translate=False)

    ray_in_local = RayIn(
        xt = coords_local[0],
        yt = coords_local[1],
        zt = coords_local[2],
        l = cosines_local[0],
        m = cosines_local[1],
        n = cosines_local[2]
    )
    return ray_in_local

def convert_ray_out_to_global(ray_out_local, pslice, transform_func):
    """Converts an output ray from local to global coordinates.
    """
    coords_local = np.array([ray_out_local.xs, ray_out_local.ys, ray_out_local.zs])
    normals_local = np.array([ray_out_local.ln, ray_out_local.mn, ray_out_local.nn])
    coords_global = transform_func(coords_local, pslice, 1, translate=True)
    normals_global = transform_func(normals_local, pslice, 1, translate=False)

    ray_out_global = RayOut(
        xs = coords_global[0],
        ys = coords_global[1],
        zs = coords_global[2],
        t = ray_out_local.t,
        ln = normals_global[0],
        mn = normals_global[1],
        nn = normals_global[2]
    )
    return ray_out_global

def slice_ray_trace(ray_in, pslice, transfer_dist_func, surface_normal_func, transform_func, normalize=True):
    # global to local
    ray_in_local = convert_ray_in_to_local(ray_in, pslice, transform_func)

    # ray trace in local coordinates
    t = transfer_dist_func(
        ray_in_local.xt, ray_in_local.yt, ray_in_local.zt,
        ray_in_local.l, ray_in_local.m, ray_in_local.n,
        pslice
    )
    xs = ray_in_local.xt + ray_in_local.l * t
    ys = ray_in_local.yt + ray_in_local.m * t
    zs = ray_in_local.zt + ray_in_local.n * t
    ln, mn, nn = surface_normal_func(
        xs, ys, pslice, normalize
    )

    ray_out_local = RayOut(
        xs = xs,
        ys = ys,
        zs = zs,
        t = t,
        ln = ln,
        mn = mn,
        nn = nn
    )

    # local to global
    ray_out_global = convert_ray_out_to_global(ray_out_local, pslice, transform_func)

    return ray_out_global

def slice_sag(x, y, pslice, transfer_dist_func, transform_func):
    # Evaluating the sag is a special case where the ray is coming in parallel 
    # to the z-axis
    ray_in = RayIn(
        xt = x,
        yt = y,
        zt = 0,
        l = 0,
        m = 0,
        n = 1
    )
    ray_out = slice_ray_trace(ray_in, pslice, transfer_dist_func,
                              no_surface_normal, transform_func, normalize=False)
    return ray_out.zs

def slice_surface_normal(x, y, pslice, transfer_dist_func, surface_normal_func, transform_func, normalize=True):
    ray_in = RayIn(
        xt = x,
        yt = y,
        zt = 0,
        l = 0,
        m = 0,
        n = 1
    )
    ray_out = slice_ray_trace(ray_in, pslice, transfer_dist_func,
                              surface_normal_func, transform_func, normalize)
    return ray_out.ln, ray_out.mn, ray_out.nn

# Solutions for 2D conic

def conic_2d_off_axis_distance(c, k, alpha, beta):
    """Off-axis distances.

    Parameters
    ----------
    c: float
        Curvature of the surface.
    alpha: float
        Angle in radians.
    beta: float
        Angle in radians.
    """
    # Put a minimum on the curvature to prevent x0, y0 from going to infinity.
    # In practice if c is close to 0, the user should be using the tilted plane
    # solutions instead.
    if (abs(c) < 1e-13): c = 1e-13

    if alpha == 0:
        y0 = 0
    else:
        tana = np.tan(alpha)
        num = ( (k-1) + np.sqrt(4 + tana*tana * (3 - k)) )
        den = tana*tana + (1 + k)
        y0 = tana / (2 * c) * num / den

    if beta == 0:
        x0 = 0
    else:
        tanb = np.tan(beta)
        num = ( (k-1) + np.sqrt(4 + tanb*tanb * (3 - k)) )
        den = tanb*tanb + (1 + k)
        x0 = tanb / (2 * c) * num / den

    # Make direction of effective rotation consistent with plane. For a reflective
    # surface where alpha and/or beta apply a rotation rather than an OAD, the
    # angles are effectively doubled.
    if c <= 0:
        x0 *= -1
    else:
        y0 *= -1

    return x0, y0

def conic_2d_off_axis_angle(x0, y0, c, k):
    if abs(c) < 1e-13:
        return 0
    sagx = c * x0 * x0 / (1 + np.sqrt(1 - (1 + k) * c * c * x0 * x0))
    denom = 1 / (2*c) - sagx
    beta = np.arctan(x0 / denom)

    sagy = c * y0 * y0 / (1 + np.sqrt(1 - (1 + k) * c * c * y0 * y0))
    denom = 1 / (2*c) - sagy
    alpha = np.arctan(y0 / denom)

    if c <= 0:
        beta *= -1
    else:
        alpha *= -1

    return alpha, beta

def conic_2d_transformation(coords, pslice, direction, translate=True):
    """Transforms the ray coordinates and direction cosines to the local
    coordinate system of the slice.

    Parameters
    ----------
    coords: array_like
        Coordinates to be transformed.
    pslice: SliceParams
        The parameters of the slice.
    direction: int
        1 for forward transformation (global to local), -1 for inverse.
    translate: bool, optional
        If True, apply the translations. If False, only apply rotations.
    """
    alpha = pslice.alpha
    beta = pslice.beta
    gamma = pslice.gamma
    zp = pslice.zp
    syx = pslice.syx
    syz = pslice.syz
    u = pslice.u

    alpha = convert_angle(alpha) * np.pi/180
    beta = convert_angle(beta) * np.pi/180
    gamma = convert_angle(gamma) * np.pi/180

    x0, y0 = conic_2d_off_axis_distance(pslice.c, pslice.k, alpha, beta)
    cosg = np.cos(gamma)
    sing = np.sin(gamma)

    # Transformation matrix
    T = np.array([[1, 0, 0, -u],
                  [0, 1, 0, 0],
                  [0, 0, 1, zp],
                  [0, 0, 0, 1]])
    
    Ry1 = np.array([[1, 0, 0, syx],
                    [0, 1, 0, 0],
                    [0, 0, 1, syz],
                    [0, 0, 0, 1]])
    Ry2 = np.array([[cosg, 0, sing, 0],
                    [0, 1, 0, 0],
                    [-sing, 0, cosg, 0],
                    [0, 0, 0, 1]])
    Ry3 = np.array([[1, 0, 0, -syx],
                    [0, 1, 0, 0],
                    [0, 0, 1, -syz],
                    [0, 0, 0, 1]])
    Ry = Ry1 @ Ry2 @ Ry3

    TOAD = np.array([[1, 0, 0, -x0],
                     [0, 1, 0, -y0],
                     [0, 0, 1, 0],
                     [0, 0, 0, 1]])
    
    Atot = T @ Ry @ TOAD
    append = 1 if translate else 0
    if direction == -1:
        Atot = np.linalg.inv(Atot)

    return (Atot @ np.append(coords, append))[:3]

def conic_2d_transfer(xt, yt, zt, l, m, n, pslice):
    cv = pslice.c
    k = pslice.k

    A = 1 + k*n*n
    B = xt*l + yt*m + zt*n*(1 + k) - n/cv
    C = xt*xt + yt*yt + zt*zt*(1 + k) - 2*zt/cv

    discrim = B*B - A*C
    if discrim < 0:
        return np.nan
    sgn = 1 if cv > 0 else -1
    
    return C / ( -B + sgn * np.sqrt(discrim) )

def conic_2d_surface_normal(x, y, pslice, normalize):
    cv = pslice.c
    k = pslice.k
    discrim = 1 - cv*cv*(1+k)*(x*x + y*y)
    if discrim < 0:
        return np.nan, np.nan, np.nan
    denom = np.sqrt(discrim)
    dervx = cv * x / denom
    dervy = cv * y / denom
    dervz = -1

    if normalize:
        norm = np.sqrt(dervx**2 + dervy**2 + dervz**2)
        return dervx / norm, dervy / norm, dervz / norm
    
    return dervx, dervy, dervz

def conic_2d_dervx(x, y, pslice):
    """Returns the partial derivative along x.
    """
    dervx, dervy, _ = slice_surface_normal(x, y, pslice, conic_2d_transfer,
                                           conic_2d_surface_normal, conic_2d_transformation,
                                           normalize=False)
    return dervx


def conic_2d_critical_xy(pslice):
    """Computes where the d/dx and d/dy of the sag equals 0.

    Check where this is undefined!
    """
    alpha = pslice.alpha
    beta = pslice.beta
    gamma = pslice.gamma
    c = pslice.c
    k = pslice.k
    zp = pslice.zp
    syx = pslice.syx
    syz = pslice.syz
    sxy = pslice.sxy
    sxz = pslice.sxz
    u = pslice.u

    # Tolerance for accepting root from secant method
    tol = 1e-13

    if abs(c) < 1e-13:
        # Curvature is very small, so this is basically a plane. In that case
        # there is no critical point because the derivative is a constant.
        return np.nan, np.nan
    
    alpha = convert_angle(alpha) * np.pi/180
    beta = convert_angle(beta) * np.pi/180
    gamma = convert_angle(gamma) * np.pi/180

    x0, y0 = conic_2d_off_axis_distance(c, k, alpha, beta)
    xc = newton(conic_2d_dervx, x0, tol=tol, args=(-y0, pslice))
    return xc, -y0

# Solutions for planar surfaces

def plane_transformation(coords, pslice, direction, translate=True):
    """Transforms the ray coordinates and direction cosines to the local
    coordinate system of the slice.
    """
    alpha = pslice.alpha
    beta = pslice.beta
    gamma = pslice.gamma
    zp = pslice.zp
    syx = pslice.syx
    syz = pslice.syz
    sxy = pslice.sxy
    sxz = pslice.sxz
    u = pslice.u

    alpha = convert_angle(alpha) * np.pi/180
    beta = convert_angle(beta) * np.pi/180
    gamma = convert_angle(gamma) * np.pi/180

    cosa = np.cos(alpha)
    sina = np.sin(alpha)
    cosbg = np.cos(beta + gamma)
    sinbg = np.sin(beta + gamma)

    # Transformation matrix
    T = np.array([[1, 0, 0, -u],
                  [0, 1, 0, 0],
                  [0, 0, 1, zp],
                  [0, 0, 0, 1]])
    
    Ry1 = np.array([[1, 0, 0, syx],
                    [0, 1, 0, 0],
                    [0, 0, 1, syz],
                    [0, 0, 0, 1]])
    Ry2 = np.array([[cosbg, 0, sinbg, 0],
                    [0, 1, 0, 0],
                    [-sinbg, 0, cosbg, 0],
                    [0, 0, 0, 1]])
    Ry3 = np.array([[1, 0, 0, -syx],
                    [0, 1, 0, 0],
                    [0, 0, 1, -syz],
                    [0, 0, 0, 1]])
    Ry = Ry1 @ Ry2 @ Ry3

    Rx1 = np.array([[1, 0, 0, 0],
                    [0, 1, 0, sxy],
                    [0, 0, 1, sxz],
                    [0, 0, 0, 1]])
    Rx2 = np.array([[1, 0, 0, 0],
                    [0, cosa, -sina, 0],
                    [0, sina, cosa, 0],
                    [0, 0, 0, 1]])
    Rx3 = np.array([[1, 0, 0, 0],
                    [0, 1, 0, -sxy],
                    [0, 0, 1, -sxz],
                    [0, 0, 0, 1]])
    Rx = Rx1 @ Rx2 @ Rx3

    Atot = T @ Ry @ Rx
    append = 1 if translate else 0
    if direction == -1:
        Atot = np.linalg.inv(Atot)

    return (Atot @ np.append(coords, append))[:3]

def plane_transfer(xt, yt, zt, l, m, n, pslice):
    """Returns the transfer distance.
    """
    if n == 0:
        return np.nan  # Ray is parallel to the plane
    return -zt / n

def plane_surface_normal(x, y, pslice, normalize):
    """Returns the surface normal vector components (the gradient).
    """
    # d/dz is always equal to -1, no need to calculate it
    # The sign in Zemax's example seems to be opposite of Shannon (1997)...
    return 0, 0, -1

def plane_critical_xy(pslice):
    """Computes where the d/dx and d/dy of the sag equals 0.
    Planes do not have critical points.
    """
    # xc, yc
    return np.nan, np.nan


# Solutions for generalized conic cylindrical surfaces. These are similar to the
# conicoid solutions.

def cylinder_transformation(coords, pslice, direction, translate=True):
    alpha = pslice.alpha
    beta = pslice.beta
    gamma = pslice.gamma
    zp = pslice.zp
    syx = pslice.syx
    syz = pslice.syz
    sxy = pslice.sxy
    sxz = pslice.sxz
    u = pslice.u

    alpha = convert_angle(alpha) * np.pi/180
    beta = convert_angle(beta) * np.pi/180
    gamma = convert_angle(gamma) * np.pi/180

    x0, y0 = conic_2d_off_axis_distance(pslice.c, pslice.k, alpha, beta)
    cosg = np.cos(gamma)
    sing = np.sin(gamma)
    cosa = np.cos(alpha)
    sina = np.sin(alpha)

    # Transformation matrix
    T = np.array([[1, 0, 0, -u],
                  [0, 1, 0, 0],
                  [0, 0, 1, zp],
                  [0, 0, 0, 1]])
    
    Ry1 = np.array([[1, 0, 0, syx],
                    [0, 1, 0, 0],
                    [0, 0, 1, syz],
                    [0, 0, 0, 1]])
    Ry2 = np.array([[cosg, 0, sing, 0],
                    [0, 1, 0, 0],
                    [-sing, 0, cosg, 0],
                    [0, 0, 0, 1]])
    Ry3 = np.array([[1, 0, 0, -syx],
                    [0, 1, 0, 0],
                    [0, 0, 1, -syz],
                    [0, 0, 0, 1]])
    Ry = Ry1 @ Ry2 @ Ry3

    Rx1 = np.array([[1, 0, 0, 0],
                    [0, 1, 0, sxy],
                    [0, 0, 1, sxz],
                    [0, 0, 0, 1]])
    Rx2 = np.array([[1, 0, 0, 0],
                    [0, cosa, -sina, 0],
                    [0, sina, cosa, 0],
                    [0, 0, 0, 1]])
    Rx3 = np.array([[1, 0, 0, 0],
                    [0, 1, 0, -sxy],
                    [0, 0, 1, -sxz],
                    [0, 0, 0, 1]])
    Rx = Rx1 @ Rx2 @ Rx3

    TOAD = np.array([[1, 0, 0, -x0],
                     [0, 1, 0, 0],
                     [0, 0, 1, 0],
                     [0, 0, 0, 1]])

    Atot = T @ Ry @ Rx @ TOAD
    append = 1 if translate else 0
    if direction == -1:
        Atot = np.linalg.inv(Atot)

    return (Atot @ np.append(coords, append))[:3]

def cylinder_transfer(xt, yt, zt, l, m, n, pslice):
    cv = pslice.c
    k = pslice.k

    A = 1 + k*n*n
    B = xt*l + zt*n*(1 + k) - n/cv
    C = xt*xt + zt*zt*(1 + k) - 2*zt/cv

    discrim = B*B - A*C
    if discrim < 0:
        return np.nan
    sgn = 1 if cv > 0 else -1
    
    return C / ( -B + sgn * np.sqrt(discrim) )

def cylinder_surface_normal(x, y, pslice, normalize):
    cv = pslice.c
    k = pslice.k
    discrim = 1 - cv*cv*(1+k)*x*x
    if discrim < 0:
        return np.nan, np.nan, np.nan
    denom = np.sqrt(discrim)
    dervx = cv * x / denom
    dervy = 0
    dervz = -1

    if normalize:
        norm = np.sqrt(dervx**2 + dervy**2 + dervz**2)
        return dervx / norm, dervy / norm, dervz / norm
    
    return dervx, dervy, dervz

def cylinder_critical_xy(pslice):
    # Strictly speaking, clyinders do not have critical points. But if the vertex
    # is within the bounds of the slice...
    return np.nan, np.nan