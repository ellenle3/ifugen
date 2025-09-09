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
    u: float


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

def calc_quadratic_derv(sgn, A, B, C, dA, dB, dC):
    """Returns the derivative of a quadratic:
        C / (-B + sgn * sqrt(B*B - A*C))

    Parameters
    ----------
    sgn: int
        Sign of the square root.
    A: float
        Coefficient of the quadratic term.
    dA: float
        Derivative of A.
    """
    discrim = B*B - A*C
    if discrim < 0:
        return np.nan  # Derivative is undefined
    sqrt_disc = np.sqrt(discrim)
    Eta = -B + sgn * sqrt_disc
    dEta = - dB + sgn * (2*B*dB - C*dA - A*dC) / 2 / sqrt_disc
    return (Eta * dC - C * dEta) / (Eta * Eta)

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
    Ry = Ry3 @ Ry2 @ Ry1

    TOAD = np.array([[1, 0, 0, -x0],
                     [0, 1, 0, -y0],
                     [0, 0, 1, 0],
                     [0, 0, 0, 1]])
    
    Atot = T @ Ry @ TOAD
    append = 1 if translate else 0
    if direction == -1:
        Atot = np.linalg.inv(Atot)

    return (Atot @ np.append(coords, append))[:3]
    
def conic_2d_sag(x, y, pslice):
    """Returns the sag of a rotationally symmetric conic.
    """
    alpha = pslice.alpha
    beta = pslice.beta
    gamma = pslice.gamma
    c = pslice.c
    k = pslice.k
    zp = pslice.zp
    syx = pslice.syx
    syz = pslice.syz
    u = pslice.u
    x-= u

    # Keep track of the angle, which determines which solution of
    # the quadratic is valid
    alpha = convert_angle(alpha) * np.pi/180
    beta = convert_angle(beta) * np.pi/180
    gamma = convert_angle(gamma) * np.pi/180
    if abs(gamma) <= np.pi/2:
        sgn = 1
    else:
        sgn = -1
    if k < -1: sgn *= -1
        
    x0, y0 = conic_2d_off_axis_distance(c, k, alpha, beta)
    
    v1 = syx + x0
    v2 = u - syx
    v3 = syz + zp

    cosg = np.cos(gamma)
    cos2g = np.cos(2*gamma)
    sing = np.sin(gamma)
    sin2g = np.sin(2*gamma)

    if k == -1:
        # Parabola handled separately
        A = c * sing**2

        B = -(
            cosg
            + c * (v2 + x) * cosg * sing
            + c * sing * (v1 + v3 * sing)
        )

        C = (
            -2 * syz
            + c * (v1**2 + (y + y0)**2)
            + 2 * (v3 + c * v1 * (v2 + x)) * cosg
            + c * (v2 + x)**2 * cosg**2
            - 2 * (v2 - c * v1 * v3 + x) * sing
            + c * v3**2 * sing**2
            + c * v3 * (v2 + x) * sin2g
        )

    else:
        A = c * (1 + k) * (2 + k + k * cos2g)

        B = - (1 + k) * (
            -2 * (-1 + c * (1 + k) * syz) * cosg
            + c * (
                (2 + k) * v3
                + k * v3 * cos2g
                + 2 * v1 * sing
                - k * (v2 + x) * sin2g
            )
        )

        C = (1 + k) * (
            -4 * syz
            + 2 * c * ((1 + k) * syz**2 + v1**2 + (y + y0)**2)
            + c * (2 + k) * (v2**2 + v3**2 + 2 * v2 * x + x**2)
            + 4 * (v3 - c * (1 + k) * syz * v3 + c * v1 * (v2 + x)) * cosg
            - c * k * (v2 - v3 + x) * (v2 + v3 + x) * cos2g
            + 4 * ((-1 + c * (1 + k) * syz) * v2 + c * v1 * v3 - x + c * (1 + k) * syz * x) * sing
            - 2 * c * k * v3 * (v2 + x) * sin2g
        )

    discrim = B*B - A*C

    # In regions where the roots are undefined, set the sag to 0 for drawing
    # purposes. We will not ray trace these regions
    return np.where( discrim < 0, 0, C /(-B + sgn*np.sqrt(discrim)))

def conic_2d_transfer(xt, yt, l, m, n, pslice):
    """Returns the transfer distance. Because the equation for t is a quadratic,
    there are two possible solutions. We almost always want the solution that
    corresponds to a smaller value of t (+ )

    Direction cosines must be normalized: l^2 + m^2 + n^2 = 1

    See Cheatham 1980.
    """
    alpha = pslice.alpha
    beta = pslice.beta
    gamma = pslice.gamma
    c = pslice.c
    k = pslice.k
    zp = pslice.zp
    syx = pslice.syx
    syz = pslice.syz
    u = pslice.u

    alpha = convert_angle(alpha) * np.pi/180
    beta = convert_angle(beta) * np.pi/180
    gamma = convert_angle(gamma) * np.pi/180
    if abs(gamma) <= np.pi/2:
        sgn = 1
    else:
        sgn = -1
    if k < -1: sgn *= -1

    x0, y0 = conic_2d_off_axis_distance(c, k, alpha, beta)

    v1 = syx + x0
    v2 = u - syx
    v3 = syz + zp

    cosg = np.cos(gamma)
    cos2g = np.cos(2*gamma)
    sing = np.sin(gamma)
    sin2g = np.sin(2*gamma)

    if k == -1:
        # Parabola handled separately
        D = c * (m**2 + (l * cosg - n * sing)**2)

        F = (
            c * m * (y0 + yt)
            + c * l * (v2 + xt) * cosg**2
            - sing * (l + c * n * v1 + c * n * v3 * sing)
            - cosg * (n - c * l * v1 + c * (-l * v3 + n * (v2 + xt)) * sing)
        )

        G = (
            -2 * syz
            + c * (v1**2 + (y0 + yt)**2)
            + 2 * (v3 + c * v1 * (v2 + xt)) * cosg
            + c * (v2 + xt)**2 * cosg**2
            - 2 * (v2 - c * v1 * v3 + xt) * sing
            + c * v3**2 * sing**2
            + c * v3 * (v2 + xt) * sin2g
        )
    
    else:
        D = 2 * c * k * (n * cosg + l * sing)**2 + 2 * c * (l**2 + m**2 + n**2)

        F = (
            c * (-(2 + k) * n * v3 + (2 + k) * l * (v2 + xt) + 2 * m * (y0 + yt))
            + 2 * (n * (-1 + c * (1 + k) * syz) + c * l * v1) * cosg
            - c * k * (n * v3 + l * (v2 + xt)) * cos2g
            + 2 * (l * (-1 + c * (1 + k) * syz) - c * n * v1) * sing
            + c * k * (-l * v3 + n * (v2 + xt)) * sin2g
        )

        G = (
            -4 * syz
            + 2 * c * ((1 + k) * syz**2 + v1**2 + (y0 + yt)**2)
            + c * (2 + k) * (v2**2 + v3**2 + 2 * v2 * xt + xt**2)
            + 4 * (v3 - c * (1 + k) * syz * v3 + c * v1 * (v2 + xt)) * cosg
            - c * k * (v2 - v3 + xt) * (v2 + v3 + xt) * cos2g
            + 4 * ((-1 + c * (1 + k) * syz) * v2 + c * v1 * v3 - xt + c * (1 + k) * syz * xt) * sing
            - 2 * c * k * v3 * (v2 + xt) * sin2g
        )
    discrim = F*F - D*G

    # Ray missed this surface
    if discrim < 0:
        return np.nan

    return G / ( -F + sgn * np.sqrt(discrim) )

def conic_2d_surface_normal(x, y, pslice, normalize):
    """Returns the surface normal vector components (the gradient).
    """

    alpha = pslice.alpha
    beta = pslice.beta
    gamma = pslice.gamma
    c = pslice.c
    k = pslice.k
    zp = pslice.zp
    syx = pslice.syx
    syz = pslice.syz
    u = pslice.u

    # Check where the derivative is undefined!!!
    
    alpha = convert_angle(alpha) * np.pi/180
    beta = convert_angle(beta) * np.pi/180
    gamma = convert_angle(gamma) * np.pi/180

    if abs(gamma) <= np.pi/2:
        sgn = 1
    else:
        sgn = -1
    if k < -1: sgn *= -1
    
    x0, y0 = conic_2d_off_axis_distance(c, k, alpha, beta)

    v1 = syx + x0
    v2 = u - syx
    v3 = syz + zp

    cosg = np.cos(gamma)
    cos2g = np.cos(2*gamma)
    sing = np.sin(gamma)
    sin2g = np.sin(2*gamma)

    if k == -1:
        # Parabola handled separately
        A = c * sing**2

        B = -(
            cosg
            + c * (v2 + x) * cosg * sing
            + c * sing * (v1 + v3 * sing)
        )

        C = (
            -2 * syz
            + c * (v1**2 + (y + y0)**2)
            + 2 * (v3 + c * v1 * (v2 + x)) * cosg
            + c * (v2 + x)**2 * cosg**2
            - 2 * (v2 - c * v1 * v3 + x) * sing
            + c * v3**2 * sing**2
            + c * v3 * (v2 + x) * sin2g
        )
        
        dAx = 0
        dBx = -c*cosg*sing
        dCx = -2 * sing + 2 * c * cosg * (v1 + (v2 + x) * cosg + v3 * sing)

        dAy = 0
        dBy = 0
        dCy = 2 * c * (y + y0)

    else: 
        A = c * (1 + k) * (2 + k + k * cos2g)

        B = - (1 + k) * (
            -2 * (-1 + c * (1 + k) * syz) * cosg
            + c * (
                (2 + k) * v3
                + k * v3 * cos2g
                + 2 * v1 * sing
                - k * (v2 + x) * sin2g
            )
        )

        C = (1 + k) * (
            -4 * syz
            + 2 * c * ((1 + k) * syz**2 + v1**2 + (y + y0)**2)
            + c * (2 + k) * (v2**2 + v3**2 + 2 * v2 * x + x**2)
            + 4 * (v3 - c * (1 + k) * syz * v3 + c * v1 * (v2 + x)) * cosg
            - c * k * (v2 - v3 + x) * (v2 + v3 + x) * cos2g
            + 4 * ((-1 + c * (1 + k) * syz) * v2 + c * v1 * v3 - x + c * (1 + k) * syz * x) * sing
            - 2 * c * k * v3 * (v2 + x) * sin2g
        )
        
        dAx = 0
        dBx = -c * k * (1 + k) * sin2g
        dCx = (1 + k) * (
            2 * c * (2 + k) * (v2 + x)
            + 4 * c * v1 * cosg
            - c * k * (v2 - v3 + x) * cos2g
            - c * k * (v2 + v3 + x) * cos2g
            - 4 * (-1 + c * (1 + k) * syz) * sing
            + 2 * c * k * v3 * sin2g
        )

        dAy = 0
        dBy = 0
        dCy = 4 * c * (1 + k) * (y + y0)
    
    dervx = calc_quadratic_derv(sgn, A, B, C, dAx, dBx, dCx)
    dervy = calc_quadratic_derv(sgn, A, B, C, dAy, dBy, dCy)

    norm = 1
    if normalize:
        norm = np.sqrt(dervx**2 + dervy**2 + 1)

    # d/dz is always equal to -1, no need to calculate it
    # The sign in Zemax's example seems to be opposite of Shannon (1997)...
    return dervx / norm, dervy / norm, -1 / norm

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

def conic_2d_dervx(x, y, pslice):
    """Returns the partial derivative along x.
    """
    dervx, dervy, _ = conic_2d_surface_normal(x, y, pslice, normalize=False)
    return dervx

# Solutions for planar surfaces

def plane_transformation(coords, dir_cosines, pslice, direction):
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
    if direction == 1:
        coords_out = Atot @ np.append(coords, 1)
        cos_out = Atot @ np.append(dir_cosines, 0)
    else:
        Atotinv = np.linalg.inv(Atot)
        coords_out = Atotinv @ np.append(coords, 1)
        cos_out = Atotinv @ np.append(dir_cosines, 0)

    return coords_out[:3], cos_out[:3]

def tilted_plane_sag(x, y, pslice):
    """
    c and k are not used, but are present so this function has the same number
    of parameters as conic_2d.
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
    tana = np.tan(alpha)
    cosbg = np.cos(beta + gamma)
    tanbg = np.tan(beta + gamma)

    # Cap to prevent these from exploding
    if abs(cosa) < 1e-13: cosa = 1e-13
    if abs(cosbg) < 1e-13: cosbg = 1e-13

    seca = 1 / cosa
    secbg = 1 / cosbg

    z = (
        secbg * (sxz - syz - sxz * seca + (y - sxy) * tana)
        - (x - syx + u) * tanbg
        + syz
        + zp
    )

    return z

def tilted_plane_transfer(xt, yt, l, m, n, pslice):
    """Returns the transfer distance.
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
    tana = np.tan(alpha)
    cosbg = np.cos(beta + gamma)
    tanbg = np.tan(beta + gamma)

    # Cap to prevent these from exploding
    if abs(cosa) < 1e-13: cosa = 1e-13
    if abs(cosbg) < 1e-13: cosbg = 1e-13

    seca = 1 / cosa
    secbg = 1 / cosbg

    arg1 = secbg * (sxz - syz - sxz * seca + (yt - sxy) * tana)
    arg2 = (xt - syx + u) * tanbg
    arg3 = syz + zp

    den = n - m * secbg * tana + l * tanbg
    if abs(den) < 1e-13:
        return np.nan

    return (arg1 - arg2 + arg3) / den

def tilted_plane_surface_normal(x, y, pslice, normalize):
    """Returns the surface normal vector components (the gradient).
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
    tana = np.tan(alpha)
    cosbg = np.cos(beta + gamma)
    tanbg = np.tan(beta + gamma)

    # Cap to prevent these from exploding
    if abs(cosa) < 1e-13: cosa = 1e-13
    if abs(cosbg) < 1e-13: cosbg = 1e-13

    secbg = 1 / cosbg

    dervx = -tanbg
    dervy = secbg * tana
    
    norm = 1
    if normalize:
        norm = np.sqrt(dervx**2 + dervy**2 + 1)

    # d/dz is always equal to -1, no need to calculate it
    # The sign in Zemax's example seems to be opposite of Shannon (1997)...
    return dervx / norm, dervy / norm, -1

def tilted_plane_critical_xy(pslice):
    """Computes where the d/dx and d/dy of the sag equals 0.
    Planes do not have critical points.
    """
    # xc, yc
    return np.nan, np.nan


# Solutions for generalized conic cylindrical surfaces

def generalized_cylinder_sag(x, y, c, k, alpha, beta, gamma):
    pass

def generalized_cylinder_transfer(xt, yt, l, m, n, c, k, alpha, beta, gamma):
    pass

def generalized_cylinder_surface_normal(x, y, c, k, alpha, beta, gamma, normalize):
    pass

def generalized_cylinder_critical_xy(c, k, alpha, beta, gamma):
    pass