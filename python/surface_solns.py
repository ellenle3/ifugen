import numpy as np
from scipy.optimize import newton


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

# Solutions for 3D conic

def conic_3d_off_axis_distance(c, alpha, beta):
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
    y0 = np.sin(alpha) / ( c * (1+np.cos(alpha)) )
    x0 = np.sin(beta) / ( c * (1+np.cos(beta)) )
    return x0, y0
    
def conic_3d_sag(x, y, c, k, alpha, beta, gamma):
    # Keep track of the angle, which determines which solution of
    # the quadratic is valid
    alpha = convert_angle(alpha) * np.pi/180
    beta = convert_angle(beta) * np.pi/180
    gamma = convert_angle(gamma) * np.pi/180
    if abs(gamma) <= np.pi/2:
        sgn = 1
    else:
        sgn = -1
        
    x0, y0 = conic_3d_off_axis_distance(c, alpha, beta)

    # Rotate about the y-axis
    cosg = np.cos(gamma)
    cosg2 = cosg*cosg
    sing = np.sin(gamma)
    sing2 = sing*sing

    asol = c*(sing2 + (k+1)*cosg2)
    bsol = 2*c*sing*(x*k*cosg - x0) - 2*cosg
    csol = c*k*x*x*sing2 - 2*x*sing + 2*c*x0*x*cosg + c*(x*x + x0*x0 + (y-y0)*(y-y0))

    # In regions where the roots are undefined, set the sag to 0 for drawing
    # purposes. We will not ray trace these regions
    return np.where(bsol**2-4*asol*csol < 0, 0, 2*csol/(-bsol + sgn*np.sqrt(bsol*bsol - 4*asol*csol)))

def conic_3d_transfer(xt, yt, l, m, n, c, k, alpha, beta, gamma):
    """Returns the transfer distance. Because the equation for t is a quadratic,
    there are two possible solutions. We almost always want the solution that
    corresponds to a smaller value of t (+ )

    Direction cosines must be normalized: l^2 + m^2 + n^2 = 1

    See Cheatham 1980.
    """
    alpha = convert_angle(alpha) * np.pi/180
    beta = convert_angle(beta) * np.pi/180
    gamma = convert_angle(gamma) * np.pi/180
    if abs(gamma) <= np.pi/2:
        sgn = 1
    else:
        sgn = -1
        
    x0, y0 = conic_3d_off_axis_distance(c, alpha, beta)

    cosg = np.cos(gamma)
    cosg2 = cosg*cosg
    sing = np.sin(gamma)
    sing2 = sing*sing

    dsol = c + c*k*(n*n*cosg2 + l*l*sing2 + 2*l*n*sing*cosg)
    fsol = c*l*(xt*(1+k*sing2)+x0*cosg) - l*sing + c*m*(yt-y0) + c*n*sing*(k*xt*cosg-x0) - n*cosg
    gsol = c*(xt*xt + x0*x0 + (yt-y0)*(yt-y0)) + 2*c*x0*xt*cosg + xt*sing*(c*k*xt*sing-2)

    # Ray missed this surface
    if fsol**2-dsol*gsol < 0:
        return np.nan
    
    print(dsol, fsol, gsol)
        
    return gsol/(-fsol + sgn*np.sqrt(fsol*fsol - dsol*gsol))

def conic_3d_surface_normal(x, y, c, k, alpha, beta, gamma, normalize):
    """Returns the surface normal vector components (the gradient).
    """
    # Check where the derivative is undefined!!!
    
    alpha = convert_angle(alpha) * np.pi/180
    beta = convert_angle(beta) * np.pi/180
    gamma = convert_angle(gamma) * np.pi/180

    if abs(gamma) <= np.pi/2:
        sgn = 1
    else:
        sgn = -1
    
    x0, y0 = conic_3d_off_axis_distance(c, alpha, beta)

    cosg = np.cos(gamma)
    cosg2 = cosg*cosg
    sing = np.sin(gamma)
    sing2 = sing*sing
    asol = c*(sing2 + (k+1)*cosg2)
    bsol = 2*c*sing*(x*k*cosg - x0) - 2*cosg
    csol = c*k*x*x*sing2 - 2*x*sing + 2*c*x0*x*cosg + c*(x*x + x0*x0 + (y-y0)*(y-y0))
    
    arg0 = bsol*bsol - 4*asol*csol
    # Partial derivatves are undefined - no surface normal
    if arg0 <= 0:
        return np.nan, np.nan, np.nan
    eta = np.sqrt(arg0)
    denom = -bsol + sgn*eta

    arg1 = 4 * (c*x*(1+k*sing2) + c*x0*cosg - sing) / denom
    arg2 = 4 * csol / (denom*denom)
    arg3 = -c*k*sing*cosg
    arg4 = c*k*bsol*sing*cosg - 2*asol*(c*x*(1+k*sing2) + c*x0*cosg - sing)
    dervx = arg1 - arg2 * (arg3 + sgn * arg4 / eta)

    arg1 = 4*c*(y-y0) / denom
    arg2 = 2*asol*csol / (eta*denom)
    dervy = arg1 * (1 + sgn*arg2)

    norm = 1
    if normalize:
        norm = np.sqrt(dervx**2 + dervy**2 + 1)

    # d/dz is always equal to -1, no need to calculate it
    # The sign in Zemax's example seems to be opposite of Shannon (1997)...
    return dervx / norm, dervy / norm, -1 / norm

def conic_3d_critical_xy(c, k, alpha, beta, gamma):
    """Computes where the d/dx and d/dy of the sag equals 0.

    Check where this is undefined!
    """
    # Tolerance for accepting root from secant method
    tol = 1e-13

    if abs(c) < 1e-13:
        # Curvature is very small, so this is basically a plane. In that case
        # there is no critical point because the derivative is a constant.
        return np.nan, np.nan
    
    alpha = convert_angle(alpha) * np.pi/180
    beta = convert_angle(beta) * np.pi/180
    gamma = convert_angle(gamma) * np.pi/180

    x0, y0 = conic_3d_off_axis_distance(c, alpha, beta)
    xc = newton(conic_3d_dervx, x0, tol=tol, args=(-y0,  c, k, alpha, beta, gamma))
    return xc, -y0

def conic_3d_dervx(x, y, c, k, alpha, beta, gamma):
    """Returns the partial derivative along x.
    """
    dervx, dervy, _ = conic_3d_surface_normal(x, y, c, k, alpha, beta, gamma, normalize=False)
    return dervx


# Solutions for planar surfaces

def tilted_plane_sag(x, y, c, k, alpha, beta, gamma):
    """
    c and k are not used, but are present so this function has the same number
    of parameters as conic_3d.
    """
    alpha = convert_angle(alpha) * np.pi/180
    beta = convert_angle(beta) * np.pi/180
    gamma = convert_angle(gamma) * np.pi/180

    sina = np.sin(alpha)
    cosa = np.cos(alpha)
    sinbg = np.sin(beta + gamma)
    cosbg = np.cos(beta + gamma)
    # Cap to prevent these from exploding
    if abs(cosa) < 1e-13: cosa = 1e-13
    if abs(cosbg) < 1e-13: cosbg = 1e-13

    return x * sinbg / (cosa*cosbg) - y * sina / cosa

def tilted_plane_critical_xy(c, k, alpha, beta, gamma):
    """Computes where the d/dx and d/dy of the sag equals 0.
    Planes do not have critical points.
    """
    # xc, yc
    return np.nan, np.nan

def tilted_plane_transfer(xt, yt, l, m, n, c, k, alpha, beta, gamma):
    """Returns the transfer distance.
    """
    alpha = convert_angle(alpha) * np.pi/180
    beta = convert_angle(beta) * np.pi/180
    gamma = convert_angle(gamma) * np.pi/180

    sina = np.sin(alpha)
    cosa = np.cos(alpha)
    sinbg = np.sin(beta + gamma)
    cosbg = np.cos(beta + gamma)
    # Cap to prevent these from exploding
    if abs(cosa) < 1e-13: cosa = 1e-13
    if abs(cosbg) < 1e-13: cosbg = 1e-13

    arg1 = xt * sinbg / (cosa * cosbg) - yt * sina / cosa
    arg2 = n - l * sinbg / (cosa * cosbg) + m * sina / cosa
    
    if abs(arg2) < 1e-13:
        return np.nanß
    return arg1 / arg2

def tilted_plane_surface_normal(x, y, c, k, alpha, beta, gamma, normalize):
    """Returns the surface normal vector components (the gradient).
    """
    alpha = convert_angle(alpha) * np.pi/180
    beta = convert_angle(beta) * np.pi/180
    gamma = convert_angle(gamma) * np.pi/180
    
    sina = np.sin(alpha)
    cosa = np.cos(alpha)
    sinbg = np.sin(beta + gamma)
    cosbg = np.cos(beta + gamma)
    # Cap to prevent these from exploding
    if abs(cosa) < 1e-13: cosa = 1e-13
    if abs(cosbg) < 1e-13: cosbg = 1e-13

    dervx = sinbg / (cosa * cosbg)
    dervy = -sina / cosa
    
    norm = 1
    if normalize:
        norm = np.sqrt(dervx**2 + dervy**2 + 1)

    # d/dz is always equal to -1, no need to calculate it
    # The sign in Zemax's example seems to be opposite of Shannon (1997)...
    return dervx / norm, dervy / norm, -1