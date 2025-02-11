import numpy as np


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
        
    # Determine the off-axis distance
    if (abs(c) < 1e-13): c = 1e-13
        
    y0 = np.sin(alpha) / ( c * (1+np.cos(alpha)) )
    x0 = np.sin(beta) / ( c * (1+np.cos(beta)) )

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

def conic_3d_critical_xy(c, k, alpha, beta, gamma):
    """Computes where the d/dx and d/dy of the sag equals 0.

    Check where this is undefined!
    """
    alpha = convert_angle(alpha) * np.pi/180
    beta = convert_angle(beta) * np.pi/180
    gamma = convert_angle(gamma) * np.pi/180

    if abs(c) < 1e-13:
        # Curvature is very small, so this is basically a plane. In that case
        # there are no critical points because the derivative is a constant.
        return None, None, None

    # Determine the off-axis distance
    if (abs(c) < 1e-13): c = 1e-13
        
    y0 = np.sin(alpha) / ( c * (1+np.cos(alpha)) )
    x0 = np.sin(beta) / ( c * (1+np.cos(beta)) )

    yc = -1*y0

    sing = np.sin(gamma)
    sin3g = np.sin(3*gamma)
    cosg = np.cos(gamma)
    cos2g = np.cos(2*gamma)

    if k == -1:
        xc1 = -x0 + ( sing*(5 + c*c*yc*yc) + sin3g*(1 + c*c*yc*yc) ) / (cosg*cosg*-8*c)
        xc2 = None

    else:
        U = -1*sing / (c*(1+k)*(-2-k+k*cos2g))
        V = -2 - k + k*cos2g
        W = cosg*k*np.sqrt( V*(-2+2*c*c*yc*yc*(1+k)) )
        #W = cosg*k*np.sqrt( V*(-2+2*c*c*(1+k)+yc*yc) )
        
        xc1 = U * (V - W) - x0
        xc2 = U * (V + W) - x0

    return xc1, xc2, yc


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
        
    # Determine the off-axis distance.
    if (abs(c) < 1e-13): c = 1e-13
        
    y0 = np.sin(alpha) / ( c * (1+np.cos(alpha)) )
    x0 = np.sin(beta) / ( c * (1+np.cos(beta)) )

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


def conic_3d_surface_normal(x, y, c, k, alpha, beta, gamma):
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
    
    # Determine the off-axis distance
    if (abs(c) < 1e-13): c = 1e-13
    y0 = np.sin(alpha) / ( c * (1+np.cos(alpha)) )
    x0 = np.sin(beta) / ( c * (1+np.cos(beta)) )

    cosg = np.cos(gamma)
    cosg2 = cosg*cosg
    sing = np.sin(gamma)
    sing2 = sing*sing
    asol = c*(sing2 + (k+1)*cosg2)
    bsol = 2*c*sing*(x*k*cosg - x0) - 2*cosg
    csol = c*k*x*x*sing2 - 2*x*sing + 2*c*x0*x*cosg + c*(x*x + x0*x0 + (y-y0)*(y-y0))
    
    arg0 = bsol*bsol - 4*asol*csol

    # Sag is undefined
    if arg0 < 0:
        return np.nan, np.nan, np.nan

    eta = np.sqrt(arg0)
    denom = -bsol + sgn*eta
    
    # Partial derivative with respect to x
    arg1 = 4 * (c*x*(1+k*sing2) + c*x0*cosg - sing) / denom
    arg2 = 4 * csol / (denom*denom)
    arg3 = -c*k*sing*cosg
    arg4 = c*k*bsol*sing*cosg - 2*asol*(c*x*(1+k*sing2) + c*x0*cosg - sing)
    dervx = arg1 - arg2 * (arg3 + sgn * arg4 / eta)

    # ...with respect to y
    arg1 = 4*c*(y-y0) / denom
    arg2 = 2*asol*csol / (eta*denom)
    dervy = arg1 * (1 + sgn*arg2)

    norm = np.sqrt(dervx**2 + dervy**2 + 1)

    # d/dz is always equal to -1, no need to calculate it
    # The sign in Zemax's example seems to be opposite of Shannon (1997)...
    return dervx / norm, dervy / norm, -1 / norm


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
    return np.nan, np.nan, np.nan

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
        return np.nan
    return arg1 / arg2

def tilted_plane_surface_normal(x, y, c, k, alpha, beta, gamma):
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
    norm = np.sqrt(dervx**2 + dervy**2 + 1)

    # d/dz is always equal to -1, no need to calculate it
    # The sign in Zemax's example seems to be opposite of Shannon (1997)...
    return dervx / norm, dervy / norm, -1