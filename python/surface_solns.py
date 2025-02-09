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
    """
    """
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
    y = y + y0
    x = x + x0

    # Rotate about the y-axis
    cosg = np.cos(gamma)
    cosg2 = cosg*cosg
    sing = np.sin(gamma)
    sing2 = sing*sing
    asol = c*(sing2 + (k+1)*cosg2)
    bsol = -2*cosg*(x*k*c*sing + 1)
    csol = 2*x*sing + c*(x*x*cosg2 + y*y + (k+1)*x*x*sing2)

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
    yt = yt + y0
    xt = xt + x0

    cosg = np.cos(gamma)
    cosg2 = cosg*cosg
    sing = np.sin(gamma)
    sing2 = sing*sing
    
    xi = cosg2 + (1+k)*sing2
    asol = c*(sing2 + (1+k)*cosg2)
    dsol = asol*n*n - 2*l*n*k*c*sing*cosg + c*l*l*xi + m*m*c
    fsol = -xt*n*k*c*sing*cosg - n*cosg + l*sing + c*xt*l*xi + yt*m*c
    gsol = 2*xt*sing + c*xt*xt*xi + c*yt*yt

    # Ray missed this surface
    if fsol**2-dsol*gsol < 0:
        return np.nan
        
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
    y = y + y0
    x = x + x0

    sq2 = np.sqrt(2)
    cosg = np.cos(gamma)
    cos2g = np.cos(2*gamma)
    sing = np.sin(gamma)
    sin2g = np.sin(2*gamma)
    arg0 = 1 - c*c*( 2*(1+k)*x*x + (2+k)*y*y ) + (1-c*c*k*y*y)*cos2g - 4*c*x*sing

    if arg0 <= 0:
        return np.nan, np.nan, np.nan

    eta = np.sqrt(arg0)
    psi = 2 * c / ( eta*sq2 + sgn*2*cosg*(1+c*k*x*sing) )**2
    
    # Partial derivative with respect to x
    arg1 = c*y*y + c*x*x*(2+k-k*cos2g)/2 + 2*x*sing
    arg2 = k*sin2g - sgn*2*sq2*( c*(1+k)*x + sing )/eta
    arg3 = 2*c*(2+k)*x - 2*c*k*x*cos2g + 4*sing
    arg4 = 2*cosg + c*k*x*sin2g + sgn*eta*sq2
    dervx = -psi * arg1 * arg2 + arg3 / arg4

    # ...with respect to y
    if sgn == -1:
        dervy = -sq2 * c * y / eta    
    else:
        arg1 = 2*eta*sq2 + 4*cosg*(1+c*k*x*sing)
        arg2 = c*sq2*(2+k+k*cos2g) / eta
        arg3 = c*y*y + c*x*x*cosg*cosg + x*sing*(2+c*(1+k)*x*sing)
        dervy = y * psi * ( arg1 + arg2 * arg3 )

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