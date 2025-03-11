import math
import numpy as np
from dataclasses import dataclass
from surface_solns import *
from custom_slicer_helpers import *

@dataclass
class ImageSlicerParams:
    """Class for storing image slicer params."""
    custom: int
    surface_type: int
    n_each: int
    n_rows: int
    n_cols: int
    angle_mode: int
    dalpha: float
    dbeta: float
    dgamma: float
    gamma_offset: float
    dzp_col: float
    dzp_row: float
    dsyx: float
    dsyz: float
    dsxy: float
    dsxz: float
    alpha_cen: float
    beta_cen: float
    gamma_cen: float
    zp_cen: float
    syx_cen: float
    syz_cen: float
    sxy_cen: float
    sxz_cen: float
    dx: float
    dy: float
    c: float
    k: float
    gx_width: float
    gx_depth: float
    gy_width: float
    gy_depth: float

@dataclass
class RayIn:
    """Class for storing input ray parameters."""
    xt: float   # x-coordinate at z = 0
    yt: float   # y-coordinate at z = 0
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

def validate_slicer_params(p):
    """Returns True if all image slicer parameters are valid. Otherwise, modifies
    p to be valid and returns False.

    Parameters
    ----------
    p: ImageSlicerParams
        Image slicer parameters.
    
    Returns
    -------
    is_valid: bool
        True if all parameters are valid. False if p needed to be modified.
    """
    is_valid = True

    # Do not touch the custom flag
    if p.surface_type not in [0, 1]:
        p.surface_type = 0
        is_valid = False

    if p.n_cols < 1:
        p.n_cols = 1
        is_valid = False
    if p.n_rows < 1:
        p.n_rows = 1
        is_valid = False
    if p.n_each < 1:
        p.n_each = 1
        is_valid = False
    if p.angle_mode not in [0, 1, 2, 3]:
        p.angle_mode = 0
        is_valid = False

    # Angles can be whatever since we will convert them to be +/- 180

    if p.dx <= 0:
        p.dx = 1
        is_valid = False
    if p.dy <= 0:
        p.dy = 1
        is_valid = False
    if p.gx_width < 0:
        p.gx_width = 0
        is_valid = False
    if p.gy_width < 0:
        p.gy_width = 0
        is_valid = False
    
    # Gap depths can also be whatever, no limitations on c or k either
    return is_valid

def get_surface_funcs(c, p):
    """Returns sag and ray tracing functions for the surface type.
    """
    if c == 0:
        return tilted_plane_sag, tilted_plane_transfer, tilted_plane_surface_normal, tilted_plane_critical_xy
    else:
        return conic_2d_sag, conic_2d_transfer, conic_2d_surface_normal, conic_2d_critical_xy
    # if p.surface_type:
    
def make_image_slicer_params_from_custom(custom_slice_params):
    """Returns an ImageSlicerParams object described by custom slice parameters
    loaded from a TXT file.

    Parameters
    ----------
    custom_slice_params: array_like
        A one-dimensional array containing parameters each slice as well as the
        configuration of the image slicer.

    Returns
    -------
    p: ImageSlicerParams
        Image slicer parameters. Some parameters are set to zero and are unused
        for custom slicers.
    """
    p = ImageSlicerParams(
        custom = 1,
        n_each = 1,
        n_rows = custom_slice_params[0],
        n_cols = custom_slice_params[1],
        surface_type = custom_slice_params[2],
        dx = custom_slice_params[3],
        dy = custom_slice_params[4],
        gx_width = custom_slice_params[5],
        gx_depth = custom_slice_params[6],
        gy_width = custom_slice_params[7],
        gy_depth = custom_slice_params[8],

        # These params are unused for custom slicers, set to 0
        angle_mode = 0,
        dalpha = 0,
        dbeta = 0,
        dgamma = 0,
        gamma_offset = 0,
        dzp_col = 0,
        dzp_row = 0,
        dsyx = 0,
        dsyz = 0,
        dsxy = 0,
        dsxz = 0,
        alpha_cen = 0,
        beta_cen = 0,
        gamma_cen = 0,
        zp_cen = 0,
        syx_cen = 0,
        syz_cen = 0,
        sxy_cen = 0,
        sxz_cen = 0,
        c = 0,
        k = 0
    )
    return p

def get_slicer_size(p):
    """Returns the x and y dimensions of the image slicer.
    """
    # Return image slicer x, y dimensions
    n_slices = p.n_each*p.n_rows
    ysize = (n_slices*p.dy + (n_slices-1)*p.gy_width)
    xsize = (p.n_cols*p.dx + (p.n_cols-1)*p.gx_width)
    return xsize, ysize

def get_slicer_index(x, y, p):
    """Returns the column, slice index for a given x, y.
    """
    xsize, ysize = get_slicer_size(p)
    # Gets column and slice indices for a given x, y
    col_num = (x + xsize/2) // (p.dx + p.gx_width)
    slice_num = (y + ysize/2) // (p.dy + p.gy_width)
    return col_num, slice_num

def is_inside_slicer_gap(x, y, p):
    """Returns whether x, y is inside of a gap between columns (x) or slices (y).
    """
    xsize, ysize = get_slicer_size(p)
    col_num, slice_num = get_slicer_index(x, y, p)
    # Check if x, y is inside of a gap
    ygap_bot = (slice_num + 1) * p.dy + slice_num*p.gy_width - ysize/2
    ygap_top = ygap_bot + p.gy_width
    in_ygap = (y > ygap_bot and y <= ygap_top)
    
    xgap_left = (col_num + 1) * p.dx + col_num*p.gx_width - xsize/2
    xgap_right = xgap_left + p.gx_width
    in_xgap = (x > xgap_left and x <= xgap_right)
    return in_xgap, in_ygap

def get_paraxial_slice_index(p, active_x, active_y):
    """Returns the indices of the paraxial slice.
    """
    col_num = p.n_cols // 2
    if (p.n_cols % 2 == 0 and active_x):
        col_num += 1
        
    n_slices = p.n_each * p.n_rows    # slices per column
    slice_num = n_slices // 2
    
    if (n_slices % 2 == 0 and active_y):
        slice_num += 1

    return col_num, slice_num

def get_slice_params(slice_num, col_num, p, custom_slice_params):
    """Returns the angles alpha and beta for the slice number. Indexing starts at
    0 from the bottom (negative y direction) of the image slicer.
    """
    if p.custom:
        return get_slice_params_custom(slice_num, col_num, custom_slice_params)
    return get_slice_params_standard(slice_num, col_num, p) 

def get_slice_params_standard(slice_num, col_num, p):
    """Returns slice parameters for standard mode.

    Parameters
    ----------
    slice_num: int
        Slice index within an column.
    col_num: int
        Column index.
    
    Returns
    -------
    out: SliceParams
        Parameters of this slice.
    """
    # Get row number as well as the subindex of the slice on that row
    row_num = slice_num // p.n_each
    slice_num_row = slice_num - row_num*p.n_each
    # Set the angle beta for each row
    alpha, beta, gamma = 0, 0, 0
    if p.n_rows % 2 == 0:
        alpha = p.alpha_cen + p.dalpha*(row_num - (p.n_rows-1)/2)
        gamma_extra = p.gamma_offset * (row_num - (p.n_rows-1)/2)
    else:
        alpha = p.alpha_cen + p.dalpha*(row_num - p.n_rows//2)
        gamma_extra = p.gamma_offset * (row_num - p.n_rows//2)

    if p.n_cols % 2 == 0:
        beta = p.beta_cen + p.dbeta*(col_num - (p.n_cols-1)/2)
    else:
        beta = p.beta_cen + p.dbeta*(col_num - p.n_cols//2)

    # Get the angles of the bottom- and top-most slices of the central row
    # If n_rows is even, there are 2 rows straddling the x=0 center line
    # Set the "central row" to the one above the x-axis (+y direction)
    if p.n_each % 2 == 0:
        gamma_bot = p.gamma_cen - p.dgamma*(p.n_each-1)/2
        gamma_top = p.gamma_cen + p.dgamma*(p.n_each-1)/2
    else:
        gamma_bot = p.gamma_cen - p.dgamma*(p.n_each//2)
        gamma_top = p.gamma_cen + p.dgamma*(p.n_each//2)

    # Determine offsets in gamma. First, check whether the mode allows the extra
    # offsets to stack or not. If no, then set gamma_extra to repeat every 2 rows.
    if (p.angle_mode == 2 or p.angle_mode == 3):
        # gamma_offset should not stack
        if row_num % 2 == 0:
            gamma_extra = -p.gamma_offset/2
        else:
            gamma_extra = p.gamma_offset/2

    if (p.angle_mode == 0 or p.angle_mode == 2):
        # If the row is even, the angles are the same as the central row
        # If odd, the top/bottom angles are flipped
        if row_num % 2 == 0:
            gamma = gamma_bot + slice_num_row*p.dgamma + gamma_extra
        else:
            gamma = gamma_top - slice_num_row*p.dgamma + gamma_extra
    # Copy the same angle pattern as the central row
    elif (p.angle_mode == 1 or p.angle_mode == 3):
        gamma = gamma_bot + slice_num_row*p.dgamma + gamma_extra
    
    zp = 0
    syx = 0
    syz = 0
    sxy = 0
    sxz = 0
    return SliceParams(alpha, beta, gamma, p.c, p.k, zp, syx, syz, sxy, sxz)

def make_image_slicer(x, y, p, custom_slice_params):
    """Returns the sag of the image slicer at a given x, y.

    Parameters
    ----------
    x: float
        x-coordinate to evaluate.
    y: float
        y-coordinate to evaluate.
    p: ImageSlicerParams
        Parameters of the image slicer.
    custom_slice_params: array_like
        Custom slice parameters. If p.custom = 0 (disabled), the parameter is
        unused.
    
    Returns
    -------
    out: float
        Sag of the image slicer at x, y.
    """
    # Get dimensions of the image slicer
    xsize, ysize = get_slicer_size(p)
    # x, y coordinate is out of bounds (x=0, y=0 is in the middle of the slicer)
    if (abs(x) >= xsize/2 or abs(y) >= ysize/2):
        return np.nan
    # Figure out which slice and column we are on to determine gaps
    col_num, slice_num = get_slicer_index(x, y, p)
    # If inside a gap, return the gap depth instead of a curved surface
    in_xgap, in_ygap = is_inside_slicer_gap(x, y, p)
    if in_xgap:
        return p.gx_depth
    if in_ygap:
        return p.gy_depth
    # Inside a slice. From the slice number, determine the angles
    slice_params = get_slice_params(slice_num, col_num, p, custom_slice_params)
    sag_func, _, _, _ = get_surface_funcs(slice_params.c, p)
    return sag_func(x, y, slice_params)

def find_bounded_extremum(x0, y0, mode, p, custom_slice_params):
    """Returns the maximum or minimum of an individual slice. This function evaluates
    the sag at each of its four corner and at a critical point if it is within bounds.

    Parameters
    ----------
    x0: float
        Initial guess for x0.
    mode: bool
        True if max, False if min.
    p: ImageSlicerParams
        Parameters of the image slicer.
    custom_slice_params: array_like
        Custom slice parameters. If p.custom = 0 (disabled), the parameter is
        unused.
    
    Returns
    -------
    out: float
        Maximum or minimum of the slice.
    """
    # If x0, y0 are inside gaps then we don't need to go further
    in_xgap, in_ygap = is_inside_slicer_gap(x0, y0, p)
    if in_xgap: return p.gx_depth
    if in_ygap: return p.gy_depth

    # Not inside a gap. We need to compute the sag at the bounds and
    # critical point(s) of the slice.
    # First, compute the bounds of the current slice
    col_num, slice_num = get_slicer_index(x0, y0, p)
    xsize, ysize = get_slicer_size(p)
    slice_params = get_slice_params(slice_num, col_num, p, custom_slice_params)
    sag_func, junk1, junk2, critical_xy_func = get_surface_funcs(slice_params.c, p)

    xlo = col_num * (p.dx + p.gx_width) - xsize/2
    xhi = xlo + p.dx
    ylo = slice_num * (p.dy + p.gy_width) - ysize/2
    yhi = ylo + p.dy

    # Compute critical point
    xc, yc = critical_xy_func(slice_params)

    # There are up to 5 points to compare depending on whether the critical
    # point is within bounds
    zsolns = np.zeros(5)
    zsolns[0] = sag_func(xlo, ylo, slice_params)
    zsolns[1] = sag_func(xlo, yhi, slice_params)
    zsolns[2] = sag_func(xhi, ylo, slice_params)
    zsolns[3] = sag_func(xhi, yhi, slice_params)

    # Number of elements to compare so far
    n_compare = 4

    # Check whether critical point is in bounds. If yes, compute the sag and
    # add to the array of points to compare.
    if (yc >= ylo and yc <= yhi and xc >= xlo and xc <= xhi):
        zsolns[4] = sag_func(xc, yc, slice_params)
        n_compare += 1

    # Compare potential solutions to get the maximum or minimum
    if mode:
        return np.max(zsolns[:n_compare])
    return np.min(zsolns[:n_compare])

def find_global_extrema_slicer(p, custom_slice_params):
    """Returns the global max and min of the image slicer.
    """
    # Determine which slice the global extrema are on. To do this, roughly sample
    # the entire image slicer.
    xsize, ysize = get_slicer_size(p)
    nx = p.n_cols * 6
    ny = p.n_rows * p.n_each * 8
    xpts = np.linspace(-xsize/2, xsize/2, nx)
    ypts = np.linspace(-ysize/2, ysize/2, ny)

    x0_max, y0_max, z0_max = 0, 0, -np.inf
    x0_min, y0_min, z0_min = 0, 0, np.inf
    for x in xpts:
        for y in ypts:
            z = make_image_slicer(x, y, p, custom_slice_params)
            # Update the (x,y) corresponding to the max and min values
            if z > z0_max: x0_max, y0_max, z0_max = x, y, z
            elif z < z0_min: x0_min, y0_min, z0_min = x, y, z

    zmin = find_bounded_extremum(x0_min, y0_min, 0, p, custom_slice_params)
    zmax = find_bounded_extremum(x0_max, y0_max, 1, p, custom_slice_params)
    
    return zmin, zmax

def transfer_equation(t, xt, yt, l, m, n, p, custom_slice_params):
    """The roots of this function give the value of t.
    """
    xs = xt + t*l
    ys = yt + t*m
    zs = t*n
    sag = make_image_slicer(xs, ys, p, custom_slice_params)
    return sag - zs

def get_ray_bounds(ray_in, zmin, zmax, p):
    """Returns the starting and ending column, slice indices of the ray.

    Parameters
    ----------
    ray_in: RayIn
        Input ray to trace.
    zmin: float
        Global minimum of the image slicer.
    zmax: float
        Global maximum of the image slicer.
    p: ImageSlicerParams
        Parameters of the image slicer.
    """
    if abs(ray_in.n) < 1e-13:
        xsize, ysize = get_slicer_size(p)
        # Ray is moving perpendicular to the z-axis! This is an unusual case.
        # Set bounds on potential xs and ys to edges of the image slicer
        xmin, xmax = -xsize/2, xsize/2
        ymin, ymax = -ysize/2, ysize/2
        if ray_in.l < 0: xmin *= -1; xmax *= -1  # propagate right to left
        if ray_in.m < 0: ymin *= -1; ymax *= -1  # propagate bottom to top
        
    else:
        # Get maximum and minimum possible values of xs and ys
        tmin = zmin / ray_in.n
        xmin, ymin = ray_in.xt + tmin*ray_in.l, ray_in.yt + tmin*ray_in.m
        tmax = zmax / ray_in.n
        xmax, ymax = ray_in.xt + tmax*ray_in.l, ray_in.yt + tmax*ray_in.m

    # Get starting and ending rows and columns
    nc_min, ns_min = get_slicer_index(xmin, ymin, p)
    nc_max, ns_max = get_slicer_index(xmax, ymax, p)

    if xmin <= xmax: sgnc = 1
    else: sgnc = -1

    if ymin <= ymax: sgns = 1
    else: sgns = -1

    return nc_min, ns_min, nc_max, ns_max, sgnc, sgns

def is_ray_in_bounds(nc_min, ns_min, nc_max, ns_max, p):
    """Returns True if the ray is in bounds at least some of the time. This is done
    by checking whether or not the col, slice indices of the ray are always out
    of bounds of the image slicer.

    Parameters
    ----------
    nc_min: int
        Column index of the array at the plane defined by the global minimum, z = zmin.
    ns_min: int
        Slice index of the array at the plane defined by the global minimum, z = zmin.
    nc_max: int
        Column index of the array at the plane defined by the global maximum, z = zmax.
    ns_max: int
        Slice index of the array at the plane defined by the global maximum, z = zmax.
    p: ImageSlicerParams
        Parameters of the image slicer.
    
    Returns
    -------
    out: bool
        True if the ray is in bounds at least some of the time. False if the ray
        is never in bounds of the image slicer.
    """
    if (nc_min < 0 and nc_max < 0) or (nc_min >= p.n_cols and nc_max >= p.n_cols):
        # x-value of the ray is too high or low
        return False
    
    n_sperc = p.n_each * p.n_rows  # number of slices per column
    if (ns_min < 0 and ns_max < 0) or (ns_min >= n_sperc and ns_max >= n_sperc):
        # y-value of the ray is too high or low
        return False
    
    return True

def ray_trace_slicer(ray_in, zmin, zmax, trace_walls, p, custom_slice_params):
    """Computes the ray trace for the image slicer.

    Parameters
    ----------
    ray_in: RayIn
        Input ray to trace.
    zmin: float
        Global minimum of the image slicer.
    zmax: float
        Global maximum of the image slicer.
    trace_walls: bool
        If True, collisions with walls and gaps will be traced.
    p: ImageSlicerParams
        Parameters of the image slicer.
    custom_slice_params: array_like
        Custom slice parameters. If p.custom = 0 (disabled), the parameter is
        unused.
    
    Returns
    -------
    ray_out: RayOut
        Output ray parameters. These are the transferred coordinates, transfer
        distance, and surface normals.
    """
    # Tolerance for accepting the transfer distance as valid
    tol = 1e-12
    ray_out = RayOut(np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)

    # Get starting and ending rows and columns
    nc_min, ns_min, nc_max, ns_max, sgnc, sgns = get_ray_bounds(ray_in, zmin, zmax, p)
    if not is_ray_in_bounds(nc_min, ns_min, nc_max, ns_max, p):
        return ray_out
    
    # Ray is in bounds at least some of the time. Start from the min col and slice
    # indices and check solutions until we hit the max
    nc_test, ns_test = nc_min, ns_min
    x_test, y_test = ray_in.xt, ray_in.yt

    dcol = abs(nc_max - nc_test)        # Number of col indices left to iterate
    dslice = abs(ns_max - ns_test)      # Number of slice (row) indices left to iterate
    
    slice_iter = 0       # Keep track of whether we're on the first slice in a column to check walls

    xsize, ysize = get_slicer_size(p)
    xt, yt = ray_in.xt, ray_in.yt
    l, m, n = ray_in.l, ray_in.m, ray_in.n
    
    while dcol >= 0:
        
        # How many slices do we need to check on this column?
           
        # No column switching needed, check all remaining slices
        if (abs(l) < 1e-13 or dcol == 0):
            n_sforc = dslice

        # Ray can potentially switch columns
        else:
            # Crossover point between current column and next one
            x_cross = nc_test * (p.dx + p.gx_width) + (1 + sgnc) * p.dx / 2
            # y-intercept between ray and the next column to check
            y_cross = yt + (x_cross - x_test) * m / l
            # Number of slices to check
            n_sforc = math.ceil(abs(y_cross - y_test) / (p.dy + p.gy_width))
            x_test, y_test = x_cross, y_cross
            
        while (dslice >= 0 and n_sforc >= 0):
            
            if (nc_test >= 0 and ns_test >= 0):
                # Check if out of bounds
                pslice = get_slice_params(ns_test, nc_test, p, custom_slice_params)
                sag_func, transfer_dist_func, surf_normal_func, _ = get_surface_funcs(pslice.c, p)

                t = transfer_dist_func(xt, yt, l, m, n, pslice)
                result = transfer_equation(t, xt, yt, l, m, n, p, custom_slice_params)
                
                # Check whether the transfer distance of the current slice is a valid
                # zero of the transfer equation.
                if abs(result) < tol:
                    # Yes - found a solution!
                    xs = xt + t*l
                    ys = yt + t*m
                    zs = t*n
                    ln, mn, nn = surf_normal_func(xs, ys, pslice, True)

                    # WAIT - Is the solution inside of a gap? If we're unlucky and zs is
                    # equal to the gap depth then this may be the case.
                    in_xgap, in_ygap = is_inside_slicer_gap(xs, ys, p)
                    if (in_xgap or in_ygap):
                        dslice = -1
                        dcol = -1
                        break

                    ray_out.xs, ray_out.ys, ray_out.zs, ray_out.t = xs, ys, zs, t
                    ray_out.ln, ray_out.mn, ray_out.nn = ln, mn, nn
                    return ray_out

                # Check if the ray is hitting a wall after this slice
                if trace_walls:

                    # Going between columns. Skip if there are no columns left to iterate
                    if (slice_iter == 0 and abs(l) > 1e-13 and dcol > 0):

                        xnear = (nc_test + (1+sgnc)/2 ) * p.dx + nc_test * p.gx_width - xsize / 2
                        tnear = (xnear - xt) / l
                        ynear = yt + tnear*m
                        znear = tnear*n

                        xfar = xnear + sgnc * p.gx_width
                        tfar = (xfar - xt) / l
                        yfar = yt + tfar*m
                        zfar = tfar*n
                        print(nc_test, p.dx, xsize /2, xnear)

                        znear_slice = sag_func(xnear, ynear, pslice)
                        pslice = get_slice_params(ns_test, nc_test + 1*sgnc, p, custom_slice_params)
                        zfar_slice = sag_func(xfar, yfar, pslice)

                        # Did it hit a near wall? This is only possible if the gap depth
                        # protrudes further in -z than the near edge
                        # THIS IS ALSO POSSIBLE IF THE RAY APPROACHES FROM THE WRONG SIDE
                        # OF THE IMAGE SLICER! But this should never happen under normal circumstances.
                        if (znear_slice > p.gx_depth and p.gx_width > 0) and (znear < znear_slice and znear > p.gx_depth):
                            ray_out.xs, ray_out.ys, ray_out.zs, ray_out.t = xnear, ynear, znear, tnear
                            ray_out.ln, ray_out.mn, ray_out.nn = -1*sgnc, 0, 0
                            return ray_out

                        # Did it hit a far wall?
                        if p.gx_width == 0: zcompare = znear_slice
                        else: zcompare = p.gx_depth
                        if ((zfar > zfar_slice and zfar < zcompare) or (zfar < zfar_slice and zfar > zcompare)):
                            ray_out.xs, ray_out.ys, ray_out.zs, ray_out.t = xfar, yfar, zfar, tfar
                            ray_out.ln, ray_out.mn, ray_out.nn = -1*sgnc, 0, 0
                            return ray_out

                
                    # Not going between columns, check walls along y-axis instead. Same as above...
                    elif (abs(m) > 1e-13 and dslice > 0):
                        ynear = (ns_test + (1+sgns)/2) * p.dy + ns_test * p.gy_width - ysize / 2
                        tnear = (ynear - yt) / m
                        xnear = xt + tnear*l
                        znear = tnear*n
                        
                        yfar = ynear + sgns * p.gy_width
                        tfar = (yfar - yt) / m
                        xfar = xt + tfar*l
                        zfar = tfar*n
                        
                        znear_slice = sag_func(xnear, ynear, pslice)
                        pslice = get_slice_params(ns_test + 1*sgns, nc_test, p, custom_slice_params)
                        zfar_slice = sag_func(xfar, yfar, pslice)

                        # Did it hit a near wall? This is only possible if the gap depth
                        # protrudes further than the near edge
                        if (znear_slice > p.gy_depth and p.gy_width > 0) and (znear < znear_slice and znear > p.gy_depth):
                            ray_out.xs, ray_out.ys, ray_out.zs, ray_out.t = xnear, ynear, znear, tnear
                            ray_out.ln, ray_out.mn, ray_out.nn = 0, -1*sgns, 0
                            return ray_out

                        # Did it hit a far wall?
                        if p.gy_width == 0: zcompare = znear_slice
                        else: zcompare = p.gy_depth
                        if ((zfar > zfar_slice and zfar < zcompare) or (zfar < zfar_slice and zfar > zcompare)):
                            ray_out.xs, ray_out.ys, ray_out.zs, ray_out.t = xfar, yfar, zfar, tfar
                            ray_out.ln, ray_out.mn, ray_out.nn = 0, -1*sgns, 0
                            return ray_out

            # Not a solution - increment slice and try again.
            slice_iter += 1
            ns_test += sgns
            dslice -= 1

        # The last slice index needs to be checked again when crossing over to the next column
        slice_iter = 0
        ns_test -= sgns
        dslice += 1
        # Go to the next column
        nc_test += sgnc
        dcol -= 1

    # If none of the above worked then we probably hit a gap
    if (not trace_walls or abs(n) < 1e-13):
        return ray_out

    # Gaps between columns take precedence over gaps between slices in rows
    t = p.gx_depth / n
    if abs(transfer_equation(t, xt, yt, l, m, n, p, custom_slice_params)) < tol:
        ray_out.xs, ray_out.ys, ray_out.zs, ray_out.t = xt + t*l, yt + t*m, p.gx_depth, t
        ray_out.ln, ray_out.mn, ray_out.nn = 0, 0, -1
        return ray_out

    t = p.gy_depth / n
    if abs(transfer_equation(t, xt, yt, l, m, n, p, custom_slice_params)) < tol:
        ray_out.xs, ray_out.ys, ray_out.zs, ray_out.t = xt + t*l, yt + t*m, p.gy_depth, t
        ray_out.ln, ray_out.mn, ray_out.nn = 0, 0, -1
        return ray_out

    # If none of that worked then this is a bizarre edge case, e.g., the direction cosine is less than 1E-13
    # but it grazed off of a wall somehow. Sweep this ray under the rug and say that it missed...
    return ray_out    # Phew!