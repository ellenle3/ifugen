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
    azps: float
    dsyx: float
    dsyz: float
    dsxy: float
    dsxz: float
    du: float
    alpha_cen: float
    beta_cen: float
    gamma_cen: float
    syx_cen: float
    syz_cen: float
    sxy_cen: float
    sxz_cen: float
    u_cen: float
    dx: float
    dy: float
    c: float
    k: float
    gx_width: float
    gx_depth: float
    gy_width: float
    gy_depth: float

@dataclass
class RayBounds:
    """Class for storing the bounds of a ray."""
    nc_min: int
    ns_min: int
    nc_max: int
    ns_max: int
    sgnc: int
    sgns: int
    xmin: float
    ymin: float
    xmax: float
    ymax: float

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

    # Angles can be whatever since we will convert them to be +/- 180. Changes
    # in axes of rotation are also unrestricted. (Please don't enter an imaginary
    # number or something else crazy.)

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

def get_surface_funcs(pslice, p):
    """Returns sag and ray tracing functions for the surface type.
    """

    if pslice.c == 0:
        return plane_transfer, plane_surface_normal, plane_critical_xy, plane_transformation
    elif p.surface_type == 1:
        return cylinder_transfer, cylinder_surface_normal, cylinder_critical_xy, cylinder_transformation
    else:
        return conic_2d_transfer, conic_2d_surface_normal, conic_2d_critical_xy, conic_2d_transformation
    
def make_image_slicer_params_from_custom(p_custom):
    """Returns an ImageSlicerParams object described by custom slice parameters
    loaded from a TXT file.

    Parameters
    ----------
    p_custom: array_like
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
        n_rows = p_custom[0],
        n_cols = p_custom[1],
        surface_type = p_custom[2],
        dx = p_custom[3],
        dy = p_custom[4],
        gx_width = p_custom[5],
        gx_depth = p_custom[6],
        gy_width = p_custom[7],
        gy_depth = p_custom[8],

        # These params are unused for custom slicers, set to 0
        angle_mode = 0,
        dalpha = 0,
        dbeta = 0,
        dgamma = 0,
        gamma_offset = 0,
        azps = 0,
        dsyx = 0,
        dsyz = 0,
        dsxy = 0,
        dsxz = 0,
        du = 0,
        alpha_cen = 0,
        beta_cen = 0,
        gamma_cen = 0,
        syx_cen = 0,
        syz_cen = 0,
        sxy_cen = 0,
        sxz_cen = 0,
        u_cen = 0,
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

def get_slicer_index(x, y, p, p_custom):
    """Returns the column, slice index for a given x, y.
    """
    xsize, ysize = get_slicer_size(p)
    # Gets column and slice indices for a given x, y
    slice_num = (y + ysize/2) // (p.dy + p.gy_width)
    row_num = slice_num // p.n_each
    # Column index depends on how much the row is shifted (u)
    u = get_u_for_row(row_num, p, p_custom)
    col_num = (x - u + xsize/2) // (p.dx + p.gx_width)
    return col_num, slice_num

def is_inside_slicer_gap(x, y, p, p_custom):
    """Returns whether x, y is inside of a gap between columns (x) or slices (y).
    """
    xsize, ysize = get_slicer_size(p)
    col_num, slice_num = get_slicer_index(x, y, p, p_custom)
    # Check if x, y is inside of a gap
    ygap_bot = (slice_num + 1) * p.dy + slice_num*p.gy_width - ysize/2
    ygap_top = ygap_bot + p.gy_width
    in_ygap = (y > ygap_bot and y <= ygap_top)
    
    u = get_u_for_row(slice_num // p.n_each, p, p_custom)
    xgap_left = (col_num + 1) * p.dx + col_num*p.gx_width - xsize/2 + u
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

def get_min_max_u(p, p_custom):
    """Returns the maximum and minimum u values of the image slicer.
    """
    # you probably want to store
    # this result so you don't have to do it every time a ray is traced
    if p.custom:
        u_all = p_custom[9: 9 + int(p.n_rows)]
        umax = np.max(u_all)
        umin = np.min(u_all)
        return umin, umax

    num_slices = p.n_each * p.n_rows
    if num_slices < 2:
        # Only one slice, so no need to compare
        return p.u_cen, p.u_cen

    if p.angle_mode in [0,1]:
        umin = get_u_for_row(0, p, p_custom)
        umax = get_u_for_row(p.n_rows - 1, p, p_custom)
    elif p.angle_mode in [2,3]:
        umin = get_u_for_row(0, p, p_custom)
        umax = get_u_for_row(1, p, p_custom)

    # Swap if alternating negative
    if p.du < 0:
        umax, umin = umin, umax

    return umin, umax

def get_slice_params(slice_num, col_num, p, p_custom):
    """Returns the angles alpha and beta for the slice number. Indexing starts at
    0 from the bottom (negative y direction) of the image slicer.
    """
    if p.custom:
        return get_slice_params_custom(slice_num, col_num, p_custom)
    return get_slice_params_standard(slice_num, col_num, p) 

def get_u_for_row(row_num, p, p_custom):
    """Returns the u value for a given row number, which is the shift along the
    x-axis. This needs to be computed seperately as it is used to calculate the
    column index of the slice.
    """
    row_num = int(row_num)
    if p.custom:
        return p_custom[9 + row_num]
    
    # angle_mode also determines whether u should stack because the way gamma is
    # constructed on the image slicer should determine how u is constructed on the
    # subpupil mirrors. This somewhat restricts which surfaces can be defined in
    # standard mode, but if this is insufficient the user can use custom mode instead.
    if (p.angle_mode == 2 or p.angle_mode == 3):
        if row_num % 2 == 0:
            u_extra = -p.du/2
        else:
            u_extra = p.du/2
    
    else:
        # Similar to get_slice_params_standard
        if p.n_rows % 2 == 0:
            u_extra = p.du * (row_num - (p.n_rows-1)/2)
        else:
            u_extra = p.du * (row_num - p.n_rows//2)

    
    return p.u_cen + u_extra

def calc_zp_from_gamma(gamma, c, k):
    """Calculates the amount to translate the slice along the z-axis (zp) to
    remove field curvature induced by the tilt gamma.
    """
    sing = np.sin( np.radians(gamma) )
    return (1/c) * sing*sing / 2

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
    u = get_u_for_row(row_num, p, [])

    row_mid = (p.n_rows - 1) / 2 if p.n_rows % 2 == 0 else p.n_rows // 2
    col_mid = (p.n_cols - 1) / 2 if p.n_cols % 2 == 0 else p.n_cols // 2
    slice_mid = (p.n_each - 1) / 2 if p.n_each % 2 == 0 else p.n_each // 2
    offset_row = row_num - row_mid
    offset_col = col_num - col_mid

    alpha = p.alpha_cen + p.dalpha * offset_row
    syx = p.syx_cen + p.dsyx * offset_row
    syz = p.syz_cen + p.dsyz * offset_row
    gamma_extra = p.gamma_offset * offset_row

    beta = p.beta_cen + p.dbeta * offset_col
    sxy = p.sxy_cen + p.dsxy * offset_col
    sxz = p.sxz_cen + p.dsxz * offset_col

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
        
    zp = p.azps * calc_zp_from_gamma(gamma, p.c, p.k)
    return SliceParams(alpha, beta, gamma, p.c, p.k, zp, syx, syz, sxy, sxz, u)

def make_image_slicer(x, y, p, p_custom):
    """Returns the sag of the image slicer at a given x, y.

    Parameters
    ----------
    x: float
        x-coordinate to evaluate.
    y: float
        y-coordinate to evaluate.
    p: ImageSlicerParams
        Parameters of the image slicer.
    p_custom: array_like
        Custom slice parameters. If p.custom = 0 (disabled), the parameter is
        unused.
    
    Returns
    -------
    out: float
        Sag of the image slicer at x, y.
    """
    # Check if out of bounds
    col_num, slice_num = get_slicer_index(x, y, p, p_custom)
    row_num = slice_num // p.n_each
    if (row_num < 0 or row_num >= p.n_rows) or (col_num < 0 or col_num >= p.n_cols):
        return np.nan

    # If inside a gap, return the gap depth unless at the outer edge (in which
    # case it should be considered out of bounds)
    in_xgap, in_ygap = is_inside_slicer_gap(x, y, p, p_custom)
    if in_xgap:
        if col_num == p.n_cols - 1:
            return np.nan
        return p.gx_depth
    if in_ygap:
        if slice_num == p.n_each * p.n_rows - 1:
            return np.nan
        return p.gy_depth
    
    # Inside a slice. From the indices, determine the slice parameters
    pslice = get_slice_params(slice_num, col_num, p, p_custom)
    transfer_dist_func, junk1, junk2, transform_func = get_surface_funcs(pslice, p)
    return slice_sag(x, y, pslice, transfer_dist_func, transform_func)

def find_bounded_extremum(x0, y0, mode, p, p_custom):
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
    p_custom: array_like
        Custom slice parameters. If p.custom = 0 (disabled), the parameter is
        unused.
    
    Returns
    -------
    out: float
        Maximum or minimum of the slice.
    """
    # If x0, y0 are inside gaps then we don't need to go further
    in_xgap, in_ygap = is_inside_slicer_gap(x0, y0, p, p_custom)
    if in_xgap: return p.gx_depth
    if in_ygap: return p.gy_depth

    # Not inside a gap. We need to compute the sag at the bounds and
    # critical point(s) of the slice.
    # First, compute the bounds of the current slice
    col_num, slice_num = get_slicer_index(x0, y0, p, p_custom)
    xsize, ysize = get_slicer_size(p)
    pslice = get_slice_params(slice_num, col_num, p, p_custom)
    transfer_dist_func, junk, critical_xy_func, transform_func = get_surface_funcs(pslice, p)

    u = get_u_for_row(slice_num // p.n_each, p, p_custom)
    xlo = col_num * (p.dx + p.gx_width) - xsize/2 + u
    xhi = xlo + p.dx
    ylo = slice_num * (p.dy + p.gy_width) - ysize/2
    yhi = ylo + p.dy

    # Compute critical point
    xc, yc = critical_xy_func(pslice)

    # There are up to 5 points to compare depending on whether the critical
    # point is within bounds
    zsolns = np.zeros(5)
    zsolns[0] = slice_sag(xlo, ylo, pslice, transfer_dist_func, transform_func)
    zsolns[1] = slice_sag(xlo, yhi, pslice, transfer_dist_func, transform_func)
    zsolns[2] = slice_sag(xhi, ylo, pslice, transfer_dist_func, transform_func)
    zsolns[3] = slice_sag(xhi, yhi, pslice, transfer_dist_func, transform_func)

    # Number of elements to compare so far
    n_compare = 4

    # Check whether critical point is in bounds. If yes, compute the sag and
    # add to the array of points to compare.
    if (yc >= ylo and yc <= yhi and xc >= xlo and xc <= xhi):
        zsolns[4] = slice_sag(xc, yc, pslice, transfer_dist_func, transform_func)
        n_compare += 1

    # Compare potential solutions to get the maximum or minimum
    if mode:
        return np.max(zsolns[:n_compare])
    return np.min(zsolns[:n_compare])

def find_global_extrema_slicer(umin, umax, p, p_custom):
    """Returns the global max and min of the image slicer.
    """
    # Determine which slice the global extrema are on. To do this, roughly sample
    # the entire image slicer.
    xsize, ysize = get_slicer_size(p)
    nx = int(p.n_cols * 6)
    ny = int(p.n_rows * p.n_each * 8)
    xpts = np.linspace(-xsize/2+umin, xsize/2+umax, nx)
    ypts = np.linspace(-ysize/2, ysize/2, ny)

    x0_max, y0_max, z0_max = 0, 0, -np.inf
    x0_min, y0_min, z0_min = 0, 0, np.inf
    for x in xpts:
        for y in ypts:
            z = make_image_slicer(x, y, p, p_custom)
            # Update the (x,y) corresponding to the max and min values
            if z > z0_max: x0_max, y0_max, z0_max = x, y, z
            elif z < z0_min: x0_min, y0_min, z0_min = x, y, z

    zmin = find_bounded_extremum(x0_min, y0_min, 0, p, p_custom)
    zmax = find_bounded_extremum(x0_max, y0_max, 1, p, p_custom)
    
    return zmin, zmax

def transfer_function(t, ray_in, p, p_custom):
    """The roots of this function give the value of t.
    """
    xs = ray_in.xt + t*ray_in.l
    ys = ray_in.yt + t*ray_in.m
    zs = t*ray_in.n
    sag = make_image_slicer(xs, ys, p, p_custom)
    return sag - zs

def get_ray_bounds(ray_in, umin, umax, zmin, zmax, p, p_custom):
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
        # Ray is moving parallel to the z-axis! This is an unusual case.
        # Set bounds on potential xs and ys to edges of the image slicer

        if abs(ray_in.l) < 1e-13:
            # Ray is moving along the y-axis only. Calculate where the ray intersects
            # the edges of the image slicer along the y-axis.
            ymin, ymax = -ysize/2, ysize/2
            if ray_in.m < 0: ymin, ymax = ymax, ymin  # propagate bottom to top
            tmin = (ymin - ray_in.yt) / ray_in.m
            xmin = ray_in.xt + tmin * ray_in.l
            tmax = (ymax - ray_in.yt) / ray_in.m
            xmax = ray_in.xt + tmax * ray_in.l
            
        else:
            # Ray is moving in the x-direction, at least partially. Repeat the
            # above steps but along the x-axis.
            xmin, xmax = -xsize/2 + umin, xsize/2 + umax
            if ray_in.l < 0: xmin, xmax = xmax, xmin  # propagate right to left
            tmin = (xmin - ray_in.xt) / ray_in.l
            ymin = ray_in.yt + tmin * ray_in.m
            tmax = (xmax - ray_in.xt) / ray_in.l
            ymax = ray_in.yt + tmax * ray_in.m
        
    else:
        # Get maximum and minimum possible values of xs and ys
        tmin = zmin / ray_in.n
        xmin, ymin = ray_in.xt + tmin*ray_in.l, ray_in.yt + tmin*ray_in.m
        tmax = zmax / ray_in.n
        xmax, ymax = ray_in.xt + tmax*ray_in.l, ray_in.yt + tmax*ray_in.m

    # Get starting and ending rows and columns
    nc_min, ns_min = get_slicer_index(xmin, ymin, p, p_custom)
    nc_max, ns_max = get_slicer_index(xmax, ymax, p, p_custom)

    if xmin <= xmax: sgnc = 1
    else: sgnc = -1

    if ymin <= ymax: sgns = 1
    else: sgns = -1

    return RayBounds(nc_min, ns_min, nc_max, ns_max, sgnc, sgns, xmin, ymin, xmax, ymax)

def is_ray_in_bounds(nc_min, ns_min, nc_max, ns_max, umin, umax, p):
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
    dcol = math.ceil( (umax - umin) / p.dx )
    col_min = -dcol
    col_max = p.n_cols + dcol
    if (nc_min < col_min and nc_max < col_min) or (nc_min >= col_max and nc_max >= col_max):
        # x-value of the ray is too low or high
        return False
    
    n_sperc = p.n_each * p.n_rows  # number of slices per column
    if (ns_min < 0 and ns_max < 0) or (ns_min >= n_sperc and ns_max >= n_sperc):
        # y-value of the ray is too low or high
        return False
    
    return True

def is_section_in_bounds(col_num, row_num, umin, umax, p):
    """Checks whether the section is within the range of x- and y-values over
    which the image slicer is defined.
    
    Note that this isn't as straightforward as checking whether nc and nr are within
    expected values since the rows may be shifted in x with respect to each other.
    For example, a ray may pass through an undefined region (nc = -1) but come back
    to a defined region after crossing through this section.

    """
    if row_num < 0 or row_num >= p.n_rows:
        return False

    # Column indices require us to consider how each row is shifted
    dcol = math.ceil( (umax - umin) / p.dx )
    if (col_num < -dcol or col_num >= p.n_cols + dcol):
        return False

    return True

def check_slice_solution(tol, ray_in, ns_test, nc_test, p, p_custom):
    """

    Returns
    -------
    ray_out: RayOut
        Output ray parameters. Values are nan if the ray missed.
    """
    ray_out = RayOut(np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)
    xt, yt = ray_in.xt, ray_in.yt
    l, m, n = ray_in.l, ray_in.m, ray_in.n

    # Check if this slice is a solution. Get params for this slice and
    # compute the transfer distance.
    pslice = get_slice_params(ns_test, nc_test, p, p_custom)
    transfer_dist_func, surf_normal_func, junk, transform_func = get_surface_funcs(pslice, p)

    ray_out_test = slice_ray_trace(ray_in, pslice, transfer_dist_func, surf_normal_func, transform_func)
    result = transfer_function(ray_out_test.t, ray_in, p, p_custom)

    # Is this a valid zero to the transfer function?
    if abs(result) < tol:
        # Yes - found a solution!

        # WAIT - Is the solution inside of a gap? If we're unlucky and zs is
        # equal to the gap depth then this may be the case.
        in_xgap, in_ygap = is_inside_slicer_gap(ray_out_test.xs, ray_out_test.ys, p, p_custom)
        if (in_xgap or in_ygap):
            return ray_out  # all nans
        
        # Not in a gap. This slice is the solution.
        ray_out = ray_out_test
        return ray_out
    
    return ray_out  # all nans

def check_ywall_collision(ray_in, ns_test, nc_test, sgns, p, p_custom):
    """Checks whether a ray has intersected a wall between the current slice and
    the next slice (between rows/slices).

    Parameters
    ----------
    ray_in: RayIn
        Input ray parameters.
    ns_test: int
        Slice index of the current slice.
    nc_test: int
        Column index of the current slice.
    p: ImageSlicerParams
        Paramters of the image slicer.

    Returns
    -------
    ray_out: RayOut
        Output ray parameters. Values are nan if the ray missed.
    """
    ray_out = RayOut(np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)
    xt, yt = ray_in.xt, ray_in.yt
    l, m, n = ray_in.l, ray_in.m, ray_in.n

    ysize = get_slicer_size(p)[1]

    if (abs(l) <= 1e-13 and abs(m) <= 1e-13 and p.gy_width > 0):
        # Ray is going straight down z-axis. If you're checking for a wall collision
        # here then it must have gone into a gap.
        t = p.gy_depth / n
        ray_out.xs, ray_out.ys, ray_out.zs, ray_out.t = xt + t*l, yt + t*m, p.gy_depth, t
        ray_out.ln, ray_out.mn, ray_out.nn = 0, 0, -1
        return ray_out

    if (abs(m) > 1e-13):

        # Ray coordinates on near wall
        ynear = ((1 + sgns) / 2) * p.dy + ns_test * (p.dy + p.gy_width) - ysize / 2
        tnear = (ynear - yt) / m
        xnear = xt + tnear * l
        znear = tnear * n
        
        # Ray coordinates on far wall
        yfar = ynear + sgns * p.gy_width
        tfar = (yfar - yt) / m
        xfar = xt + tfar * l
        zfar = tfar * n
        
        # Sag of near and far slices
        pslice_near = get_slice_params(ns_test, nc_test, p, p_custom)
        transfer_dist_func, junk1, junk2, transform_func = get_surface_funcs(pslice_near, p)
        znear_slice = slice_sag(xnear, ynear, pslice_near, transfer_dist_func, transform_func)

        pslice_far = get_slice_params(ns_test + 1*sgns, nc_test, p, p_custom)
        transfer_dist_func, junk1, junk2, transform_func = get_surface_funcs(pslice_far, p)
        zfar_slice = slice_sag(xfar, yfar, pslice_far, transfer_dist_func, transform_func)

        # Did it hit a near wall? This is only possible if the gap depth
        # protrudes further than the near edge
        if (znear_slice > p.gy_depth and p.gy_width > 0) and (znear <= znear_slice and znear >= p.gy_depth):
            ray_out.xs, ray_out.ys, ray_out.zs, ray_out.t = xnear, ynear, znear, tnear
            ray_out.ln, ray_out.mn, ray_out.nn = 0, -1*sgns, 0
            return ray_out

        # Did it hit a far wall?
        if p.gy_width == 0: zcompare = znear_slice
        else: zcompare = p.gy_depth
        if ((zfar >= zfar_slice and zfar <= zcompare) or (zfar <= zfar_slice and zfar >= zcompare)):
            ray_out.xs, ray_out.ys, ray_out.zs, ray_out.t = xfar, yfar, zfar, tfar
            ray_out.ln, ray_out.mn, ray_out.nn = 0, -1*sgns, 0
            return ray_out
        
        # Neither - did it hit a gap?
        if (zfar >= p.gy_depth and abs(n) > 1e-13 and p.gy_width > 0):
            t = p.gy_depth / n
            ray_out.xs, ray_out.ys, ray_out.zs, ray_out.t = xt + t*l, yt + t*m, p.gy_depth, t
            ray_out.ln, ray_out.mn, ray_out.nn = 0, 0, -1
            return ray_out
    
    # Ray missed...
    return ray_out

def check_xwall_collision(ray_in, ns_test, nc_test, sgnc, p, p_custom):
    """Analagous to check_ywall_collision, but for the x walls (between columns).
    """
    ray_out = RayOut(np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)
    xt, yt = ray_in.xt, ray_in.yt
    l, m, n = ray_in.l, ray_in.m, ray_in.n

    xsize = get_slicer_size(p)[0]

    if (abs(l) <= 1e-13 and abs(m) <= 1e-13 and p.gx_width > 0):
        # Ray is going straight down z-axis. If you're checking for a wall collision
        # here then it must have gone into a gap.
        t = p.gx_depth / n
        ray_out.xs, ray_out.ys, ray_out.zs, ray_out.t = xt + t*l, yt + t*m, p.gx_depth, t
        ray_out.ln, ray_out.mn, ray_out.nn = 0, 0, -1
        return ray_out

    if abs(l) > 1e-13:
        
        u = get_u_for_row(ns_test // p.n_each, p, p_custom)
        # Ray coordinates on near wall
        xnear = ((1 + sgnc) / 2) * p.dx + nc_test * (p.dx + p.gx_width) - xsize / 2 + u
        tnear = (xnear - xt) / l
        ynear = yt + tnear * m
        znear = tnear * n

        # Ray coordinates on far wall
        xfar = xnear + sgnc * p.gx_width
        tfar = (xfar - xt) / l
        yfar = yt + tfar * m
        zfar = tfar * n

        # Sag of near and far slices
        pslice_near = get_slice_params(ns_test, nc_test, p, p_custom)
        transfer_dist_func, junk1, junk2, transform_func = get_surface_funcs(pslice_near, p)
        znear_slice = slice_sag(xnear, ynear, pslice_near, transfer_dist_func, transform_func)

        pslice_far = get_slice_params(ns_test, nc_test + 1*sgnc, p, p_custom)
        transfer_dist_func, junk1, junk2, transform_func = get_surface_funcs(pslice_far, p)
        zfar_slice = slice_sag(xfar, yfar, pslice_far, transfer_dist_func, transform_func)

        # Did it hit a near wall? This is only possible if the gap depth
        # protrudes further in -z than the near edge
        # THIS IS ALSO POSSIBLE IF THE RAY APPROACHES FROM THE WRONG SIDE
        # OF THE IMAGE SLICER! But this should never happen under normal circumstances.
        if (znear_slice > p.gx_depth and p.gx_width > 0) and (znear <= znear_slice and znear >= p.gx_depth):
            ray_out.xs, ray_out.ys, ray_out.zs, ray_out.t = xnear, ynear, znear, tnear
            ray_out.ln, ray_out.mn, ray_out.nn = -1*sgnc, 0, 0
            return ray_out

        # Did it hit a far wall?
        if p.gx_width == 0: zcompare = znear_slice
        else: zcompare = p.gx_depth
        if ((zfar >= zfar_slice and zfar <= zcompare) or (zfar <= zfar_slice and zfar >= zcompare)):
            ray_out.xs, ray_out.ys, ray_out.zs, ray_out.t = xfar, yfar, zfar, tfar
            ray_out.ln, ray_out.mn, ray_out.nn = -1*sgnc, 0, 0
            return ray_out
        
        # Neither - did it hit a gap?
        if (zfar >= p.gx_depth and abs(n) > 1e-13 and p.gx_width > 0):
            t = p.gx_depth / n
            ray_out.xs, ray_out.ys, ray_out.zs, ray_out.t = xt + t*l, yt + t*m, p.gx_depth, t
            ray_out.ln, ray_out.mn, ray_out.nn = 0, 0, -1
            return ray_out
                
    # Ray missed...
    return ray_out

def calc_next_coords(ray_in, sgnc, sgns, nc_test, nr_test, x_test, y_test, xmax, ymax, p, p_custom):
    """

    Returns
    -------
    x_next, y_next: float
        Coordinates to go next.
    code: int
        0 if next coordinate corresponds to zmax, 1 if next column, 2 if next row.
    """
    l, m = ray_in.l, ray_in.m

    # Ray in coming in straight-on in the z-direction. No need to change sections.
    if (abs(l) < 1e-13 and abs(m) < 1e-13):
        return xmax, ymax, 0
    
    # Avoid division by zero. Doesn't matter what exactly we set this to because
    # the non-infinitesimal value of l or m will win out.
    if abs(l) < 1e-13:
        l = 1e-13
    if abs(m) < 1e-13:
        m = 1e-13
    xsize, ysize = get_slicer_size(p)
    
    # x- and y-coordinates of crossover points to the next section, incremented
    # either column- or row-wise.
    u = get_u_for_row(nr_test, p, p_custom)
    x_nextcol = (nc_test + (1 + sgnc) / 2) * (p.dx + p.gx_width) - xsize / 2 + u
    y_nextcol = y_test + (x_nextcol - x_test) * m / l

    y_nextrow = (nr_test + (1 + sgns) / 2) * p.n_each * (p.dy + p.gy_width) - ysize / 2
    x_nextrow = x_test + (y_nextrow - y_test) * l / m

    # Whether the next section is incremented column- or row-wise depends on which
    # distance is smaller
    dtocol = np.sqrt((x_nextcol - x_test)**2 + (y_nextcol - y_test)**2)
    dtorow = np.sqrt((x_nextrow - x_test)**2 + (y_nextrow - y_test)**2)
    dtomax = np.sqrt((xmax - x_test)**2 + (ymax - y_test)**2)

    if dtomax <= dtocol and dtomax <= dtorow:
        # Stay in this section
        return xmax, ymax, 0
    elif dtocol <= dtorow and dtocol <= dtomax:
        # Go to next column
        return x_nextcol, y_nextcol, 1
    else:
        # Go to next row
        return x_nextrow, y_nextrow, 2
    
def calc_num_slices_to_check(sgnc, sgns, nc_test, nr_test, x_test, y_test, x_next, y_next, code, p, p_custom):

    # Number of slices to check can be found by taking the difference between
    # starting and ending slice indices. For nc2 and ns2, subtract an infinitesimal
    # value from the coordinates to ensure that they are in the same section as (x_test, y_test).
    ns1 = get_slicer_index(x_test, y_test, p, p_custom)[1]
    ns2 = get_slicer_index(x_next - 1e-15*sgnc, y_next - 1e-15*sgns, p, p_custom)[1]
    n_stocheck = int( abs(ns2 - ns1) + 1 )
    # Cap the number of slices to the number of slices within a single section
    n_stocheck = min(n_stocheck, p.n_each)

    match code:
        case 0:
            # Stay in this section
            nc_new = nc_test
            nr_new = nr_test
        case 1:
            # Go to next column
            nc_new = nc_test + sgnc
            nr_new = nr_test
        case 2:
            # Go to next row
            nc_new = nc_test
            nr_new = nr_test + sgns
            
    return n_stocheck, nc_new, nr_new

def is_last_slice_in_section(ns_test, sgns, p):
    """ Returns True if ns_test is the last slice in the current section depending
    on the direction of iteration.
    """
    if sgns > 0:
        return ns_test % p.n_each == p.n_each - 1
    elif sgns < 0:
        return ns_test % p.n_each == 0
    return False


def ray_trace_slicer(ray_in, zmin, zmax, umin, umax, trace_walls, p, p_custom):
    """
    """
    # Tolerance for accepting the transfer distance as valid
    tol = 1e-12
    ray_out = RayOut(np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)
    # Get starting and ending rows and columns
    bounds = get_ray_bounds(ray_in, umin, umax, zmin, zmax, p, p_custom)
    nc_min, ns_min, nc_max, ns_max = bounds.nc_min, bounds.ns_min, bounds.nc_max, bounds.ns_max
    sgnc, sgns = bounds.sgnc, bounds.sgns
    xmin, ymin = bounds.xmin, bounds.ymin
    xmax, ymax = bounds.xmax, bounds.ymax

    if not is_ray_in_bounds(nc_min, ns_min, nc_max, ns_max, umin, umax, p):
        return ray_out
    # Ray is in bounds at least some of the time. Start from the min col and slice
    # indices and check solutions until we hit the max
    nc_test, ns_test = nc_min, ns_min
    nr_test = ns_min // p.n_each
    x_test, y_test = xmin, ymin
    
    # If starting index isn't in bounds, set it to where it will be in bounds
    while not is_section_in_bounds(nc_test, nr_test, umin, umax, p):
        x_next, y_next, code = calc_next_coords(ray_in, sgnc, sgns, nc_test, nr_test, x_test, y_test, xmax, ymax, p, p_custom)
        n_stocheck, nc_test, nr_test = calc_num_slices_to_check(sgnc, sgns, nc_test, nr_test, x_test, y_test, x_next, y_next, code, p, p_custom)  # ugh
        x_test, y_test = x_next, y_next

    ns_test = get_slicer_index(x_test, y_test, p, p_custom)[1]
    nc_new, nr_new = -1, -1
    is_same_section = False

    while is_section_in_bounds(nc_test, nr_test, umin, umax, p) and not is_same_section:
        # Number of slices to check in this section
        x_next, y_next, code = calc_next_coords(ray_in, sgnc, sgns, nc_test, nr_test, x_test, y_test, xmax, ymax, p, p_custom)
        n_stocheck, nc_new, nr_new = calc_num_slices_to_check(sgnc, sgns, nc_test, nr_test, x_test, y_test, x_next, y_next, code, p, p_custom)
        x_test, y_test = x_next, y_next

        # Iterate slices within this section
        for n in range(n_stocheck):

            # Check whether each slice is a solution to the transfer equation
            ray_out = check_slice_solution(tol, ray_in, ns_test, nc_test, p, p_custom)
            if not np.isnan(ray_out.t):
                # Found a solution within a slice.
                return ray_out
            
            # Before going to the next slice, check if there is a collision with
            # a wall. Skip if this is the last slice in the section.
            if trace_walls and not is_last_slice_in_section(ns_test, sgns, p): #(ns_test % p.n_each != 0 or ns_test % p.n_each != p.n_each - 1):
                ray_out = check_ywall_collision(ray_in, ns_test, nc_test, sgns, p, p_custom)
                if not np.isnan(ray_out.t):
                    # Found a wall collision.
                    return ray_out
            
            # Increment the slice index
            ns_test += sgns

        # If crossing columns, need to double check the last slice index in the
        # next section
        if nc_new != nc_test:
            ns_test -= sgns

        # Before continuing, check walls between sections.
        is_same_section = (nc_new == nc_test and nr_new == nr_test)

        if ( trace_walls and is_section_in_bounds(nc_new, nr_new, umin, umax, p) ):

            # Section changed columns, check collision with x-wall
            if nc_new != nc_test:
                ray_out = check_xwall_collision(ray_in, ns_test, nc_test, sgnc, p, p_custom)
                if not np.isnan(ray_out.t):
                    # Found a wall collision.
                    return ray_out
                
            # Must have changed rows instead. Check collision with y-wall
            elif nr_new != nr_test:
                ray_out = check_ywall_collision(ray_in, ns_test, nc_test, sgns, p, p_custom)
                if not np.isnan(ray_out.t):
                    # Found a wall collision.
                    return ray_out
                
            # Staying on the same section because this is the last one. Check
            # both x- and y-walls...
            else:
                # Need to subtract off an infinitesimal amount if xmax and ymax
                # correspond exactly to the boundary
                in_xgap, in_ygap = is_inside_slicer_gap(xmax-1e-15*sgnc, y_test-1e-15*sgns, p, p_custom)
                if in_xgap:
                    ray_out = check_xwall_collision(ray_in, ns_test, nc_test, sgnc, p, p_custom)
                    if not np.isnan(ray_out.t):
                        # Found a wall collision.
                        return ray_out
                elif in_ygap:
                    ray_out = check_ywall_collision(ray_in, ns_test, nc_test, sgns, p, p_custom)
                    if not np.isnan(ray_out.t):
                        # Found a wall collision.
                        return ray_out
                # If neither in_xgap or in_ygap are True, then a slice intersection should
                # have already been found so we don't need to worry about this case.
            
        nc_test = nc_new
        nr_test = nr_new

    # If none of that worked then this is a bizarre edge case, e.g., the direction
    # cosine is less than 1E-13 but it grazed off of a wall somehow. Sweep this ray
    # under the rug and say that it missed... ray_out will be filled with NaN values
    # at this point.
    return ray_out    # Phew!

def paraxial_ray_trace_slicer(ray_in, n1, n2, active_x, active_y, p, p_custom):

    nc, ns = get_paraxial_slice_index(p, active_x, active_y)
    pslice = get_slice_params(ns, nc, p, p_custom)
    transfer_dist_func, surf_normal_func, junk, transform_func = get_surface_funcs(pslice, p)

    # Transform the ray into the local coordinates of the slice
    junk1, junk2, surf_normal_func, transform_func = get_surface_funcs(pslice, p)
    ray_in_local = convert_ray_in_to_local(ray_in, pslice, transform_func)
    coords = np.array([ray_in.xt, ray_in.yt, 0])
    cosines = np.array([ray_in.l, ray_in.m, ray_in.n])
    coords_local = transform_func(coords, pslice, -1, True)
    cosines_local = transform_func(cosines, pslice, -1, False)

    power = (n2 - n1) * pslice.c
    if p.surface_type == 0:
        powerx = power
        powery = power
    elif p.surface_type == 1: # cylinder never has power along y-direction
        powerx = power
        powery = 0.0

    x, y, z = coords_local
    l, m, n = cosines_local

    if abs(n) < 1e-13:
        # Ray is traveling parallel to the surface. No refraction.
        return RayOut(np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)
    
    # Convert to slopes
    l /= n
    m /= n

    l = (n1 * l - x * powerx) / n2
    m = (n1 * m - y * powery) / n2

    # Convert back to direction cosines and normalize
    n = np.sqrt(1/(1 + l*l + m*m))
    l *= n
    m *= n

    cosines_local = np.array([l, m, n])

    # Transform back to global coordinates
    coords_back = transform_func(coords_local, pslice, 1, True)
    cosines_back = transform_func(cosines_local, pslice, 1, False)

    # Calculate surface normals normally. haha
    ln, mn, nn = slice_surface_normal(x, y, pslice, transfer_dist_func, surf_normal_func, transform_func)

    ray_out = RayOut(
        xs=coords_back[0],
        ys=coords_back[1],
        zs=coords_back[2],
        t=0,
        # above should be same as input
        ln=ln,
        mn=mn,
        nn=nn
    )

    return ray_out, cosines_back[0], cosines_back[1], cosines_back[2]