import math
import numpy as np
from dataclasses import dataclass
from surface_solns import *
from custom_array_helpers import *

@dataclass
class ImageSlicerParams:
    """Class for storing image slicer params."""
    custom: int
    cylinder: int
    n_each: int
    n_rows: int
    n_cols: int
    angle_mode: int
    dalpha: float
    dbeta: float
    dgamma: float
    gamma_offset: float
    alpha_cen: float
    beta_cen: float
    gamma_cen: float
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
    xt: float
    yt: float
    l: float
    m: float
    n: float

@dataclass
class RayOut:
    xs: float
    ys: float
    zs: float
    t: float
    ln: float
    mn: float
    nn: float

def check_slicer_params(p):
    """Returns True if all image slicer parameters are valid.
    """
    pass

def get_surface_funcs(p):
    if p.c == 0:
        return tilted_plane_sag,

def get_slicer_size(p):
    # Return image slicer x, y dimensions
    n_slices = p.n_each*p.n_rows
    ysize = (n_slices*p.dy + (n_slices-1)*p.gy_width)
    xsize = (p.n_cols*p.dx + (p.n_cols-1)*p.gx_width)
    return xsize, ysize

def get_slicer_index(x, y, p):
    xsize, ysize = get_slicer_size(p)
    # Gets column and slice indices for a given x, y
    col_num = (x + xsize/2) // (p.dx + p.gx_width)
    slice_num = (y + ysize/2) // (p.dy + p.gy_width)
    return col_num, slice_num

def is_inside_slicer_gap(x, y, p):
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

def get_slice_params(slice_num, col_num, p):
    """Returns the angles alpha and beta for the slice number. Indexing starts at
    0 from the bottom (negative y direction) of the image slicer.
    """
    if p.custom:
        return get_slice_params_custom(slice_num, col_num, p)
    return get_slice_params_standard(slice_num, col_num, p) 

def get_slice_params_standard(slice_num, col_num, p):
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
    return alpha, beta, gamma, p.c, p.k

def make_image_slicer(x, y, p, sag_func):
    """Returns the sag of the image slicer.

    Parameters
    ---------
    p: ImageSlicerParams
    
    Returns
    -------
    z: nd_array
        Sag
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
    alpha, beta, gamma, c, k = get_slice_params(slice_num, col_num, p)
    return sag_func(x, y, c, k, alpha, beta, gamma)

def find_bounded_extremum(x0, y0, p, mode, sag_func, critical_xy_func):
    """
    Parameters
    ----------
    x0: float
        Initial guess for x0.
    mode: bool
        True if max, False if min.
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
    alpha, beta, gamma, c, k = get_slice_params(slice_num, col_num, p)

    xlo = col_num * (p.dx + p.gx_width) - xsize/2
    xhi = xlo + p.dx
    ylo = slice_num * (p.dy + p.gy_width) - ysize/2
    yhi = ylo + p.dy

    # Compute critical point
    xc, yc = critical_xy_func(c, k, alpha, beta, gamma)

    # There are up to 5 points to compare depending on whether the critical
    # point is within bounds
    zsolns = np.zeros(5)
    zsolns[0] = sag_func(xlo, ylo, c, k, alpha, beta, gamma)
    zsolns[1] = sag_func(xlo, yhi, c, k, alpha, beta, gamma)
    zsolns[2] = sag_func(xhi, ylo, c, k, alpha, beta, gamma)
    zsolns[3] = sag_func(xhi, yhi, c, k, alpha, beta, gamma)

    # Number of elements to compare so far
    n_compare = 4

    # Check whether critical point is in bounds. If yes, compute the sag and
    # add to the array of points to compare.
    if (yc >= ylo and yc <= yhi and xc >= xlo and xc <= xhi):
        zsolns[4] = sag_func(xc, yc, c, k, alpha, beta, gamma)
        n_compare += 1

    # Compare potential solutions to get the maximum or minimum
    if mode:
        return np.max(zsolns[:n_compare])
    return np.min(zsolns[:n_compare])

def find_global_extrema_slicer(p, sag_func, critical_xy_func):
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
            z = make_image_slicer(x, y, p, sag_func)
            # Update the (x,y) corresponding to the max and min values
            if z > z0_max: x0_max, y0_max, z0_max = x, y, z
            elif z < z0_min: x0_min, y0_min, z0_min = x, y, z

    zmin = find_bounded_extremum(x0_min, y0_min, p, False, sag_func, critical_xy_func)
    zmax = find_bounded_extremum(x0_max, y0_max, p, True, sag_func, critical_xy_func)
    
    return zmin, zmax

def transfer_equation(t, xt, yt, l, m, n, p, sag_func):
    """The roots of this function give the value of t.
    """
    xs = xt + t*l
    ys = yt + t*m
    zs = t*n
    sag = make_image_slicer(xs, ys, p, sag_func)
    return sag - zs

def get_ray_bounds(ray_in, zmin, zmax, p):

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

    return nc_min, ns_min, nc_max, ns_max

def is_ray_in_bounds(nc_min, ns_min, nc_max, ns_max, p):
    """Returns True if the ray is in bounds at least some of the time.
    """
    if (nc_min < 0 and nc_max < 0) or (nc_min >= p.n_cols and nc_max >= p.n_cols):
        # x-value of the ray is too high or low
        return False
    
    n_sperc = p.n_each * p.n_rows  # number of slices per column
    if (ns_min < 0 and ns_max < 0) or (ns_min >= n_sperc and ns_max >= n_sperc):
        # y-value of the ray is too high or low
        return False
    
    return True

def ray_trace_slicer(ray_in, zmin, zmax, p, trace_walls, sag_func, transfer_dist_func, surf_normal_func):
    """Computes ray trace.
    
    Returns
    -------
    """
    # Tolerance for accepting the transfer distance as valid
    tol = 1e-11
    ray_out = RayOut(np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)

    # Get starting and ending rows and columns
    nc_min, ns_min, nc_max, ns_max = get_ray_bounds(ray_in, zmin, zmax, p)
    if not is_ray_in_bounds(nc_min, ns_min, nc_max, ns_max, p):
        return ray_out
    print("bounds: " + str((nc_min, ns_min, nc_max, ns_max)))
    
    # Ray is in bounds at least some of the time. Start from the min col and slice
    # indices and check solutions until we hit the max
    nc_test, ns_test = nc_min, ns_min
    x_test, y_test = ray_in.xt, ray_in.yt

    dcol = abs(nc_max - nc_test)        # Number of col indices left to iterate
    dslice = abs(ns_max - ns_test)      # Number of slice (row) indices left to iterate

    # Keep track of signs - which way to iterate
    sgnc = 0
    sgns = 0
    if dcol > 0:
        sgnc = (nc_max - nc_min) / dcol     
    if dslice > 0:
        sgns = (ns_max - ns_min) / dslice 
    
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
                alpha, beta, gamma, c, k = get_slice_params(ns_test, nc_test, p)
                t = transfer_dist_func(xt, yt, l, m, n, c, k, alpha, beta, gamma)
                result = transfer_equation(t, xt, yt, l, m, n, p, sag_func)

                print("checking " + str((nc_test, ns_test)))
                print("t is " + str(t))
                print("result is " + str(result))
                
                # Check whether the transfer distance of the current slice is a valid
                # zero of the transfer equation.
                if abs(result) < tol:
                    # Yes - found a solution!
                    xs = xt + t*l
                    ys = yt + t*m
                    zs = t*n
                    ln, mn, nn = surf_normal_func(xs, ys, c, k, alpha, beta, gamma, True)

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
                        xnear = (nc_test + 1) * p.dx - xsize / 2
                        tnear = (xnear - xt) / l
                        ynear = yt + tnear*m
                        znear = tnear*n
                        
                        xfar = (nc_test + 1) * p.dx + sgnc * p.gx_width - xsize / 2
                        tfar = (xfar - xt) / l
                        yfar = yt + tfar*m
                        zfar = tfar*n

                        znear_slice = sag_func(xnear, ynear, c, k, alpha, beta, gamma)
                        alpha, beta, gamma, c, k = get_slice_params(ns_test, nc_test + 1*sgnc, p)
                        zfar_slice = sag_func(xfar, yfar, c, k, alpha, beta, gamma)

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
                        ynear = (ns_test + 1) * p.dy - ysize / 2
                        tnear = (ynear - yt) / m
                        xnear = xt + tnear*l
                        znear = tnear*n
                        
                        yfar = (ns_test + 1) * p.dy + sgns * p.gy_width - ysize / 2
                        tfar = (yfar - yt) / m
                        xfar = xt + tfar*l
                        zfar = tfar*n
                        
                        znear_slice = sag_func(xnear, ynear, c, k, alpha, beta, gamma)
                        alpha, beta, gamma, c, k = get_slice_params(ns_test + 1*sgns, nc_test, p)
                        zfar_slice = sag_func(xfar, yfar, c, k, alpha, beta, gamma)

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
    if abs(transfer_equation(t, xt, yt, l, m, n, p, sag_func)) < tol:
        ray_out.xs, ray_out.ys, ray_out.zs, ray_out.t = xt + t*l, yt + t*m, p.gx_depth, t
        ray_out.ln, ray_out.mn, ray_out.nn = 0, 0, -1
        return ray_out

    t = p.gy_depth / n
    if abs(transfer_equation(t, xt, yt, l, m, n, p, sag_func)) < tol:
        ray_out.xs, ray_out.ys, ray_out.zs, ray_out.t = xt + t*l, yt + t*m, p.gy_depth, t
        ray_out.ln, ray_out.mn, ray_out.nn = 0, 0, -1
        return ray_out

    # If none of that worked then this is a bizarre edge case, e.g., the direction cosine is less than 1E-13
    # but it grazed off of a wall somehow. Sweep this ray under the rug and say that it missed...
    return ray_out    # Phew!