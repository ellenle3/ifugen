import numpy as np
from slicer_generation import ImageSlicerParams
from dataclasses import dataclass

def load_slice_params(file_num, trace_walls, active_x, active_y):
    """Load parameters from a txt file. The two lines of the file *must* be
    as such:

    n_slices_per_col, n_cols
    dx, dy, gx_width, gx_depth, gy_width, gy_depth

    All entries in the first line must be integers and all entries in the second
    line must be floats. From then on, the entires must be as follows:

    alpha, beta, gamma, c, k, solution

    The number of entires *must* be equal to n_slices_per_col * n_cols. The row
    at slice_num * col_num + 2 corresponds to the entries for that slice and column
    index (plus 2 to account for the first two rows of data). For entries that are
    missing, a values are replaced by zeros.

    If there are repeat entries, they will be overwritten so that the last line
    determines the parameters for the slice at the given col_num, slice_num.
    
    Returns
    -------
    arr: nd.array of shape (n_slices_per_col * n_cols, 6)
        Array of parameters for each slice.
    p: ImageSlicerParams
        Contains data from the first two rows.
    """
    # Read in row by row
    data = np.loadtxt("custom_mirror_array_params" + file_num + ".txt")
    n_slices_per_col, n_cols = data[0]
    dx, dy, gx_width, gx_depth, gy_width, gy_depth = data[1]

    slice_data = data[2:]
    arr = np.zeros((n_slices_per_col * n_cols, 6))
    arr[:slice_data.shape[0]] = slice_data

    p = ImageSlicerParams(
        n_each = n_slices_per_col,
        n_rows = 1,
        n_cols = n_cols,
        mode = -1,
        trace_walls = False,
        active_x = False,
        active_y = False,
        dalpha = -1,
        dbeta = -1,
        dgamma = -1,
        gamma_offset = -1,
        alpha_cen = -1,
        beta_cen = -1,
        gamma_cen = -1,
        dx = dx,
        dy = dy,
        c = -1,
        k = -1,
        gx_width = gx_width,
        gx_depth = gx_depth,
        gy_width = gy_width,
        gy_depth = gy_depth,
        custom = True
    )

def make_sample_txt(p, fn_out):
    """Make a sample txt file for loading slice parameters.
    
    Parameters
    ----------
    p: ImageSliceParams
        Parameter file to use to create the sample txt file.
    fn_out: str
        Filename to write the sample txt file to.

    Returns
    -------
    arr:

    """
    pass