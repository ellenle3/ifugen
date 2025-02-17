import numpy as np

def load_slice_params_file(file_num):
    """Load parameters from a txt file. The two lines of the file *must* be
    as such:

    n_slices_per_col, n_cols, surface_type
    dx, dy, gx_width, gx_depth, gy_width, gy_depth

    All entries in the first line must be integers and all entries in the second
    line must be floats. solution refers to what type of surface to use: 0 for
    2D conic, 1 for cylinder. ofFrom then on, the entires must
    be as follows:

    alpha, beta, gamma, c, k

    The number of entires *must* be equal to n_slices_per_col * n_cols. The row
    at slice_num * col_num + 2 corresponds to the entries for that slice and column
    index (plus 2 to account for the first two rows of data). For entries that are
    missing, a values are replaced by zeros.

    If there are repeat entries, they will be overwritten so that the last line
    determines the parameters for the slice at the given col_num, slice_num.
    
    Returns
    -------
    slice_param_arr: nd.array of shape (n_slices_per_col * n_cols, 6)
        Array of parameters for each slice.
    p: ImageSlicerParams
        Contains data from the first two rows.
    """
    # Read in row by row
    data = np.loadtxt("custom_mirror_array_params" + file_num + ".txt")
    n_slices_per_col, n_cols, surface_type = data[:3]

    slice_param_arr = np.zeros((n_slices_per_col * n_cols + 9, 5))
    slice_param_arr[:3] = n_slices_per_col, n_cols, surface_type

    slice_param_arr[3:9] = data[1]
    slice_data = data[2:]
    slice_param_arr[9:slice_data.shape[0] + 9] = slice_data
    
    return slice_param_arr

def get_slice_params_custom(col_num, slice_num, slice_param_arr):
    """Returns parameters that define a slice at a given col, slice index.
    """
    alpha = slice_param_arr[slice_num * col_num + 9, 0]
    beta = slice_param_arr[slice_num * col_num + 9, 1]
    gamma = slice_param_arr[slice_num * col_num + 9, 2]
    c = slice_param_arr[slice_num * col_num + 9, 3]
    k = slice_param_arr[slice_num * col_num + 9, 4]
    return alpha, beta, gamma, c, k

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