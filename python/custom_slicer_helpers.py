import numpy as np

def load_slice_params_file(file_num):
    """Load parameters from a txt file. The two lines of the file *must* be
    as such:

    n_slices_per_col n_cols surface_type
    dx dy gx_width gx_depth gy_width gy_depth

    Entires are space delimited. All entries in the first line must be integers
    and all remaining entries must be floats. solution refers to what type of
    surface to use: 0 for 2D conic, 1 for cylinder. From then on, the entires must
    be as follows:

    alpha beta gamma c k

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
    fname = "custom_mirror_array_params_" + str(file_num) + ".txt"
    file = open(fname, 'r')

    data = []
    for line in file.readlines():
        split_line = line.strip().split(' ')
        data.append([float(num) for num in split_line])

    # Set the first 9 entries to the parameters listed in the first two lines
    n_slices_per_col, n_cols, _ = data[0]
    n_slices_per_col = int(n_slices_per_col)
    n_cols = int(n_cols)
    slice_param_arr = np.zeros(n_slices_per_col * n_cols * 5 + 9)

    slice_param_arr[:3] = data[0]
    slice_param_arr[3:9] = data[1]

    # Set the remaining entries to the parameters for each slice. This isn't the
    # most efficient way to do this in Python, but it's more readable and is
    # somewhat closer to how it is implemented in C.

    # Iterate through rows in the file, skipping the first two
    for i in range(len(data) - 2):

        if i > n_slices_per_col * n_cols * 5 - 1:
            # Too many slices defined in the TXT file, truncate here
            break

        slice_param_arr[9 + 5*i: 9 + 5*(i+1)] = data[i + 2]
    
    return slice_param_arr

def get_slice_params_custom(slice_num, col_num, custom_slice_params):
    """Returns parameters that define a slice at a given col, slice index.
    """
    slice_num = int(slice_num)
    col_num = int(col_num)
    n_slices_per_col = int(custom_slice_params[0])
    start_idx = 9 + 5 * (col_num * n_slices_per_col + slice_num)

    alpha = custom_slice_params[start_idx]
    beta = custom_slice_params[start_idx + 1]
    gamma = custom_slice_params[start_idx + 2]
    c = custom_slice_params[start_idx + 3]
    k = custom_slice_params[start_idx + 4]

    return alpha, beta, gamma, c, k