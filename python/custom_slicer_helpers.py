import numpy as np
from surface_solns import SliceParams

def load_slice_params_file(file_num):
    """Load parameters from a txt file. The three lines of the file *must* be
    as such:

    n_rows n_cols surface_type
    dx dy gx_width gx_depth gy_width gy_depth
    u0 u1 u2 ... uN

    Entires are space delimited. All entries in the first line must be integers
    and all remaining entries must be floats. solution refers to what type of
    surface to use: 0 for 2D conic, 1 for cylinder. From then on, the entires must
    be as follows:

    alpha beta gamma c k zp syx syz sxy sxz

    The number of entires *must* be equal to n_rows * n_cols. The row
    at slice_num * col_num + 2 corresponds to the entries for that slice and column
    index (plus 2 to account for the first two rows of data). For entries that are
    missing, all values are replaced by zeros. Excess parameters are truncated.
    
    Returns
    -------
    slice_param_arr: nd.array of shape (n_rows * n_cols, 6)
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
    n_rows, n_cols, _ = data[0]
    n_rows = int(n_rows)
    n_cols = int(n_cols)
    num_slices = n_rows * n_cols

    # Next row of entries are the shifts for each row
    custom_slice_params = np.zeros(5 * num_slices + 9 + n_rows)

    custom_slice_params[:3] = data[0]
    custom_slice_params[3:9] = data[1]
    custom_slice_params[9: 9 + n_rows] = data[2]

    # Set the remaining entries to the parameters for each slice. This isn't the
    # most efficient way to do this in Python, but it's more readable and is
    # somewhat closer to how it is implemented in C.

    # Iterate through rows in the file, skipping the first two
    for i in range(num_slices):

        # Stop if not enough rows in file or if fewer values than expected in a row
        if (i > len(data) - 2 or len(data[i + 2]) < 5):
            break

        base_idx = 9 + 5 * i
        custom_slice_params[base_idx: base_idx + 5] = data[i + 2]
    
    return custom_slice_params

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

    zp = 0
    syx = 0
    syz = 0
    sxy = 0
    sxz = 0
    u = 0
    return SliceParams(alpha, beta, gamma, c, k, zp, syx, syz, sxy, sxz, u)