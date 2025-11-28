import numpy as np
from surface_solns import SliceParams

NUM_PARAMS_PER_SLICE = 10

def load_slice_params_file(file_num):
    """Load parameters from a txt file. The three lines of the file *must* be
    as such:

    n_rows n_cols surface_type
    dx dy gx_width gx_depth gy_width gy_depth
    u0 u1 u2 ... uN

    (N is the number of rows.) Entries are space delimited. All entries in the
    first line must be integers and all remaining entries must be floats. 
    From then on, the entries must be as follows:

    alpha beta gamma c k zp syx syz sxy sxz

    The number of entries *must* be equal to n_rows * n_cols. The row
    at slice_num * col_num + 2 corresponds to the entries for that slice and column
    index (plus 2 to account for the first two rows of data). If entries are missing,
    all values following the missing entry will be truncated and replaced with
    zeroes. Excess parameters are also truncated.
    
    Returns
    -------
    p_custom: nd.array of shape (9 + n_rows + n_params * num_slices,)
        Array of parameters for each slice.
    """
    # Read in row by row
    fname = "p_custom_" + str(file_num) + ".txt"

    data = []
    with open(fname, 'r') as file:
        for line in file.readlines():
            split_line = line.strip().split()
            data.append([float(num) for num in split_line])

    # Set the first 9 entries to the parameters listed in the first two lines
    n_rows, n_cols, _ = data[0]
    n_rows = int(n_rows)
    n_cols = int(n_cols)
    num_slices = n_rows * n_cols

    # Next row of entries are the shifts for each row
    p_custom = np.zeros(NUM_PARAMS_PER_SLICE * num_slices + 9 + n_rows)

    p_custom[:3] = data[0]
    p_custom[3:9] = data[1]
    p_custom[9: 9 + n_rows] = data[2]

    # Set the remaining entries to the parameters for each slice. This isn't the
    # most efficient way to do this in Python, but it's more readable and is
    # somewhat closer to how it is implemented in C.

    # Iterate through rows in the file, skipping the first three
    for i in range(num_slices):

        # Stop if not enough rows in file or if fewer values than expected in a row
        if (i + 3 >= len(data) or len(data[i + 3]) < NUM_PARAMS_PER_SLICE):
            break

        base_idx = 9 + n_rows + NUM_PARAMS_PER_SLICE * i
        p_custom[base_idx: base_idx + NUM_PARAMS_PER_SLICE] = data[i + 3]

    return p_custom

def get_slice_params_custom(row_num, col_num, p_custom):
    """Returns parameters that define a slice at a given col, slice index.
    """
    row_num = int(row_num)
    col_num = int(col_num)
    
    n_rows = int(p_custom[0])
    n_cols = int(p_custom[1])

    # Safety check
    if not (0 <= row_num < n_rows and 0 <= col_num < n_cols):
        raise IndexError(f"Requested slice ({row_num}, {col_num}) out of bounds ({n_rows}, {n_cols})")

    u_start_idx = 9
    slice_params_start_idx = 9 + n_rows

    # Compute the 1D index of the slice in the flattened slice array
    row_idx = row_num + col_num * n_rows

    # Get the u value corresponding to this row
    u = p_custom[u_start_idx + row_num]

    # Get the 10 slice parameters
    base_idx = slice_params_start_idx + NUM_PARAMS_PER_SLICE * row_idx
    params = p_custom[base_idx : base_idx + NUM_PARAMS_PER_SLICE]

    # Return as a SliceParams object
    return SliceParams(*params, u)