import numpy as np
import math
from surface_solns import SliceParams, conic_2d_off_axis_angle

NUM_BASE_PARAMS = 10
NUM_PARAMS_PER_SLICE = 13

def load_slice_params_file(fname):
    """
    Python equivalent of LoadCustomParamsFromFile with linear-to-angular conversion
    """
    data = []
    with open(fname, "r") as f:
        for line in f:
            if line.strip():
                data.append([float(x) for x in line.split()])

    # ---- First line ----
    if len(data[0]) != 4:
        raise ValueError("First line must have 4 entries")

    n_rows, n_cols, surface_type, is_linear = map(int, data[0])
    num_slices = n_rows * n_cols

    # ---- Second line ----
    if len(data[1]) != 6:
        raise ValueError("Second line must have 6 entries")

    dx, dy, gx_width, gx_depth, gy_width, gy_depth = data[1]

    # ---- Allocate array ----
    array_size = NUM_BASE_PARAMS + n_rows + NUM_PARAMS_PER_SLICE * num_slices
    p_custom = np.zeros(array_size, dtype=float)

    # ---- Base params (must match C order) ----
    p_custom[0] = n_rows
    p_custom[1] = n_cols
    p_custom[2] = 1.0              # n_each (hard-coded in C)
    p_custom[3] = surface_type
    p_custom[4] = dx
    p_custom[5] = dy
    p_custom[6] = gx_width
    p_custom[7] = gx_depth
    p_custom[8] = gy_width
    p_custom[9] = gy_depth

    # ---- Third line: u values ----
    u_vals = data[2] if len(data) > 2 else []
    for i in range(n_rows):
        if i < len(u_vals):
            p_custom[NUM_BASE_PARAMS + i] = u_vals[i]
        else:
            p_custom[NUM_BASE_PARAMS + i] = 0.0

    # ---- Slice parameters ----
    slice_start = NUM_BASE_PARAMS + n_rows
    for slice_idx in range(num_slices):
        file_row = 3 + slice_idx
        if file_row >= len(data):
            break
        if len(data[file_row]) < NUM_PARAMS_PER_SLICE:
            break

        alpha, beta, gamma, cv, k, zp, syx, syz, sxy, sxz, theta, szx, szy = data[file_row][:NUM_PARAMS_PER_SLICE]

        # Linear-to-angular conversion if is_linear is true (f = 1)
        if is_linear:
            # alpha and beta are already OADs; convert via Conic2DOffAxisAngle
            alpha, beta = conic_2d_off_axis_angle(beta, alpha, cv, k)
            alpha = math.degrees(alpha)
            beta  = math.degrees(beta)
            gamma = math.degrees(math.atan(gamma))   # gamma = d / f, f = 1
            theta = math.degrees(math.atan(theta))

        base_idx = slice_start + NUM_PARAMS_PER_SLICE * slice_idx
        p_custom[base_idx:base_idx + NUM_PARAMS_PER_SLICE] = [
            alpha, beta, gamma, cv, k, zp,
            syx, syz, sxy, sxz, theta, szx, szy
        ]

    return p_custom


def get_slice_params_custom(slice_num, col_num, p_custom):
    """
    Python equivalent of GetSliceParams
    """
    slice_num = int(slice_num)
    col_num = int(col_num)

    n_rows = int(p_custom[0])
    n_cols = int(p_custom[1])
    n_each = int(p_custom[2])

    # Safety check
    if (slice_num < 0 or
        slice_num >= n_rows * n_each or
        col_num < 0 or
        col_num >= n_cols):
        raise IndexError(
            f"Requested slice ({slice_num}, {col_num}) out of bounds "
            f"({n_rows * n_each}, {n_cols})"
        )

    u_start_idx = NUM_BASE_PARAMS
    slice_params_start_idx = NUM_BASE_PARAMS + n_rows

    # Global slice index
    slice_idx = slice_num + col_num * n_rows * n_each
    row_num = int(math.floor(slice_num / n_each))

    # u value
    u = p_custom[u_start_idx + row_num]

    # Slice parameters (13)
    base_idx = slice_params_start_idx + NUM_PARAMS_PER_SLICE * slice_idx
    params = p_custom[base_idx:base_idx + NUM_PARAMS_PER_SLICE]

    return SliceParams(*params, u)
