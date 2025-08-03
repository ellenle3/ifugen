#ifndef CUSTOM_SLICER_HELPERS_H
#define CUSTOM_SLICER_HELPERS_H

void LoadCustomParamsFromFile(double *custom_slice_params, int file_num, char params_dir[], int max_elements);

SLICE_PARAMS GetSliceParamsCustom(int slice_num, int col_num, double custom_slice_params[]);

#endif