#ifndef CUSTOM_SLICER_HELPERS_H
#define CUSTOM_SLICER_HELPERS_H

void LoadCustomParamsFromFile(double *custom_slice_params, int file_num, char params_dir[], int max_elements, int max_path_length);

void GetSliceParamsCustom(double* alpha, double* beta, double* gamma, double* cv, double* k,
int slice_num, int col_num, double custom_slice_params[]);

#endif