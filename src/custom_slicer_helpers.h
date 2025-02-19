#ifndef CUSTOM_SLICER_HELPERS_H
#define CUSTOM_SLICER_HELPERS_H

double* LoadCustomParamsFromFile(int file_num, int *array_size);

void GetSliceParamsCustom(double* alpha, double* beta, double* gamma, double* cv, double* k,
int slice_num, int col_num, double custom_slice_params[]);

#endif