#ifndef CUSTOM_SLICER_HELPERS_H
#define CUSTOM_SLICER_HELPERS_H
#include "surface_solns.h"

void LoadCustomParamsFromFile(double p_custom[], int file_num, char params_dir[], int max_elements);

SLICE_PARAMS GetSliceParamsCustom(int row_num, int col_num, double p_custom[]);

#endif