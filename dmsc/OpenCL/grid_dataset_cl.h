#ifndef GRID_DATASET_CL_H
#define GRID_DATASET_CL_H

#include <grid.h>

void assign_max_facet_opencl
(grid::rect_t ext_rect,grid::cell_fn_t *func, grid::cell_flag_t *flag);

void init_opencl();

#endif
