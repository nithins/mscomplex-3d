#ifndef GRID_DATASET_CL_H
#define GRID_DATASET_CL_H

#include <grid.h>

void init_opencl();

void assign_gradient_opencl(grid::rect_t rct,grid::rect_t ext,grid::rect_t dom,
                            grid::cell_fn_t *func, grid::cell_flag_t *flag);

void check_assign_gradient_opencl
  (grid::dataset_ptr_t ds,int dim,grid::rect_t check_rect,grid::cell_flag_t mask = 0xff);


#endif
