#ifndef GRID_DATASET_CL_H
#define GRID_DATASET_CL_H

#include <grid.h>

namespace grid
{
  namespace opencl
  {

    void init();

    void assign_gradient(dataset_ptr_t ds, mscomplex_ptr_t msc);

//    void check_assign_gradient(dataset_ptr_t ds, int dim,rect_t check_rect,cell_flag_t mask);
  }
}


#endif
