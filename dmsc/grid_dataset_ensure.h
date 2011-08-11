#ifndef GRID_DATASET_ENSURE_H_INCLUDED
#define GRID_DATASET_ENSURE_H_INCLUDED

#include <grid_dataset.h>

// bunch of predicates that throw when I suspect something could be logically
// wrong with the state of the MS complex.. to be disabled in release builds

#ifndef NDEBUG
#define USE_ENSURE_PREDICATES
#endif

namespace grid
{

  inline cell_coord_t get_cell_dim(const cellid_t &c)
  {
    cell_coord_t d = 0;

    for(cell_coord_t i = cellid_t::static_size ; i != 0 ; )
      d += c[--i]%2;

    return d;
  }

  inline void ensure_cell_incidence(cellid_t c1,cellid_t c2)
  {
#ifdef USE_ENSURE_PREDICATES
    uint num_coords_diff = 0;

    for(uint i = 0 ;i < gc_grid_dim;++i)
    {
      num_coords_diff += (c1[i] != c2[i])?(1):(0);

      if(std::abs(c1[i]-c2[i]) != 1 &&
         std::abs(c1[i]-c2[i]) != 0)
        throw std::logic_error("failed to ensure cell incidences");
    }

    if(num_coords_diff != 1)
      throw std::logic_error("failed to ensure cell incidences");
#endif
  }

  inline void order_pair(cellid_t &c,cellid_t& p)
  {
    if(get_cell_dim(c) > get_cell_dim(p))
      std::swap(c,p);
  }

  inline void ensure_pairable(const dataset_t *ds, cellid_t c, cellid_t p)
  {
#ifdef USE_ENSURE_PREDICATES
    order_pair(c,p);

    ASSERT(ds->getCellDim(c)+1 == ds->getCellDim(p));

    ensure_cell_incidence(c,p);

    ASSERT(!ds->isCellPaired(c));

    ASSERT(!ds->isCellPaired(p));
#endif
  }


}

#endif
