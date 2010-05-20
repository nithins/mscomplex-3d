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

  inline void ensure_cell_max_facet_known(const dataset_t *ds , cellid_t c)
  {
#ifdef USE_ENSURE_PREDICATES

    if((*ds->m_cell_mxfct)(c) == dataset_t::CELLADJDIR_UNKNOWN)
    {
      std::stringstream ss;

      ss<<"failed to ensure cell max fct is known\n";
      ss<<"cellid = "<<c;

      throw std::logic_error(ss.str());
    }
#endif
  }

  inline void ensure_cell_dim(const dataset_t *ds,cellid_t c,uint dim)
  {
#ifdef USE_ENSURE_PREDICATES
    if(ds->getCellDim(c) != dim)
      throw std::logic_error("falied to ensure cell dim");
#endif
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

  inline void ensure_never_reached(const std::string &m)
  {
#ifdef USE_ENSURE_PREDICATES
    throw std::logic_error(std::string("Code should not reach here\n")+m);
#endif

  }

  inline void ensure_cell_not_paired(dataset_t * ds,cellid_t c)
  {
#ifdef USE_ENSURE_PREDICATES
    if(ds->isCellPaired(c) == true)
      throw std::logic_error("failed to ensure that cell is not paired");
#endif
  }

  inline void ensure_cell_paired(const dataset_t *ds , cellid_t c)
  {
#ifdef USE_ENSURE_PREDICATES
    if(ds->isCellPaired(c) == false)
      throw std::logic_error("failed to ensure that cell is paired");
#endif
  }

  inline void ensure_valid_pair(const dataset_t *ds , cellid_t c,cellid_t p)
  {
#ifdef USE_ENSURE_PREDICATES
    if(ds->getCellPairId(c) != p ||
       ds->getCellPairId(p) != c)
      throw std::logic_error("failed to ensure valid pairing");
#endif
  }

  inline void ensure_cell_marked(dataset_t * ds,cellid_t c)
  {
#ifdef USE_ENSURE_PREDICATES
    if(ds->isCellMarked(c) == false)
      throw std::logic_error("failed to ensure that cell is marked");
#endif
  }
}

#endif
