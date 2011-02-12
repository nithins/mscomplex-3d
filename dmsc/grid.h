#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

#include <vector>
#include <cpputils.h>
#include <aabb.h>

#include <boost/static_assert.hpp>
#define static_assert BOOST_STATIC_ASSERT

namespace grid
{
  const uint gc_grid_dim = 3;

  typedef int16_t                                         cell_coord_t;
  typedef float                                           cell_fn_t;
  typedef aabb::aabb_t<cell_coord_t,gc_grid_dim>          rect_t;
  typedef aabb::aabb_t<cell_coord_t,gc_grid_dim>::point_t cellid_t;
  typedef aabb::aabb_t<cell_coord_t,gc_grid_dim>::point_t rect_point_t;
  typedef aabb::aabb_t<cell_coord_t,gc_grid_dim>::point_t rect_size_t;
  typedef aabb::aabb_t<cell_coord_t,gc_grid_dim>::range_t rect_range_t;
  typedef std::vector<cellid_t>                           cellid_list_t;

  enum eGradientDirection
  {
    GRADDIR_DESCENDING,
    GRADDIR_ASCENDING,
    GRADDIR_COUNT,
  };

}

#endif
