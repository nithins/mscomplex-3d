/***************************************************************************
 *   Copyright (C) 2009 by nithin,,,   *
 *   nithin@gauss   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifndef __GRID_DATASET_H_INCLUDED_
#define __GRID_DATASET_H_INCLUDED_

#include <vector>

#include <boost/multi_array.hpp>

#include <grid.h>

namespace grid
{

  class mscomplex_t;

  class dataset_t
  {

  public:

    enum eCellFlags
    {
      CELLFLAG_UNKNOWN = 0,
      CELLFLAG_PAIRED  = 1,
      CELLFLAG_CRITCAL = 2,
    };

    typedef int8_t                                          cell_flag_t;
    typedef boost::multi_array<cell_fn_t,gc_grid_dim>       varray_t;
    typedef boost::multi_array<cellid_t,gc_grid_dim>        cellpair_array_t;
    typedef boost::multi_array<cell_flag_t,gc_grid_dim>     cellflag_array_t;
    typedef boost::multi_array_ref<cell_fn_t,gc_grid_dim>   varray_ref_t;

  public:

    class pt_comp_t
    {
      dataset_t *pOwn;
    public:
      pt_comp_t(dataset_t *o):pOwn(o){}

      bool operator()(cellid_t c1,cellid_t c2)
      {
        return pOwn->ptLt(c1,c2);
      }
    };


    rect_t             m_rect;
    rect_t             m_ext_rect;

    varray_ref_t      *m_vert_fns_ref;

    cellpair_array_t  *m_cell_pairs;
    cellflag_array_t  *m_cell_flags;
    cellid_list_t      m_critical_cells;

    pt_comp_t          m_ptcomp;

  public:

    // initialization of the dataset

    dataset_t ( const rect_t &r,const rect_t &e );

    dataset_t ();

    ~dataset_t ();

    void  init();

    void  clear();

    void  init_fnref(cell_fn_t * pData);

    void  clear_fnref();

    // actual algorithm work
  public:

    void  work();

    void  writeout_connectivity(mscomplex_t *msgraph);

    void  assignGradients();

    void  collateCriticalPoints();

    void  postMergeFillDiscs(mscomplex_t *msgraph);


    // dataset interface
  public:

    cellid_t   getCellPairId ( cellid_t ) const;

    inline bool   ptLt ( cellid_t c1,cellid_t c2) const
    {
      static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

      cell_fn_t f1 = (*m_vert_fns_ref)[c1[0]>>1][c1[1]>>1][c1[2]>>1];
      cell_fn_t f2 = (*m_vert_fns_ref)[c2[0]>>1][c2[1]>>1][c2[2]>>1];

      if (f1 != f2)
        return f1 < f2;

      return c1<c2;
    }

    bool   compareCells( cellid_t ,cellid_t ) const;

    uint   getCellPoints ( cellid_t ,cellid_t  * ) const;

    uint   getCellFacets ( cellid_t ,cellid_t * ) const;

    inline uint   getCellIncCells( cellid_t ,cellid_t * ) const;

    uint   getCellCofacets ( cellid_t ,cellid_t * ) const;

    bool   isPairOrientationCorrect ( cellid_t c, cellid_t p ) const;

    bool   isCellMarked ( cellid_t c ) const;

    bool   isCellCritical ( cellid_t c ) const;

    bool   isCellPaired ( cellid_t c ) const;

    void   pairCells ( cellid_t c,cellid_t p );

    void   markCellCritical ( cellid_t c );

    inline uint getCellDim ( cellid_t c ) const;

    bool   isTrueBoundryCell ( cellid_t c ) const;

    bool   isFakeBoundryCell ( cellid_t c ) const;

    bool   isCellExterior ( cellid_t c ) const;

    std::string  getCellFunctionDescription ( cellid_t pt ) const;

    std::string getCellDescription ( cellid_t cellid ) const;

    // misc functions
  public:
    inline static uint s_getCellDim ( cellid_t c )
    {
      uint dim = 0;

      for (size_t i = 0 ; i < gc_grid_dim;++i)
        dim += ( c[i]&0x01 );

      return (dim);
    }

    inline cell_fn_t get_cell_fn(cellid_t c)
    {
      cellid_t f[20];

      uint f_ct = getCellPoints(c,f);

      cell_fn_t ret = 0.0;

      for(uint i = 0 ;i < f_ct;++i)
      {
        ret += (*m_vert_fns_ref)(f[i]/2);
      }

      return ret/f_ct;
    }

    inline rect_t get_rect()
    {
      return m_rect;
    }

    inline rect_t get_ext_rect()
    {
      return m_ext_rect;
    }

    void log_flags();

    void log_pairs();
  };

  inline uint dataset_t::getCellDim ( cellid_t c ) const
  {
    return ( s_getCellDim(c));
  }

}

namespace boost
{
  namespace serialization
  {
    template<class Archive>
    void serialize(Archive & ar, grid::dataset_t & d, const unsigned int );

  } // namespace serialization
}
#endif
