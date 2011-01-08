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

#include <boost/shared_ptr.hpp>

#include <grid.h>

namespace grid
{

  class mscomplex_t;

  class dataset_t
  {

  public:

    // used as a bit mask.. cells can be critical and paired..in theory they all are
    enum eCellFlags
    {
      CELLFLAG_UNKNOWN  = 0,
      CELLFLAG_PAIRED   = 1,
      CELLFLAG_CRITICAL = 2,
      CELLFLAG_MASK     = 0x03,
    };

    enum eCellAdjDirection
    {
      CELLADJDIR_UNKNOWN   = (0),
      CELLADJDIR_LEFT      = (1),
      CELLADJDIR_RIGHT     = (2),
      CELLADJDIR_DOWN      = (3),
      CELLADJDIR_UP        = (4),
      CELLADJDIR_BACK      = (5),
      CELLADJDIR_FRONT     = (6),
    };


    typedef u_int8_t                                        cell_flag_t;
    typedef boost::multi_array<cell_flag_t,gc_grid_dim>     cellflag_array_t;
    typedef boost::multi_array_ref<cell_fn_t,gc_grid_dim>   varray_ref_t;
    typedef boost::shared_ptr<varray_ref_t>                 varray_ref_ptr_t;

  public:

    rect_t             m_rect;
    rect_t             m_ext_rect;

    varray_ref_ptr_t   m_vert_fns_ref;

    cellflag_array_t   m_cell_flags;
    cellflag_array_t   m_cell_pairs;
    cellflag_array_t   m_cell_mxfct;
    cellflag_array_t   m_cell_efdim_a;
    cellflag_array_t   m_cell_efdim_d;
    cellid_list_t      m_critical_cells;

  public:

    // initialization of the dataset

    dataset_t ( const rect_t &r,const rect_t &e );

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

    void  pairCellsWithinEst();

    void  assignMaxFacets();

    void  collateCriticalPoints();

    void  postMergeFillDiscs(mscomplex_t *msgraph);

    void  aggregateEffCellDim();


    // dataset interface
  public:

    cellid_t   getCellPairId ( cellid_t ) const;

    cellid_t   getCellMaxFacetId ( cellid_t ) const;

    cellid_t   getCellSecondMaxFacetId ( cellid_t ) const;

    inline bool   ptLt ( cellid_t c1,cellid_t c2) const
    {
      cell_fn_t f1 = (*m_vert_fns_ref)(c1/2);
      cell_fn_t f2 = (*m_vert_fns_ref)(c2/2);

      if (f1 != f2)
        return f1 < f2;

      return c1<c2;
    }

    bool   compareCells( cellid_t ,cellid_t ) const;

    uint   getCellPoints ( cellid_t ,cellid_t  * ) const;

    uint   getCellFacets ( cellid_t ,cellid_t * ) const;

    inline uint   getCellIncCells( cellid_t ,cellid_t * ) const;

    uint   getCellCofacets ( cellid_t ,cellid_t * ) const;

    uint   getCellCofaces ( cellid_t ,cellid_t * ) const;

    uint   getCellEst (cellid_t,cellid_t*) const;

    bool   isPairOrientationCorrect ( cellid_t c, cellid_t p ) const;

    bool   isCellMarked ( cellid_t c ) const;

    bool   isCellCritical ( cellid_t c ) const;

    bool   isCellPaired ( cellid_t c ) const;

    bool   areCellsIncident(cellid_t c1,cellid_t c2) const;

    void   pairCells ( cellid_t c,cellid_t p );

    void   unpairCells ( cellid_t c,cellid_t p );

    void   setCellMaxFacet (cellid_t c,cellid_t f);

    void   markCellCritical ( cellid_t c );

    inline uint getCellDim ( cellid_t c ) const;

    bool   isTrueBoundryCell ( cellid_t c ) const;

    bool   isFakeBoundryCell ( cellid_t c ) const;

    bool   isCellExterior ( cellid_t c ) const;

    // misc functions
  public:

    inline cell_fn_t get_cell_fn(cellid_t c)
    {
      cellid_t f[20];

      uint f_ct = getCellPoints(c,f);

      cell_fn_t ret = (*m_vert_fns_ref)(f[0]/2);

      for(uint i = 1 ;i < f_ct;++i)
      {
        ret = std::max(ret,(*m_vert_fns_ref)(f[i]/2));
      }

      return ret;
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

    void log_eff_dim();
  };

  inline uint dataset_t::getCellDim ( cellid_t c ) const
  {
    uint dim = 0;

    for (size_t i = 0 ; i < gc_grid_dim;++i)
      dim += ( c[i]&0x01 );

    return (dim);
  }

}

//namespace boost
//{
//  namespace serialization
//  {
//    template<class Archive>
//    void serialize(Archive & ar, grid::dataset_t & d, const unsigned int );
//
//  } // namespace serialization
//}
#endif
