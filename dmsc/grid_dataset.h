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
#include <boost/function.hpp>

#include <grid.h>

namespace grid
{

  class mscomplex_t;

  typedef u_int8_t cell_flag_t;

  class dataset_t:public boost::enable_shared_from_this<dataset_t>
  {

  public:

    // used as a bit mask.. cells can be critical and paired..in theory they all are
    enum eCellFlags
    {
      CELLFLAG_VISITED   = 0x80,
      CELLFLAG_CRITICAL = 0x40,
      CELLFLAG_MASK     = 0xc0,
    };

    // bits [0,3) max facet of a cell
    // bits [3,6) pair of a cell
    // bit 6 ..  mark bit used by bfs to say visted or not
    // bit 7 .. is cell critical or not.


    typedef boost::multi_array<cellid_t,gc_grid_dim>        cellid_array_t;
    typedef boost::multi_array<cell_flag_t,gc_grid_dim>     cellflag_array_t;
    typedef boost::multi_array<cell_fn_t,gc_grid_dim>       varray_t;


  public:

    rect_t             m_rect;
    rect_t             m_ext_rect;
    rect_t             m_domain_rect;

    varray_t           m_vert_fns;

    cellflag_array_t   m_cell_flags;
    cellid_list_t      m_critical_cells[g_num_threads];

    boost::function<bool (cellid_t,cellid_t)> cmp_ftor;
    boost::function<bool (cellid_t,cellid_t)> cmp_ftors[2];


  public:

    // initialization of the dataset
    dataset_t ( const rect_t &r,const rect_t &e,const rect_t &d );
    ~dataset_t ();

    void  init(const std::string &filename);
    void  clear();

  // the actual work routines
  public:
    void  assignGradient();

    void  markBoundryCritical(const rect_t &b);

    void  computeMsGraph(mscomplex_ptr_t msgraph);

    void  saveManifolds(mscomplex_const_ptr_t msgraph,const std::string &);

  // subroutines to main functions
  public:

    void  assignMaxFacets_thd(int tid,int dim);

    void  pairCellsWithinEst_thd(int tid);


    void  markBoundryCritical_thd(const rect_t &b,int tid);

    void  saddle_visit(mscomplex_ptr_t msgraph);

    void  extrema_connect_thd(mscomplex_ptr_t msgraph,cp_producer_ptr_t p);

    void  saddle_connect_thd(mscomplex_ptr_t msgraph,cp_producer_ptr_t p);

    template<typename mfold_t>
    void  get_mfold(mfold_t * mfold,mscomplex_const_ptr_t msc,int i,int dir) const;

  // dataset interface
  public:

    cellid_t   getCellPairId ( cellid_t ) const;

    cellid_t   getCellMaxFacetId ( cellid_t ) const;

    cellid_t   getCellSecondMaxFacetId ( cellid_t ) const;

    inline bool   ptLt ( cellid_t c1,cellid_t c2) const
    {
      cell_fn_t f1 = m_vert_fns(c1/2);
      cell_fn_t f2 = m_vert_fns(c2/2);

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

    bool   isCellCritical ( cellid_t c ) const;

    bool   isCellPaired ( cellid_t c ) const;

    bool   isCellVisited ( cellid_t c ) const;

    bool   areCellsIncident(cellid_t c1,cellid_t c2) const;

    void   pairCells ( cellid_t c,cellid_t p );

    void   visitCell( cellid_t c);

    void   setCellMaxFacet (cellid_t c,cellid_t f);

    void   markCellCritical ( cellid_t c );

    inline int getCellDim ( cellid_t c ) const;

    bool   isTrueBoundryCell ( cellid_t c ) const;

    bool   isFakeBoundryCell ( cellid_t c ) const;

    bool   isCellExterior ( cellid_t c ) const;

    // misc functions
  public:

    inline cell_fn_t get_cell_fn(cellid_t c)
    {
      cellid_t f[20];

      uint f_ct = getCellPoints(c,f);

      cell_fn_t ret = m_vert_fns(f[0]/2);

      for(uint i = 1 ;i < f_ct;++i)
      {
        ret = std::max(ret,m_vert_fns(f[i]/2));
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

    void log_pairs(std::ostream &os = std::cout);
    void log_pairs(const std::string &s);

    void log_visits(std::ostream &os = std::cout);
    void log_visits(const std::string &s);

    void log_pair_visits(std::ostream &os = std::cout);
    void log_pair_visits(const std::string &s);


    void log_max_facets();

    void extract_vdata_subarray(rect_t r,const std::string &filename);
  };

  inline int dataset_t::getCellDim ( cellid_t c ) const
  {
    int dim = 0;

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
