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

#ifndef __GRID_MSCOMPLEX_H_INCLUDED_
#define __GRID_MSCOMPLEX_H_INCLUDED_

#include <set>
#include <map>

#include <grid.h>

namespace grid
{
  struct critpt_t;

  typedef std::vector<uint>          critpt_idx_list_t;
  typedef std::vector<cell_fn_t>     cp_fn_list_t;
  typedef n_vector_t<uint,2>         uint_pair_t;
  typedef std::vector<uint_pair_t>   uint_pair_list_t;
  typedef std::multiset<uint>        conn_set_t;
  typedef std::vector<cellid_t>      disc_t;
  typedef std::vector<uint>          disc_contrib_t;
  typedef std::map<cellid_t,uint>    id_cp_map_t;
  typedef std::vector<critpt_t *>    critpt_list_t;

  struct critpt_t
  {
    cellid_t       cellid;
    cellid_t       vert_cell;
    char           index;
    cell_fn_t      fn;

    bool           is_cancelled;
    bool           is_paired;
    int            pair_idx;

    disc_contrib_t contrib[GRADDIR_COUNT];
    disc_t         disc[GRADDIR_COUNT] ;
    conn_set_t     conn[GRADDIR_COUNT];

    critpt_t(cellid_t c,uchar i,cell_fn_t f, cellid_t v);
    critpt_t(const critpt_t &);
  };


  class mscomplex_t
  {
  public:

    critpt_list_t m_cps;
    id_cp_map_t   m_id_cp_map;

    rect_t        m_rect;
    rect_t        m_ext_rect;

  public:

    mscomplex_t(rect_t r,rect_t e);
    ~mscomplex_t();

    int  add_critpt(cellid_t c,uchar i,cell_fn_t f,cellid_t vert_cell);
    int  add_critpt(const critpt_t &);

    void connect_cps(cellid_t c1,cellid_t c2);
    void connect_cps(uint_pair_t p);

    void dir_connect_cps(cellid_t c1,cellid_t c2);
    void dir_connect_cps(uint_pair_t p);

    void pair_cps(cellid_t c1,cellid_t c2);
    void pair_cps(uint_pair_t p);

  public:

    void simplify(uint_pair_list_t &,double simplification_treshold);

    void un_simplify(const uint_pair_list_t &);

    void simplify_un_simplify(double simplification_treshold );

    void add_disc_tracking_seed_cps();

    void clear();

    void merge_up  (const mscomplex_t& ,const mscomplex_t& ,const rect_t&);
    void merge_down(mscomplex_t& ,mscomplex_t& ,const rect_t&);

    void write_manifolds(std::ostream &os);
    void write_graph(std::ostream & os) const;

    void write_manifolds(const std::string &fn);
    void write_graph(const std::string & fn) const;
  };

  typedef conn_set_t::iterator       conn_iter_t;
  typedef conn_set_t::const_iterator const_conn_iter_t;

  inline void order_pr_by_cp_index(mscomplex_t *msc,uint_pair_t &e)
  {
    if(msc->m_cps[e[0]]->index < msc->m_cps[e[1]]->index)
      std::swap(e[0],e[1]);
  }
}


//#include <boost/serialization/array.hpp>
//#include <boost/serialization/base_object.hpp>
//
//namespace boost
//{
//  namespace serialization
//  {
//    template<class Archive>
//    void serialize(Archive & ar, grid::mscomplex_t & g, const unsigned int );
//
//  } // namespace serialization
//}
#endif
