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
  typedef unsigned int                         critpt_idx_t;
  typedef std::vector<critpt_idx_t>            critpt_idx_list_t;
  typedef std::vector<cell_fn_t>               cp_fn_list_t;
  typedef std::pair<critpt_idx_t,critpt_idx_t> crit_idx_pair_t;
  typedef std::vector<crit_idx_pair_t>         crit_idx_pair_list_t;

  struct critpt_t
  {
    typedef std::multiset<critpt_idx_t>     conn_t;
    typedef std::vector<cellid_t>           disc_t;

    cellid_t     cellid;
    critpt_idx_t pair_idx;
    cell_fn_t    fn;

    bool isCancelled;
    bool isOnStrangulationPath;
    bool isBoundryCancelable;

    critpt_t()
    {
      isCancelled           = false;
      isOnStrangulationPath = false;
      isBoundryCancelable   = false;
      pair_idx              = -1;
    }

    disc_t asc_disc;
    disc_t des_disc;

    conn_t asc;
    conn_t des;
  };


  class mscomplex_t
  {
  public:

    typedef std::map<cellid_t,critpt_idx_t>  id_cp_map_t;
    typedef std::vector<critpt_t *>          critpt_list_t;

    critpt_list_t m_cps;
    id_cp_map_t   m_id_cp_map;

    rect_t        m_rect;
    rect_t        m_ext_rect;

    void add_critpt(cellid_t);

    void simplify(crit_idx_pair_list_t &,double simplification_treshold);

    void un_simplify(const crit_idx_pair_list_t &);

    void simplify_un_simplify(double simplification_treshold );


    void clear();

    static mscomplex_t * merge_up(const mscomplex_t& msc1,
                                  const mscomplex_t& msc2);

    void merge_down(mscomplex_t& msc1,
                    mscomplex_t& msc2);

    mscomplex_t(rect_t r,rect_t e):m_rect(r),m_ext_rect(e){}

    mscomplex_t(){}

    ~mscomplex_t();

    void write_discs(const std::string &fn_prefix);

    void print_connections(std::ostream & os);
  };

  typedef critpt_t::conn_t                 conn_t;
  typedef critpt_t::disc_t                 critpt_disc_t;
  typedef critpt_t::conn_t::iterator       conn_iter_t;
  typedef critpt_t::conn_t::const_iterator const_conn_iter_t;

}


#include <boost/serialization/array.hpp>
#include <boost/serialization/base_object.hpp>

namespace boost
{
  namespace serialization
  {
    template<class Archive>
    void serialize(Archive & ar, grid::mscomplex_t & g, const unsigned int );

  } // namespace serialization
}
#endif
