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

  typedef std::vector<cell_fn_t>     cp_fn_list_t;
  typedef n_vector_t<int,2>          int_pair_t;
  typedef std::vector<int_pair_t>    int_pair_list_t;
  typedef std::map<cellid_t,uint>    id_cp_map_t;
  typedef std::vector<critpt_t >     critpt_list_t;

  typedef std::multiset<uint>                 conn_t;
  typedef std::multiset<uint>::iterator       conn_iter_t;
  typedef std::multiset<uint>::const_iterator const_conn_iter_t;
  typedef std::vector<conn_t>                 conn_list_t;

  struct critpt_t
  {
    cellid_t       m_cellid;
    cellid_t       m_vertid;
    char           m_index;
    cell_fn_t      m_fn;

    bool           m_is_cancelled;
    bool           m_is_paired;
    int            m_pair_idx;

    critpt_t(cellid_t c,uchar i,cell_fn_t f, cellid_t v);
    critpt_t(const critpt_t &);
  };


  class mscomplex_t
  {
  public:

    critpt_list_t m_cps;
    id_cp_map_t   m_id_cp_map;

    conn_list_t   m_conn[GDIR_CT];
    conn_list_t  &m_des_conn;
    conn_list_t  &m_asc_conn;

    rect_t        m_rect;
    rect_t        m_ext_rect;

  public:

    mscomplex_t(rect_t r,rect_t e);
    ~mscomplex_t();

    int  add_critpt(cellid_t c,uchar i,cell_fn_t f,cellid_t vert_cell);
    int  add_critpt(const critpt_t &);

    void connect_cps(cellid_t c1,cellid_t c2);
    void connect_cps(int p, int q);

    void dir_connect_cps(cellid_t c1,cellid_t c2);
    void dir_connect_cps(int p , int q);

    void pair_cps(cellid_t c1,cellid_t c2);
    void pair_cps(int p , int q);

    inline char& index(int i)
    {
      try{ASSERT(is_in_range(i,0,m_cps.size()));}
      catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

      return m_cps[i].m_index;
    }

    inline const char& index(int i) const
    {
      try{ASSERT(is_in_range(i,0,m_cps.size()));}
      catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

      return m_cps[i].m_index;
    }

    inline int& pair_idx(int i)
    {
      try{ASSERT(is_in_range(i,0,m_cps.size()));}
      catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

      return m_cps[i].m_pair_idx;
    }

    inline const int& pair_idx(int i) const
    {
      try{ASSERT(is_in_range(i,0,m_cps.size()));}
      catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

      return m_cps[i].m_pair_idx;
    }

    inline bool& is_paired(int i)
    {
      try{ASSERT(is_in_range(i,0,m_cps.size()));}
      catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

      return m_cps[i].m_is_paired;
    }

    inline const bool& is_paired(int i) const
    {
      try{ASSERT(is_in_range(i,0,m_cps.size()));}
      catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

      return m_cps[i].m_is_paired;
    }

    inline bool& is_canceled(int i)
    {
      try{ASSERT(is_in_range(i,0,m_cps.size()));}
      catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

      return m_cps[i].m_is_cancelled;
    }

    inline const bool& is_canceled(int i) const
    {
      try{ASSERT(is_in_range(i,0,m_cps.size()));}
      catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

      return m_cps[i].m_is_cancelled;
    }

    inline cellid_t& cellid(int i)
    {
      try{ASSERT(is_in_range(i,0,m_cps.size()));}
      catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

      return m_cps[i].m_cellid;
    }

    inline const cellid_t& cellid(int i) const
    {
      try{ASSERT(is_in_range(i,0,m_cps.size()));}
      catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

      return m_cps[i].m_cellid;
    }

    inline cellid_t& vertid(int i)
    {
      try{ASSERT(is_in_range(i,0,m_cps.size()));}
      catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

      return m_cps[i].m_vertid;
    }

    inline const cellid_t& vertid(int i) const
    {
      try{ASSERT(is_in_range(i,0,m_cps.size()));}
      catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

      return m_cps[i].m_vertid;
    }

    inline cell_fn_t& fn(int i)
    {
      try{ASSERT(is_in_range(i,0,m_cps.size()));}
      catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

      return m_cps[i].m_fn;
    }

    inline const cell_fn_t& fn(int i) const
    {
      try{ASSERT(is_in_range(i,0,m_cps.size()));}
      catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

      return m_cps[i].m_fn;
    }


  public:

    void simplify(int_pair_list_t &,double simplification_treshold);

    void un_simplify(const int_pair_list_t &);

    void simplify_un_simplify(double simplification_treshold );

    void invert_for_collection();

    void cancel_pair(int p, int q);
    void uncancel_pair( int p, int q);

    void clear();

    void merge_up  (const mscomplex_t& ,const mscomplex_t& ,const rect_t&);
    void merge_down(mscomplex_t& ,mscomplex_t& ,const rect_t&);

    void write_graph(std::ostream & os) const;
    void write_graph(const std::string & fn) const;


    inline std::string cp_info (int cp_no) const
    {
      std::stringstream ss;

      ss<<std::endl;
      ss<<"cp_no        ::"<<cp_no<<std::endl;
      ss<<"cellid       ::"<<cellid(cp_no)<<std::endl;
      ss<<"vert cell    ::"<<vertid(cp_no)<<std::endl;
      ss<<"index        ::"<<(int)index(cp_no)<<std::endl;
//      ss<<"fn           ::"<<m_cps[cp_no].fn<<std::endl;
      ss<<"is_cancelled ::"<<is_canceled(cp_no)<<std::endl;
      ss<<"is_paired    ::"<<is_paired(cp_no)<<std::endl;
      ss<<"pair_idx     ::"<<pair_idx(cp_no)<<std::endl;
      return ss.str();
    }

  };

  inline void order_pr_by_cp_index(mscomplex_t &msc,int &p,int &q)
  {
    if(msc.index(p) < msc.index(q))
      std::swap(p,q);

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
