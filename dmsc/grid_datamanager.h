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

#ifndef GRID_DATAMANAGER_H_INCLUDED_
#define GRID_DATAMANAGER_H_INCLUDED_

#include <fstream>
#include <vector>

#include <grid.h>
#include <boost/shared_ptr.hpp>

namespace grid
{
  class dataset_t;
  class mscomplex_t;

  typedef boost::shared_ptr<dataset_t>   dataset_ptr_t;
  typedef boost::shared_ptr<mscomplex_t> mscomplex_ptr_t;

  struct octtree_piece_t
  {
    dataset_ptr_t   m_dataset;
    mscomplex_ptr_t m_msgraph;

//    rect_t                         m_rct,m_ext;
    int                            m_level;

    octtree_piece_t(rect_t r,rect_t e,rect_t d,int l);

  };

  class data_manager_t
  {

  public:
    typedef boost::shared_ptr<octtree_piece_t> piece_ptr_t;
    typedef std::vector<piece_ptr_t>           piece_ptr_list_t;

  public:

    piece_ptr_list_t             m_pieces;

    cellid_t                     m_size;
    std::string                  m_filename;
    double                       m_simp_tresh;
    int                          m_max_levels;
    int                          m_num_pieces;

    data_manager_t
        ( std::string  filename,
          cellid_t     size,
          int          max_levels,
          double       simp_tresh
          );

    virtual ~data_manager_t ();

    void work();

    void createPieces();

    void split_dataset();

    void compute_subdomain_msgraphs ();

    void merge_up_subdomain_msgraphs ();

    void merge_down_subdomain_msgraphs ();

    void save_results ();

    void destoryPieces();
  };

  void compute_mscomplex_basic(std::string filename, cellid_t size, double simp_tresh);

}

#endif
