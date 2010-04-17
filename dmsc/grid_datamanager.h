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

#include <grid_dataset.h>

namespace grid
{
  struct octtree_piece
  {
    dataset_t   *dataset;
    mscomplex_t *msgraph;

    uint level;

    uint m_pieceno;

    octtree_piece (uint l);

    std::string label();
  };

  class data_manager_t
  {
    typedef std::vector<octtree_piece *> pieces_list_t;

  public:

    pieces_list_t                m_pieces;

    cellid_t                     m_size;
    std::string                  m_filename;
    double                       m_simp_tresh;
    bool                         m_use_ocl;
    cell_fn_t                   *m_pData;

    data_manager_t
        ( std::string  filename,
          cellid_t     size,
          bool         use_ocl,
          double       simp_tresh
          );

    void work();

    virtual ~data_manager_t ();

    void createDataPieces();

    void destoryDataPieces();

    void computeMsGraph ( octtree_piece  * );

    void collectManifold( octtree_piece  * );

    void readDataToMem();

    void logAllConnections(const std::string &prefix);

  };
}

#endif
