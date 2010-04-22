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

#include <iostream>
#include <fstream>
#include <timer.h>

#include <logutil.h>

#include <grid_dataset.h>
#include <grid_mscomplex.h>
#include <grid_datamanager.h>

using namespace std;

namespace grid
{

  void data_manager_t::createDataPieces ()
  {
    rect_t r(cellid_t::zero,(m_size-cellid_t::one)*2);
    rect_t e(cellid_t::zero,(m_size-cellid_t::one)*2);

    octtree_piece *dp = new octtree_piece(m_pieces.size());
    dp->dataset = new dataset_t(r,e);
    dp->msgraph = new mscomplex_t(r,e);
    m_pieces.push_back(dp);

    return;
  }

  void data_manager_t::destoryDataPieces()
  {
    for(uint i = 0 ; i<m_pieces.size();++i)
    {
      octtree_piece *dp = m_pieces[i];

      if (dp->msgraph != NULL)
        delete dp->msgraph;

      if(dp->dataset != NULL)
      {
        dp->dataset->clear_fnref();
        dp->dataset->clear();

        delete dp->dataset;
      }
    }

    m_pieces.clear();
  }

  void data_manager_t::computeMsGraph( octtree_piece *dp )
  {
    //  if(m_use_ocl != true)
    {
      dp->dataset->work();

      dp->dataset->writeout_connectivity(dp->msgraph);
    }
    //  else
    //  {
    //    dp->dataset->work_ocl();
    //
    //    dp->dataset->writeout_connectivity_ocl(dp->msgraph);
    //  }
  }

  void data_manager_t::collectManifold( octtree_piece  * dp)
  {
    dp->dataset->postMergeFillDiscs(dp->msgraph);

    dp->msgraph->write_discs("dp_disc_");
  }

  void data_manager_t::readDataToMem()
  {
    m_pData = new cell_fn_t[m_size[0]*m_size[1]*m_size[2]];

    ifstream ifs(m_filename.c_str(),std::ios::in|std::ios::binary);

    ifs.read((char*)(void*)m_pData,sizeof(cell_fn_t)*m_size[0]*m_size[1]*m_size[2]);
  }

  data_manager_t::data_manager_t
      ( std::string filename,
        cellid_t     size,
        bool         use_ocl,
        double       simp_tresh):
      m_filename(filename),
      m_size(size),
      m_use_ocl(use_ocl),
      m_simp_tresh(simp_tresh),
      m_pData(NULL)
  {
    //  if(m_use_ocl)
    //    GridDataset::init_opencl();

    createDataPieces();
  }

  data_manager_t::~data_manager_t ()
  {

    destoryDataPieces();

    //  if(m_use_ocl)
    //    GridDataset::stop_opencl();

    if(m_pData != NULL)
      delete m_pData;

    m_pData = NULL;
  }

  void data_manager_t::work()
  {

    Timer t;

    _LOG ( "==========================" );
    _LOG ( "Starting Processing       " );
    _LOG ( "--------------------------" );

    t.start();

    readDataToMem();

    _LOG("timer_time = "<<t.getElapsedTimeInMilliSec());

    m_pieces[0]->dataset->init();
    m_pieces[0]->dataset->init_fnref(m_pData);

    computeMsGraph(m_pieces[0]);

    _LOG("timer_time = "<<t.getElapsedTimeInMilliSec());

    if(m_simp_tresh > 0.0)
      m_pieces[0]->msgraph->simplify_un_simplify(m_simp_tresh);

    _LOG("timer_time = "<<t.getElapsedTimeInMilliSec());

    collectManifold( m_pieces[0]);

    _LOG("timer_time = "<<t.getElapsedTimeInMilliSec());

    _LOG ( "--------------------------" );
    _LOG ( "Finished Processing       " );
    _LOG ( "==========================" );

  }

  void data_manager_t::logAllConnections(const std::string &prefix)
  {

    for(uint i = 0 ; i <m_pieces.size();++i)
    {
      octtree_piece *dp = m_pieces[i];

      std::string filename(prefix+dp->label()+string(".txt"));

      ofstream outfile;
      outfile.open(filename.c_str(),  std::ios::out|std::ios::trunc);

      if(outfile.is_open() == false )
      {
        _LOG("failed to open log file");
        break;
      }

      std::stringstream ss;

      dp->msgraph->print_connections( (ostream&)ss);

      outfile<<ss.str();
    }

  }

  octtree_piece::octtree_piece (uint pno):
      dataset(NULL),
      msgraph(NULL),
      level(0),
      m_pieceno(pno)
  {
  }

  std::string octtree_piece::label()
  {
    std::stringstream ss;
    ss<<m_pieceno;

    return ss.str();
  }
}
