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

#include <boost/format.hpp>
#include <boost/lambda/lambda.hpp>

#include <logutil.h>
#include <timer.h>
#include <cpputils.h>

#include <grid_dataset.h>
#include <grid_mscomplex.h>
#include <grid_datamanager.h>

using namespace std;

namespace grid
{
  rect_range_t split_range(rect_range_t r_in, int i,int n)
  {
    ASSERT(is_in_range(i,0,n));

    int l = r_in[0] + divide_ri(r_in.span()*i,n);
    int u = r_in[0] + divide_ri(r_in.span()*(i+1),n);

    if( i+1 != n ) u = min<int>(u+1,r_in[1]);

    return rect_range_t(l,u);
  }

  rect_range_t expand_range(rect_range_t r_in, int i,int n)
  {
    ASSERT(is_in_range(i,0,n));

    if(i != 0 )  r_in[0] = r_in[0]-1;
    if(i+1 != n) r_in[1] = r_in[1]+1;

    return r_in;
  }

  void data_manager_t::createPieces ()
  {
    rect_t d(cellid_t::zero,(m_size-cellid_t::one)*2);

    for(int lev = 0 ;lev<m_max_levels+1;++lev)
    {
      int num_lev_pieces = two_power(lev);

      for(int i = 0 ;i<num_lev_pieces; ++i)
      {
        rect_t pr(cellid_t::zero,m_size);
        pr[2] = split_range(pr[2],i,num_lev_pieces);

        rect_t pe(pr);
        pe[2] = expand_range(pe[2],i,num_lev_pieces);

        rect_t r = rect_t(pr.lower_corner()*2,(pr.upper_corner()-cellid_t::one)*2);
        rect_t e = rect_t(pe.lower_corner()*2,(pe.upper_corner()-cellid_t::one)*2);

        piece_ptr_t p(new octtree_piece_t(r,e,d,m_max_levels-lev));
        m_pieces.push_back(p);
      }
    }

    return;
  }

  void data_manager_t::split_dataset()
  {
    std::ifstream ifs(m_filename.c_str(),std::ios::in|std::ios::binary);

    ensure(ifs.is_open(),"unable to open file");

    int pData_size = (divide_ri(m_size[2],m_num_pieces)+3)*m_size[0]*m_size[1];

    cell_fn_t * pData = new cell_fn_t[pData_size];

    for(int i = 0 ; i <m_num_pieces;++i)
    {
      rect_range_t z_rng(0,m_size[2]);
      z_rng = split_range(z_rng,i,m_num_pieces);
      z_rng = expand_range(z_rng,i,m_num_pieces);

      int beg = z_rng[0] * m_size[1] *m_size[0];
      int end = z_rng[1] * m_size[1] *m_size[0];

      ensure((end-beg) <= pData_size,"miscalculated max pData size");

      ifs.seekg(beg*sizeof(cell_fn_t),std::ios::beg);
      ensure(ifs.fail()==false,"failed to seek");

      ifs.read((char*)(void*)pData,sizeof(cell_fn_t)*(end-beg));
      ensure(ifs.fail()==false,"failed to read some data");

      std::string filename = m_filename+"."+to_string(i);

      std::ofstream ofs(filename.c_str(),std::ios::out|std::ios::binary);

      ofs.write((char*)(void*)pData,sizeof(cell_fn_t)*(end-beg));

      ofs.close();
    }

    delete []pData;

    ifs.close();
  }

  void data_manager_t::compute_subdomain_msgraphs ()
  {
    using namespace boost::lambda;

    int pData_size = (divide_ri(m_size[2],m_num_pieces)+3)*m_size[0]*m_size[1];

    cell_fn_t * pData = new cell_fn_t[pData_size];

    for(int i = 0 ; i <m_num_pieces;++i)
    {
      piece_ptr_t dp = m_pieces[m_num_pieces-1+i];
      rect_t e = dp->m_dataset->get_ext_rect();
      int num_pts = rect_t(e.lower_corner()/2,e.upper_corner()/2+cellid_t::one).volume();

      string filename = m_filename+"."+to_string(i);
      ifstream ifs(filename.c_str(),ios::in|ios::binary);
      ensure(ifs.is_open(),"unable to open file");

      ifs.read((char*)(void*)pData,sizeof(cell_fn_t)*num_pts);
      ensure(ifs.fail()==false,"failed to read some data");
      ifs.seekg(0,ios::end);
      ensure(ifs.tellg()==num_pts*sizeof(cell_fn_t),"file/piece size mismatch");

      dp->m_dataset->init(pData);
      dp->m_dataset->assignGradient();

      if( i != 0 )
      {
        rect_t bnd = e;
        bnd[2] = rect_range_t(e[2][0]+2,e[2][0]+2);
        dp->m_dataset->markBoundryCritical(bnd);
      }

      if( i+1 != m_num_pieces )
      {
        rect_t bnd = e;
        bnd[2] = rect_range_t(e[2][1]-2,e[2][1]-2);
        dp->m_dataset->markBoundryCritical(bnd);
      }

//      ofstream ofs((filename+".pairs").c_str(),ios::out);
//      dp->m_dataset->log_pairs(ofs);
//      ofs.close();

      dp->m_dataset->computeMsGraph(dp->m_msgraph.get());
      dp->m_dataset->clear();
    }

    delete []pData;
  }

  void data_manager_t::merge_subdomain_msgraphs ()
  {
    for(int lev = m_max_levels-1 ;lev >= 0 ;--lev)
    {
      int n = two_power(lev);

      for(int i = 0 ;i < n; ++i)
      {
        mscomplex_ptr_t msc  = m_pieces[n+i-1]->m_msgraph;
        mscomplex_ptr_t msc1 = m_pieces[(n+i)*2 -1]->m_msgraph;
        mscomplex_ptr_t msc2 = m_pieces[(n+i)*2]->m_msgraph;

        rect_t bnd = msc1->m_rect.intersection(msc2->m_rect);

        msc->merge_up(*msc1,*msc2,bnd);
      }
    }
  }

  void data_manager_t::save_results()
  {
    for(int i = 0 ;i < 2*m_num_pieces-1; ++i)
    {
      mscomplex_ptr_t msc  = m_pieces[i]->m_msgraph;

      msc->write_graph(string("msc_graph.txt.")+to_string(i));
    }
  }

  void data_manager_t::destoryPieces()
  {
    m_pieces.clear();
  }

  void data_manager_t::work()
  {
    split_dataset();

    createPieces();

    compute_subdomain_msgraphs();

    merge_subdomain_msgraphs();

    save_results();

    destoryPieces();
  }

  data_manager_t::data_manager_t
      ( std::string filename,
        cellid_t     size,
        int          max_levels,
        double       simp_tresh):
      m_filename(filename),
      m_size(size),
      m_max_levels(max_levels),
      m_simp_tresh(simp_tresh),
      m_num_pieces(two_power(max_levels))
  {
  }

  data_manager_t::~data_manager_t ()
  {
  }


  void compute_mscomplex_basic(std::string filename, cellid_t size, double simp_tresh)
  {
    cout<<"===================================="<<endl;
    cout<<"         Starting Processing        "<<endl;
    cout<<"------------------------------------"<<endl;

    Timer t;
    t.start();

    rect_t d(cellid_t::zero,(size-cellid_t::one)*2);
    boost::shared_ptr<dataset_t>   dataset(new dataset_t(d,d,d));
    boost::shared_ptr<mscomplex_t> msgraph(new mscomplex_t(d,d));

    int num_pts = size[0]*size[1]*size[2];
    std::vector<cell_fn_t> pt_data(num_pts);

    ifstream ifs(filename.c_str(),ios::in|ios::binary);
    ensure(ifs.is_open(),"unable to open file");

    ifs.read((char*)(void*)pt_data.data(),sizeof(cell_fn_t)*num_pts);
    ensure(ifs.fail()==false,"failed to read some data");

    ifs.seekg(0,ios::end);
    ensure(ifs.tellg()==num_pts*sizeof(cell_fn_t),"file/piece size mismatch");
    cout<<"data read ---------------- "<<t.getElapsedTimeInMilliSec()<<endl;

    dataset->init(pt_data.data());
    dataset->assignGradient();
    cout<<"gradient done ------------ "<<t.getElapsedTimeInMilliSec()<<endl;

    dataset->computeMsGraph(msgraph.get());
    cout<<"msgraph done ------------- "<<t.getElapsedTimeInMilliSec()<<endl;

    msgraph->simplify_un_simplify(simp_tresh);
    cout<<"simplification done ------ "<<t.getElapsedTimeInMilliSec()<<endl;

    dataset->collectManifolds(msgraph.get());
    cout<<"collect manifolds done --- "<<t.getElapsedTimeInMilliSec()<<endl;

    msgraph->write_manifolds("msc_manifolds.txt");
    msgraph->write_graph("msc_graph.txt");
    cout<<"results write done ------- "<<t.getElapsedTimeInMilliSec()<<endl;

    cout<<"------------------------------------"<<endl;
    cout<<"        Finished Processing         "<<endl;
    cout<<"===================================="<<endl;
  }

  octtree_piece_t::octtree_piece_t (rect_t r,rect_t e,rect_t d,int l):
      m_level(l),
//      m_rct(r),
//      m_ext(e),
      m_msgraph(new mscomplex_t(r,e))
  {
    if( l == 0 )
      m_dataset.reset(new dataset_t(r,e,d));
  }


}
