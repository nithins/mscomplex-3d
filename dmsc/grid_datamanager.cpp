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
#include <boost/multi_array.hpp>

#include <timer.h>
#include <cpputils.h>

#include <grid_dataset.h>
#include <grid_mscomplex.h>
#include <grid_datamanager.h>

using namespace std;

namespace grid
{
  void split_rect(rect_t r, rect_t &r0,rect_t &r1 ,int dir)
  {
    r0 = r;
    r1 = r;

    r0[dir][1] = r[dir][0] + divide_ri(r[dir].span(),2)+1;
    r1[dir][0] = r[dir][0] + divide_ri(r[dir].span(),2);
  }

  void data_manager_t::createPieces ()
  {
    int level = 0;

    rect_t d(cellid_t::zero,m_size);

    piece_ptr_t p(new octtree_piece_t(d,d,m_level_ct-level));
    m_pieces.push_back(p);

    for ( int s_dir = gc_grid_dim-1 ; s_dir >= 0; s_dir--)
    {
      for( int s_level = 0 ; s_level < m_levels[s_dir]; ++s_level)
      {
        int s = two_power(level)-1;
        int e = two_power(level+1)-1;

        ++level;

        for(int i = s ; i < e ; ++i)
        {
          piece_ptr_t dp = m_pieces[i];

          rect_t pr1,pr2;

          split_rect(dp->m_prct,pr1,pr2,s_dir);

          piece_ptr_t p1(new octtree_piece_t(pr1,d,m_level_ct-level));
          m_pieces.push_back(p1);

          piece_ptr_t p2(new octtree_piece_t(pr2,d,m_level_ct-level));
          m_pieces.push_back(p2);

        }
      }
    }

    return;
  }

  void data_manager_t::split_dataset()
  {
    typedef boost::multi_array_ref<cell_fn_t, 2>   zslice_ref_t;
    typedef boost::multi_array<cell_fn_t, 2>       zslice_t;
    typedef boost::multi_array_types::index_range  range_t;
    typedef boost::shared_ptr<ofstream>            ofstream_ptr_t;



    std::ifstream ifs(m_filename.c_str(),std::ios::in|std::ios::binary);
    ensure(ifs.is_open(),"unable to open file");

    int xy_size = m_size[0]*m_size[1];
    int xy_ct   = two_power(m_levels[0])*two_power(m_levels[1]);

    cell_fn_ptr_t pData(new cell_fn_t[xy_size]);

    zslice_ref_t zslice(pData.get(),boost::extents[m_size[0]][m_size[1]],
                    boost::fortran_storage_order());

    int slc_i = two_power(m_levels[2])-1;
    int slc_e = two_power(m_levels[2]+1)-1;

    int bnd_z = -1;

    if(slc_i +1 != slc_e)
      bnd_z = m_pieces[slc_i+1]->m_prct[2][0];

    vector< ofstream_ptr_t> xy_files;

    int pc_beg = two_power(m_level_ct)-1;
    int pc_end = pc_beg + xy_ct;

    for( int pc_i = pc_beg ; pc_i < pc_end;++pc_i)
    {
      string filename = m_pieces[pc_i]->get_basename(m_basename)+".raw";
      ofstream_ptr_t p(new ofstream(filename.c_str(),ios::binary));
      ensure(p->is_open(),"unable to open piece file");
      xy_files.push_back(p);
    }

    for( int z = 0; z < m_size[2]; ++z)
    {
      if (bnd_z >=0 && z == bnd_z - 1)
      {
        pc_end += xy_ct;

        for( int pc_i = pc_beg+xy_ct ; pc_i < pc_end;++pc_i)
        {
          string filename = m_pieces[pc_i]->get_basename(m_basename)+".raw";
          ofstream_ptr_t p(new ofstream(filename.c_str(),ios::binary));
          ensure(p->is_open(),"unable to open piece file");
          xy_files.push_back(p);
        }
      }

      ifs.read((char*)(void*)pData.get(),xy_size*sizeof(cell_fn_t));

      for(int pc_i = pc_beg; pc_i != pc_end; ++pc_i)
      {
        rect_range_t x_rng = m_pieces[pc_i]->m_ext_prct[0];
        rect_range_t y_rng = m_pieces[pc_i]->m_ext_prct[1];

        zslice_t pc_slice(zslice[boost::indices[range_t(x_rng[0],x_rng[1])]
                                 [range_t(y_rng[0],y_rng[1])]],boost::fortran_storage_order());

        int sz = x_rng.span()*y_rng.span()*sizeof(cell_fn_t);

        xy_files[pc_i-pc_beg]->write((const char *)(const void*)pc_slice.data(),sz);
      }

      if(bnd_z >=0 && z  == bnd_z + 1)
      {
        pc_beg += xy_ct;

        xy_files.erase(xy_files.begin(),xy_files.begin()+xy_ct);

        if(slc_i +1 != slc_e)
          ++slc_i;

        if(slc_i +1 != slc_e)
          bnd_z = m_pieces[slc_i+1]->m_prct[2][0];
        else
          bnd_z = -1;
      }
    }

    xy_files.erase(xy_files.begin(),xy_files.begin()+xy_ct);

    ifs.close();
  }

  void data_manager_t::compute_subdomain_msgraphs ()
  {
    using namespace boost::lambda;

    int pData_size = (divide_ri(m_size[2],two_power(m_levels[2]))+3)*
                     (divide_ri(m_size[1],two_power(m_levels[1]))+3)*
                     (divide_ri(m_size[0],two_power(m_levels[0]))+3);

    cell_fn_ptr_t pData(new cell_fn_t[pData_size]);

    int pc_beg = two_power(m_level_ct)-1;
    int pc_end = pc_beg + two_power(m_level_ct);

    for(int pc_i = pc_beg ; pc_i < pc_end;++pc_i)
    {
      piece_ptr_t dp = m_pieces[pc_i];

      dp->m_dataset->init(dp->get_basename(m_basename)+".raw");
      dp->m_dataset->computeMsGraph(dp->m_msgraph);
    }
  }

  int get_aaplane_normal_dir(const rect_t & plane)
  {
    try
    {

      ASSERT(plane.eff_dim() ==gc_grid_dim-1);

      for( int d = 0 ; d < gc_grid_dim;++d)
      {
        if(plane[d][0] == plane[d][1])
          return d;
      }

      ASSERT(false);
    }
    catch(assertion_error e)
    {
      e.push(_FFL);
      e.push(SVAR(plane));

      throw;
    }
  }

  void data_manager_t::merge_up_subdomain_msgraphs ()
  {
    for(int lev = m_level_ct-1 ;lev >= 0 ;--lev)
    {
      int n = two_power(lev);

      for(int i = 0 ;i < n; ++i)
      {
        mscomplex_ptr_t msc  = m_pieces[n+i-1]->m_msgraph;
        mscomplex_ptr_t msc1 = m_pieces[(n+i)*2 -1]->m_msgraph;
        mscomplex_ptr_t msc2 = m_pieces[(n+i)*2]->m_msgraph;

        rect_t bnd = msc1->m_ext_rect.intersection(msc2->m_ext_rect);
        rect_t ixn = msc1->m_rect.intersection(msc2->m_rect);

        int d = get_aaplane_normal_dir(ixn);

        bnd[d] = ixn[d];

        msc->merge_up(*msc1,*msc2,bnd);
      }
    }
  }

  void data_manager_t::merge_down_subdomain_msgraphs ()
  {
    for(int lev = 0 ;lev < m_level_ct ;++lev)
    {
      int n = two_power(lev);

      for(int i = 0 ;i < n; ++i)
      {
        mscomplex_ptr_t msc  = m_pieces[n+i-1]->m_msgraph;
        mscomplex_ptr_t msc1 = m_pieces[(n+i)*2 -1]->m_msgraph;
        mscomplex_ptr_t msc2 = m_pieces[(n+i)*2]->m_msgraph;

        rect_t bnd = msc1->m_ext_rect.intersection(msc2->m_ext_rect);
        rect_t ixn = msc1->m_rect.intersection(msc2->m_rect);

        int d = get_aaplane_normal_dir(ixn);

        bnd[d] = ixn[d];

        msc->merge_down(*msc1,*msc2,bnd);
      }
    }
  }

  void data_manager_t::save_graphs()
  {
    for(int i = 0 ;i < m_pieces.size(); ++i)
    {
      piece_ptr_t dp = m_pieces[i];
      dp->m_msgraph->write_graph(dp->get_basename(m_basename)+".graph.txt");
    }
  }

  void data_manager_t::save_mfolds()
  {
    for(int i = two_power(m_level_ct)-1 ;i < m_pieces.size(); ++i)
    {
      piece_ptr_t dp = m_pieces[i];
      dp->m_msgraph->invert_for_collection();
      dp->m_dataset->saveManifolds(dp->m_msgraph,dp->get_basename(m_basename));
    }
  }

  void data_manager_t::destoryPieces()
  {
    m_pieces.clear();
  }

  void data_manager_t::work()
  {
    g_timer.start();

    cout<<"===================================="<<endl;
    cout<<"         Starting Processing        "<<endl;
    cout<<"------------------------------------"<<endl;

    createPieces();
    cout<<"create pieces done ------- "<<g_timer.getElapsedTimeInMilliSec()<<endl;

    split_dataset();
    cout<<"split dataset done ------- "<<g_timer.getElapsedTimeInMilliSec()<<endl;

    compute_subdomain_msgraphs();
    cout<<"subdomain msgraph done --- "<<g_timer.getElapsedTimeInMilliSec()<<endl;

    merge_up_subdomain_msgraphs();
    cout<<"merge up done ------------ "<<g_timer.getElapsedTimeInMilliSec()<<endl;

    merge_down_subdomain_msgraphs();
    cout<<"merge down done ---------- "<<g_timer.getElapsedTimeInMilliSec()<<endl;

    save_graphs();
    cout<<"save graphs done --------- "<<g_timer.getElapsedTimeInMilliSec()<<endl;

    save_mfolds();
    cout<<"save mfolds done --------- "<<g_timer.getElapsedTimeInMilliSec()<<endl;

    destoryPieces();
    cout<<"destroy pieces done ------ "<<g_timer.getElapsedTimeInMilliSec()<<endl;

    cout<<"------------------------------------"<<endl;
    cout<<"        Finished Processing         "<<endl;
    cout<<"===================================="<<endl;
  }

  data_manager_t::data_manager_t
      ( std::string filename,
        cellid_t     size,
        cellid_t     levels,
        double       simp_tresh):
      m_filename(filename),
      m_basename(basename(filename.c_str())),
      m_size(size),
      m_levels(levels),
      m_simp_tresh(simp_tresh)
  {
    m_level_ct   = m_levels[0] + m_levels[1] + m_levels[2];

    int ext_pos = m_basename.size() -4;

    if(ext_pos >=0 && m_basename.substr(ext_pos,4) == ".raw")
      m_basename = m_basename.substr(0,ext_pos);
  }

  data_manager_t::~data_manager_t ()
  {
  }


  void compute_mscomplex_basic(std::string filename, cellid_t size, double simp_tresh)
  {
    g_timer.start();

    cout<<"===================================="<<endl;
    cout<<"         Starting Processing        "<<endl;
    cout<<"------------------------------------"<<endl;

    rect_t d(cellid_t::zero,(size-cellid_t::one)*2);
    dataset_ptr_t   dataset(new dataset_t(d,d,d));
    mscomplex_ptr_t msgraph(new mscomplex_t(d,d));

    dataset->init(filename);
    cout<<"data read ---------------- "<<g_timer.getElapsedTimeInMilliSec()<<endl;

    dataset->computeMsGraph(msgraph);
    cout<<"msgraph done ------------- "<<g_timer.getElapsedTimeInMilliSec()<<endl;

    msgraph->simplify_un_simplify(simp_tresh);
    cout<<"simplification done ------ "<<g_timer.getElapsedTimeInMilliSec()<<endl;

    msgraph->write_graph("graph.txt");
    cout<<"write msgraph done ------- "<<g_timer.getElapsedTimeInMilliSec()<<endl;

    msgraph->invert_for_collection();
    dataset->saveManifolds(msgraph,"mfolds.bin");
    cout<<"write msmfolds done ------- "<<g_timer.getElapsedTimeInMilliSec()<<endl;

    cout<<"------------------------------------"<<endl;
    cout<<"        Finished Processing         "<<endl;
    cout<<"===================================="<<endl;
  }

  octtree_piece_t::octtree_piece_t (rect_t p,rect_t pd,int l):
      m_level(l),
      m_prct(p),
      m_ext_prct(pd.intersection(rect_t(p.lc()-1,p.uc()+1)))
  {
    rect_t r = rect_t( p.lc()*2,( p.uc()-1)*2);
    rect_t d = rect_t(pd.lc()*2,(pd.uc()-1)*2);

    rect_t e = d.intersection(rect_t((p.lc()-1)*2,p.uc()*2));

    ASSERT(r.intersection(d) == r);

    m_msgraph.reset(new mscomplex_t(r,e));

    if(l == 0)
    {
      m_dataset.reset(new dataset_t(r,e,d));
    }
  }

  std::string octtree_piece_t::get_basename(const std::string& basename)
  {
    return basename+"."+utls::to_string(m_ext_prct);
  }
}
