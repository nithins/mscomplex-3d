#include <sstream>
#include <cstring>
#include <iostream>

#include <boost/algorithm/string_regex.hpp>
#include <boost/lambda/lambda.hpp>

#include <GL/glew.h>

#include <glutils.h>
#include <GLSLProgram.h>
#include <logutil.h>

#include <grid_viewer.h>
#include <grid_datamanager.h>
#include <grid_mscomplex.h>
#include <grid_mscomplex_ensure.h>
#include <grid_dataset.h>

#include <shadersources.h>

GLSLProgram * s_cell_shaders[grid::GRADDIR_COUNT] = {NULL,NULL};

glutils::color_t g_grid_cp_colors[] =
{
  glutils::color_t(0.0,0.0,1.0),
  glutils::color_t(0.0,1.0,0.0),
  glutils::color_t(0.0,1.0,1.0),
  glutils::color_t(1.0,0.0,0.0),
};

glutils::color_t g_grid_grad_colors[] =
{
  glutils::color_t(0.0,0.5,0.5 ),
  glutils::color_t(0.5,0.0,0.5 ),
  glutils::color_t(0.5,0.5,0.0 ),
};

glutils::color_t g_disc_colors[grid::GRADDIR_COUNT][grid::gc_grid_dim+1] =
{
  {
    glutils::color_t(0.15,0.45,0.35 ),
    glutils::color_t(0.25,0.15,0.75 ),
    glutils::color_t(0.65,0.95,0.35 ),
    glutils::color_t(0.0,0.0,0.0 ),
  },

  {
    glutils::color_t(0.0,0.0,0.0 ),
    glutils::color_t(0.35,0.25,0.65 ),
    glutils::color_t(0.65,0.25,0.15 ),
    glutils::color_t(0.15,0.25,0.75 ),
  },
};

glutils::color_t g_grid_cp_conn_colors[] =
{
  glutils::color_t(0.0,0.5,0.5 ),
  glutils::color_t(0.5,0.0,0.5 ),
  glutils::color_t(0.5,0.5,0.0 ),
};

glutils::color_t g_roiaabb_color = glutils::color_t(0.85,0.75,0.65);

const char * shader_consts[grid::GRADDIR_COUNT]
    = {"const float even_sz = 0.1;"\
       "const float odd_sz  = 1.0;",
       "const float even_sz = 1.0;"\
       "const float odd_sz  = 0.1;"};

namespace grid
{
  void disc_rendata_t::init()
  {
    for(uint i = 0 ;i < GRADDIR_COUNT;++i)
    {

      if(s_cell_shaders[i] != NULL )
        continue;


      std::string geom_glsl(cell_shader_geom_glsl);

      boost::replace_regex
          ( geom_glsl,
            boost::regex("//HEADER_REPLACE_BEGIN(.*)//HEADER_REPLACE_END"),
            std::string(shader_consts[i]) );


      s_cell_shaders[i] = GLSLProgram::createFromSourceStrings
                          (cell_shader_vert_glsl,
                           geom_glsl,
                           std::string(),
                           GL_POINTS,GL_TRIANGLES);

      std::string log;

      s_cell_shaders[i]->GetProgramLog ( log );

      if(log.size() !=0 )
        std::cout<<"shader log ::\n"<<log<<"\n";
    }

  }

  void disc_rendata_t::cleanup()
  {
    for(uint i = 0 ;i < GRADDIR_COUNT;++i)
    {
      if(s_cell_shaders[i] != NULL )
        continue;

      delete s_cell_shaders[i];

      s_cell_shaders[i] = NULL;
    }
  }

  grid_viewer_t::grid_viewer_t
      (data_manager_t * gdm):
      m_size(gdm->m_size),m_scale_factor(0),
      m_bRebuildRens(true),m_bShowRoiBB(false),m_bCenterToRoi(false),
      m_gdm(gdm)

  {
    m_roi = rect_t(cellid_t::zero,(m_size-cellid_t::one)*2);

    for(uint i = 0 ;i < m_gdm->m_pieces.size();++i)
      m_grid_piece_rens.push_back(new octtree_piece_rendata(m_gdm->m_pieces.at(i)));

    m_roi_base_pt  = ((m_roi.upper_corner() +  m_roi.lower_corner())/2);

  }

  grid_viewer_t::~grid_viewer_t()
  {
    for ( uint i = 0 ; i < m_grid_piece_rens.size();i++ )
      delete m_grid_piece_rens[i];

    m_grid_piece_rens.clear();

    disc_rendata_t::cleanup();

    delete m_gdm;
  }

  void grid_viewer_t::set_roi_dim_range_nrm(double l,double u,int dim)
  {
    if(!(l<u && 0.0 <= l && u <=1.0 && 0<=dim && dim < 3))
      return;

    rect_t roi = rect_t(cellid_t::zero,(m_size-cellid_t::one)*2);

    double span = roi[dim].span();

    m_roi[dim][0]  = (uint)(l*span);
    m_roi[dim][1]  = (uint)(u*span);

    m_roi_base_pt  = ((m_roi.upper_corner() +  m_roi.lower_corner())/2);
  }

  int grid_viewer_t::render()
  {
    if(m_bRebuildRens)
    {
      build_rens();

      m_bRebuildRens = false;
    }

    glPushAttrib(GL_ENABLE_BIT);

    glEnable(GL_NORMALIZE);

    glScalef(m_scale_factor,
             m_scale_factor,
             m_scale_factor);

    if(m_bCenterToRoi)
      glTranslatef(-m_roi_base_pt[0],
                   -m_roi_base_pt[1],
                   -m_roi_base_pt[2]);
    else
      glTranslatef(std::min(-m_size[0]+1,-1),
                   std::min(-m_size[1]+1,-1),
                   std::min(-m_size[2]+1,-1));

    if(m_bShowRoiBB)
    {
      glPushAttrib(GL_ENABLE_BIT);

      glDisable(GL_LIGHTING);

      glColor3dv(g_roiaabb_color.data());

      glutils::draw_aabb_line(m_roi.lower_corner(),m_roi.upper_corner());

      glPopAttrib();
    }

    for ( uint i = 0 ; i < m_grid_piece_rens.size();i++ )
    {
      m_grid_piece_rens[i]->render();
    }

    glPopAttrib();
  }

  void grid_viewer_t::build_rens()
  {
    for ( uint i = 0 ; i < m_grid_piece_rens.size();i++ )
    {
      m_grid_piece_rens[i]->create_cp_rens(m_roi);
      m_grid_piece_rens[i]->create_grad_rens(m_roi);
      m_grid_piece_rens[i]->create_structure_rens(m_roi);
    }
  }

  void grid_viewer_t::init()
  {
    glutils::init();

    disc_rendata_t::init();

    m_scale_factor = 0.5/ std::max((double) *std::max_element
                                   (m_size.begin(),m_size.end())-1.0,1.0);

    /*turn back face culling off */
    glEnable ( GL_CULL_FACE );

    /*cull backface */
    glCullFace ( GL_BACK );

    /*polymode */
    glPolygonMode ( GL_FRONT, GL_FILL );

    glPolygonMode ( GL_BACK, GL_LINE );

    for ( uint i = 0 ; i < m_grid_piece_rens.size();i++ )
    {
      m_grid_piece_rens[i]->create_cp_loc_bo();
      m_grid_piece_rens[i]->create_disc_rds();
    }
  }

  configurable_t::data_index_t grid_viewer_t::dim()
  {
    return data_index_t(13,m_grid_piece_rens.size());
  }
  bool grid_viewer_t::exchange_field(const data_index_t &i,boost::any &v)
  {
    octtree_piece_rendata * dprd = m_grid_piece_rens[i[1]];

    switch(i[0])
    {
    case 0: return s_exchange_data_ro(dprd->dp->label(),v);
    case 1: return s_exchange_data_rw(dprd->m_bShowAllCps,v);
    case 2: return s_exchange_data_rw(dprd->m_bShowCps[0],v);
    case 3: return s_exchange_data_rw(dprd->m_bShowCps[1],v);
    case 4: return s_exchange_data_rw(dprd->m_bShowCps[2],v);
    case 5: return s_exchange_data_rw(dprd->m_bShowCps[3],v);
    case 6: return s_exchange_data_rw(dprd->m_bShowCpLabels,v);
    case 7: return s_exchange_data_rw(dprd->m_bShowMsGraph,v);
    case 8: return s_exchange_data_rw(dprd->m_bShowGrad,v);
    case 9: return s_exchange_data_rw(dprd->m_bShowCancCps,v);
    case 10: return s_exchange_data_rw(dprd->m_bShowCancMsGraph,v);
    case 11: return s_exchange_data_rw(dprd->m_bShowStructure[1],v);
    case 12: return s_exchange_data_rw(dprd->m_bShowStructure[2],v);

    }

    throw std::logic_error("unknown index");
  }
  configurable_t::eFieldType grid_viewer_t::exchange_header
      (const int &i,boost::any &v)
  {
    switch(i)
    {

    case 0: v =  std::string("oct tree piece"); return EFT_DATA_RO;
    case 1: v =  std::string("all cps");return EFT_DATA_RW;
    case 2: v =  std::string("minima");return EFT_DATA_RW;
    case 3: v =  std::string("1 saddle");return EFT_DATA_RW;
    case 4: v =  std::string("2 saddle");return EFT_DATA_RW;
    case 5: v =  std::string("maxima");return EFT_DATA_RW;
    case 6: v =  std::string("cp labels");return EFT_DATA_RW;
    case 7: v =  std::string("msgraph");return EFT_DATA_RW;
    case 8: v =  std::string("gradient");return EFT_DATA_RW;
    case 9: v =  std::string("cancelled cps");return EFT_DATA_RW;
    case 10: v =  std::string("cancelled cp msgraph");return EFT_DATA_RW;
    case 11: v =  std::string("structure 1");return EFT_DATA_RW;
    case 12: v =  std::string("structure 2");return EFT_DATA_RW;

    }

    throw std::logic_error("unknown index");
  }

  octtree_piece_rendata::octtree_piece_rendata (octtree_piece * _dp):
      m_bShowAllCps(false),
      m_bShowCpLabels ( false ),
      m_bShowMsGraph ( false ),
      m_bShowGrad ( false ),
      m_bShowCancCps(false),
      m_bShowCancMsGraph(false),
      m_bNeedUpdateDiscRens(false),
      dp(_dp)
  {
    using namespace boost::lambda;

    std::for_each(m_bShowCps,m_bShowCps+gc_grid_dim+1,_1 = false);
    std::for_each(m_bShowStructure,m_bShowStructure+gc_grid_dim+1,_1 = false);

  }

  void octtree_piece_rendata::create_cp_loc_bo()
  {
    if(dp->msgraph == NULL)
      return;

    std::vector<glutils::vertex_t>  cp_loc;

    for(uint i = 0; i < dp->msgraph->m_cps.size(); ++i)
    {
      cellid_t c = (dp->msgraph->m_cps[i]->cellid);
      cp_loc.push_back(glutils::vertex_t(c[0],c[1],c[2]));
    }

    cp_loc_bo = glutils::make_buf_obj(cp_loc);
  }

  void  octtree_piece_rendata::create_cp_rens(const rect_t & roi)
  {
    if(dp->msgraph == NULL)
      return;

    glutils::string_list_t              crit_labels[gc_grid_dim+1];
    glutils::vertex_list_t              crit_label_locations[gc_grid_dim+1];
    glutils::point_idx_list_t           crit_pt_idxs[gc_grid_dim+1];
    glutils::line_idx_list_t            crit_conn_idxs[gc_grid_dim];

    for(uint i = 0; i < dp->msgraph->m_cps.size(); ++i)
    {
      if(dp->msgraph->m_cps[i]->is_paired)
        continue;

      cellid_t c = (dp->msgraph->m_cps[i]->cellid);

      if(!roi.contains(c))
        continue;

      uint index = dp->msgraph->m_cps[i]->index;

      std::stringstream ss;

      ((std::ostream&)ss)<<c;

      if(!dp->msgraph->m_cps[i]->is_paired)
      {
        crit_labels[index].push_back(ss.str());
        crit_label_locations[index].push_back(glutils::vertex_t(c[0],c[1],c[2]) );
        crit_pt_idxs[index].push_back(i);
      }
    }

    for(uint i = 0 ; i < gc_grid_dim+1; ++i)
    {
      ren_cp_labels[i].reset(glutils::create_buffered_text_ren
                             (crit_labels[i],crit_label_locations[i]));

      ren_cp[i].reset(glutils::create_buffered_points_ren
                      (cp_loc_bo,glutils::make_buf_obj(crit_pt_idxs[i])));
    }

    for(uint i = 0 ; i < dp->msgraph->m_cps.size(); ++i)
    {
      if(dp->msgraph->m_cps[i]->isCancelled)
        continue;

      if(dp->msgraph->m_cps[i]->is_paired)
        continue;

      cellid_t c = (dp->msgraph->m_cps[i]->cellid);

      if(!roi.contains(c))
        continue;

      uint index = dp->msgraph->m_cps[i]->index;

      for(conn_iter_t it  = dp->msgraph->m_cps[i]->conn[0].begin();
      it != dp->msgraph->m_cps[i]->conn[0].end(); ++it)
      {
        if(!roi.contains(dp->msgraph->m_cps[*it]->cellid))
          continue;

        crit_conn_idxs[index-1].push_back
            (glutils::line_idx_t(i,*it));
      }
    }

    for(uint i = 0 ; i < gc_grid_dim; ++i)
    {
      ren_cp_conns[i].reset(glutils::create_buffered_lines_ren
                            (cp_loc_bo,
                             glutils::make_buf_obj(crit_conn_idxs[i])));
    }

  }

  void octtree_piece_rendata::create_structure_rens(const rect_t & roi)
  {
    if(dp->dataset == NULL)
      return;

    rect_t r;
    if(!dp->dataset->get_ext_rect().intersection(roi,r))
      return;

    std::vector<glutils::vertex_t>  cell_locations[gc_grid_dim+1];

    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    cellid_t c;

    for(c[2] = r[2][0] ; c[2] <= r[2][1]; ++c[2])
    {
      for(c[1] = r[1][0] ; c[1] <= r[1][1]; ++c[1])
      {
        for(c[0] = r[0][0] ; c[0] <= r[0][1]; ++c[0])
        {
          cell_locations[(*dp->dataset->m_cell_eff_dim)(c)].push_back
              (glutils::vertex_t(c[0],c[1],c[2]));
        }
      }
    }

    for(uint i = 0 ; i < gc_grid_dim+1; ++i)
    {
      ren_structure[i].reset(glutils::create_buffered_points_ren
                        (glutils::make_buf_obj(cell_locations[i])));
    }

  }

  void octtree_piece_rendata::create_grad_rens(const rect_t & roi)
  {
    if(dp->dataset == NULL)
      return;

    rect_t r;
    if(!dp->dataset->get_ext_rect().intersection(roi,r))
      return;

    std::vector<glutils::vertex_t>      cell_locations;
    std::vector<glutils::line_idx_t>    pair_idxs[gc_grid_dim];

    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    cellid_t c;

    for(c[2] = r[2][0] ; c[2] <= r[2][1]; ++c[2])
    {
      for(c[1] = r[1][0] ; c[1] <= r[1][1]; ++c[1])
      {
        for(c[0] = r[0][0] ; c[0] <= r[0][1]; ++c[0])
        {
          uint dim = dp->dataset->getCellDim(c);

          if(dp->dataset->isCellPaired(c))
          {
            cellid_t p = dp->dataset->getCellPairId(c);

            if(dp->dataset->isPairOrientationCorrect(c,p))
            {
              cell_locations.push_back(glutils::vertex_t(c[0],c[1],c[2]) );

              cell_locations.push_back(glutils::vertex_t(p[0],p[1],p[2]) );

              pair_idxs[dim].push_back
                  (glutils::line_idx_t(cell_locations.size()-2,
                                       cell_locations.size()-1));
            }
          }
        }
      }
    }

    glutils::bufobj_ptr_t cell_bo= glutils::make_buf_obj(cell_locations);

    for(uint i = 0 ; i < gc_grid_dim; ++i)

    {
      ren_grad[i].reset(glutils::create_buffered_lines_ren
                        (cell_bo,
                         glutils::make_buf_obj(pair_idxs[i])));
    }
  }


  void octtree_piece_rendata::create_disc_rds()
  {
    if(dp->msgraph == NULL)
      return;

    boost::shared_ptr<disc_rendata_t> sptr;

    for(uint i = 0 ; i < dp->msgraph->m_cps.size();++i)
    {
      critpt_t * cp = dp->msgraph->m_cps[i];

      if(cp->is_paired) continue;

      sptr.reset(new disc_rendata_t(cp->cellid,cp->index));

      disc_rds.push_back(sptr);
    }
  }

  void octtree_piece_rendata::update_active_disc_rens()
  {
    if(dp->msgraph == NULL)
      return;

    for(uint i = 0 ; i < disc_rds.size();++i)
    {
      if(disc_rds[i]->update(dp->msgraph))
      {
        if(active_disc_rens.count(disc_rds[i]) == 0)
        {
          active_disc_rens.insert(disc_rds[i]);
        }
        else
        {
          active_disc_rens.erase(disc_rds[i]);
        }
      }
    }
  }

  void octtree_piece_rendata::render()
  {
    if(m_bNeedUpdateDiscRens)
    {
      update_active_disc_rens();
      m_bNeedUpdateDiscRens = false;
    }

    glPushMatrix();
    glPushAttrib ( GL_ENABLE_BIT );

    glDisable ( GL_LIGHTING );

    if(m_bShowGrad)
    {
      for(uint i = 0 ; i < gc_grid_dim; ++i)
      {
        if(ren_grad[i])
        {
          glColor3dv ( g_grid_grad_colors[i].data() );

          ren_grad[i]->render();
        }
      }
    }

    for(uint i = 0 ; i < gc_grid_dim+1; ++i)
    {
      s_cell_shaders[grid::GRADDIR_ASCENDING]->use();

      if(ren_grad[i] && m_bShowStructure[i])
      {
        glColor3dv ( g_disc_colors[grid::GRADDIR_ASCENDING][i].data() );

        ren_structure[i]->render();
      }

      s_cell_shaders[grid::GRADDIR_ASCENDING]->disable();
    }

    glPointSize ( 4.0 );

    for(uint i = 0 ; i < gc_grid_dim+1;++i)
    {
      if(ren_cp[i]&& (m_bShowCps[i]||m_bShowAllCps))
      {
        glColor3dv(g_grid_cp_colors[i].data());

        ren_cp[i]->render();

        if(ren_cp_labels[i] && m_bShowCpLabels)
          ren_cp_labels[i]->render();
      }
    }


    if ( m_bShowCancCps)
    {
      for(uint i = 0 ; i < gc_grid_dim;++i)
      {
        if(ren_canc_cp[i])
        {
          glColor3dv(g_grid_cp_colors[i].data());

          ren_canc_cp[i]->render();

          if(ren_canc_cp_labels[i] &&
             m_bShowCpLabels)
            ren_canc_cp_labels[i]->render();
        }
      }
    }

    if (m_bShowMsGraph)
    {
      for(uint i = 0 ; i < gc_grid_dim;++i)
      {
        if(ren_cp_conns[i])
        {
          glColor3dv(g_grid_cp_conn_colors[i].data());

          ren_cp_conns[i]->render();
        }
      }
    }

    if (m_bShowCancMsGraph)
    {
      for(uint i = 0 ; i < gc_grid_dim;++i)
      {
        if(ren_canc_cp_conns[i])
        {
          glColor3dv(g_grid_cp_conn_colors[i].data());

          ren_canc_cp_conns[i]->render();
        }
      }
    }

    for(uint i = 0 ; i < disc_rds.size();++i)
    {
      disc_rds[i]->render();
    }


    glPopAttrib();
    glPopMatrix();
  }

  configurable_t::data_index_t octtree_piece_rendata::dim()
  {
    return data_index_t(6,disc_rds.size());
  }

  bool octtree_piece_rendata::exchange_field(const data_index_t &i,boost::any &v)
  {
    boost::shared_ptr<disc_rendata_t> drd = disc_rds[i[1]];

    switch(i[0])
    {
    case 0:
      return s_exchange_data_ro(drd->cellid.to_string(),v);
    case 1:
      return s_exchange_data_ro((int)drd->index,v);
    case 2:
    case 3:
      {
        bool need_update = false;

        bool is_read     = v.empty();

        need_update =  s_exchange_data_rw(drd->show[i[0]%2],v);

        if(need_update && is_read == false )
          m_bNeedUpdateDiscRens = true;

        return need_update;
      }
    case 4:
    case 5:
      return s_exchange_data_rw(drd->color[i[0]%2],v);
    };

    throw std::logic_error("invalid index");
  }

  configurable_t::eFieldType octtree_piece_rendata::exchange_header
      (const int &i,boost::any &v)
  {
    switch(i)
    {
    case 0: v = std::string("cellid"); return EFT_DATA_RO;
    case 1: v = std::string("index"); return EFT_DATA_RO;
    case 2: v = std::string("des mfold"); return EFT_DATA_RW;
    case 3: v = std::string("asc mfold"); return EFT_DATA_RW;
    case 4: v = std::string("des mfold color"); return EFT_DATA_RW;
    case 5: v = std::string("asc mfold color"); return EFT_DATA_RW;
    }
    throw std::logic_error("invalid index");
  }

  disc_rendata_t::disc_rendata_t(cellid_t c,uint i):cellid(c),index(i)
  {
    color[0] = g_disc_colors[1][index];
    color[1] = g_disc_colors[0][index];

    show[0] =false; ren[0] =NULL;
    show[1] =false; ren[1] =NULL;
  }

  disc_rendata_t::~disc_rendata_t()
  {
    show[0] =false;
    show[1] =false;

    update(NULL);

  }

  void disc_rendata_t::render()
  {
    for(uint dir = 0 ; dir<2;++dir)
    {

      s_cell_shaders[dir]->use();

      if(show[dir])
      {
        glColor3dv(g_grid_cp_colors[index].data());

        glBegin(GL_POINTS);
        glVertex3sv(cellid.data());
        glEnd();

        glColor3dv(color[dir].data());

        ren[dir]->render();
      }

      s_cell_shaders[dir]->disable();
    }
  }

  bool disc_rendata_t::update(mscomplex_t *msc )
  {
    uint ret = false;

    for(uint dir = 0 ; dir<2;++dir)
    {
      if(show[dir] && this->ren[dir] == NULL && msc)
      {
        ensure_cellid_critical(msc,cellid);

        critpt_t *cp = msc->m_cps[msc->m_id_cp_map[cellid]];

        std::set<cellid_t> vset;

        for(uint j = 0 ; j < cp->contrib[dir].size();++j)
        {

          critpt_t *cp_contrib = msc->m_cps[cp->contrib[dir][j]];


          for(uint i = 0; i < cp_contrib->disc[dir].size(); ++i)
          {
            cellid_t c = cp_contrib->disc[dir][i];

            if(vset.count(c) == 0)
              vset.insert(c);
          }
        }

        std::vector<glutils::vertex_t> vlist(vset.size());
        std::copy(vset.begin(),vset.end(),vlist.begin());

        ren[dir] = glutils::create_buffered_points_ren
                   (glutils::make_buf_obj(vlist));

        ret = true;
      }

      if(!show[dir] && this->ren[dir] != NULL )
      {
        delete ren[dir];

        ren[dir] = NULL;

        ret = true;

      }
    }
    return ret;
  }

}
