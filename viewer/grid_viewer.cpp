#include <sstream>

#include <GL/glew.h>

#include <glutils.h>

#include <grid_viewer.h>

namespace grid
{

  glviewer_t::glviewer_t
      (std::vector<octtree_piece *> * p ,uint size_x,uint size_y):
      m_size_x(size_x),
      m_size_y(size_y)
  {
    for(uint i = 0 ;i < p->size();++i)
      m_grid_piece_rens.push_back(new octtree_piece_rendata(p->at(i)));
  }

  glviewer_t::~glviewer_t()
  {
    for ( uint i = 0 ; i < m_grid_piece_rens.size();i++ )
      delete m_grid_piece_rens[i];

    m_grid_piece_rens.clear();
  }

  void glviewer_t::draw()
  {

    glPushAttrib(GL_ENABLE_BIT);

    glEnable(GL_NORMALIZE);

    glTranslatef(-1,0,-1);

    glScalef(0.5/(double) m_size_x,0.1,0.5/(double) m_size_y);

    for ( uint i = 0 ; i < m_grid_piece_rens.size();i++ )
    {
      m_grid_piece_rens[i]->render();
    }

    glPopAttrib();
  }

  void glviewer_t::init()
  {
    glutils::init();

    // Restore previous viewer state.
    restoreStateFromFile();

    for ( uint i = 0 ; i < m_grid_piece_rens.size();i++ )
    {
      m_grid_piece_rens[i]->create_cp_rens();
      m_grid_piece_rens[i]->create_grad_rens();
      m_grid_piece_rens[i]->create_surf_ren();

    }

  }

  QString glviewer_t::helpString() const
  {
    QString text("<h2>MS Complex Viewer</h2>");
    return text;
  }

  octtree_piece_rendata::octtree_piece_rendata (octtree_piece * _dp):
      m_bShowSurface ( false ),
      m_bShowCps ( false ),
      m_bShowCpLabels ( false ),
      m_bShowMsGraph ( false ),
      m_bShowGrad ( false ),
      m_bShowCancCps(false),
      m_bShowCancMsGraph(false),
      dp(_dp)
  {
    memset(ren_cp,0,sizeof(ren_cp));
    memset(&ren_surf,0,sizeof(ren_surf));
    memset(ren_grad,0,sizeof(ren_grad));
    memset(ren_cp_labels,0,sizeof(ren_cp_labels));
    memset(ren_cp_conns,0,sizeof(ren_cp_conns));
    memset(ren_canc_cp,0,sizeof(ren_canc_cp));
    memset(ren_canc_cp_labels,0,sizeof(ren_canc_cp_labels));
    memset(ren_canc_cp_conns,0,sizeof(ren_canc_cp_conns));
  }

  void  octtree_piece_rendata::create_cp_rens()
  {
    if(dp->msgraph == NULL)
      return;

    std::vector<std::string>            crit_labels[gc_grid_dim+1];
    std::vector<glutils::vertex_t>      crit_label_locations[gc_grid_dim+1];
    std::vector<glutils::point_idx_t>   crit_pt_idxs[gc_grid_dim+1];
    std::vector<glutils::line_idx_t>    crit_conn_idxs[gc_grid_dim];


    std::vector<std::string>            crit_canc_labels[gc_grid_dim+1];
    std::vector<glutils::vertex_t>      crit_canc_label_locations[gc_grid_dim+1];
    std::vector<glutils::point_idx_t>   crit_canc_pt_idxs[gc_grid_dim+1];
    std::vector<glutils::line_idx_t>    crit_canc_conn_idxs[gc_grid_dim];

    std::vector<glutils::vertex_t>      crit_locations;
    std::map<uint,uint>                 crit_ms_idx_ren_idx_map;

    for(uint i = 0; i < dp->msgraph->m_cps.size(); ++i)
    {
      if(dp->msgraph->m_cps[i]->isCancelled)
        continue;

      cellid_t c = (dp->msgraph->m_cps[i]->cellid);

      uint dim = dataset_t::s_getCellDim(c);

      std::stringstream ss;

      ((std::ostream&)ss)<<c;

      if(dp->msgraph->m_cps[i]->isBoundryCancelable)
      {
        crit_canc_labels[dim].push_back(ss.str());
        crit_canc_label_locations[dim].push_back(glutils::vertex_t(c[0],c[1],c[2]) );
        crit_canc_pt_idxs[dim].push_back(glutils::point_idx_t(crit_locations.size()));
      }
      else
      {
        crit_labels[dim].push_back(ss.str());
        crit_label_locations[dim].push_back(glutils::vertex_t(c[0],c[1],c[2]) );
        crit_pt_idxs[dim].push_back(glutils::point_idx_t(crit_locations.size()));
      }

      crit_ms_idx_ren_idx_map[i] = crit_locations.size();
      crit_locations.push_back(glutils::vertex_t(c[0],c[1],c[2]));
    }

    glutils::bufobj_ptr_t crit_loc_bo = glutils::make_buf_obj(crit_locations);

    for(uint i = 0 ; i < gc_grid_dim+1; ++i)
    {
      ren_cp_labels[i] =
          glutils::create_buffered_text_ren
          (crit_labels[i],crit_label_locations[i]);

      ren_cp[i] =glutils::create_buffered_points_ren
                 (crit_loc_bo,
                  glutils::make_buf_obj(crit_pt_idxs[i]),
                  glutils::make_buf_obj());

      ren_canc_cp_labels[i] =
          glutils::create_buffered_text_ren
          (crit_canc_labels[i],crit_canc_label_locations[i]);

      ren_canc_cp[i] =glutils::create_buffered_points_ren
                      (crit_loc_bo,
                       glutils::make_buf_obj(crit_canc_pt_idxs[i]),
                       glutils::make_buf_obj());
    }

    for(uint i = 0 ; i < dp->msgraph->m_cps.size(); ++i)
    {
      if(dp->msgraph->m_cps[i]->isCancelled)
        continue;

      conn_t *cp_acdc[] = {&dp->msgraph->m_cps[i]->des,&dp->msgraph->m_cps[i]->asc};

      uint acdc_ct = 1;

      if(dp->msgraph->m_cps[i]->isBoundryCancelable)
        acdc_ct = 2;

      uint cp_ren_idx = crit_ms_idx_ren_idx_map[i];

      uint dim = dataset_t::s_getCellDim
                 (dp->msgraph->m_cps[i]->cellid);

      for (uint j = 0 ; j < acdc_ct; ++j)
      {
        for(conn_t::iterator it = cp_acdc[j]->begin();
        it != cp_acdc[j]->end(); ++it)
        {
          if(dp->msgraph->m_cps[*it]->isCancelled)
            throw std::logic_error("this cancelled cp should not be present here");

          if(dp->msgraph->m_cps[*it]->isBoundryCancelable)
            throw std::logic_error("a true cp should not be connected to a bc cp");

          uint conn_cp_ren_idx = crit_ms_idx_ren_idx_map[*it];

          if(dp->msgraph->m_cps[i]->isBoundryCancelable)
          {
            crit_canc_conn_idxs[dim-1+j].push_back
                (glutils::line_idx_t(cp_ren_idx,conn_cp_ren_idx));
          }
          else
          {
            crit_conn_idxs[dim-1+j].push_back
                (glutils::line_idx_t(cp_ren_idx,conn_cp_ren_idx));
          }
        }
      }
    }

    for(uint i = 0 ; i < gc_grid_dim; ++i)
    {
      ren_cp_conns[i] = glutils::create_buffered_lines_ren
                        (crit_loc_bo,
                         glutils::make_buf_obj(crit_conn_idxs[i]),
                         glutils::make_buf_obj());

      ren_canc_cp_conns[i] = glutils::create_buffered_lines_ren
                             (crit_loc_bo,
                              glutils::make_buf_obj(crit_canc_conn_idxs[i]),
                              glutils::make_buf_obj());
    }

  }

  void octtree_piece_rendata::create_grad_rens()
  {
    if(dp->dataset == NULL)
      return;

    rect_t r = dp->dataset->get_ext_rect();

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
          uint dim = dataset_t::s_getCellDim(c);

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

    ren_grad[0] = glutils::create_buffered_lines_ren
                  (cell_bo,
                   glutils::make_buf_obj(pair_idxs[0]),
                   glutils::make_buf_obj());

    ren_grad[1] = glutils::create_buffered_lines_ren
                  (cell_bo,
                   glutils::make_buf_obj(pair_idxs[1]),
                   glutils::make_buf_obj());


  }

  void octtree_piece_rendata::create_surf_ren()
  {
    if(dp->dataset == NULL)
      return;

#warning "create surf ren not implemented"

  }

  glutils::color_t g_grid_cp_colors[] =
  {
    glutils::color_t(0.0,0.0,1.0),
    glutils::color_t(0.0,1.0,0.0),
    glutils::color_t(1.0,0.0,0.0),
  };


  void octtree_piece_rendata::render() const
  {
    glPushMatrix();
    glPushAttrib ( GL_ENABLE_BIT );

    glScalef ( 2.0,2.0,2.0 );
    glTranslatef ( -0.5,0.0,-0.5 );

    if ( m_bShowSurface && ren_surf)
    {
      glColor3f ( 0.75,0.75,0.75 );
      ren_surf->render();
    }

    glDisable ( GL_LIGHTING );

    glTranslatef ( 0.0,0.02,0.0 );

    if ( m_bShowGrad && ren_grad[0] && ren_grad[1])
    {
      glColor3f ( 0.5,0.0,0.5 );
      ren_grad[0]->render();

      glColor3f ( 0.0,0.5,0.5 );
      ren_grad[1]->render();

    }

    glPointSize ( 4.0 );

    if ( m_bShowCps)
    {
      for(uint i = 0 ; i < 3;++i)
      {
        if(ren_cp[i])
        {
          glColor3dv(g_grid_cp_colors[i].data());

          ren_cp[i]->render();

          if(ren_cp_labels[i] && m_bShowCpLabels)
            ren_cp_labels[i]->render();
        }
      }
    }

    if ( m_bShowCancCps)
    {
      for(uint i = 0 ; i < 3;++i)
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

    if ( m_bShowMsGraph && ren_cp_conns[0] && ren_cp_conns[1])
    {
      glColor3f ( 0.0,0.5,1.0 );
      ren_cp_conns[0]->render();

      glColor3f ( 1.0,0.5,0.0 );
      ren_cp_conns[1]->render();

    }

    if ( m_bShowCancMsGraph&& ren_canc_cp_conns[0] && ren_canc_cp_conns[1])
    {
      glColor3f ( 0.0,0.5,1.0 );
      ren_canc_cp_conns[0]->render();

      glColor3f ( 1.0,0.5,0.0 );
      ren_canc_cp_conns[1]->render();

    }

    glPopAttrib();
    glPopMatrix();
  }
}
