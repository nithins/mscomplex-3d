#include <sstream>

#include <GL/glew.h>

#include <glutils.h>

#include <grid_viewer.h>
#include <grid_datamanager.h>
#include <grid_mscomplex.h>
#include <grid_dataset.h>
#include <cstring>

glutils::color_t g_grid_cp_colors[] =
{
  glutils::color_t(0.0,0.0,1.0),
  glutils::color_t(0.0,1.0,0.0),
  glutils::color_t(1.0,0.0,0.0),
  glutils::color_t(1.0,0.0,1.0),
};

glutils::color_t g_grid_grad_colors[] =
{
  glutils::color_t(0.0,0.5,0.5 ),
  glutils::color_t(0.5,0.0,0.5 ),
  glutils::color_t(0.5,0.5,0.0 ),
};

glutils::color_t g_grid_cp_conn_colors[] =
{
  glutils::color_t(0.0,0.5,0.5 ),
  glutils::color_t(0.5,0.0,0.5 ),
  glutils::color_t(0.5,0.5,0.0 ),
};

namespace grid
{

  grid_viewer_t::grid_viewer_t
      (data_manager_t * gdm,const rect_t &roi):
      m_size(gdm->m_size)
  {
    m_roi = rect_t(cellid_t::zero,(m_size-cellid_t::one)*2);

    rect_t s_roi(roi.lower_corner()*2,roi.upper_corner()*2);

    if(m_roi.intersects(s_roi))
      m_roi.intersection(s_roi,m_roi);

    for(uint i = 0 ;i < gdm->m_pieces.size();++i)
      m_grid_piece_rens.push_back(new octtree_piece_rendata(gdm->m_pieces.at(i)));

  }

  grid_viewer_t::~grid_viewer_t()
  {
    for ( uint i = 0 ; i < m_grid_piece_rens.size();i++ )
      delete m_grid_piece_rens[i];

    m_grid_piece_rens.clear();
  }

  int grid_viewer_t::render()
  {

    glPushAttrib(GL_ENABLE_BIT);

    glEnable(GL_NORMALIZE);

    glTranslatef(-0.5,-0.5,-0.5);

    glScalef(0.5/std::max(1.0,(double)(m_size[0]-1)),
             0.5/std::max(1.0,(double)(m_size[1]-1)),
             0.5/std::max(1.0,(double)(m_size[2]-1)));

    for ( uint i = 0 ; i < m_grid_piece_rens.size();i++ )
    {
      m_grid_piece_rens[i]->render();
    }

    glPopAttrib();
  }

  void grid_viewer_t::init()
  {
    glutils::init();

    for ( uint i = 0 ; i < m_grid_piece_rens.size();i++ )
    {
      m_grid_piece_rens[i]->create_cp_loc_bo();

      m_grid_piece_rens[i]->create_cp_rens(m_roi);
      m_grid_piece_rens[i]->create_grad_rens(m_roi);
      m_grid_piece_rens[i]->create_surf_ren(m_roi);
      m_grid_piece_rens[i]->create_disc_rds();
    }

  }

  int grid_viewer_t::rows()
  {
    return m_grid_piece_rens.size();

  }
  int grid_viewer_t::columns()
  {
    return 8;
  }
  bool grid_viewer_t::exchange_data(const data_index_t &idx,
                                    boost::any &v,
                                    const eExchangeMode &m)
  {

    switch(idx[0])
    {
    case 0: return s_exchange_read_only(m_grid_piece_rens[idx[1]]->dp->label(),v,m);
    case 1: return s_exchange_read_write(m_grid_piece_rens[idx[1]]->m_bShowSurface,v,m);
    case 2: return s_exchange_read_write(m_grid_piece_rens[idx[1]]->m_bShowCps,v,m);
    case 3: return s_exchange_read_write(m_grid_piece_rens[idx[1]]->m_bShowCpLabels,v,m);
    case 4: return s_exchange_read_write(m_grid_piece_rens[idx[1]]->m_bShowMsGraph,v,m);
    case 5: return s_exchange_read_write(m_grid_piece_rens[idx[1]]->m_bShowGrad,v,m);
    case 6: return s_exchange_read_write(m_grid_piece_rens[idx[1]]->m_bShowCancCps,v,m);
    case 7: return s_exchange_read_write(m_grid_piece_rens[idx[1]]->m_bShowCancMsGraph,v,m);
    }

    throw std::logic_error("unknown index");
  }
  std::string grid_viewer_t::get_header(int i)
  {
    switch(i)
    {

    case 0: return "oct tree piece";
    case 1: return "surface";
    case 2: return "cps";
    case 3: return "cp labels";
    case 4: return "msgraph";
    case 5: return "gradient";
    case 6: return "cancelled cps";
    case 7: return "cancelled cp msgraph";
    }

    throw std::logic_error("unknown index");
  }


  octtree_piece_rendata::octtree_piece_rendata (octtree_piece * _dp):
      m_bShowSurface ( false ),
      m_bShowCps ( false ),
      m_bShowCpLabels ( false ),
      m_bShowMsGraph ( false ),
      m_bShowGrad ( false ),
      m_bShowCancCps(false),
      m_bShowCancMsGraph(false),
      m_bNeedUpdateDiscRens(false),
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

    std::vector<std::string>            crit_labels[gc_grid_dim+1];
    std::vector<glutils::vertex_t>      crit_label_locations[gc_grid_dim+1];
    std::vector<glutils::point_idx_t>   crit_pt_idxs[gc_grid_dim+1];
    std::vector<glutils::line_idx_t>    crit_conn_idxs[gc_grid_dim];

    for(uint i = 0; i < dp->msgraph->m_cps.size(); ++i)
    {
      if(dp->msgraph->m_cps[i]->isCancelled)
        continue;

      cellid_t c = (dp->msgraph->m_cps[i]->cellid);

      if(!roi.contains(c))
        continue;

      uint dim = dataset_t::s_getCellDim(c);

      std::stringstream ss;

      ((std::ostream&)ss)<<c;

      if(!dp->msgraph->m_cps[i]->is_paired)
      {
        crit_labels[dim].push_back(ss.str());
        crit_label_locations[dim].push_back(glutils::vertex_t(c[0],c[1],c[2]) );
        crit_pt_idxs[dim].push_back(i);
      }
    }

    for(uint i = 0 ; i < gc_grid_dim+1; ++i)
    {
      ren_cp_labels[i] =
          glutils::create_buffered_text_ren
          (crit_labels[i],crit_label_locations[i]);

      ren_cp[i] =glutils::create_buffered_points_ren
                 (cp_loc_bo,
                  glutils::make_buf_obj(crit_pt_idxs[i]),
                  glutils::make_buf_obj());
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

      uint dim = dataset_t::s_getCellDim(c);

      for(conn_t::iterator it  = dp->msgraph->m_cps[i]->des.begin();
      it != dp->msgraph->m_cps[i]->des.end(); ++it)
      {
        if(!roi.contains(dp->msgraph->m_cps[*it]->cellid))
          continue;

        crit_conn_idxs[dim-1].push_back
            (glutils::line_idx_t(i,*it));
      }
    }

    for(uint i = 0 ; i < gc_grid_dim; ++i)
    {
      ren_cp_conns[i] = glutils::create_buffered_lines_ren
                        (cp_loc_bo,
                         glutils::make_buf_obj(crit_conn_idxs[i]),
                         glutils::make_buf_obj());
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

    for(uint i = 0 ; i < gc_grid_dim; ++i)

    {
      ren_grad[i] = glutils::create_buffered_lines_ren
                    (cell_bo,
                     glutils::make_buf_obj(pair_idxs[i]),
                     glutils::make_buf_obj());
    }
  }

  void octtree_piece_rendata::create_surf_ren(const rect_t & roi)
  {
    if(dp->dataset == NULL)
      return;

#warning "create surf ren not implemented"

  }

  void octtree_piece_rendata::create_disc_rds()
  {
    if(dp->msgraph == NULL)
      return;

    boost::shared_ptr<disc_rendata_t> sptr;

    for(uint i = 0 ; i < dp->msgraph->m_cps.size();++i)
    {
      if(dp->msgraph->m_cps[i]->is_paired) continue;

      sptr.reset(new disc_rendata_t(dp->msgraph->m_cps[i]->cellid));
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

    if ( m_bShowSurface && ren_surf)
    {
      glColor3f ( 0.75,0.75,0.75 );
      ren_surf->render();
    }

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

    glPointSize ( 4.0 );

    if ( m_bShowCps)
    {
      for(uint i = 0 ; i < gc_grid_dim+1;++i)
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

  int octtree_piece_rendata::rows()
  {
    return disc_rds.size();
  }
  int octtree_piece_rendata::columns()
  {
    return 5;
  }
  bool octtree_piece_rendata::exchange_data(const data_index_t &idx,
                                            boost::any &v,
                                            const eExchangeMode &m)
  {
    switch(idx[0])
    {
    case 0:return s_exchange_read_only(disc_rds[idx[1]]->cellid.to_string(),v,m);
    case 1:return s_exchange_read_write(disc_rds[idx[1]]->m_bShowAsc,v,m);
    case 2:return s_exchange_read_write(disc_rds[idx[1]]->asc_color,v,m);
    case 3:return s_exchange_read_write(disc_rds[idx[1]]->m_bShowDes,v,m);
    case 4:return s_exchange_read_write(disc_rds[idx[1]]->des_color,v,m);
    };
    throw std::logic_error("invalid index");
  }

  std::string octtree_piece_rendata::get_header(int i)
  {
    switch(i)
    {
    case 0: return "cellid";
    case 1: return "asc disc";
    case 2: return "asc disc color";
    case 3: return "des disc";
    case 4: return "des disc color";
    }
    throw std::logic_error("invalid index");
  }

  disc_rendata_t::disc_rendata_t(cellid_t c):
      m_bShowAsc(false),
      m_bShowDes(false),
      asc_ren(NULL),
      des_ren(NULL),
      cellid(c){}

  disc_rendata_t::~disc_rendata_t()
  {
    m_bShowAsc = false;
    m_bShowDes = false;

    update(NULL);

  }

  void disc_rendata_t::render()
  {
    uint dim = dataset_t::s_getCellDim(cellid);

    if(m_bShowAsc || m_bShowDes)
    {
//      glColor3dv(g_grid_cp_colors[dim].data());
      glColor3dv(asc_color.data());

      glBegin(GL_POINTS);
      glVertex3sv(cellid.data());
      glEnd();
    }

    bool *bShow[] = {&m_bShowAsc,&m_bShowDes};
    glutils::renderable_t *ren[] = {asc_ren,des_ren};

    for(uint i = 0 ;i <2 ;++i )
    {
      if(*bShow[i] && ren[i])
      {
        ren[i]->render();
      }
    }
  }
  bool disc_rendata_t::update(mscomplex_t *msc )
  {
    uint ret = false;

    bool *bShow[] = {&m_bShowAsc,&m_bShowDes};
    glutils::renderable_t **ren[] = {&asc_ren,&des_ren};

    for(uint i = 0 ;i <2 ;++i )
    {
      if(bShow[i] && *ren[i] == NULL && msc != NULL)
      {
        if(msc->m_id_cp_map.count(cellid) == 0)
          throw std::logic_error("cannot find cellid in id_cp_map while updating ren");

        uint cp_idx = msc->m_id_cp_map[cellid];

        std::vector<glutils::vertex_t> vlist;

        *ren[i] = glutils::create_buffered_points_ren
                  (glutils::make_buf_obj(vlist),
                   glutils::make_buf_obj(),
                   glutils::make_buf_obj());

        ret = true;
      }

      if(!bShow[i] && *ren[i] != NULL)
      {
        delete *ren[i];
        *ren[i] = NULL;
        ret = true;
      }
    }

    return ret;
  }

}
