#include <QMenu>
#include <QTreeView>
#include <QColorDialog>
#include <QDebug>

#include <grid_viewer.h>
#include <grid_viewer_mainwindow.h>
#include <grid_datamanager.h>

namespace grid
{
  QColor glutils_to_q_color(const glutils::color_t &c)
  {
    return  QColor::fromRgbF(c[0],c[1],c[2],1.0);
  }

  glutils::color_t q_to_glutils_color(const QColor &qc)
  {
    return  glutils::color_t(qc.redF(),qc.greenF(),qc.blueF());
  }

  typedef QVariant (*get_val_for_action_fnt )
      (const uint &,          // action
       viewer_mainwindow * ,
       const uint &);         // index

  typedef void (*set_val_for_action_fnt )
      (const uint &,          // action
       viewer_mainwindow * ,
       const uint &,          // index
       const QVariant &);     // value

  struct ctx_menu_data_t
  {
    const char * text;
    get_val_for_action_fnt get_fn;
    set_val_for_action_fnt set_fn;
  };


  enum eContextMenuAction
  {
    RD_MA_BEGIN               = 0,
    OTPRD_MA_BEGIN            = RD_MA_BEGIN,

    OTPRD_MA_SURF             = OTPRD_MA_BEGIN,
    OTPRD_MA_CPS,
    OTPRD_MA_CPLABELS,
    OTPRD_MA_GRAPH,
    OTPRD_MA_GRAD,

    OTPRD_MA_END,
    CPDRD_MA_BEGIN            = OTPRD_MA_END,

    CPDRD_MA_SHOW_ASC_DISC    = CPDRD_MA_BEGIN,
    CPDRD_MA_ASC_DISC_COLOR,
    CPDRD_MA_SHOW_DES_DISC,
    CPDRD_MA_DES_DISC_COLOR,

    CPDRD_MA_END,
    RD_MA_END                 = CPDRD_MA_END,
  };

  QVariant get_val_for_otprd_ma
      (const uint &action,
       viewer_mainwindow * mw ,
       const uint &idx)
  {
    octtree_piece_rendata * otp_rd = NULL;

    otp_rd = mw->m_viewer->m_grid_piece_rens[idx];

    switch (action)
    {
    case OTPRD_MA_SURF: return otp_rd->m_bShowSurface;
    case OTPRD_MA_CPS: return otp_rd->m_bShowCps;
    case OTPRD_MA_CPLABELS: return otp_rd->m_bShowCpLabels;
    case OTPRD_MA_GRAPH:return otp_rd->m_bShowMsGraph;
    case OTPRD_MA_GRAD:return otp_rd->m_bShowGrad;
    default:throw std::logic_error("cannot deal with this action");
    }
    return QVariant();
  }

  void set_val_for_otprd_ma
      (const uint &action,
       viewer_mainwindow * mw ,
       const uint &idx,
       const QVariant &val)
  {
    octtree_piece_rendata * otp_rd = mw->m_viewer->m_grid_piece_rens[idx];

    switch (action)
    {
    case OTPRD_MA_SURF: otp_rd->m_bShowSurface = val.value<bool>();break;
    case OTPRD_MA_CPS: otp_rd->m_bShowCps = val.value<bool>();break;
    case OTPRD_MA_CPLABELS: otp_rd->m_bShowCpLabels = val.value<bool>();break;
    case OTPRD_MA_GRAPH: otp_rd->m_bShowMsGraph = val.value<bool>();break;
    case OTPRD_MA_GRAD: otp_rd->m_bShowGrad = val.value<bool>();break;
    default:throw std::logic_error("cannot deal with this action");
    }
  }

  QVariant get_val_for_cpdrd_ma
      (const uint &action,
       viewer_mainwindow * mw ,
       const uint &idx)
  {
    octtree_piece_rendata * otp_rd =
        mw->m_viewer->m_grid_piece_rens[mw->m_active_otp_idx];

    boost::shared_ptr<disc_rendata_t> cpd_rd = otp_rd->disc_rds[idx];

    switch (action)
    {
    case CPDRD_MA_ASC_DISC_COLOR:return glutils_to_q_color(cpd_rd->asc_color);
    case CPDRD_MA_SHOW_ASC_DISC:return cpd_rd->m_bShowAsc;
    case CPDRD_MA_DES_DISC_COLOR:return glutils_to_q_color(cpd_rd->des_color);
    case CPDRD_MA_SHOW_DES_DISC:return cpd_rd->m_bShowDes;
    default:throw std::logic_error("cannot deal with this action");
    }
    return QVariant();
  }

  void set_val_for_cpdrd_ma
      (const uint &action,
       viewer_mainwindow * mw ,
       const uint &idx,
       const QVariant &val)
  {
    octtree_piece_rendata * otp_rd =
        mw->m_viewer->m_grid_piece_rens[mw->m_active_otp_idx];

    boost::shared_ptr<disc_rendata_t> cpd_rd = otp_rd->disc_rds[idx];

    switch (action)
    {
    case CPDRD_MA_ASC_DISC_COLOR:cpd_rd->asc_color  = q_to_glutils_color(val.value<QColor>());break;
    case CPDRD_MA_SHOW_ASC_DISC:cpd_rd->m_bShowAsc  = val.value<bool>();break;
    case CPDRD_MA_DES_DISC_COLOR:cpd_rd->des_color  = q_to_glutils_color(val.value<QColor>());break;
    case CPDRD_MA_SHOW_DES_DISC:cpd_rd->m_bShowDes  = val.value<bool>();break;
    default:
      throw std::logic_error("cannot deal with this action");
    }
  }

  ctx_menu_data_t s_ctx_menu_data [RD_MA_END - RD_MA_BEGIN] =
  {
    {"show surface",get_val_for_otprd_ma,set_val_for_otprd_ma},
    {"show critical points",get_val_for_otprd_ma,set_val_for_otprd_ma},
    {"show critical point labels",get_val_for_otprd_ma,set_val_for_otprd_ma},
    {"show graph",get_val_for_otprd_ma,set_val_for_otprd_ma},
    {"show grad",get_val_for_otprd_ma,set_val_for_otprd_ma},

    {"show asc disc",get_val_for_cpdrd_ma,set_val_for_cpdrd_ma},
    {"set asc disc color",get_val_for_cpdrd_ma,set_val_for_cpdrd_ma},
    {"show des disc",get_val_for_cpdrd_ma,set_val_for_cpdrd_ma},
    {"set des disc color",get_val_for_cpdrd_ma,set_val_for_cpdrd_ma},
  };

  QAbstractItemView * get_action_view(uint action,viewer_mainwindow *mw)
  {
    if(OTPRD_MA_BEGIN <= action && action < OTPRD_MA_END)
      return mw->datapiece_view;

    if(CPDRD_MA_BEGIN <= action && action < CPDRD_MA_END)
      return mw->critpt_view;

    throw std::logic_error("unknown action type");

    return NULL;
  }

  void toggled_signal_retransmitter::toggled(bool state)
  {
    QModelIndexList l
        = get_action_view(m_act,m_pMw)->selectionModel()->selectedIndexes();

    QVariant updated_value;

    if(m_val.canConvert<bool>())
      updated_value = state;

    if(m_val.canConvert<QColor>())
    {
      QColor ic = m_val.value<QColor>();

      if (ic.isValid() == false)  ic = Qt::white;

      QColor c = QColorDialog::getColor(ic,m_pMw);

      if(c.isValid())
        updated_value=c;
    }

    if(updated_value.isValid())
    {
      for(QModelIndexList::iterator it = l.begin(); it != l.end();++it)
      {
        (*s_ctx_menu_data[m_act].set_fn)(m_act,m_pMw,(*it).row(),updated_value);
      }
    }

    m_pMw->m_viewer->updateGL();
  }


  QAction * add_menu_action(QMenu *m,viewer_mainwindow *mw,const uint & act)
  {
    QAction * action  = m->addAction ( s_ctx_menu_data[act].text);

    QModelIndex i = get_action_view(act,mw)->selectionModel()->currentIndex();

    QVariant cur_val = (*s_ctx_menu_data[act].get_fn)(act,mw,i.row());

    if(cur_val.canConvert<bool>())
    {
      action->setCheckable(true);
      action->setChecked(cur_val.value<bool>());
    }

    toggled_signal_retransmitter * re_trns =
        new toggled_signal_retransmitter(mw,act,cur_val,m);

    re_trns->connect(action,SIGNAL ( triggered ( bool ) ),re_trns,SLOT(toggled ( bool )));

    return action;
  }

  void viewer_mainwindow::on_datapiece_view_customContextMenuRequested  ( const QPoint &pos )
  {
    QMenu m;

    for(uint i = OTPRD_MA_BEGIN; i < OTPRD_MA_END; ++i)
      add_menu_action(&m,this,i);

    m.exec ( datapiece_view->mapToGlobal ( pos ) );
  }

  void viewer_mainwindow::on_critpt_view_customContextMenuRequested ( const QPoint &p )
  {
    QMenu m;

    for(uint i = CPDRD_MA_BEGIN; i < CPDRD_MA_END; ++i)
      add_menu_action(&m,this,i);

    m.exec ( critpt_view->mapToGlobal ( p ) );
  }

  void viewer_mainwindow::on_datapiece_view_activated ( const QModelIndex & index  )
  {
    if(m_active_otp_idx == index.row())
      return;

    m_active_otp_idx = index.row();

    critpt_item_model * cp_model =
        dynamic_cast<critpt_item_model *>(critpt_view->model());

    if(cp_model)
      cp_model->active_otp_changed();
  }

  viewer_mainwindow::viewer_mainwindow
      (data_manager_t * gdm,const rect_t & roi):m_gdm(gdm),
      m_active_otp_idx(NULL)
  {
    setupUi (this);

    m_viewer = new glviewer_t(&gdm->m_pieces,gdm->m_size,roi);

    m_viewer->setParent(glviewer);

    m_viewer->resize(glviewer->size());

    octtree_piece_item_model *otp_model = new octtree_piece_item_model ( this);

    datapiece_view->setModel ( otp_model );

    critpt_item_model *cp_model = new critpt_item_model ( this);

    critpt_view->setModel ( cp_model );
  }

  viewer_mainwindow::~viewer_mainwindow()
  {
    delete m_gdm;
  }

  QVariant octtree_piece_item_model::data ( const QModelIndex &index, int role ) const
  {
    if ( !index.isValid() )
      return QVariant();

    if ( role != Qt::DisplayRole )
      return QVariant();

    std::string s = m_mw->m_viewer->m_grid_piece_rens[index.row()]->dp->label();

    return QString(s.c_str());
  }

  QVariant octtree_piece_item_model::headerData
      ( int /*section*/, Qt::Orientation orientation,int role ) const
  {
    if ( orientation == Qt::Horizontal && role == Qt::DisplayRole )
      return "Data Pieces";

    return QVariant();
  }

  int octtree_piece_item_model::rowCount ( const QModelIndex &parent ) const
  {
    return m_mw->m_viewer->m_grid_piece_rens.size();
  }

  QVariant critpt_item_model::data ( const QModelIndex &index, int role ) const
  {
    if ( !index.isValid() )
      return QVariant();

    if ( role != Qt::DisplayRole )
      return QVariant();

    int active_otp = m_mw->m_active_otp_idx;

    cellid_t c = m_mw->m_viewer->m_grid_piece_rens[active_otp]->disc_rds[index.row()]->cellid;

    return QString(c.to_string().c_str());
  }

  QVariant critpt_item_model::headerData
      ( int /*section*/, Qt::Orientation orientation,int role ) const
  {
    if ( orientation == Qt::Horizontal && role == Qt::DisplayRole )
      return "Cript of piece no:";

    return QVariant();
  }

  int critpt_item_model::rowCount ( const QModelIndex &parent ) const
  {
    int active_otp = m_mw->m_active_otp_idx;

    return m_mw->m_viewer->m_grid_piece_rens[active_otp]->disc_rds.size();
  }

}
