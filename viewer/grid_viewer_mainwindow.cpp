#include <QMenu>
#include <QTreeView>


#include <grid_viewer.h>
#include <grid_viewer_mainwindow.h>
#include <grid_datamanager.h>
#include <grid_mscomplex.h>

namespace grid
{
  bool &bool_menuaction_ref
      (const viewer_mainwindow::eBoolMenuAction &action ,
       octtree_piece_rendata *gp_rd )
  {
    switch (action)
    {
    case viewer_mainwindow::VA_SURF: return gp_rd->m_bShowSurface;
    case viewer_mainwindow::VA_CPS: return gp_rd->m_bShowCps;
    case viewer_mainwindow::VA_CPLABELS: return gp_rd->m_bShowCpLabels;
    case viewer_mainwindow::VA_GRAPH:return gp_rd->m_bShowMsGraph;
    case viewer_mainwindow::VA_GRAD:return gp_rd->m_bShowGrad;
    case viewer_mainwindow::VA_CANC_CPS:return gp_rd->m_bShowCancCps;
    case viewer_mainwindow::VA_CANC_GRAPH:return gp_rd->m_bShowCancMsGraph;
    }

    throw std::invalid_argument("undefined tva action");

    return gp_rd->m_bShowSurface;
  }

  bool get_bool_menuaction_value
      ( viewer_mainwindow * mw,
        const viewer_mainwindow::eBoolMenuAction &action )
  {

    uint num_checked_items = 0;
    uint num_unchecked_items = 0;

    QModelIndexList indexes
        = mw->datapiece_view->selectionModel()->selectedIndexes();

    for ( QModelIndexList::iterator ind_it = indexes.begin();
    ind_it != indexes.end(); ++ind_it )
    {
      octtree_piece_rendata * gp_rd = mw->m_viewer->m_grid_piece_rens[(*ind_it).row()];

      if ( bool_menuaction_ref(action,gp_rd) )
        ++num_checked_items;
      else
        ++num_unchecked_items;
    }

    if ( num_checked_items > num_unchecked_items )
      return true;
    else
      return false;
  }

  void update_bool_menuaction_value
      ( viewer_mainwindow * mw,const viewer_mainwindow::eBoolMenuAction &action,
        const bool & state )
  {

    QModelIndexList indexes = mw->datapiece_view->selectionModel()->selectedIndexes();

    bool need_update = false;

    for ( QModelIndexList::iterator ind_it = indexes.begin();
    ind_it != indexes.end(); ++ind_it )
    {
      octtree_piece_rendata * gp_rd = mw->m_viewer->m_grid_piece_rens[(*ind_it).row()];

      if(state != bool_menuaction_ref(action,gp_rd))
      {
        bool_menuaction_ref(action,gp_rd) = state;
        need_update = true;
      }
    }

    if(need_update )
    {
      mw->m_viewer->updateGL();
    }
  }

  void toggled_signal_retransmitter::toggled(bool state)
  {
    update_bool_menuaction_value(m_pMw,m_act,state);
  }

  QAction * add_bool_menu_action(QMenu *m,viewer_mainwindow *mw,
                                 const QString &str,viewer_mainwindow::eBoolMenuAction act)
  {
    QAction * action  = m->addAction ( str );
    action->setCheckable ( true );
    action->setChecked ( get_bool_menuaction_value(mw,act));

    toggled_signal_retransmitter * re_trns =
        new toggled_signal_retransmitter(mw,act,m);

    re_trns->connect(action,SIGNAL ( toggled ( bool ) ),re_trns,SLOT(toggled ( bool )));

    return action;
  }

  void viewer_mainwindow::on_datapiece_view_customContextMenuRequested  ( const QPoint &pos )
  {
    QMenu m;

    QAction * a = NULL;

    a = add_bool_menu_action(&m,this,"show surface",VA_SURF);
    a = add_bool_menu_action(&m,this,"show cps ",VA_CPS);
    a = add_bool_menu_action(&m,this,"show cp labels",VA_CPLABELS);
    a = add_bool_menu_action(&m,this,"show graph",VA_GRAPH);
    a = add_bool_menu_action(&m,this,"show grad",VA_GRAD);

    m.exec ( datapiece_view->mapToGlobal ( pos ) );
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

    mscomplex_t * msc = m_mw->m_viewer->m_grid_piece_rens[active_otp]->dp->msgraph;

    std::string s = msc->m_cps[index.row()]->cellid.to_string();

    return QString(s.c_str());
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

    mscomplex_t * msc = m_mw->m_viewer->m_grid_piece_rens[active_otp]->dp->msgraph;

    return msc->m_cps.size();
  }

}
