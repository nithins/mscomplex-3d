#include <QMenu>
#include <QTreeView>
#include <QColorDialog>
#include <QDebug>

#include <grid_viewer.h>
#include <grid_viewer_mainwindow.h>
#include <grid_datamanager.h>

namespace grid
{
  void viewer_mainwindow::on_datapiece_view_customContextMenuRequested  ( const QPoint &p )
  {
    QModelIndexList l =  datapiece_view->selectionModel()->selectedIndexes();

    std::vector<configureable_t *> c_list;

    for(uint i = 0 ; i < l.size();++i)
      c_list.push_back(m_viewer->m_grid_piece_rens[l[i].row()]);

    configure_ctx_menu(c_list,datapiece_view->mapToGlobal(p));

    m_viewer->updateGL();

  }

  void viewer_mainwindow::on_critpt_view_customContextMenuRequested ( const QPoint &p )
  {
    QModelIndexList l =  critpt_view->selectionModel()->selectedIndexes();

    std::vector<configureable_t *> c_list;

    for(uint i = 0 ; i < l.size();++i)
      c_list.push_back(m_viewer->m_grid_piece_rens[m_active_otp_idx]->disc_rds[l[i].row()].get());

    configure_ctx_menu(c_list,critpt_view->mapToGlobal(p));

    m_viewer->updateGL();
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

  void configure_ctx_menu(const std::vector<configureable_t *> &l, const QPoint &p)
  {
    if(l.size() == 0)
      return;

    configureable_t * c = l[0];

    QMenu m;

    for(uint i = 0 ; i < c->get_num_items();++i)
    {
      QAction * action  = m.addAction ( c->get_description(i).c_str());

      boost::any val;

      c->update_item(i,val,false);

      if(val.type() == typeid(bool))
      {
        action->setCheckable(true);
        action->setChecked(boost::any_cast<bool>(val));
      }

      configure_ctx_menu_sig_collector * coll =
          new configure_ctx_menu_sig_collector(l,val,i,&m);

      coll->connect(action,SIGNAL ( triggered ( bool ) ),coll,SLOT(triggered ( bool )));
    }

    m.exec(p);
  }

  void configure_ctx_menu_sig_collector::triggered(bool state)
  {
    boost::any out_val;

    if(m_val.type() == typeid(bool))
      out_val = boost::any(state);

    if(m_val.type() == typeid(glutils::color_t))
    {

      glutils::color_t c = boost::any_cast<glutils::color_t>(m_val);

      QColor ic = QColor::fromRgbF(c[0],c[1],c[2],1.0);

      QColor qc = QColorDialog::getColor(ic);

      if(qc.isValid())
        out_val = glutils::color_t(qc.redF(),qc.greenF(),qc.blueF());
    }

    if(out_val.empty())
      return;

    for(uint i = 0 ; i < m_list.size();++i)
    {
      configureable_t * c = m_list[i];
      c->update_item(m_i,out_val,true);
    }
  }

}
