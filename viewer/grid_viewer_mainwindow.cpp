#include <QMenu>
#include <QTreeView>


#include <grid_viewer.h>
#include <grid_viewer_mainwindow.h>
#include <grid_datamanager.h>

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
      octtree_piece_item_model::tree_item *item =
          static_cast<octtree_piece_item_model::tree_item*> ( ( *ind_it ).internalPointer());

      if ( bool_menuaction_ref(action,item->node) )
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
      octtree_piece_item_model::tree_item *item =
          static_cast<octtree_piece_item_model::tree_item*> ( ( *ind_it ).internalPointer());

      if(state != bool_menuaction_ref(action,item->node))
      {
        bool_menuaction_ref(action,item->node) = state;
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



  viewer_mainwindow::viewer_mainwindow
      (data_manager_t * gdm,const rect_t & roi):m_gdm(gdm)
  {
    setupUi (this);

    m_viewer = new glviewer_t(&gdm->m_pieces,gdm->m_size,roi);

    m_viewer->setParent(glviewer);

    m_viewer->resize(glviewer->size());

    octtree_piece_item_model *model = new octtree_piece_item_model ( &m_viewer->m_grid_piece_rens);

    datapiece_view->setModel ( model );
  }

  viewer_mainwindow::~viewer_mainwindow()
  {
    delete m_gdm;
  }
}


octtree_piece_item_model::octtree_piece_item_model ( std::vector<grid::octtree_piece_rendata *> * dpList, QObject *parent )
  : QAbstractItemModel ( parent )
{
  setupModelData ( dpList );
}

octtree_piece_item_model::~octtree_piece_item_model()
{
  delete m_tree;
  m_tree = NULL;
}


int octtree_piece_item_model::columnCount ( const QModelIndex &/*parent*/ ) const
{
  return 1;
}

QVariant octtree_piece_item_model::data ( const QModelIndex &index, int role ) const
{
  if ( !index.isValid() )
    return QVariant();

  if ( role != Qt::DisplayRole )
    return QVariant();

  tree_item *item = static_cast<tree_item*> ( index.internalPointer() );

  return QString(item->node->dp->label().c_str());
}

Qt::ItemFlags octtree_piece_item_model::flags ( const QModelIndex &index ) const
{
  if ( !index.isValid() )
    return 0;

  return Qt::ItemIsEnabled | Qt::ItemIsSelectable;
}

QVariant octtree_piece_item_model::headerData ( int /*section*/, Qt::Orientation orientation,
                                     int role ) const
{
  if ( orientation == Qt::Horizontal && role == Qt::DisplayRole )
    return "Data Pieces";

  return QVariant();
}

QModelIndex octtree_piece_item_model::index ( int row, int column, const QModelIndex &parent ) const
{
  if ( !hasIndex ( row, column, parent ) )
    return QModelIndex();

  tree_item *parentItem;

  if ( !parent.isValid() )
    parentItem = m_tree;
  else
    parentItem = static_cast<tree_item*> ( parent.internalPointer() );

  if ( row < ( int ) parentItem->children.size() )
    return createIndex ( row, column, parentItem->children[row] );
  else
    return QModelIndex();

}

QModelIndex octtree_piece_item_model::parent ( const QModelIndex &index ) const
{
  if ( !index.isValid() )
    return QModelIndex();

  tree_item *childItem  = static_cast<tree_item*> ( index.internalPointer() );

  tree_item *parentItem = childItem->parent;

  if ( parentItem == m_tree )
    return QModelIndex();

  return createIndex ( parentItem->row(), 0, parentItem );

}

int octtree_piece_item_model::rowCount ( const QModelIndex &parent ) const
{
  tree_item *parentItem;

  if ( parent.column() > 0 )
    return 0;

  if ( !parent.isValid() )
    parentItem = m_tree;
  else
    parentItem = static_cast<tree_item*> ( parent.internalPointer() );

  return parentItem->children.size();

}

void octtree_piece_item_model::setupModelData( std::vector<grid::octtree_piece_rendata *> * dpList)
{
  m_tree = new tree_item();


  for ( std::vector<grid::octtree_piece_rendata *>::iterator dp_it =  dpList->begin();
  dp_it != dpList->end(); ++dp_it )
  {
    grid::octtree_piece_rendata *dp = *dp_it;

    tree_item *parentItem = m_tree;

    tree_item * dpItem = new tree_item ( dp, parentItem );

    parentItem->children.push_back ( dpItem );

  }
}

octtree_piece_item_model::tree_item::tree_item
    ( grid::octtree_piece_rendata * _node , tree_item * par ) :
    node ( _node ),parent ( par )
{
}

octtree_piece_item_model::tree_item::tree_item()
{
  node = NULL;
  parent = NULL;
}

int octtree_piece_item_model::tree_item::row()
{
  return std::find ( parent->children.begin(),parent->children.end(),this )
      - parent->children.begin();
}
