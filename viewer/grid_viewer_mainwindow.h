#ifndef GRID_VIEWER_MAINWINDOW_INCLUDED
#define GRID_VIEWER_MAINWINDOW_INCLUDED

#include <QDialog>
#include <QAbstractItemModel>
#include <QModelIndex>
#include <QVariant>
#include <QItemSelectionModel>
#include <QTreeView>

#include <ui_grid_viewer_mainwindow.h>
#include <grid.h>
#include <iostream>

namespace grid
{

  class glviewer_t;

  class data_manager_t;

  class octtree_piece_rendata;

  class viewer_mainwindow:
      public QDialog,
      public Ui::grid_viewer_mainwindow_Dialog
  {

  public:

    glviewer_t      *m_viewer;
    data_manager_t  *m_gdm;

    viewer_mainwindow
        (data_manager_t *gdm,const rect_t & roi);

    ~viewer_mainwindow();

    Q_OBJECT

  public:

    enum eBoolMenuAction
    {
      VA_SURF,
      VA_CPS,
      VA_CPLABELS,
      VA_GRAPH,
      VA_GRAD,
      VA_CANC_CPS,
      VA_CANC_GRAPH,
    };

  private slots:
    void on_datapiece_view_customContextMenuRequested ( const QPoint &p );

  };

  class toggled_signal_retransmitter:public QObject
  {
    Q_OBJECT

  public:

    viewer_mainwindow *m_pMw;
    viewer_mainwindow::eBoolMenuAction m_act;

    toggled_signal_retransmitter
        (viewer_mainwindow *pMw,viewer_mainwindow::eBoolMenuAction act,
         QObject *par):m_pMw(pMw),m_act(act)
    {setParent(par);}

  private slots:
    void toggled(bool state);

  };
}

class octtree_piece_item_model : public QAbstractItemModel
{
  Q_OBJECT

public:

  octtree_piece_item_model ( std::vector<grid::octtree_piece_rendata *> *, QObject *parent = 0 );
  ~octtree_piece_item_model();

  QVariant data ( const QModelIndex &index, int role ) const;

  Qt::ItemFlags flags ( const QModelIndex &index ) const;

  QVariant headerData ( int section, Qt::Orientation orientation,
                        int role = Qt::DisplayRole ) const;

  QModelIndex index ( int row, int column,
                      const QModelIndex &parent = QModelIndex() ) const;

  QModelIndex parent ( const QModelIndex &index ) const;

  int rowCount ( const QModelIndex &parent = QModelIndex() ) const;

  int columnCount ( const QModelIndex &parent = QModelIndex() ) const;

  struct tree_item
  {
    std::vector<tree_item *>               children;
    grid::octtree_piece_rendata          * node;
    tree_item                            * parent;

    tree_item ( grid::octtree_piece_rendata * _node , tree_item * par );
    tree_item();
    int row();
  };


private:
  void setupModelData
      ( std::vector<grid::octtree_piece_rendata *> *);

  tree_item *m_tree;
};

#endif
