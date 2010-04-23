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

    glviewer_t              *m_viewer;
    data_manager_t          *m_gdm;

    uint                     m_active_otp_idx;

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

    void on_datapiece_view_activated ( const QModelIndex & index  );
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

  class octtree_piece_item_model : public QAbstractListModel
  {
    Q_OBJECT

  public:

    octtree_piece_item_model ( viewer_mainwindow *mw,QObject *parent = 0 ):
        QAbstractListModel ( parent ),m_mw(mw){}

    QVariant data ( const QModelIndex &index, int role ) const;

    QVariant headerData ( int section, Qt::Orientation orientation,
                          int role = Qt::DisplayRole ) const;

    int rowCount ( const QModelIndex &parent = QModelIndex() ) const;

  private:

    viewer_mainwindow * m_mw;

  };

  class critpt_item_model : public QAbstractListModel
  {
    Q_OBJECT

  public:

    critpt_item_model ( viewer_mainwindow *mw,QObject *parent = 0 ):
        QAbstractListModel ( parent ),m_mw(mw){}

    QVariant data ( const QModelIndex &index, int role ) const;

    QVariant headerData ( int section, Qt::Orientation orientation,
                          int role = Qt::DisplayRole ) const;

    int rowCount ( const QModelIndex &parent = QModelIndex() ) const;

    void active_otp_changed(){reset();}

  private:

    viewer_mainwindow * m_mw;

  };
}



#endif
