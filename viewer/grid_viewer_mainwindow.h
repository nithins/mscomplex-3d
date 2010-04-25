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
#include <boost/any.hpp>

namespace grid
{

  class glviewer_t;

  class data_manager_t;

  class octtree_piece_rendata;

  class viewer_mainwindow:
      public QDialog,
      public Ui::grid_viewer_mainwindow_Dialog
  {

  Q_OBJECT

  public:

    glviewer_t              *m_viewer;
    data_manager_t          *m_gdm;

    uint                     m_active_otp_idx;

    viewer_mainwindow
        (data_manager_t *gdm,const rect_t & roi);

    ~viewer_mainwindow();

  private slots:
    void on_datapiece_view_customContextMenuRequested ( const QPoint &p );

    void on_critpt_view_customContextMenuRequested ( const QPoint &p );

    void on_datapiece_view_activated ( const QModelIndex & index  );
  };

  class configureable_t;

  static void configure_ctx_menu(const std::vector<configureable_t *> &,const QPoint &p );

  class configure_ctx_menu_sig_collector:public QObject
  {
    Q_OBJECT

  public:

    boost::any m_val;

    int m_i;

    const std::vector<configureable_t *> & m_list;

    configure_ctx_menu_sig_collector
        (const std::vector<configureable_t *> & l,
         const boost::any & val,
         const int & i,
         QObject *par):m_list(l),m_val(val),m_i(i)
    {setParent(par);}

  private slots:
    void triggered(bool state);
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
