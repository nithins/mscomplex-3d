#ifndef GRID_VIEWER_MAINWINDOW_INCLUDED
#define GRID_VIEWER_MAINWINDOW_INCLUDED

#include <grid.h>

#include <QDialog>
#include <QAbstractItemModel>
#include <QModelIndex>
#include <QVariant>
#include <QItemSelectionModel>
#include <QTreeView>
#include <QGLViewer/qglviewer.h>

#include <ui_grid_viewer_mainwindow.h>
#include <boost/any.hpp>

class configurable_t;

namespace grid
{


  class grid_viewer_t;

  class data_manager_t;

  class glviewer_t : public QGLViewer
  {

  public:

    grid_viewer_t *m_ren;

    glviewer_t(data_manager_t * p ,const rect_t &r);

    ~glviewer_t();

  protected:

    virtual void draw();
    virtual void init();
    virtual QString helpString() const;
  };

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

  void configurable_ctx_menu
      (configurable_t *c,
       const QModelIndexList & l,
       const QPoint &p);

  class configurable_ctx_menu_sig_collector:public QObject
  {
    Q_OBJECT

  public:

    boost::any              m_val;
    int                     m_col;
    configurable_t *        m_conf;
    const QModelIndexList & m_rows;

    configurable_ctx_menu_sig_collector
        (configurable_t * conf,
         const boost::any & val,
         const int & col,
         const QModelIndexList & rows,
         QObject *par):
        m_conf(conf),
        m_val(val),
        m_col(col),
        m_rows(rows)
    {setParent(par);}

  private slots:
    void triggered(bool state);
  };

  class configurable_item_model : public QAbstractListModel
  {
    Q_OBJECT

  public:

    configurable_item_model ( configurable_t *conf,QObject *parent = 0 ):
        QAbstractListModel ( parent ),m_conf(conf){}

    QVariant data ( const QModelIndex &index, int role ) const;

    QVariant headerData ( int section, Qt::Orientation orientation,
                          int role = Qt::DisplayRole ) const;

    int rowCount ( const QModelIndex &parent = QModelIndex() ) const;

    void reset_configurable(configurable_t *conf);

  private:

    configurable_t * m_conf;

  };
}



#endif
