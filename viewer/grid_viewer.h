#ifndef GRID_VIEWER_H_INCLUDED
#define GRID_VIEWER_H_INCLUDED

#include <QGLViewer/qglviewer.h>
#include <grid_datamanager.h>

typedef unsigned char uchar;
typedef unsigned int  uint;

namespace glutils
{
  class renderable_t;
}

namespace grid
{
  class octtree_piece ;

  class octtree_piece_rendata
  {
  public:

    octtree_piece * dp;

    // set externally to control what is rendered
    bool m_bShowSurface;
    bool m_bShowCps;
    bool m_bShowCpLabels;
    bool m_bShowMsGraph;
    bool m_bShowGrad;
    bool m_bShowCancCps;
    bool m_bShowCancMsGraph;

    glutils::renderable_t  *ren_surf;
    glutils::renderable_t  *ren_grad[gc_grid_dim];
    glutils::renderable_t  *ren_cp_labels[gc_grid_dim+1];
    glutils::renderable_t  *ren_cp[gc_grid_dim+1];
    glutils::renderable_t  *ren_cp_conns[gc_grid_dim];
    glutils::renderable_t  *ren_canc_cp_labels[gc_grid_dim+1];
    glutils::renderable_t  *ren_canc_cp[gc_grid_dim+1];
    glutils::renderable_t  *ren_canc_cp_conns[gc_grid_dim];

    void create_cp_rens();
    void create_grad_rens();
    void create_surf_ren();
    void render() const ;

    octtree_piece_rendata(octtree_piece *);
  };

  class glviewer_t : public QGLViewer
  {
  public:

    std::vector<octtree_piece_rendata * >  m_grid_piece_rens;

    cellid_t                               m_size;


  public:
    glviewer_t(std::vector<octtree_piece *> * p ,
               cellid_t size);
    ~glviewer_t();

  protected :

      virtual void draw();
  virtual void init();
  virtual QString helpString() const;
};
}
#endif //VIEWER_H_INCLUDED
