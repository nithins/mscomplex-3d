#ifndef GRID_VIEWER_H_INCLUDED
#define GRID_VIEWER_H_INCLUDED

#include <QGLViewer/qglviewer.h>
#include <grid.h>

#include <glutils.h>
#include <set>

#include <boost/any.hpp>

namespace grid
{
  class octtree_piece ;

  class mscomplex_t;

  class configureable_t
  {
  public:
    virtual int         get_num_items() = 0 ;
    virtual bool        update_item(const int &,boost::any &,const bool &) = 0;
    virtual std::string get_description(int i) = 0;


    template <typename T>
        static bool s_update_item(T &p_val,boost::any &c_val,const bool &update_self)
    {
      bool ret = false;

      if(update_self)
      {
        ret   = (p_val != boost::any_cast<T>(c_val));
        p_val = boost::any_cast<T>(c_val);
      }
      else
        c_val = boost::any(p_val);

      return ret;
    }
  };

  class disc_rendata_t:public configureable_t
  {

  public:

    glutils::renderable_t *asc_ren;
    glutils::renderable_t *des_ren;
    cellid_t               cellid;

    bool                   m_bShowAsc;
    bool                   m_bShowDes;

    glutils::color_t       asc_color;
    glutils::color_t       des_color;


    disc_rendata_t(cellid_t c);
    ~disc_rendata_t();

    void render();
    bool update(mscomplex_t *);

    virtual int         get_num_items();
    virtual bool        update_item(const int & ,boost::any &,const bool &);
    virtual std::string get_description(int i);
  };

  class octtree_piece_rendata:public configureable_t
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

    // set externally .. cleared by render
    bool m_bNeedUpdateDiscRens;

    glutils::renderable_t  *ren_surf;
    glutils::renderable_t  *ren_grad[gc_grid_dim];
    glutils::renderable_t  *ren_cp_labels[gc_grid_dim+1];
    glutils::renderable_t  *ren_cp[gc_grid_dim+1];
    glutils::renderable_t  *ren_cp_conns[gc_grid_dim];
    glutils::renderable_t  *ren_canc_cp_labels[gc_grid_dim+1];
    glutils::renderable_t  *ren_canc_cp[gc_grid_dim+1];
    glutils::renderable_t  *ren_canc_cp_conns[gc_grid_dim];

    glutils::bufobj_ptr_t   cp_loc_bo;

    std::vector<boost::shared_ptr<disc_rendata_t> > disc_rds;

    std::set<boost::shared_ptr<disc_rendata_t> >    active_disc_rens;


    void create_disc_rds();
    void update_active_disc_rens();

    void create_cp_loc_bo();

    void create_cp_rens(const rect_t &roi);
    void create_grad_rens(const rect_t &roi);
    void create_surf_ren(const rect_t &roi);
    void render() ;

    octtree_piece_rendata(octtree_piece *);

    virtual int         get_num_items();
    virtual bool        update_item(const int &,boost::any &,const bool &);
    virtual std::string get_description(int i);

  };

  class glviewer_t : public QGLViewer
  {
  public:

    std::vector<octtree_piece_rendata * >  m_grid_piece_rens;
    cellid_t                               m_size;
    rect_t                                 m_roi;


  public:
    glviewer_t(std::vector<octtree_piece *> * p ,
               cellid_t size,const rect_t &roi);

    ~glviewer_t();

  protected:

    virtual void draw();
    virtual void init();
    virtual QString helpString() const;
  };
}
#endif //VIEWER_H_INCLUDED
