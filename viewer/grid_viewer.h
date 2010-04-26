#ifndef GRID_VIEWER_H_INCLUDED
#define GRID_VIEWER_H_INCLUDED
#include <grid.h>

#include <glutils.h>
#include <set>

#include <boost/any.hpp>

namespace grid
{
  class octtree_piece ;

  class mscomplex_t;

  class configurable_t
  {
  public:

    typedef two_tuple_t<int> data_index_t;

    enum eExchangeMode {EXCHANGE_READ,EXCHANGE_WRITE};

    virtual int         rows()    = 0 ;

    virtual int         columns() = 0 ;

    virtual bool        exchange_data(const data_index_t &,
                                      boost::any &,
                                      const eExchangeMode &) = 0;

    virtual std::string get_header(int i)
    {
      std::stringstream ss;
      ss<<i;
      return ss.str();
    }

    template <typename T>
        static bool s_exchange_read_write
        (T &p_val,boost::any &c_val,const eExchangeMode & mode)
    {
      bool ret = true;

      switch(mode)
      {
      case EXCHANGE_WRITE:
        ret   = (p_val != boost::any_cast<T>(c_val));
        p_val = boost::any_cast<T>(c_val);
        break;
      case EXCHANGE_READ:
        c_val = boost::any(p_val);
        break;
      }
      return ret;
    }

    template <typename T>
        static bool s_exchange_read_only
        (const T &p_val,boost::any &c_val,const eExchangeMode & mode)
    {
      switch(mode)
      {
      case EXCHANGE_WRITE:
        throw std::logic_error("read only property cannot write");
        break;
      case EXCHANGE_READ:
        c_val = boost::any(p_val);
      }
      return false;
    }
  };

  class disc_rendata_t
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
  };

  class octtree_piece_rendata:public configurable_t
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


    // configurable_t interface
  public:
    int rows();
    int columns();
    bool exchange_data(const data_index_t &,boost::any &,const eExchangeMode &);
    std::string get_header(int i);
  };

  class data_manager_t;

  class grid_viewer_t:
      public glutils::renderable_t,
      public configurable_t
  {
  public:
    std::vector<octtree_piece_rendata * >  m_grid_piece_rens;
    cellid_t                               m_size;
    rect_t                                 m_roi;

  public:

    grid_viewer_t(data_manager_t * p ,const rect_t &roi);

    ~grid_viewer_t();


    void init();

    // renderable_t interface
  public:
    int  render();

    // configurable_t interface
  public:
    int rows();
    int columns();
    bool exchange_data(const data_index_t &,boost::any &,const eExchangeMode &);
    std::string get_header(int i);
  };
}
#endif //VIEWER_H_INCLUDED
