#include <stack>
#include <queue>
#include <list>
#include <set>

#include <fstream>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/function.hpp>

#include <logutil.h>

#include <grid_dataset.h>
#include <grid_mscomplex.h>

namespace bl = boost::lambda;

using namespace std;

#define TRACEM(msg) (cout<<msg<<endl)
#define TRACEV(v)   (cout<<SVAR(v)<<endl)

namespace grid
{
  static uint ( dataset_t::*getcets[2] ) ( cellid_t,cellid_t * ) const =
  {
    &dataset_t::getCellFacets,
    &dataset_t::getCellCofacets
  };

  inline uint   dataset_t::getCellIncCells( cellid_t c,cellid_t * inc) const
  {
    for(uint i = 0; i < gc_grid_dim; ++i)
    {
      inc[i*2+0] = c;
      inc[i*2+1] = c;

      inc[i*2+0][i] -= 1;
      inc[i*2+0][i] += 1;
    }
    return gc_grid_dim*2;
  }

  dataset_t::dataset_t (const rect_t &r,const rect_t &e,const rect_t &d) :
      m_rect (r),
      m_ext_rect (e),
      m_domain_rect(d),
      m_cell_flags(cellid_t::zero,boost::fortran_storage_order()),
      m_vert_fns(cellid_t::zero,boost::fortran_storage_order())

  {
    // TODO: assert that the given rect is of even size..
    //       since each vertex is in the even positions

    using namespace boost::lambda;

    cmp_ftor     = bind(&dataset_t::compareCells,this,_1,_2);
    cmp_ftors[0] = bind(&dataset_t::compareCells,this,_1,_2);
    cmp_ftors[1] = bind(&dataset_t::compareCells,this,_2,_1);
  }

  dataset_t::~dataset_t ()
  {
    clear();
  }

  void dataset_t::init(const std::string &filename)
  {
    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    rect_size_t   span   = m_ext_rect.span() + 1;
    rect_size_t  pt_span = (m_ext_rect.span()/2)+1;
    rect_point_t bl = m_ext_rect.lower_corner();

    m_cell_flags.resize(span);
    m_vert_fns.resize(pt_span);

    uint num_cells = span[0]*span[1]*span[2];
    uint num_pts   = pt_span[0]*pt_span[1]*pt_span[2];

    m_cell_flags.reindex(bl);
    m_vert_fns.reindex(bl/2);

    std::fill_n(m_cell_flags.data(),num_cells,0);

    ifstream ifs(filename.c_str(),ios::in|ios::binary);
    ensure(ifs.is_open(),"unable to open file");

    ifs.read((char*)(void*)m_vert_fns.data(),sizeof(cell_fn_t)*num_pts);
    ensure(ifs.fail()==false,"failed to read some data");

    ifs.seekg(0,ios::end);
    ensure(ifs.tellg()==num_pts*sizeof(cell_fn_t),"file/piece size mismatch");

    ifs.close();
  }

  void  dataset_t::clear()
  {
    m_critical_cells.clear();

    m_cell_flags.resize(cellid_t::zero);

    m_vert_fns.resize(cellid_t::zero);
  }

  inline cellid_t flag_to_mxfct(cellid_t c,cell_flag_t f)
  {
    cell_flag_t d = f&0x07;
    ASSERT(is_in_range(d,1,7));
    c[(d-1)>>1] += (d&1)?(-1):(+1);
    return c;
  }

  inline cell_flag_t mxfct_to_flag(cellid_t c,cellid_t fct)
  {
    ASSERT(euclid_norm2(c-fct) == 1);

    int d = 0;

    if(c[1] != fct[1])
      d = 1;
    else if(c[2] != fct[2])
      d = 2;

    if(c[d] > fct[d])
      return (1 + d*2 + 0);
    else
      return (1 + d*2 + 1);
  }

  inline cellid_t flag_to_pair(cellid_t c,cell_flag_t f)
  {
    return flag_to_mxfct(c,(f>>3)&0x07);
  }

  inline cell_flag_t pair_to_flag(cellid_t c,cellid_t p)
  {
    return (mxfct_to_flag(c,p)<<3);
  }

  bool dataset_t::areCellsIncident(cellid_t c1,cellid_t c2) const
  {
    return ( euclid_norm2(c1-c2) == 1);
  }

  cellid_t dataset_t::getCellPairId (cellid_t c) const
  {
    ASSERT(isCellPaired(c));

    return flag_to_pair(c,m_cell_flags(c));
  }

  cellid_t dataset_t::getCellMaxFacetId (cellid_t c) const
  {
    return flag_to_mxfct(c,m_cell_flags(c));
  }

  cellid_t dataset_t::getCellSecondMaxFacetId (cellid_t c) const
  {
    return (2*c - flag_to_mxfct(c,m_cell_flags(c)));
  }

  inline bool compareCells_dim
      (const dataset_t * ds,
       const cellid_t &c1,
       const cellid_t &c2,
       const int & dim)
  {
    if(dim == 0)
      return ds->ptLt(c1,c2);

    cellid_t f1 = ds->getCellMaxFacetId(c1);
    cellid_t f2 = ds->getCellMaxFacetId(c2);

    if(f1 != f2)
      return compareCells_dim(ds,f1,f2,dim-1);

    f1 = ds->getCellSecondMaxFacetId(c1);
    f2 = ds->getCellSecondMaxFacetId(c2);

    int boundry_ct1 = ds->m_domain_rect.boundryCount(f1);
    int boundry_ct2 = ds->m_domain_rect.boundryCount(f2);

    if(boundry_ct1 != boundry_ct2)
      return (boundry_ct2 < boundry_ct1);

    return compareCells_dim(ds,f1,f2,dim-1);
  }

  bool dataset_t::compareCells( cellid_t c1,cellid_t  c2 ) const
  {
    int c1_dim = getCellDim(c1);
    int c2_dim = getCellDim(c2);

    cellid_t &p = (c1_dim > c2_dim)?(c1):(c2);

    for(int i = 0 ; i < std::abs(c1_dim-c2_dim);++i)
      p = getCellMaxFacetId(p);

    if(c1 == c2)
      return (c1_dim < c2_dim);

    return compareCells_dim(this,c1,c2,std::min(c1_dim,c2_dim));
  }

  uint dataset_t::getCellPoints (cellid_t c,cellid_t  *p) const
  {
    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    static_assert(((cell_coord_t)-1) < 0 && "coord_t needs to support -1 ");

    uint pos = 0;

    cellid_t i;

    for(i[2] = -(c[2]&1) ; i[2] <= (c[2]&1) ;i[2]+=2)
    {
      for(i[1] = -(c[1]&1) ; i[1] <= (c[1]&1) ;i[1]+=2)
      {
        for(i[0] = -(c[0]&1) ; i[0] <= (c[0]&1) ;i[0]+=2)
        {
          p[pos++] = c+i;
        }
      }
    }

    return (1<<getCellDim (c));
  }

  uint dataset_t::getCellFacets (cellid_t c,cellid_t *f) const
  {
    uint pos = 0;

    for(uint d = 0; d< gc_grid_dim; ++d)
    {
      for(uint i = 0 ; i < (c[d]&1);++i)
      {
        f[pos] = c; f[pos++][d] += 1;
        f[pos] = c; f[pos++][d] -= 1;
      }
    }
    return getCellDim (c)*2;
  }

  uint dataset_t::getCellCofacets (cellid_t c,cellid_t *cf) const
  {
    uint cf_ct = (gc_grid_dim - getCellDim (c))*2 ;

    uint pos = 0;

    for(uint d = 0; d< gc_grid_dim; ++d)
    {
      for(uint i = 0 ; i < ((c[d]+1)&1);++i)
      {

        cf[pos] = c; cf[pos++][d] += 1;
        cf[pos] = c; cf[pos++][d] -= 1;
      }
    }

    uint cf_nv_pos = 0;

    for (uint i = 0 ;i < cf_ct;++i)
      if (m_ext_rect.contains (cf[i]))
        cf[cf_nv_pos++] = cf[i];

    return cf_nv_pos;

  }

  uint dataset_t::getCellCofaces (cellid_t c,cellid_t *cf) const
  {
    uint cf_ct = std::pow(3,(gc_grid_dim - getCellDim (c))) ;

    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    static_assert(((cell_coord_t)-1) < 0 && "coord_t needs to support -1 ");

    uint pos = 0;

    cellid_t i,l = (c+cellid_t::one)&(cellid_t::one);

    for(i[2] = -(l[2]) ; i[2] <= (l[2]) ;i[2]+=1)
    {
      for(i[1] = -(l[1]) ; i[1] <= (l[1]) ;i[1]+=1)
      {
        for(i[0] = -(l[0]) ; i[0] <= (l[0]) ;i[0]+=1)
        {
          cf[pos++] = c+i;
        }
      }
    }

    uint cf_nv_pos = 0;

    for (uint i = 0 ;i < cf_ct;++i)
      if (m_ext_rect.contains (cf[i]))
        cf[cf_nv_pos++] = cf[i];

    return cf_nv_pos;
  }

  uint dataset_t::getCellEst (cellid_t c,cellid_t* est)  const
  {
    cellid_t cfs[40];

    int cfs_ct = getCellCofaces(c,cfs);

    ASSERT(is_in_range(cfs_ct,0,40));

    uint pos = 0;

    for(int i = 0 ; i< cfs_ct;++i)
    {
      cellid_t cf = cfs[i];

      ASSERT(m_ext_rect.contains(cf));

      uint c_dim  = getCellDim(c);
      uint cf_dim = getCellDim(cf);

      for(uint j = c_dim; j < cf_dim ; ++j)
        cf = getCellMaxFacetId(cf);

      if(cf == c)
        est[pos++] = cfs[i];
    }

    return pos;
  }

  bool dataset_t::isPairOrientationCorrect (cellid_t c, cellid_t p) const
  {
    return (getCellDim (c) <getCellDim (p));
  }

//  bool dataset_t::isCellMarked (cellid_t c) const
//  {
//    return ! (m_cell_flags (c) == CELLFLAG_UNKNOWN);
//  }

  bool dataset_t::isCellCritical (cellid_t c) const
  {
    return (m_cell_flags (c) & CELLFLAG_CRITICAL);
  }

  bool dataset_t::isCellPaired (cellid_t c) const
  {
    return (((m_cell_flags(c)>>3) & 0x07) !=0);
  }

  bool dataset_t::isCellVisited (cellid_t c) const
  {
    return (m_cell_flags (c) & CELLFLAG_VISITED);
  }


  void dataset_t::pairCells (cellid_t c,cellid_t p)
  {
    ASSERT(isCellPaired(c) == false);
    ASSERT(isCellPaired(p) == false);

    m_cell_flags (c) |= (pair_to_flag(c,p));
    m_cell_flags (p) |= (pair_to_flag(p,c));

    ASSERT(getCellPairId(c) == p);
    ASSERT(getCellPairId(p) == c);
  }

  void dataset_t::visitCell(cellid_t c)
  {
    ASSERT(isCellVisited(c) == false);
    m_cell_flags (c) |= CELLFLAG_VISITED;
    ASSERT(isCellVisited(c) == true);
  }

//  void  dataset_t::unpairCells ( cellid_t c,cellid_t p )
//  {
//    ASSERT(getCellPairId(c) == p);
//    ASSERT(getCellPairId(p) == c);

//    m_cell_pairs (c) = CELLADJDIR_UNKNOWN;
//    m_cell_pairs (p) = CELLADJDIR_UNKNOWN;

//    m_cell_flags (c) &= (CELLFLAG_MASK^CELLFLAG_PAIRED);
//    m_cell_flags (p) &= (CELLFLAG_MASK^CELLFLAG_PAIRED);

//    ASSERT(isCellPaired(c) == false);
//    ASSERT(isCellPaired(p) == false);
//  }

  void dataset_t::setCellMaxFacet (cellid_t c,cellid_t f)
  {
    ASSERT(getCellDim(c) == getCellDim(f)+1);
    m_cell_flags (c) |= mxfct_to_flag(c,f);
    ASSERT(getCellMaxFacetId(c) == f);
  }

  void dataset_t::markCellCritical (cellid_t c)
  {
    ASSERT(isCellCritical(c) == false);
    m_cell_flags (c) |= CELLFLAG_CRITICAL;
    ASSERT(isCellCritical(c) == true);
  }

  bool dataset_t::isTrueBoundryCell (cellid_t c) const
  {
    return (m_domain_rect.isOnBoundry (c));
  }

  bool dataset_t::isFakeBoundryCell (cellid_t c) const
  {
    return (m_rect.isOnBoundry (c) && (!m_domain_rect.isOnBoundry (c)));
  }

  bool dataset_t::isCellExterior (cellid_t c) const
  {
    return (!m_rect.contains (c));
  }

  void dataset_t::assignMaxFacets()
  {
    using namespace boost::lambda;

    cellid_t f[20],c;

    for(uint dim = 1 ; dim <= gc_grid_dim; ++dim)
    {
      cellid_t   seq(cellid_t::zero);

      std::for_each(seq.begin(),seq.begin()+dim,_1=1);

      do
      {
        static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

        for(c[2] = m_ext_rect[2][0] + seq[2] ; c[2] <= m_ext_rect[2][1]; c[2] += 2)
        {
          for(c[1] = m_ext_rect[1][0] + seq[1]; c[1] <= m_ext_rect[1][1]; c[1] += 2)
          {
            for(c[0] = m_ext_rect[0][0] + seq[0]; c[0] <= m_ext_rect[0][1]; c[0] += 2)
            {
              int f_ct = getCellFacets(c,f);

              setCellMaxFacet(c,*std::max_element(f,f+f_ct,cmp_ftor));
            }
          }
        }
      }
      while(std::next_permutation(seq.rbegin(),seq.rend()));
    }
  }

  void  dataset_t::pairCellsWithinEst()
  {
    using namespace boost::lambda;

    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

//    cellid_array_t est_vert;

//    cellid_t s = m_ext_rect.span() + cellid_t::one;
//    int num_cells = s[0]*s[1]*s[2];
//    cellid_t bl = m_ext_rect.lower_corner();

//    est_vert.resize(s);
//    std::fill_n(est_vert.data(),num_cells,cellid_t(-1,-1,-1));
//    est_vert.reindex(bl);

    for(cellid_t c = m_rect.lower_corner() ; c[2] <= m_rect[2][1]; c[2] += 2)
    {
      for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; c[1] += 2)
      {
        for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; c[0] += 2)
        {
          cellid_t est_arr[40];

          uint est_ct = getCellEst(c,est_arr);

          std::sort(est_arr,est_arr+est_ct,cmp_ftor);

//          for(int i = 0; i < est_ct; ++i)
//          {
//            if(est_vert(est_arr[i]) != cellid_t(-1,-1,-1))
//            {
//              cout<<c<<endl;
//              cout<<est_arr[i];
//              cout<<est_vert(est_arr[i]);

//              throw runtime_error("gotcha");
//            }

//            est_vert(est_arr[i]) = c;
//          }

          for(int j = 0; j < est_ct; ++j)
          {
            cellid_t c = est_arr[j];

            if(j+1 < est_ct)
            {
              cellid_t p = est_arr[j+1];

              bool is_adj        = areCellsIncident(c,p);
              bool is_same_bndry = (isTrueBoundryCell(c) == isTrueBoundryCell(p));

              if(is_adj && is_same_bndry)
              {
                pairCells(c,p);++j;
                continue;
              }
            }

            if(m_rect.contains(c))
            {
              markCellCritical(c);
              m_critical_cells.push_back(c);
            }
          }


//          for(cellid_llist_t::iterator j = est_list.begin(); j != est_list.end();)
//          {
//            cellid_llist_t::iterator c_it = j,p_it = j;
//            ++j;
//            while( p_it != est_list.begin())
//            {
//              p_it--;
//              bool is_adj        = areCellsIncident(*c_it,*p_it);
//              bool is_same_bndry = (isTrueBoundryCell(*c_it) == isTrueBoundryCell(*p_it));
//              if(is_adj && is_same_bndry)
//              {
//                ++j;
//                pairCells(*c_it,*p_it);
//                est_list.erase(c_it);
//                est_list.erase(p_it);
//              }
//            }
//          }

        }
      }
    }

//    for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; c[2] += 1)
//    {
//      for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; c[1] += 1)
//      {
//        for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; c[0] += 1)
//        {
//          if(est_vert(c) == cellid_t(-1,-1,-1))
//            throw runtime_error("gotcha");

//          if(!isCellPaired(c) && !isCellCritical(c))
//            throw runtime_error("gotcha");
//        }
//      }
//    }
  }

  void dataset_t::assignGradient()
  {
    assignMaxFacets();
    pairCellsWithinEst();
  }

  void dataset_t::markBoundryCritical(const rect_t &bnd)
  {
    try
    {
      ASSERT(m_ext_rect.intersection(bnd) == bnd);
    }
    catch(assertion_error e)
    {
      e.push(_FFL);
      e.push(SVAR(bnd));
      e.push(SVAR(m_rect));
      e.push(SVAR(m_ext_rect));

      throw;
    }

    rect_t ixn = bnd.intersection(m_rect);

    for(cellid_t c = ixn.lower_corner() ; c[2] <= ixn[2][1]; c[2] += 1)
    {
      for(c[1] = ixn[1][0] ; c[1] <= ixn[1][1]; c[1] += 1)
      {
        for(c[0] = ixn[0][0] ; c[0] <= ixn[0][1]; c[0] += 1)
        {
          if(!isCellPaired(c))
            continue;

          cellid_t p = getCellPairId(c);

          if(!isPairOrientationCorrect(c,p))
            continue;

          if(bnd.contains(p))
            continue;

          try
          {
            markCellCritical(c);
            m_critical_cells.push_back(c);

            if(!m_rect.contains(p))
              continue;

            markCellCritical(p);
            m_critical_cells.push_back(p);
          }
          catch(assertion_error e)
          {
            e.push(_FFL);
            e.push(SVAR(c));
            e.push(SVAR(p));
            e.push(SVAR(ixn));
            e.push(SVAR(m_rect));

            throw;
          }
        }
      }
    }
  }

  namespace dfs
  {
    using namespace grid;
    using namespace boost::lambda;

    typedef std::stack<cellid_t>                              stack_t;
    typedef boost::function<bool (cellid_t,const stack_t &)>  can_visit_ftor_t;
    typedef boost::function<void (cellid_t,const stack_t &)>  visit_ftor_t;
    typedef boost::function<void (cellid_t,const stack_t &)>  cp_visit_ftor_t;

    void do_dfs
        (dataset_ptr_t ds,
         cellid_t start_cell,
         eGradientDirection dir,
         can_visit_ftor_t can_visit_ftor,
         visit_ftor_t visit_ftor,
         cp_visit_ftor_t  cp_visit_ftor)
    {
      stack_t cell_stack;

      cell_stack.push ( start_cell );

      uint dim = ds->getCellDim(start_cell);

      while ( !cell_stack.empty() )
      {
        cellid_t top_cell = cell_stack.top();

        try
        {
          ASSERT(ds->m_rect.contains(top_cell));

          visit_ftor(top_cell,cell_stack);

          cell_stack.pop();

          cellid_t      cets[20];

          uint cet_ct = ( ds.get()->*getcets[dir] ) ( top_cell,cets );

          for ( uint i = 0 ; i < cet_ct ; i++ )
          {
            try
            {
              if ( ds->isCellExterior ( cets[i] ) )
                continue;

              if ( ds->isCellCritical ( cets[i] ) )
              {
                cp_visit_ftor(cets[i],cell_stack);
                continue;
              }

              cellid_t next_cell = ds->getCellPairId ( cets[i] );

              ASSERT(ds->m_rect.contains(next_cell));

              if(can_visit_ftor(next_cell,cell_stack) == false)
                continue;

              bool is_dim        = (dim  == ds->getCellDim ( next_cell ));
              bool is_vnext      = ds->cmp_ftors[dir](next_cell,top_cell);

              if (is_dim && is_vnext)
              {
                cell_stack.push ( next_cell );
              }
            }
            catch(assertion_error e)
            {
              e.push(_FFL);
              e.push(SVAR(cets[i]));
              e.push(SVAR(ds->isCellPaired(cets[i])));
              e.push(SVAR(ds->isCellCritical(cets[i])));
              throw;
            }
          }
        }
        catch(assertion_error e)
        {
          e.push(_FFL);
          e.push(SVAR(dim));
          e.push(SVAR(top_cell));
          e.push(SVAR(start_cell));

          throw;
        }
      }
    }

    typedef std::stack<int> age_stack_t;
    typedef boost::function<bool (cellid_t,const stack_t &,const stack_t &)>  path_can_visit_ftor_t;
    typedef boost::function<void (cellid_t,const stack_t &,const stack_t &)>  path_visit_ftor_t;
    typedef boost::function<void (cellid_t,const stack_t &,const stack_t &)>  path_cp_visit_ftor_t;

    void visit_ftor_wrapper
        (path_visit_ftor_t path_visit_ftor,
         cellid_t c,
         const stack_t &cell_stack,
         stack_t *path_stack,
         age_stack_t *age_stack)
    {
      ASSERT(cell_stack.top() == c);

      while(age_stack->top() > cell_stack.size())
      {
        age_stack->pop();
        path_stack->pop();
      }

      path_stack->push(c);
      age_stack->push(cell_stack.size());

      path_visit_ftor(c,cell_stack,*path_stack);
    }

    void do_path_dfs
        (dataset_ptr_t ds,
         cellid_t c,
         eGradientDirection dir,
         path_visit_ftor_t    path_visit_ftor,
         path_cp_visit_ftor_t path_cp_visit_ftor,
         path_can_visit_ftor_t path_can_visit_ftor)
    {
      stack_t     path_stack;
      age_stack_t age_stack;

      do_dfs
          (ds,c,dir,
           bind(path_can_visit_ftor,_1,_2,path_stack),
           bind(visit_ftor_wrapper,path_visit_ftor,_1,_2,&path_stack,&age_stack),
           bind(path_cp_visit_ftor,_1,_2,path_stack));
    }

    void pass_visit(cellid_t ,const stack_t &)
    {
      return;
    }

    bool pass_can_visit(cellid_t ,const stack_t &)
    {
      return true;
    }

    void connect_cps(mscomplex_ptr_t msc,cellid_t c1,cellid_t c2,const stack_t &)
    {
      msc->connect_cps(c1,c2);
    }

    void do_dfs_connect
        (dataset_ptr_t ds,
         mscomplex_ptr_t msc,
         cellid_t c,
         eGradientDirection dir)
    {
      do_dfs(ds,c,dir,pass_can_visit,pass_visit,bind(connect_cps,msc,c,_1,_2));
    }

    void add_to_disc(std::set<cellid_t> *mfold,cellid_t c,const stack_t &)
    {
      mfold->insert(c);
    }

    void do_dfs_collect_manifolds
        (dataset_ptr_t ds,
         std::set<cellid_t> *mfold,
         cellid_t c,
         eGradientDirection dir)
    {
      do_dfs(ds,c,dir,pass_can_visit,bind(add_to_disc,mfold,_1,_2),pass_visit);
    }

    bool visit_if_not_visited(dataset_ptr_t ds,cellid_t c,const stack_t &)
    {
      return (ds->isCellVisited(c) == false);
    }

    void mark_visit(dataset_ptr_t ds,cellid_t c,const stack_t &)
    {
      ds->visitCell(c);
    }

    void do_dfs_mark_visit(dataset_ptr_t ds,cellid_t c,eGradientDirection dir)
    {
      do_dfs(ds,c,dir,bind(visit_if_not_visited,ds,_1,_2),
             bind(mark_visit,ds,_1,_2),pass_visit);
    }

    bool visit_if_pair_visited(dataset_ptr_t ds,cellid_t c,const stack_t &)
    {
      return (ds->isCellVisited(c) && ds->isCellVisited(ds->getCellPairId(c)));
    }

    void do_dfs_connect_thru_visted_pairs
        (dataset_ptr_t ds,
         mscomplex_ptr_t msc,
         cellid_t c,
         eGradientDirection dir)
    {
      do_dfs(ds,c,dir,bind(visit_if_pair_visited,ds,_1,_2),
             pass_visit,bind(connect_cps,msc,c,_1,_2));
    }
  };

  void  dataset_t::computeMsGraph(mscomplex_ptr_t msgraph)
  {
    using namespace boost::lambda;

    for(uint i = 0 ; i < m_critical_cells.size();++i)
    {
      cellid_t c = m_critical_cells[i];

      cellid_t v = c;

      for( int j = 0 ; j < getCellDim(c);++j)
        v=getCellMaxFacetId(v);

      try
      {
        msgraph->add_critpt(c,getCellDim(c),get_cell_fn(c),v);
      }
      catch (assertion_error e)
      {
        e.push(_FFL);
        e.push(SVAR(c));
        e.push(SVAR(isCellPaired(c)));
        e.push(SVAR(isCellCritical(c)));

        if(isCellPaired(c))
          e.push(SVAR(isCellPaired(c)));

        throw;
      }
    }

    for(uint i = 0 ; i < m_critical_cells.size();++i)
    {
      cellid_t c = m_critical_cells[i];

      if(!isCellPaired(c))
        continue;

      cellid_t p = getCellPairId(c);

      if(!isPairOrientationCorrect(c,p))
        continue;

      if(!m_rect.contains(p))
        continue;

      msgraph->pair_cps(c,p);
    }

    try
    {
      for(uint i = 0 ; i < m_critical_cells.size();++i)
      {
        cellid_t c = m_critical_cells[i];

        switch(getCellDim(c))
        {
        case 0:
          dfs::do_dfs_connect(shared_from_this(),msgraph,c,GRADDIR_ASCENDING);
          break;
        case 1:
          dfs::do_dfs_mark_visit(shared_from_this(),c,GRADDIR_ASCENDING);
          break;
        case 2:
          dfs::do_dfs_mark_visit(shared_from_this(),c,GRADDIR_DESCENDING);
          break;
        case 3:
          dfs::do_dfs_connect(shared_from_this(),msgraph,c,GRADDIR_DESCENDING);
          break;
        }
      }

      for(uint i = 0 ; i < m_critical_cells.size();++i)
      {
        cellid_t c = m_critical_cells[i];

        if(getCellDim(c) == 2)
          dfs::do_dfs_connect_thru_visted_pairs
              (shared_from_this(),msgraph,c,GRADDIR_DESCENDING);
      }

    }
    catch(assertion_error e)
    {
      e.push(_FFL);
      e.push(SVAR(m_rect));
      e.push(SVAR(m_ext_rect));
      e.push(SVAR(m_domain_rect));

      throw;
    }
  }

  template<typename T>
  void bin_write(std::ostream & os,const T & d)
  {
    os.write((const char*)(const void*)&d,sizeof(T));
//    os<<d<<endl;
  }

  int dataset_t::saveManifolds(mscomplex_ptr_t msc,std::ostream &os,int i,int dir)
  {
    critpt_t * cp = msc->m_cps[i];

    set<cellid_t> mfold;

    try
    {
      ASSERT(cp->is_paired == false);

      if(m_rect.contains(cp->cellid))
        dfs::do_dfs_collect_manifolds
            (shared_from_this(),&mfold,cp->cellid,(eGradientDirection)dir);

      for( conn_iter_t it = cp->conn[dir].begin(); it != cp->conn[dir].end();++it)
      {
        critpt_t * ccp    = msc->m_cps[*it];
        critpt_t * ccp_pr = msc->m_cps[ccp->pair_idx];

        ASSERT(ccp->is_paired == true);
        ASSERT(cp->index == ccp_pr->index);

        dfs::do_dfs_collect_manifolds
            (shared_from_this(),&mfold,ccp_pr->cellid,(eGradientDirection)dir);
      }

      for(set<cellid_t>::iterator it = mfold.begin();it != mfold.end(); ++it)
        bin_write(os,*it);
    }
    catch(assertion_error e)
    {
      e.push(_FFL);
      e.push(SVAR(m_rect));
      e.push(SVAR(m_ext_rect));

      throw;
    }

    return mfold.size();
  }


  void  dataset_t::saveManifolds(mscomplex_ptr_t msc,std::ostream &os)
  {
    bin_write(os,m_rect);
    bin_write(os,m_ext_rect);
    bin_write(os,m_domain_rect);

    int cp_ct_wpos = os.tellp();
    os.seekp(sizeof(int),ios::cur);
    int num_cps = 0;

    for(int i = 0 ; i < msc->m_cps.size();++i)
    {
      critpt_t * cp = msc->m_cps[i];

      if(cp->is_paired)
        continue;

      bin_write<cellid_t>(os,cp->cellid);
      num_cps++;
    }

    os.seekp(cp_ct_wpos,ios::beg);
    bin_write<int>(os,num_cps);
    os.seekp(0,ios::end);

    std::vector<int> mfold_offsets;
    mfold_offsets.reserve(2*num_cps+1);
    int mfold_off_wpos= os.tellp();
    os.seekp(sizeof(int)*(2*num_cps+1),ios::cur);
    int mfold_offset  = 0;
    mfold_offsets.push_back(mfold_offset);

    try
    {
      for(int i = 0 ; i < msc->m_cps.size();++i)
      {
        if(msc->m_cps[i]->is_paired)
          continue;

        for(int d = 0 ; d < 2;++d)
        {
          mfold_offset += saveManifolds(msc,os,i,d);
          mfold_offsets.push_back(mfold_offset);
        }
      }
      ASSERT(mfold_offsets.size() == 2*num_cps+1);
    }
    catch(assertion_error e)
    {
      e.push(_FFL);
      e.push(SVAR(num_cps));
      e.push(SVAR(mfold_offset));
      e.push(SVAR(mfold_offsets.size()));

      throw;
    }

    os.seekp(mfold_off_wpos,ios::beg);
    for(int i = 0 ; i < mfold_offsets.size();++i)
      bin_write<int>(os,mfold_offsets[i]);
    os.seekp(0,ios::end);
  }

  void  dataset_t::saveManifolds(mscomplex_ptr_t msc,const std::string &s)
  {
    std::ofstream fs(s.c_str());
    ensure(fs.is_open(),"unable to open file");
    saveManifolds(msc,fs);
    fs.close();
  }

  void dataset_t::log_flags()
  {
    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    cellid_t c;

    for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; ++c[2])
    {
      for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; ++c[1])
      {
        for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; ++c[0])
        {
          std::cout<<m_cell_flags(c)<<" ";
        }
        std::cout<<std::endl;
      }
      std::cout<<std::endl;
    }
  }

  char get_dir_txt(cellid_t c,cellid_t p)
  {
    int dir = 0;

    if (c[1] != p[1])
      dir = 1;
    else if( c[2] != p[2])
      dir = 2;

    if(dir == 2)
    {
      if( c[dir] > p[dir])
        return 'd';
      else
        return 'u';
    }

    if(c[dir]&1)
    {
      if(c[dir] > p[dir])
      {
        switch(dir)
        {
        case 0: return '>';
        case 1: return 'v';
        }
      }
      else
      {
        switch(dir)
        {
        case 0: return '<';
        case 1: return '^';
        }
      }
    }
    else
    {
      switch(dir)
      {
      case 0: return '-';
      case 1: return '|';
      }
    }

    return '#';
  }

  void dataset_t::log_pairs(std::ostream &os)
  {
    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    cellid_t c;

    for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; ++c[2])
    {
      os<<"sheet no:: "<<c[2]<<std::endl;
      for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; ++c[1])
      {
        for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; ++c[0])
        {
          if(isCellCritical(c))
          {
            if(isCellPaired(c))
              os<<get_dir_txt(c,getCellPairId(c))<<"c";
            else
              os<<"C ";
          }
          else if(isCellPaired(c))
            os<<get_dir_txt(c,getCellPairId(c))<<" ";
          else
            os<<"? ";
        }
        os<<std::endl;
      }

    }
  }

  void dataset_t::log_visits(std::ostream &os)
  {
    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    cellid_t c;

    for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; ++c[2])
    {
      os<<"sheet no:: "<<c[2]<<std::endl;
      for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; ++c[1])
      {
        for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; ++c[0])
        {
          if(isCellVisited(c))
            os<<"x ";
          else
            os<<". ";
        }
        os<<std::endl;
      }
    }
  }

  void dataset_t::log_pair_visits(std::ostream &os)
  {
    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    cellid_t c;

    for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; ++c[2])
    {
      os<<"sheet no:: "<<c[2]<<std::endl;
      for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; ++c[1])
      {
        for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; ++c[0])
        {
          if(isCellVisited(c)&&isCellPaired(c)&&isCellVisited(getCellPairId(c)))
            os<<"x ";
          else
            os<<". ";
        }
        os<<std::endl;
      }
    }
  }

  void dataset_t::log_pairs(const std::string &s)
  {
    std::ofstream fs(s.c_str());
    ensure(fs.is_open(),"unable to open file");
    log_pairs(fs);
    fs.close();
  }

  void dataset_t::log_pair_visits(const std::string &s)
  {
    std::ofstream fs(s.c_str());
    ensure(fs.is_open(),"unable to open file");
    log_pair_visits(fs);
    fs.close();
  }

  void dataset_t::log_visits(const std::string &s)
  {
    std::ofstream fs(s.c_str());
    ensure(fs.is_open(),"unable to open file");
    log_visits(fs);
    fs.close();
  }


  void dataset_t::log_max_facets()
  {
    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    cellid_t c;

    for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; ++c[2])
    {
      for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; ++c[1])
      {
        for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; ++c[0])
        {
          if(getCellDim(c) != 0 )
            std::cout<< getCellMaxFacetId(c);
          else
            std::cout<< "(.,.,.)";
        }
        std::cout<<std::endl;
      }
      std::cout<<std::endl;
    }
  }

  void dataset_t::extract_vdata_subarray(rect_t r,const std::string &filename)
  {
    if(r.lower_corner()%2 != cellid_t::zero ||
       r.upper_corner()%2 != cellid_t::zero )
    {
      throw std::runtime_error("r must specify an aabb with vertex end pts");
    }

    std::ofstream ofs(filename.c_str(),std::ios::out|std::ios::binary);

    if(ofs.is_open() == false)
      throw std::runtime_error("unable to open file");

    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    cellid_t c;

    for(c[2] = r[2][0] ; c[2] <= r[2][1]; c[2] +=2)
    {
      for(c[1] = r[1][0] ; c[1] <= r[1][1]; c[1] +=2)
      {
        for(c[0] = r[0][0] ; c[0] <= r[0][1]; c[0] +=2)
        {
          cell_fn_t fn = get_cell_fn(c);

          ofs.write((char*)(void*)&fn,sizeof(cell_fn_t));
        }
      }
    }
  }
}
