#include <stack>
#include <queue>
#include <list>
#include <set>

#include <fstream>

#include <boost/typeof/typeof.hpp>
#include <boost/function.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/static_assert.hpp>

#define static_assert BOOST_STATIC_ASSERT

#include <config.h>

#include <grid_dataset.h>
#include <grid_mscomplex.h>

#ifdef BUILD_EXEC_OPENCL
#include <grid_dataset_cl.h>
#endif

using namespace std;

#define MUTEX_CALL(__c) \
  {static boost::mutex __mutex;\
   boost::mutex::scoped_lock scoped_lock(__mutex);\
   (__c);}




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
    for(int i = 0 ; i < g_num_threads; ++i)
      m_critical_cells[i].clear();

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

  void dataset_t::assignMaxFacets_thd(int thid,int dim)
  {
    cellid_t f[20],c,s(0,0,0),stride(2,2,2);

    for(int i = 0 ; i < dim; ++i)
      s[i] = 1;

    while(true)
    {
      rect_t rect  = rect_t(m_ext_rect.lc()+s,m_ext_rect.uc()-s);

      int n = c_to_i(rect.uc(),rect,stride) + 1;

      for( int i = thid; i < n; i += g_num_threads)
      {
        c = i_to_c(i,rect,stride);

        ASSERT(getCellDim(c) == dim);

        int f_ct   = getCellFacets(c,f);

        setCellMaxFacet(c,*std::max_element(f,f+f_ct,cmp_ftor));
      }

      if(!next_permutation(s.rbegin(),s.rend()))
        break;
    }
  }

  void  dataset_t::pairCellsWithinEst_thd(int thid)
  {
    cellid_t stride(2,2,2);

    int n = c_to_i(m_rect.uc(),m_rect,stride)+1;

    for(int i = thid; i < n; i += g_num_threads)
    {
      cellid_t c = i_to_c(i,m_rect,stride);

      cellid_t est_arr[40];

      uint est_ct = getCellEst(c,est_arr);

      std::sort(est_arr,est_arr+est_ct,cmp_ftor);

      for(int j = 0; j < est_ct; ++j)
      {
        cellid_t p = est_arr[j];

        if(j+1 < est_ct)
        {
          cellid_t q = est_arr[j+1];

          bool is_adj        = areCellsIncident(p,q);
          bool is_same_bndry = (isTrueBoundryCell(p) == isTrueBoundryCell(q));

          if(is_adj && is_same_bndry)
          {
            if(m_rect.contains(p))
            {
              pairCells(p,q);

              if(m_rect.contains(q) == false)
              {
                markCellCritical(p);
                m_critical_cells[thid].push_back(p);
              }
              else if(m_rect.boundryCount(p) != m_rect.boundryCount(q))
              {
                markCellCritical(p);
                m_critical_cells[thid].push_back(p);

                markCellCritical(q);
                m_critical_cells[thid].push_back(q);
              }
            }

            ++j;
            continue;
          }
        }

        if(m_rect.contains(p))
        {
          markCellCritical(p);
          m_critical_cells[thid].push_back(p);
        }
      }
    }
  }

  namespace dfs
  {
    using namespace grid;

    typedef std::stack<cellid_t>                              stack_t;
    typedef boost::function<bool (cellid_t,const stack_t &)>  can_visit_ftor_t;
    typedef boost::function<void (cellid_t,const stack_t &)>  visit_ftor_t;
    typedef boost::function<void (cellid_t,const stack_t &)>  cp_visit_ftor_t;

    void do_dfs
        (dataset_const_ptr_t ds,
         cellid_t start_cell,
         eGDIR dir,
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
         eGDIR dir,
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

    void connect_cps(dataset_ptr_t ds, mscomplex_ptr_t msc,cellid_t p,cellid_t q,const stack_t &)
    {
      if(ds->isCellPaired(p) && ds->getCellDim(ds->getCellPairId(p)) != ds->getCellDim(q))
        return;

      if(ds->isCellPaired(q) && ds->getCellDim(ds->getCellPairId(q)) != ds->getCellDim(p))
        return;

      MUTEX_CALL(msc->connect_cps(p,q));
    }

    void do_dfs_connect
        (dataset_ptr_t ds,
         mscomplex_ptr_t msc,
         cellid_t c,
         eGDIR dir)
    {
      try
      {
        do_dfs(ds,c,dir,pass_can_visit,pass_visit,bind(connect_cps,ds,msc,c,_1,_2));
      }
      catch(assertion_error e)
      {
        e.push(_FFL).push(SVAR(c))
         .push(SVAR(ds->m_rect))
         .push(SVAR(msc->m_rect));

        throw;
      }
    }

    template<typename MFOLD_T>
    void add_to_mfold(MFOLD_T *mfold,cellid_t c,const stack_t &);

    template<>
    void add_to_mfold(std::set<cellid_t> *mfold,cellid_t c,const stack_t &)
    {
      mfold->insert(c);
    }

    template<>
    void add_to_mfold(cellid_list_t *mfold,cellid_t c,const stack_t &)
    {
      mfold->push_back(c);
    }

    template<typename mfold_t>
    void do_dfs_collect_manifolds
        (dataset_const_ptr_t ds,
         mfold_t *mfold,
         cellid_t c,
         eGDIR dir)
    {
      do_dfs(ds,c,dir,pass_can_visit,bind(add_to_mfold<mfold_t>,mfold,_1,_2),pass_visit);
    }

    bool visit_if_not_visited(dataset_ptr_t ds,cellid_t c,const stack_t &)
    {
      return (ds->isCellVisited(c) == false);
    }

    void mark_visit(dataset_ptr_t ds,cellid_t c,const stack_t &)
    {
      ds->visitCell(c);
    }

    void do_dfs_mark_visit(dataset_ptr_t ds,cellid_t c,eGDIR dir)
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
         eGDIR dir)
    {
      do_dfs(ds,c,dir,bind(visit_if_pair_visited,ds,_1,_2),
             pass_visit,bind(connect_cps,ds,msc,c,_1,_2));
    }
  }

  void  dataset_t::extrema_connect_thd
      (mscomplex_ptr_t msgraph, cp_producer_ptr_t prd)
  {
    for(int i ; prd->next(i);)
    {
      cellid_t c = msgraph->cellid(i);

      eGDIR dir = (msgraph->index(i) == 3)?(GDIR_DES):(GDIR_ASC);

      dfs::do_dfs_connect(shared_from_this(),msgraph,c,dir);
    }
  }

  void  dataset_t::saddle_visit(mscomplex_ptr_t msgraph)
  {
    for(int i = 0 ; i < msgraph->get_num_critpts();++i)
    {
      if(msgraph->index(i) == 0 ||msgraph->index(i) == gc_grid_dim )
        continue;

      eGDIR dir = (msgraph->index(i) == 2)?(GDIR_DES):(GDIR_ASC);

      dfs::do_dfs_mark_visit(shared_from_this(),msgraph->cellid(i),dir);
    }
  }

  void  dataset_t::saddle_connect_thd
      (mscomplex_ptr_t msgraph, cp_producer_ptr_t prd)
  {
    for(int i ; prd->next(i);)
    {
      cellid_t c = msgraph->cellid(i);

      dfs::do_dfs_connect_thru_visted_pairs(shared_from_this(),msgraph,c,GDIR_DES);
    }
  }

  void  dataset_t::setupCriticalPoint_thd(int tid, mscomplex_ptr_t msgraph, int cp_offset)
  {
    for(int i = 0 ; i < m_critical_cells[tid].size();++i)
    {
      cellid_t c = m_critical_cells[tid][i];

      cellid_t v = c;

      for( int j = 0 ; j < getCellDim(c);++j)
        v=getCellMaxFacetId(v);

      msgraph->set_critpt(cp_offset++,c,getCellDim(c),get_cell_fn(c),v);

      if(!isCellPaired(c))
        continue;

      cellid_t p = getCellPairId(c);

      if(isPairOrientationCorrect(c,p))
        continue;

      if(!m_rect.contains(p))
        continue;

      ASSERT(m_critical_cells[tid][i-1] == p);

      msgraph->pair_cps(cp_offset-1,cp_offset-2);
    }
  }

  void  dataset_t::computeMsGraph(mscomplex_ptr_t msgraph)
  {
#ifdef BUILD_EXEC_OPENCL
    opencl::assign_gradient(shared_from_this(),msgraph);
#else

    for(int dim = 1 ; dim <= gc_grid_dim; ++dim)
    {
      boost::thread_group group;

      for(int tid = 0 ; tid < g_num_threads; ++tid)
        group.create_thread(bind(&dataset_t::assignMaxFacets_thd,this,tid,dim));

      group.join_all();
    }

    {
      boost::thread_group group;

      for(int tid = 0 ; tid < g_num_threads; ++tid)
        group.create_thread(bind(&dataset_t::pairCellsWithinEst_thd,this,tid));

      group.join_all();
    }

    {
      int cp_offset[g_num_threads+1];

      cp_offset[0] = 0;

      for(int tid = 0 ; tid < g_num_threads; ++tid)
        cp_offset[tid + 1] = cp_offset[tid] + m_critical_cells[tid].size();

      msgraph->resize(cp_offset[g_num_threads]);

      boost::thread_group group;

      for(int tid = 0 ; tid < g_num_threads; ++tid)
        group.create_thread(bind(&dataset_t::setupCriticalPoint_thd,this,tid,msgraph,cp_offset[tid]));

      group.join_all();
    }

#endif

    msgraph->build_id_cp_map();

    {
      boost::thread_group group;

      cp_producer_ptr_t prd(new cp_producer_t(msgraph,cp_producer_t::extrema_filter));

      for(int tid = 0 ; tid < g_num_threads-1; ++tid)
        group.create_thread(bind(&dataset_t::extrema_connect_thd,this,msgraph,prd));

      saddle_visit(msgraph);

      group.join_all();
    }

    {
      boost::thread_group group;

      cp_producer_ptr_t prd(new cp_producer_t(msgraph,cp_producer_t::twosaddle_filter));

      for(int tid = 0 ; tid < g_num_threads; ++tid)
        group.create_thread(bind(&dataset_t::saddle_connect_thd,this,msgraph,prd));

      group.join_all();
    }
  }

  template<typename mfold_t>
  void  dataset_t::get_mfold
      (mfold_t *mfold, mscomplex_const_ptr_t msc, int i, int dir) const
  {
    try
    {
      ASSERT(msc->is_paired(i) == false);

      if(m_rect.contains(msc->cellid(i)))
        dfs::do_dfs_collect_manifolds
            (shared_from_this(),mfold,msc->cellid(i),(eGDIR)dir);

      for( conn_iter_t j  = msc->m_conn[dir][i].begin();
                       j != msc->m_conn[dir][i].end();++j)
      {
        ASSERT(msc->is_paired(*j) == true);
        ASSERT(msc->index(i) == msc->index(msc->pair_idx(*j)));

        dfs::do_dfs_collect_manifolds
            (shared_from_this(),mfold,msc->cellid(msc->pair_idx(*j)),(eGDIR)dir);
      }
    }
    catch(assertion_error e)
    {
      e.push(_FFL);
      e.push(SVAR(m_rect));
      e.push(SVAR(m_ext_rect));

      throw;
    }
  }

  template<typename T>
  class producer_consumer_t:boost::noncopyable
  {
    queue<T>                      m_queue;
    boost::mutex                  m_mutex;
    boost::condition_variable_any m_cond;

  public:
    void put(const T & t)
    {
      boost::mutex::scoped_lock scoped_lock(m_mutex);

      m_queue.push(t);

      m_cond.notify_one();
    }

    T get()
    {
      boost::mutex::scoped_lock scoped_lock(m_mutex);

      while (m_queue.empty())
        m_cond.wait(m_mutex);

      T t = m_queue.front();

      m_queue.pop();

      return t;
    }
  };

  namespace save_mfolds
  {
    typedef cellid_list_t                                     mfold_t;
    typedef boost::shared_ptr<mfold_t>                        mfold_ptr_t;
    typedef boost::tuples::tuple<int,mfold_ptr_t,mfold_ptr_t> cp_no_mfolds_t;
    typedef producer_consumer_t<cp_no_mfolds_t>               mfolds_queue_t;
    typedef boost::shared_ptr<mfolds_queue_t>                 mfolds_queue_ptr_t;

    void process_mfold
        (dataset_const_ptr_t  ds,
         mscomplex_const_ptr_t msc,
         cp_producer_ptr_t prd,
         mfolds_queue_ptr_t mque)
    {
      using namespace boost::tuples;

      for(int i; prd->next(i) ;)
      {
        mfold_ptr_t des_mfold(new mfold_t);
        mfold_ptr_t asc_mfold(new mfold_t);

        ds->get_mfold<mfold_t>(des_mfold.get(),msc,i,GDIR_DES);
        ds->get_mfold<mfold_t>(asc_mfold.get(),msc,i,GDIR_ASC);

        mque->put(make_tuple(i,des_mfold,asc_mfold));
      }
    }

    void write_mfold
        (mfolds_queue_ptr_t mque,
         std::ostream & os,
         int_list_t & cp_order,
         int_list_t & mfold_offsets,
         int num_cps)
    {
      using boost::tuples::get;

      int offset = 0;

      for(int i = 0 ; i < num_cps; ++i)
      {
        cp_no_mfolds_t cp_no_mfold = mque->get();

        int cp_no = get<0>(cp_no_mfold);

        cp_order.push_back(cp_no);
        mfold_ptr_t des_mfold = get<1>(cp_no_mfold);
        mfold_ptr_t asc_mfold = get<2>(cp_no_mfold);

        os.write((char*)(void*)des_mfold->data(),des_mfold->size()*sizeof(cellid_t));
        os.write((char*)(void*)asc_mfold->data(),asc_mfold->size()*sizeof(cellid_t));

        mfold_offsets.push_back(offset); offset += des_mfold->size();
        mfold_offsets.push_back(offset); offset += asc_mfold->size();
      }

      mfold_offsets.push_back(offset);
    }

    int get_header_size(int num_cps)
    {
      return sizeof(rect_t)*3         + // rects
             sizeof(int)              + // num_cps
             sizeof(cellid_t)*num_cps + // cellids
             sizeof(int)*(2*num_cps+1); // offsets
    }

    template<typename T>
    void bin_write(std::ostream & os,const T & d)
    {
      os.write((const char*)(const void*)&d,sizeof(T));
    }

    void write_header(dataset_const_ptr_t ds,
                      mscomplex_const_ptr_t msc,
                      const int_list_t & cp_order,
                      const int_list_t & offsets,
                      std::ostream & os)
    {
      os.seekp(0,ios::beg);

      bin_write(os,ds->m_rect);
      bin_write(os,ds->m_ext_rect);
      bin_write(os,ds->m_domain_rect);

      bin_write(os,(int)cp_order.size());

      for( int i = 0 ; i < cp_order.size(); ++i)
        bin_write(os,msc->cellid(cp_order[i]));

      os.write((char*)(void*)offsets.data(),offsets.size()*sizeof(int));
    }


    void save(std::ostream & os,
              dataset_const_ptr_t ds,
              mscomplex_const_ptr_t msc)
    {
      cp_producer_ptr_t prd
          (new cp_producer_t(msc,cp_producer_t::unpaired_cp_filter));

      int num_cps = prd->count();

      os.seekp(get_header_size(num_cps),ios::beg);

      mfolds_queue_ptr_t mque(new mfolds_queue_t);

      int_list_t cp_order,offsets;

      boost::thread_group group;

      for(int tid = 0 ; tid < g_num_threads; ++tid)
        group.create_thread(bind(process_mfold,ds,msc,prd,mque));

      write_mfold(mque,os,cp_order,offsets,num_cps);

      write_header(ds,msc,cp_order,offsets,os);
    }
  }

  void  dataset_t::saveManifolds(mscomplex_const_ptr_t msc,const std::string &s)
  {
    std::ofstream fs(s.c_str());
    ensure(fs.is_open(),"unable to open file");
    save_mfolds::save(fs,shared_from_this(),msc);
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
