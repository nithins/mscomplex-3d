#include <stack>
#include <queue>
#include <list>

#include <fstream>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/typeof/typeof.hpp>

#include <logutil.h>

#include <grid_dataset.h>
#include <grid_dataset_ensure.h>
#include <grid_mscomplex.h>

namespace bl = boost::lambda;

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

  inline void pair_with_lowest_cofacet (dataset_t *dataset,cellid_t c)
  {
    cellid_t p;

    cellid_t cf[20];

    uint cf_count = dataset->getCellCofacets ( c,cf );

    bool p_usable = false;

    bool is_c_boundry = dataset->isTrueBoundryCell(c);

    for ( uint i =0 ; i < cf_count;i++ )
    {
      bool is_p_boundry = dataset->isTrueBoundryCell(cf[i]);

      if(dataset->getCellMaxFacetId(cf[i]) == c)
      {
        if((p_usable == false ||dataset->compareCells ( cf[i],p )) &&
           (is_c_boundry == is_p_boundry))
        {
          p_usable = true;
          p = cf[i];
        }
      }
    }

    if(p_usable)
      dataset->pairCells(c,p);

  }

  void track_gradient_tree_bfs
      (dataset_t *dataset,
       mscomplex_t *msgraph,
       cellid_t start_cellId,
       eGradientDirection gradient_dir
       )
  {
    std::queue<cellid_t> cell_queue;

    cell_queue.push ( start_cellId );

    while ( !cell_queue.empty() )
    {
      cellid_t top_cell = cell_queue.front();

      cell_queue.pop();

      cellid_t      cets[20];

      uint cet_ct = ( dataset->*getcets[gradient_dir] ) ( top_cell,cets );

      for ( uint i = 0 ; i < cet_ct ; i++ )
      {
        if ( dataset->isCellCritical ( cets[i] ) )
        {
          msgraph->connect_cps(start_cellId,cets[i]);
        }
        else
        {
          if ( !dataset->isCellExterior ( cets[i] ) )
          {
            cellid_t next_cell = dataset->getCellPairId ( cets[i] );

            if ( dataset->getCellDim ( top_cell ) ==
                 dataset->getCellDim ( next_cell ) &&
                 next_cell != top_cell )
            {
              cell_queue.push ( next_cell );
            }
          }
        }
      }
    }
  }

  void compute_disc_bfs
      (dataset_t *dataset,
       critpt_disc_t *disc,
       cellid_t start_cellId,
       eGradientDirection gradient_dir
       )
  {
    typedef cellid_t id_type;

    std::stack<id_type> cell_stack;

    cell_stack.push ( start_cellId );

    while ( !cell_stack.empty() )
    {
      id_type top_cell = cell_stack.top();

      cell_stack.pop();

      disc->push_back(top_cell);

      id_type cets[20];

      uint cet_ct = ( dataset->*getcets[gradient_dir] ) ( top_cell,cets );

      for ( uint i = 0 ; i < cet_ct ; i++ )
      {
        if ( !dataset->isCellCritical ( cets[i] ) )
        {
          if ( !dataset->isCellExterior ( cets[i] ) )
          {
            id_type next_cell = dataset->getCellPairId ( cets[i] );

            if ( dataset->getCellDim ( top_cell ) ==
                 dataset->getCellDim ( next_cell ) &&
                 next_cell != top_cell )
            {
              cell_stack.push ( next_cell );
            }
          }
        }
      }
    }
  }

  dataset_t::dataset_t (const rect_t &r,const rect_t &e) :
      m_rect (r),
      m_ext_rect (e),
      m_cell_flags(cellid_t::zero,boost::fortran_storage_order()),
      m_cell_pairs(cellid_t::zero,boost::fortran_storage_order()),
      m_cell_mxfct(cellid_t::zero,boost::fortran_storage_order()),
      m_cell_efdim_a(cellid_t::zero,boost::fortran_storage_order()),
      m_cell_efdim_d(cellid_t::zero,boost::fortran_storage_order())
  {
    // TODO: assert that the given rect is of even size..
    //       since each vertex is in the even positions
  }

  dataset_t::~dataset_t ()
  {
    clear();

    clear_fnref();
  }

  void dataset_t::init()
  {
    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    rect_size_t   s = m_ext_rect.size() + cellid_t::one;

    m_cell_flags.resize(s);
    m_cell_pairs.resize(s);
    m_cell_mxfct.resize(s);
    m_cell_efdim_a.resize(s);
    m_cell_efdim_d.resize(s);

    uint num_cells = s[0]*s[1]*s[2];

    std::fill_n(m_cell_flags.data(),num_cells,CELLFLAG_UNKNOWN);
    std::fill_n(m_cell_pairs.data(),num_cells,CELLADJDIR_UNKNOWN);
    std::fill_n(m_cell_mxfct.data(),num_cells,CELLADJDIR_UNKNOWN);

    rect_point_t bl = m_ext_rect.lower_corner();

    m_cell_flags.reindex (bl);
    m_cell_pairs.reindex (bl);
    m_cell_mxfct.reindex (bl);
    m_cell_efdim_a.reindex (bl);
    m_cell_efdim_d.reindex (bl);
  }

  void  dataset_t::clear()
  {
    m_critical_cells.clear();

    m_cell_flags.resize(cellid_t::zero);
    m_cell_pairs.resize(cellid_t::zero);
    m_cell_mxfct.resize(cellid_t::zero);
    m_cell_efdim_a.resize(cellid_t::zero);
    m_cell_efdim_d.resize(cellid_t::zero);
  }

  void dataset_t::init_fnref(cell_fn_t * pData)
  {
    rect_size_t   s = m_ext_rect.size();

    rect_point_t bl = m_ext_rect.lower_corner();

    if(pData != NULL)
    {
      m_vert_fns_ref.reset
          (new varray_ref_t(pData,boost::extents[1+s[0]/2][1+s[1]/2][1+s[2]/2],
                           boost::fortran_storage_order()));

      (*m_vert_fns_ref).reindex (bl/2);
    }
  }

  void dataset_t::clear_fnref()
  {
    m_vert_fns_ref.reset();
  }

  inline cellid_t get_adj_cell(cellid_t c,uint d)
  {
    c[(d-1)>>1] += (d&1)?(-1):(+1);
    return c;
  }

  inline uint get_cell_adj_dir(cellid_t c,cellid_t p)
  {
    ensure_cell_incidence(c,p);

    for(uint i = 0 ;i < gc_grid_dim;++i)
      if(c[i] != p[i])
        return (1+i*2+(((c[i]-p[i])==1)?(0):(1)));

    ensure_never_reached("cell coords are the same");

    return (uint)-1;
  }

  bool   dataset_t::areCellsIncident(cellid_t c1,cellid_t c2) const
  {
    uint dim_diff;

    for(uint i = 0 ;i < gc_grid_dim;++i)
      dim_diff += (c1[i] >c2[i])?(c1[i] - c2[i]):(c2[i]  - c1[i]);

    return (dim_diff == 1);
  }

  cellid_t dataset_t::getCellPairId (cellid_t c) const
  {
    ensure_cell_paired(this,c);

    return get_adj_cell(c,m_cell_pairs(c));
  }

  cellid_t dataset_t::getCellMaxFacetId (cellid_t c) const
  {
    ensure_cell_max_facet_known(this,c);

    return get_adj_cell(c,m_cell_mxfct(c));
  }

  cellid_t dataset_t::getCellSecondMaxFacetId (cellid_t c) const
  {
    ensure_cell_max_facet_known(this,c);

    uint mxfct = m_cell_mxfct(c);
    return get_adj_cell(c,(mxfct&1)?(mxfct+1):(mxfct-1));
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

    bool is_boundry_c1 = ds->isTrueBoundryCell(c1);
    bool is_boundry_c2 = ds->isTrueBoundryCell(c2);

    if(is_boundry_c1 != is_boundry_c2)
      return is_boundry_c1;

    f1 = ds->getCellSecondMaxFacetId(c1);
    f2 = ds->getCellSecondMaxFacetId(c2);

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
    cellid_list_t cofaces(40);

    cofaces.resize(getCellCofaces(c,cofaces.data()));

    uint pos = 0;

    for(uint i = 0 ; i< cofaces.size();++i)
    {
      cellid_t cf = cofaces[i];

      uint c_dim = getCellDim(c);
      uint cf_dim = getCellDim(cf);

      for(uint j = c_dim; j < cf_dim ; ++j)
        cf = getCellMaxFacetId(cf);

      if(cf == c)
        est[pos++] = cofaces[i];
    }

    return pos;
  }

  bool dataset_t::isPairOrientationCorrect (cellid_t c, cellid_t p) const
  {
    return (getCellDim (c) <getCellDim (p));
  }

  bool dataset_t::isCellMarked (cellid_t c) const
  {
    return ! (m_cell_flags (c) == CELLFLAG_UNKNOWN);
  }

  bool dataset_t::isCellCritical (cellid_t c) const
  {
    return (m_cell_flags (c) & CELLFLAG_CRITICAL);
  }

  bool dataset_t::isCellPaired (cellid_t c) const
  {
    return (m_cell_flags (c) & CELLFLAG_PAIRED);
  }

  void dataset_t::pairCells (cellid_t c,cellid_t p)
  {
    ensure_pairable(this,c,p);

    m_cell_pairs (c) = get_cell_adj_dir(c,p);
    m_cell_pairs (p) = get_cell_adj_dir(p,c);

    m_cell_flags (c) |= CELLFLAG_PAIRED;
    m_cell_flags (p) |= CELLFLAG_PAIRED;

    ensure_valid_pair(this,p,c);
  }

  void  dataset_t::unpairCells ( cellid_t c,cellid_t p )
  {
    ensure_valid_pair(this,p,c);

    m_cell_pairs (c) = CELLADJDIR_UNKNOWN;
    m_cell_pairs (p) = CELLADJDIR_UNKNOWN;

    m_cell_flags (c) &= (CELLFLAG_MASK^CELLFLAG_PAIRED);
    m_cell_flags (p) &= (CELLFLAG_MASK^CELLFLAG_PAIRED);

    ensure_cell_not_paired(this,c);
    ensure_cell_not_paired(this,p);
  }

  void dataset_t::setCellMaxFacet (cellid_t c,cellid_t f)
  {
    ensure_cell_dim(f,getCellDim(c)-1);

    m_cell_mxfct (c) = get_cell_adj_dir(c,f);
  }

  void dataset_t::markCellCritical (cellid_t c)
  {
    m_cell_flags (c) |= CELLFLAG_CRITICAL;
  }

  bool dataset_t::isTrueBoundryCell (cellid_t c) const
  {
    return (m_ext_rect.isOnBoundry (c));
  }

  bool dataset_t::isFakeBoundryCell (cellid_t c) const
  {
    return (m_rect.isOnBoundry (c) && (!m_ext_rect.isOnBoundry (c)));
  }

  bool dataset_t::isCellExterior (cellid_t c) const
  {
    return (!m_rect.contains (c) && m_ext_rect.contains (c));
  }

  void dataset_t::work()
  {
    assignMaxFacets();

    pairCellsWithinEst();

    collateCriticalPoints();

//    aggregateEffCellDim();
  }

  void dataset_t::assignMaxFacets()
  {
    using namespace boost::lambda;

    BOOST_AUTO(cmp,bind(&dataset_t::compareCells,this,_1,_2));

    cellid_t f[20],c;

    for(uint dim = 1 ; dim <= gc_grid_dim; ++dim)
    {
      cellid_t   seq(cellid_t::zero);

      std::for_each(seq.begin(),seq.begin()+dim,_1=1);

      do
      {
        static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

        for(c[2] = m_rect[2][0] + seq[2] ; c[2] <= m_rect[2][1]; c[2] += 2)
        {
          for(c[1] = m_rect[1][0] + seq[1]; c[1] <= m_rect[1][1]; c[1] += 2)
          {
            for(c[0] = m_rect[0][0] + seq[0]; c[0] <= m_rect[0][1]; c[0] += 2)
            {
              int f_ct = getCellFacets(c,f);

              setCellMaxFacet(c,*std::max_element(f,f+f_ct,cmp));

            }
          }
        }
      }
      while(std::next_permutation(seq.rbegin(),seq.rend()));
    }
  }

  void  dataset_t::assignGradients()
  {
    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    cellid_t c;

    for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; ++c[2])
    {
      for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; ++c[1])
      {
        for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; ++c[0])
        {
          pair_with_lowest_cofacet(this,c);
        }
      }
    }
  }

  void  dataset_t::pairCellsWithinEst()
  {
    using namespace boost::lambda;

    typedef std::list<cellid_t> cellid_llist_t;

    BOOST_AUTO(cmp,bind(&dataset_t::compareCells,this,_1,_2));

    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    cellid_t c;

    for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; c[2] += 2)
    {
      for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; c[1] += 2)
      {
        for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; c[0] += 2)
        {
          cellid_t est_arr[40];

          uint est_ct = getCellEst(c,est_arr);

          std::sort(est_arr,est_arr+est_ct,cmp);

          cellid_llist_t est(est_arr,est_arr+est_ct);

          bool bCanPair = true;

          int num_rounds = gc_grid_dim;

          while(bCanPair && num_rounds-- != 0)
          {
            bCanPair = false;

            for(cellid_llist_t::iterator est_it = est.begin();
                est_it != est.end() && ++est_it != est.end();)
            {
              cellid_t p = *est_it,c = *--est_it;

              if(getCellDim(c)+1      == getCellDim(p)  &&
                 isTrueBoundryCell(c) == isTrueBoundryCell(p) &&
                 areCellsIncident(c,p))
              {
                pairCells(c,p);

                est.erase(est_it++);
                est.erase(est_it++);

                bCanPair = true;
              }
              else
                ++est_it;
            }
          }
        }
      }
    }
  }

  void  dataset_t::collateCriticalPoints()
  {
    static_assert(gc_grid_dim == 3&&"defined for 3-manifolds only");

    cellid_t c;

    for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; ++c[2])
    {
      for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; ++c[1])
      {
        for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; ++c[0])
        {
          if (!isCellMarked (c) )
          {
            markCellCritical(c);

            m_critical_cells.push_back(c);
          }
        }
      }
    }
  }

  void  dataset_t::aggregateEffCellDim()
  {
    using namespace boost::lambda;

    BOOST_AUTO(cmp,bind(&dataset_t::compareCells,this,_1,_2));

    cellid_list_t c_list;

    cellid_t c,p;

    static_assert(gc_grid_dim == 3&&"defined for 3-manifolds only");

    for(c[2] = m_ext_rect[2][0] ; c[2] <= m_ext_rect[2][1]; ++c[2])
    {
      for(c[1] = m_ext_rect[1][0] ; c[1] <= m_ext_rect[1][1]; ++c[1])
      {
        for(c[0] = m_ext_rect[0][0] ; c[0] <= m_ext_rect[0][1]; ++c[0])
        {
          m_cell_efdim_a(c) = getCellDim(c);

          m_cell_efdim_d(c) = gc_grid_dim - getCellDim(c);

          c_list.push_back(c);
        }
      }
    }

    std::sort(c_list.begin(),c_list.end(),cmp);

    cellid_t fcts[20];

    uint fcts_eff_dim[20];

    for(int i = 0 ; i< c_list.size();++i)
    {
      c = c_list[i];

      if(isCellCritical(c))
        continue;

      p = getCellPairId(c);

      if(getCellDim(c) < getCellDim(p))
        continue;

      uint fct_ct = getCellFacets(c,fcts);

      for(uint j = 0; j< fct_ct; ++j)
      {
        if(fcts[j] == p)
          fcts_eff_dim[j] = 0;
        else
          fcts_eff_dim[j] = m_cell_efdim_a(fcts[j]);
      }

      uint d = *std::max_element(fcts_eff_dim,fcts_eff_dim + fct_ct);

      m_cell_efdim_a(c) = d;
      m_cell_efdim_a(p) = d;
    }

    cellid_t cofcts[20];

    uint cofcts_eff_dim[20];

    for(int i = c_list.size() ; i >0 ;--i)
    {
      c = c_list[i-1];

      if(isCellCritical(c))
        continue;

      p = getCellPairId(c);

      if(getCellDim(c) > getCellDim(p))
        continue;

      uint cofct_ct = getCellCofacets(c,cofcts);

      for(uint j = 0; j< cofct_ct; ++j)
      {
        if(cofcts[j] == p)
          cofcts_eff_dim[j] = 0;
        else
          cofcts_eff_dim[j] = m_cell_efdim_d(cofcts[j]);
      }

      uint d = *std::max_element(cofcts_eff_dim,cofcts_eff_dim + cofct_ct);

      m_cell_efdim_d(c) = d;
      m_cell_efdim_d(p) = d;
    }
  }

  void  dataset_t::writeout_connectivity(mscomplex_t *msgraph)
  {
    for(uint i = 0 ; i < m_critical_cells.size();++i)
    {
      cellid_t c = m_critical_cells[i];
      msgraph->add_critpt(c,getCellDim(c),get_cell_fn(c));
    }

    for(uint i = 0 ; i < m_critical_cells.size();++i)
    {
      std::cout<<i<<" of "<<m_critical_cells.size()<<" c ="<<m_critical_cells[i]<<"\n";

      track_gradient_tree_bfs
          (this,msgraph,m_critical_cells[i],GRADDIR_DESCENDING);
    }

  }

  void dataset_t::postMergeFillDiscs(mscomplex_t *msgraph)
  {
    msgraph->add_disc_tracking_seed_cps();

    for(uint i = 0 ; i < msgraph->m_cps.size() ; ++i)
    {
      critpt_t * cp = msgraph->m_cps[i];

      for(uint dir = 0 ; dir < GRADDIR_COUNT;++dir)
      {
        if(cp->disc[dir].size() == 1)
        {
          cp->disc[dir].clear();
          compute_disc_bfs(this,&cp->disc[dir],cp->cellid,(eGradientDirection)dir);
        }
      }
    }
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

  void dataset_t::log_pairs()
  {
    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    cellid_t c;

    for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; ++c[2])
    {
      for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; ++c[1])
      {
        for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; ++c[0])
        {
          std::cout<<getCellPairId(c)<<" ";
        }
        std::cout<<std::endl;
      }
      std::cout<<std::endl;
    }
  }

  void dataset_t::log_eff_dim()
  {
    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    cellid_t c;

    for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; ++c[2])
    {
      for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; ++c[1])
      {
        for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; ++c[0])
        {
          std::cout<<(int)m_cell_efdim_a(c)<<(int)m_cell_efdim_d(c)<<" ";
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
