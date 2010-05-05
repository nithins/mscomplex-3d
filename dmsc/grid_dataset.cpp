#include <queue>

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
  cellid_t get_cp_cellid(mscomplex_t *msgraph,uint idx)
  {
    return msgraph->m_cps[idx]->cellid;
  }

  static uint ( dataset_t::*getcets[2] ) ( cellid_t,cellid_t * ) const =
  {
    &dataset_t::getCellFacets,
    &dataset_t::getCellCofacets
  };

  inline uint   dataset_t::getCellIncCells( cellid_t c,cellid_t * inc) const
  {
    for(uint i = 0; i < cellid_t::base_t::static_size; ++i)
    {
      inc[i*2+0] = c;
      inc[i*2+1] = c;

      inc[i*2+0][i] -= 1;
      inc[i*2+0][i] += 1;
    }
    return cellid_t::base_t::static_size*2;
  }

  inline bool lowestPairableCoFacet
      (dataset_t *dataset,
       cellid_t cellId,
       cellid_t& pairid)
  {
    cellid_t cofacets[20];
    bool    cofacet_usable[20];

    uint cofacet_count = dataset->getCellCofacets ( cellId,cofacets );

    bool isTrueBoundryCell = dataset->isTrueBoundryCell ( cellId ) ;

    // for each co facet
    for ( uint i = 0 ; i < cofacet_count ; i++ )
    {
      cofacet_usable[i] = false;

      if ( isTrueBoundryCell &&
           !dataset->isTrueBoundryCell ( cofacets[i] ) )
        continue;

      if(dataset->getCellMaxFacetId(cofacets[i]) == cellId)
        cofacet_usable[i] = true;
    }

    bool pairid_usable = false;

    for ( uint i =0 ; i < cofacet_count;i++ )
    {
      if ( cofacet_usable[i] == false )
        continue;

      if(pairid_usable == false)
      {
        pairid_usable = true;
        pairid = cofacets[i];
        continue;
      }

      if ( dataset->compareCells ( cofacets[i],pairid ) )
        pairid = cofacets[i];

    }
    return pairid_usable;
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

    std::queue<id_type> cell_queue;

    cell_queue.push ( start_cellId );

    while ( !cell_queue.empty() )
    {
      id_type top_cell = cell_queue.front();

      cell_queue.pop();

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
              cell_queue.push ( next_cell );
            }
          }
        }
      }
    }
  }

  dataset_t::dataset_t (const rect_t &r,const rect_t &e) :
      m_rect (r),m_ext_rect (e)
  {

    // TODO: assert that the given rect is of even size..
    //       since each vertex is in the even positions
    //

    m_vert_fns_ref = NULL;
    m_cell_flags   = NULL;
    m_cell_pairs   = NULL;
    m_cell_mxfct   = NULL;

  }

  dataset_t::dataset_t ()
  {
    m_vert_fns_ref = NULL;
    m_cell_flags   = NULL;
    m_cell_pairs   = NULL;
    m_cell_mxfct   = NULL;

  }

  dataset_t::~dataset_t ()
  {
    clear();

    clear_fnref();
  }

  void dataset_t::init()
  {

    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    rect_size_t   s = m_ext_rect.size();

    m_cell_flags = new cellflag_array_t(boost::extents[1+s[0]][1+s[1]][1+s[2]],
                                        boost::fortran_storage_order());
    m_cell_pairs = new cellflag_array_t(boost::extents[1+s[0]][1+s[1]][1+s[2]],
                                        boost::fortran_storage_order());
    m_cell_mxfct = new cellflag_array_t(boost::extents[1+s[0]][1+s[1]][1+s[2]],
                                        boost::fortran_storage_order());

    cellid_t c;

    for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; ++c[2])
    {
      for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; ++c[1])
      {
        for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; ++c[0])
        {
          (*m_cell_flags)(c) = CELLFLAG_UNKNOWN;
          (*m_cell_pairs)(c) = CELLADJDIR_UNKNOWN;
          (*m_cell_mxfct)(c) = CELLADJDIR_UNKNOWN;
        }
      }
    }


    rect_point_t bl = m_ext_rect.lower_corner();

    (*m_cell_flags).reindex (bl);
    (*m_cell_pairs).reindex (bl);
    (*m_cell_mxfct).reindex (bl);
  }

  void  dataset_t::clear()
  {
    if(m_cell_flags != NULL)
      delete m_cell_flags;

    if(m_cell_pairs != NULL)
      delete m_cell_pairs;

    if(m_cell_mxfct != NULL)
      delete m_cell_mxfct;

    m_critical_cells.clear();

    m_cell_flags   = NULL;
    m_cell_pairs   = NULL;
    m_cell_mxfct   = NULL;

  }

  void dataset_t::init_fnref(cell_fn_t * pData)
  {
    rect_size_t   s = m_ext_rect.size();

    if(pData != NULL)
      m_vert_fns_ref =
          new varray_ref_t(pData,boost::extents[1+s[0]/2][1+s[1]/2][1+s[2]/2],
                           boost::fortran_storage_order());

    rect_point_t bl = m_ext_rect.lower_corner();

    if(pData != NULL)
      (*m_vert_fns_ref).reindex (bl/2);

  }

  void dataset_t::clear_fnref()
  {
    if(m_vert_fns_ref != NULL)
      delete m_vert_fns_ref;

    m_vert_fns_ref = NULL;
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

  cellid_t dataset_t::getCellPairId (cellid_t c) const
  {
    ensure_cell_paired(this,c);

    return get_adj_cell(c,(*m_cell_pairs)(c));
  }

  cellid_t dataset_t::getCellMaxFacetId (cellid_t c) const
  {
    ensure_cell_max_facet_known(this,c);

    return get_adj_cell(c,(*m_cell_mxfct)(c));
  }

  cellid_t dataset_t::getCellSecondMaxFacetId (cellid_t c) const
  {
    ensure_cell_max_facet_known(this,c);

    uint mxfct = (*m_cell_mxfct)(c);
    return get_adj_cell(c,(mxfct&1)?(mxfct+1):(mxfct-1));
  }

  bool dataset_t::compareCells( cellid_t c1,cellid_t  c2 ) const
  {
    ensure_cell_dim(this,c1,getCellDim(c2));

    if(c1 == c2)
      return false;

    if(getCellDim(c1) == 0)
      return ptLt(c1,c2);

    cellid_t f1 = getCellMaxFacetId(c1);
    cellid_t f2 = getCellMaxFacetId(c2);


    if(f1 == f2)
    {
      f1 = getCellSecondMaxFacetId(c1);
      f2 = getCellSecondMaxFacetId(c2);
    }

    return compareCells(f1,f2);
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

  bool dataset_t::isPairOrientationCorrect (cellid_t c, cellid_t p) const
  {
    return (getCellDim (c) <getCellDim (p));
  }

  bool dataset_t::isCellMarked (cellid_t c) const
  {
    return ! ((*m_cell_flags) (c) == CELLFLAG_UNKNOWN);
  }

  bool dataset_t::isCellCritical (cellid_t c) const
  {
    return ((*m_cell_flags) (c) & CELLFLAG_CRITCAL);
  }

  bool dataset_t::isCellPaired (cellid_t c) const
  {
    return ((*m_cell_flags) (c) & CELLFLAG_PAIRED);
  }

  void dataset_t::pairCells (cellid_t c,cellid_t p)
  {
    (*m_cell_pairs)(c) = get_cell_adj_dir(c,p);
    (*m_cell_pairs)(p) = get_cell_adj_dir(p,c);

    (*m_cell_flags) (c) |= CELLFLAG_PAIRED;
    (*m_cell_flags) (p) |= CELLFLAG_PAIRED;
  }

  void dataset_t::setCellMaxFacet (cellid_t c,cellid_t f)
  {
    ensure_cell_dim(this,f,getCellDim(c)-1);

    (*m_cell_mxfct)(c) = get_cell_adj_dir(c,f);
  }

  void dataset_t::markCellCritical (cellid_t c)
  {
    (*m_cell_flags) (c) |= CELLFLAG_CRITICAL;
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

    assignGradients();

    collateCriticalPoints();
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

    cellid_t c,p;

    for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; ++c[2])
    {
      for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; ++c[1])
      {
        for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; ++c[0])
        {
          if (lowestPairableCoFacet (this,c,p))
          {
            ensure_cell_not_marked(this,c);

            ensure_cell_not_marked(this,p);

            pairCells (c,p);
          }
        }
      }
    }
#warning "havent implemnted mark boundry critpts yet"
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

  void  dataset_t::writeout_connectivity(mscomplex_t *msgraph)
  {
    for(uint i = 0 ; i < m_critical_cells.size();++i)
    {
      cellid_t c = m_critical_cells[i];
      msgraph->add_critpt(c,getCellDim(c),get_cell_fn(c));
    }

    for(uint i = 0 ; i < m_critical_cells.size();++i)
      track_gradient_tree_bfs
          (this,msgraph,m_critical_cells[i],GRADDIR_DESCENDING);

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
          std::cout<<(*m_cell_flags)(c)<<" ";
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
}
