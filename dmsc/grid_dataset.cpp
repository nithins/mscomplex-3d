#include <queue>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include <logutil.h>

#include <grid_dataset.h>
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
       cellid_t& pairid
       )
  {
    typedef cellid_t id_type;

    id_type cofacets[20];
    bool    cofacet_usable[20];

    uint cofacet_count = dataset->getCellCofacets ( cellId,cofacets );

    bool isTrueBoundryCell = dataset->isTrueBoundryCell ( cellId ) ;

    // for each co facet
    for ( uint i = 0 ; i < cofacet_count ; i++ )
    {
      id_type facets[20];
      uint facet_count = dataset->getCellFacets ( cofacets[i],facets );

      cofacet_usable[i] = true;

      if ( isTrueBoundryCell &&
           !dataset->isTrueBoundryCell ( cofacets[i] ) )
      {
        cofacet_usable[i] = false;
        continue;
      }

      for ( uint j = 0 ; j < facet_count ; j++ )
      {
        if ( dataset->compareCells ( cellId,facets[j] ))
        {
          cofacet_usable[i] = false;
          break;
        }
      }
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
       eDirection gradient_dir
       )
  {
    typedef cellid_t id_type;

    std::queue<id_type> cell_queue;

    cell_queue.push ( start_cellId );

    while ( !cell_queue.empty() )
    {
      id_type top_cell = cell_queue.front();

      cell_queue.pop();

      id_type      cets[20];

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

  void compute_disc_bfs
      (dataset_t *dataset,
       critpt_disc_t *disc,
       cellid_t start_cellId,
       eDirection gradient_dir
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
      m_rect (r),m_ext_rect (e),m_ptcomp(this),
      m_vert_fns_ref(NULL)
  {

    // TODO: assert that the given rect is of even size..
    //       since each vertex is in the even positions
    //

    m_vert_fns_ref = NULL;
    m_cell_flags   = NULL;
    m_cell_pairs   = NULL;

  }

  dataset_t::dataset_t () :
      m_ptcomp(this)
  {
    m_vert_fns_ref = NULL;
    m_cell_flags   = NULL;
    m_cell_pairs   = NULL;
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

    m_cell_flags = new cellflag_array_t(boost::extents[1+s[0]][1+s[1]][1+s[2]],boost::fortran_storage_order());
    m_cell_pairs = new cellpair_array_t(boost::extents[1+s[0]][1+s[1]][1+s[2]],boost::fortran_storage_order());

    cellid_t c;

    for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; ++c[2])
      for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; ++c[1])
        for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; ++c[0])
          (*m_cell_flags)(c) = CELLFLAG_UNKNOWN;


    rect_point_t bl = m_ext_rect.lower_corner();

    (*m_cell_flags).reindex (bl);
    (*m_cell_pairs).reindex (bl);
  }

  void  dataset_t::clear()
  {
    if(m_cell_flags != NULL)
      delete m_cell_flags;

    if(m_cell_pairs != NULL)
      delete m_cell_pairs;

    m_critical_cells.clear();

    m_cell_flags   = NULL;
    m_cell_pairs   = NULL;
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

  cellid_t   dataset_t::getCellPairId (cellid_t c) const
  {
    if ((*m_cell_flags) (c) &CELLFLAG_PAIRED == 0)
      throw std::logic_error ("invalid pair requested");

    return (*m_cell_pairs) (c);
  }

  bool dataset_t::compareCells( cellid_t c1,cellid_t  c2 ) const
  {
    if(getCellDim(c1) == 0)
      return ptLt(c1,c2);

    cellid_t pts1[20];
    cellid_t pts2[20];

    uint pts1_ct = getCellPoints ( c1,pts1);
    uint pts2_ct = getCellPoints ( c2,pts2);

    std::sort ( pts1,pts1+pts1_ct,m_ptcomp );
    std::sort ( pts2,pts2+pts2_ct,m_ptcomp);

    return std::lexicographical_compare
        ( pts1,pts1+pts1_ct,pts2,pts2+pts2_ct,
          m_ptcomp );
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
    (*m_cell_pairs) (c) = p;
    (*m_cell_pairs) (p) = c;

    (*m_cell_flags) (c) = (*m_cell_flags) (c) |CELLFLAG_PAIRED;
    (*m_cell_flags) (p) = (*m_cell_flags) (p) |CELLFLAG_PAIRED;
  }

  void dataset_t::markCellCritical (cellid_t c)
  {
    (*m_cell_flags) (c) = (*m_cell_flags) (c) |CELLFLAG_CRITCAL;
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
    assignGradients();

    collateCriticalPoints();
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
          if (isCellMarked (c))
            continue;

          if (lowestPairableCoFacet (this,c,p))
            pairCells (c,p);
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
          (this,msgraph,m_critical_cells[i],DIRECTION_DESCENDING);

  }

  void dataset_t::postMergeFillDiscs(mscomplex_t *msgraph)
  {
    msgraph->add_disc_tracking_seed_cps();

    for(uint i = 0 ; i < msgraph->m_cps.size() ; ++i)
    {
      critpt_t * cp = msgraph->m_cps[i];

      for(uint dir = 0 ; dir < DIRECTION_COUNT;++dir)
      {
        if(cp->disc[dir].size() == 1)
        {
          cp->disc[dir].clear();
          compute_disc_bfs(this,&cp->disc[dir],cp->cellid,(eDirection)dir);
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
          std::cout<<(*m_cell_pairs)(c)<<" ";
        }
        std::cout<<std::endl;
      }
      std::cout<<std::endl;
    }
  }
}
