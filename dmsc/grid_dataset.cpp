#include <grid_dataset.h>

#include <discreteMorseAlgorithm.h>
#include <vector>

#include <timer.h>

#include <QFile>
#include <prefix_scan.h>
#include <bitonic_sort.h>

typedef GridDataset::cellid_t cellid_t;

cellid_t get_cp_cellid(GridDataset::mscomplex_t *msgraph,uint idx)
{
  return msgraph->m_cps[idx]->cellid;
}

static uint ( GridDataset::*getcets[2] ) ( cellid_t,cellid_t * ) const =
{
  &GridDataset::getCellFacets,
  &GridDataset::getCellCofacets
};

inline uint   GridDataset::getCellIncCells( cellid_t c,cellid_t * inc) const
{
  inc[0] = cellid_t (c[0]  ,c[1]+1);
  inc[1] = cellid_t (c[0]  ,c[1]-1);
  inc[2] = cellid_t (c[0]-1,c[1]);
  inc[3] = cellid_t (c[0]+1,c[1]);
  return 4;
}

int GridDataset::postMergeFillDiscs(mscomplex_t *msgraph)
{

  std::map<cellid_t,critpt_disc_t *> surv_crit_asc_disc_map;
  std::map<cellid_t,critpt_disc_t *> surv_crit_des_disc_map;
  std::vector<uint> surv_saddle_idxs;

  // update local copy of owner extrema with the cancllation info
  for(uint i = 0 ; i < msgraph->m_cps.size();++i)
  {
    critpt_t * cp = msgraph->m_cps[i];

    int cp_dim = getCellDim(cp->cellid);

    if(cp->isBoundryCancelable == true)
    {
      critpt_t * pair_cp = msgraph->m_cps[cp->pair_idx];

      if( cp_dim == 1)
      {
        critpt_conn_t * cp_conn = (getCellDim(pair_cp->cellid) == 0)? (&cp->des):(&cp->asc);

        if(cp_conn->size() != 1)
        {
          log_range(cp_conn->begin(),cp_conn->end(),boost::bind(&get_cp_cellid,msgraph,_1),"conns");

          throw std::logic_error("I should be connected to exactly one surv extrema");
        }

        critpt_t * extrema_cp =msgraph->m_cps[*cp_conn->begin()];

        if(pair_cp->pair_idx != i)
          throw std::logic_error("My pair does not know Im paired with him");

        if(m_rect.contains(cp->cellid))
          (*m_cell_own)(cp->cellid) = extrema_cp->cellid;

        if(m_rect.contains(pair_cp->cellid))
          (*m_cell_own)(pair_cp->cellid) = extrema_cp->cellid;
      }
      else if(cp_dim == 2)
      {
        for(critpt_conn_t::iterator it = cp->des.begin() ;it!= cp->des.end();++it)
        {
          critpt_t * conn_saddle_cp  = msgraph->m_cps[*it];

          if(m_rect.contains(pair_cp->cellid))
            conn_saddle_cp->asc_disc.push_back(pair_cp->cellid);
        }
      }
      else if(cp_dim == 0)
      {
        for(critpt_conn_t::iterator it = cp->asc.begin() ;it!= cp->asc.end();++it)
        {
          critpt_t * conn_saddle_cp  = msgraph->m_cps[*it];

          if(m_rect.contains(pair_cp->cellid))
            conn_saddle_cp->des_disc.push_back(pair_cp->cellid);
        }
      }
    }
    else// not boundry cancellabele
    {
      switch(cp_dim)
      {
      case 0 :
        surv_crit_asc_disc_map.insert(std::make_pair(cp->cellid,&cp->asc_disc));
        break;
      case 1:
        {
          cellid_t inc_cells[4];

          getCellIncCells(cp->cellid,inc_cells);

          for(uint j = 0 ; j < 4 ; ++j)
          {
            if(m_rect.contains(inc_cells[j])
              && isCellPaired(inc_cells[j])
              && !isCellCritical(inc_cells[j])
              )
            {
              critpt_disc_t * disc =
                  (getCellDim(inc_cells[j]) ==0)?(&cp->des_disc):(&cp->asc_disc);

              disc->push_back(getCellPairId(inc_cells[j]));
            }
          }
          surv_saddle_idxs.push_back(i);

          break;
        }

      case 2:
        surv_crit_des_disc_map.insert(std::make_pair(cp->cellid,&cp->des_disc));
      }
    }
  }

  // all 0 d cells are owned by some minima
  for (cell_coord_t y = m_rect.bottom(); ;y += 2)
  {
    for (cell_coord_t x = m_rect.left(); ;x += 2)
    {
      cellid_t c (x,y);

      cellid_t o = (*m_cell_own)(c);

      if(m_rect.contains(o))
      {
        o = (*m_cell_own)(o);
      }

      if(o[0] != -1 && o[1] != -1)
      {
        surv_crit_asc_disc_map[o]->push_back(c);
      }

      if(x == m_rect.right()) break;
    }
    if(y == m_rect.top()) break;
  }

  // all 2 d cells are owned by some maxima
  for (cell_coord_t y = m_rect.bottom()+1;;y += 2)
  {
    for (cell_coord_t x = m_rect.left()+1;;x += 2)
    {
      cellid_t c (x,y);

      cellid_t o = (*m_cell_own)(c);

      if(m_rect.contains(o))
      {
        o = (*m_cell_own)(o);
      }

      if(o[0] != -1 && o[1] != -1)
      {
        surv_crit_des_disc_map[o]->push_back(c);
      }
      if(x + 1 == m_rect.right()) break;
    }
    if(y + 1== m_rect.top()) break;
  }
  // all surv saddles must now contain seed points to track their 1 manifold

  for(uint i = 0 ;i < surv_saddle_idxs.size();++i)
  {
    critpt_t * cp = msgraph->m_cps[surv_saddle_idxs[i]];

    critpt_disc_t * disc[] = {&cp->asc_disc,&cp->des_disc};

    for(uint j = 0 ; j <2;++j)
    {
      uint path_cell_idx = 0;

      while(path_cell_idx != disc[j]->size())
      {
        cellid_t path_cell = (*disc[j])[path_cell_idx];

        cellid_t cets[2];

        uint cet_ct = ( this->*getcets[(j+1)%2] )(path_cell,cets);

        for(uint k = 0 ; k < cet_ct;++k)
        {
          if(isCellCritical(cets[k]))
            continue;

          if(m_rect.contains(cets[k]) == false)
            continue;

          cellid_t p = getCellPairId(cets[k]);

          if( p== path_cell || m_rect.contains(cets[k]) == false)
            continue;

          disc[j]->push_back(p);
        }

        path_cell_idx++;
      }
    }

  }

  return 0;
}

void connectCps (GridDataset::mscomplex_t *msgraph,
                 uint cp1_ind,
                 uint cp2_ind)
{
  if(GridDataset::s_getCellDim(msgraph->m_cps[cp1_ind]->cellid) <
     GridDataset::s_getCellDim(msgraph->m_cps[cp2_ind]->cellid))
    std::swap(cp1_ind,cp2_ind);

  GridDataset::critpt_t *cp1 = msgraph->m_cps[cp1_ind];

  GridDataset::critpt_t *cp2 = msgraph->m_cps[cp2_ind];

  cp1->des.insert (cp2_ind);

  cp2->asc.insert (cp1_ind);
}

void connectCps (GridDataset::mscomplex_t *msgraph,
                 GridDataset::cellid_t c1,
                 GridDataset::cellid_t c2)
{
  if (GridDataset::s_getCellDim (c1) <GridDataset::s_getCellDim (c2))
    std::swap (c1,c2);

  if (GridDataset::s_getCellDim (c1) != GridDataset::s_getCellDim (c2) +1)
    throw std::logic_error ("must connect i,i+1 cp (or vice versa)");

  if (msgraph->m_id_cp_map.find (c1) == msgraph->m_id_cp_map.end())
    throw std::logic_error (_SSTR ("cell not in id_cp_map c1="<<c1));

  if (msgraph->m_id_cp_map.find (c2) == msgraph->m_id_cp_map.end())
    throw std::logic_error (_SSTR ("cell not in id_cp_map c2="<<c2));

  uint cp1_ind = msgraph->m_id_cp_map[c1];

  uint cp2_ind = msgraph->m_id_cp_map[c2];

  GridDataset::critpt_t *cp1 = msgraph->m_cps[cp1_ind];

  GridDataset::critpt_t *cp2 = msgraph->m_cps[cp2_ind];

  cp1->des.insert (cp2_ind);

  cp2->asc.insert (cp1_ind);
}

inline bool lowestPairableCoFacet
    (GridDataset *dataset,
     GridDataset::cellid_t cellId,
     GridDataset::cellid_t& pairid
     )
{
  typedef GridDataset::cellid_t id_type;

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
    (GridDataset *dataset,
     GridDataset::cellid_t start_cellId,
     eGradDirection gradient_dir
     )
{
  typedef GridDataset::cellid_t id_type;

  std::queue<id_type> cell_queue;

  // mark here that that cellid has no parent.

  cell_queue.push ( start_cellId );

  while ( !cell_queue.empty() )
  {
    id_type top_cell = cell_queue.front();

    cell_queue.pop();

    (*dataset->m_cell_own)(top_cell) = start_cellId;

    id_type      cets[20];

    uint cet_ct = ( dataset->*getcets[gradient_dir] ) ( top_cell,cets );

    for ( uint i = 0 ; i < cet_ct ; i++ )
    {
      if ( dataset->isCellCritical ( cets[i] ) )
      {
//        connectCps(msgraph,start_cellId,cets[i]);
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
            (*dataset->m_cell_own)(cets[i]) = start_cellId;

            // mark here that the parent of next cell is top_cell
            cell_queue.push ( next_cell );
          }
        }
      }
    }
  }
}

GridDataset::GridDataset (const rect_t &r,const rect_t &e) :
    m_rect (r),m_ext_rect (e),m_ptcomp(this),
    m_vert_fns_ref(NULL)
{

  // TODO: assert that the given rect is of even size..
  //       since each vertex is in the even positions
  //

  m_vert_fns_ref = NULL;
  m_cell_flags   = NULL;
  m_cell_pairs   = NULL;
  m_cell_own     = NULL;

}

GridDataset::GridDataset () :
    m_ptcomp(this)
{
  m_vert_fns_ref = NULL;
  m_cell_flags   = NULL;
  m_cell_pairs   = NULL;
  m_cell_own     = NULL;
}

GridDataset::~GridDataset ()
{
  clear();

  clear_fnref();
}

void GridDataset::init()
{
  rect_size_t   s = m_ext_rect.size();

  m_cell_flags = new cellflag_array_t( (boost::extents[1+s[0]][1+s[1]]));
  m_cell_pairs = new cellpair_array_t( (boost::extents[1+s[0]][1+s[1]]));
  m_cell_own   = new cellpair_array_t( (boost::extents[1+s[0]][1+s[1]]));

  for (int y = 0 ; y<=s[1];++y)
    for (int x = 0 ; x<=s[0];++x)
      (*m_cell_flags)[x][y] = CELLFLAG_UNKNOWN;

  rect_point_t bl = m_ext_rect.bottom_left();

  (*m_cell_flags).reindex (bl);
  (*m_cell_pairs).reindex (bl);
  (*m_cell_own).reindex (bl);
}

void  GridDataset::clear()
{
  if(m_cell_flags != NULL)
    delete m_cell_flags;

  if(m_cell_pairs != NULL)
    delete m_cell_pairs;

  if(m_cell_own != NULL)
    delete m_cell_own;

  m_critical_cells.clear();

  m_cell_flags   = NULL;
  m_cell_pairs   = NULL;
  m_cell_own     = NULL;
}

void GridDataset::init_fnref(cell_fn_t * pData)
{
  rect_size_t   s = m_ext_rect.size();

  if(pData != NULL)
    m_vert_fns_ref =
        new varray_ref_t(pData,boost::extents[1+s[0]/2][1+s[1]/2],
                         boost::fortran_storage_order());

  rect_point_t bl = m_ext_rect.bottom_left();

  if(pData != NULL)
    (*m_vert_fns_ref).reindex (bl/2);

}

void GridDataset::clear_fnref()
{
  if(m_vert_fns_ref != NULL)
    delete m_vert_fns_ref;

  m_vert_fns_ref = NULL;
}

GridDataset::cellid_t   GridDataset::getCellPairId (cellid_t c) const
{
  if ((*m_cell_flags) (c) &CELLFLAG_PAIRED == 0)
    throw std::logic_error ("invalid pair requested");

  return (*m_cell_pairs) (c);
}

bool GridDataset::compareCells( cellid_t c1,cellid_t  c2 ) const
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

GridDataset::cell_fn_t GridDataset::get_cell_fn (cellid_t c) const
{

  if(m_vert_fns_ref == NULL)
    return 0.0;

  cell_fn_t  fn = 0.0;

  cellid_t pts[20];

  uint pts_ct = getCellPoints (c,pts);

  for (int j = 0 ; j <pts_ct ;++j)
    fn += (*m_vert_fns_ref) (pts[j]/2);

  fn /= pts_ct;

  return fn;
}

void GridDataset::set_cell_fn (cellid_t c,cell_fn_t f)
{
  if (getCellDim (c) != 0)
    throw std::logic_error ("values only for vertices are specified");

  c[0] /=2;

  c[1] /=2;

  (*m_vert_fns_ref) (c) = f;
}

uint GridDataset::getCellPoints (cellid_t c,cellid_t  *p) const
{
  switch (getCellDim (c))
  {
  case 0:
    p[0] = c;
    return 1;
  case 1:
    {
      cell_coord_t d0 = c[0]&0x01,d1 = c[1]&0x01;
      p[0] = cellid_t (c[0]+d0,c[1]+d1);
      p[1] = cellid_t (c[0]-d0,c[1]-d1);
    }

    return 2;
  case 2:
    p[0] = cellid_t (c[0]+1,c[1]+1);
    p[1] = cellid_t (c[0]+1,c[1]-1);
    p[2] = cellid_t (c[0]-1,c[1]-1);
    p[3] = cellid_t (c[0]-1,c[1]+1);
    return 4;
  default:
    throw std::logic_error ("impossible dim");
    return 0;
  }
}

uint GridDataset::getCellFacets (cellid_t c,cellid_t *f) const
{
  switch (getCellDim (c))
  {
  case 0:
    return 0;
  case 1:
    {
      cell_coord_t d0 = c[0]&0x01,d1 = c[1]&0x01;
      f[0] = cellid_t (c[0]+d0,c[1]+d1);
      f[1] = cellid_t (c[0]-d0,c[1]-d1);
    }

    return 2;
  case 2:
    f[0] = cellid_t (c[0]  ,c[1]+1);
    f[1] = cellid_t (c[0]  ,c[1]-1);
    f[2] = cellid_t (c[0]-1,c[1]);
    f[3] = cellid_t (c[0]+1,c[1]);
    return 4;
  default:
    throw std::logic_error ("impossible dim");
    return 0;
  }
}

uint GridDataset::getCellCofacets (cellid_t c,cellid_t *cf) const
{
  uint cf_ct = 0;

  switch (getCellDim (c))
  {
  case 0:
    cf[0] = cellid_t (c[0]  ,c[1]+1);
    cf[1] = cellid_t (c[0]  ,c[1]-1);
    cf[2] = cellid_t (c[0]-1,c[1]);
    cf[3] = cellid_t (c[0]+1,c[1]);
    cf_ct =  4;
    break;
  case 1:
    {
      cell_coord_t d0 = c[0]&0x01,d1 = c[1]&0x01;
      cf[0] = cellid_t (c[0]+d1,c[1]+d0);
      cf[1] = cellid_t (c[0]-d1,c[1]-d0);
      cf_ct =  2;
    }

    break;
  case 2:
    return 0;
  default:
    throw std::logic_error ("impossible dim");
    return 0;
  }

  // position in cf[] where the next valid cf should be placed
  uint cf_nv_pos = 0;

  for (uint i = 0 ;i < cf_ct;++i)
    if (m_ext_rect.contains (cf[i]))
      cf[cf_nv_pos++] = cf[i];

  return cf_nv_pos;

}

bool GridDataset::isPairOrientationCorrect (cellid_t c, cellid_t p) const
{
  return (getCellDim (c) <getCellDim (p));
}

bool GridDataset::isCellMarked (cellid_t c) const
{
  return ! ((*m_cell_flags) (c) == CELLFLAG_UNKNOWN);
}

bool GridDataset::isCellCritical (cellid_t c) const
{
  return ((*m_cell_flags) (c) & CELLFLAG_CRITCAL);
}

bool GridDataset::isCellPaired (cellid_t c) const
{
  return ((*m_cell_flags) (c) & CELLFLAG_PAIRED);
}

void GridDataset::pairCells (cellid_t c,cellid_t p)
{
  (*m_cell_pairs) (c) = p;
  (*m_cell_pairs) (p) = c;

  (*m_cell_flags) (c) = (*m_cell_flags) (c) |CELLFLAG_PAIRED;
  (*m_cell_flags) (p) = (*m_cell_flags) (p) |CELLFLAG_PAIRED;
}

void GridDataset::markCellCritical (cellid_t c)
{
  (*m_cell_flags) (c) = (*m_cell_flags) (c) |CELLFLAG_CRITCAL;
}

bool GridDataset::isTrueBoundryCell (cellid_t c) const
{
  return (m_ext_rect.isOnBoundry (c));
}

bool GridDataset::isFakeBoundryCell (cellid_t c) const
{
  return (m_rect.isOnBoundry (c) && (!m_ext_rect.isOnBoundry (c)));
}

bool GridDataset::isCellExterior (cellid_t c) const
{
  return (!m_rect.contains (c) && m_ext_rect.contains (c));
}

std::string GridDataset::getCellFunctionDescription (cellid_t c) const
{
  std::stringstream ss;

  ( (std::ostream &) ss) <<c;

  return ss.str();

}

std::string GridDataset::getCellDescription (cellid_t c) const
{

  std::stringstream ss;

  ( (std::ostream &) ss) <<c;

  return ss.str();

}

void GridDataset::work()
{
  assignGradients();

  collateCriticalPoints();

  assignCellOwnerExtrema();
}

void  GridDataset::assignGradients()
{

  // determine all the pairings of all cells in m_rect
  for (cell_coord_t y = m_rect.bottom(); y <= m_rect.top();y += 1)
    for (cell_coord_t x = m_rect.left(); x <= m_rect.right();x += 1)
    {
    cellid_t c (x,y),p;

    if (isCellMarked (c))
      continue;

    if (lowestPairableCoFacet (this,c,p))
      pairCells (c,p);
  }

  for (cell_coord_t y = m_rect.bottom(); y <= m_rect.top();y += 1)
    for (cell_coord_t x = m_rect.left(); x <= m_rect.right();x += 1)
    {
    cellid_t c (x,y);

    if (!isCellMarked (c)) markCellCritical (c);
  }

  // mark artificial boundry as critical

  for (cell_coord_t x = m_rect.left(); x <= m_rect.right();x += 1)
  {
    cellid_t bcs[] = {cellid_t (x,m_rect.bottom()),cellid_t (x,m_rect.top()) };

    for (uint i = 0 ; i <sizeof (bcs) /sizeof (cellid_t);++i)
    {
      cellid_t &c = bcs[i];

      if (isCellCritical (c)) continue;

      cellid_t cf[20];

      u_int cf_ct =  getCellCofacets (c,cf);

      for (u_int j = 0 ; j <cf_ct;++j)
      {
        if (isCellExterior (cf[j]))
        {
          markCellCritical (c);
          markCellCritical (getCellPairId (c));
          break;
        }
      }
    }
  }

  for (cell_coord_t y = m_rect.bottom() +1; y < m_rect.top();y += 1)
  {
    cellid_t bcs[] = {cellid_t (m_rect.left(),y),cellid_t (m_rect.right(),y) };

    for (uint i = 0 ; i <sizeof (bcs) /sizeof (cellid_t);++i)
    {
      cellid_t &c = bcs[i];

      if (isCellCritical (c)) continue;

      cellid_t cf[20];

      u_int cf_ct =  getCellCofacets (c,cf);

      for (u_int j = 0 ; j <cf_ct;++j)
      {
        if (isCellExterior (cf[j]))
        {
          markCellCritical (c);
          markCellCritical (getCellPairId (c));
          break;
        }
      }
    }
  }
}

void  GridDataset::collateCriticalPoints()
{
  for (cell_coord_t y = m_ext_rect.bottom(); y <= m_ext_rect.top();y += 1)
    for (cell_coord_t x = m_ext_rect.left(); x <= m_ext_rect.right();x += 1)
    {
    cellid_t c (x,y);

    if (isCellCritical (c))
      m_critical_cells.push_back(c);
  }
}


void  GridDataset::assignCellOwnerExtrema()
{
  for (cellid_list_t::iterator it = m_critical_cells.begin() ;
  it != m_critical_cells.end();++it)
  {

    (*m_cell_own)(*it) = *it;

    switch (getCellDim (*it))
    {
    case 0:
      track_gradient_tree_bfs(this,*it,GRADIENT_DIR_UPWARD);
      break;
    case 2:
      track_gradient_tree_bfs(this,*it,GRADIENT_DIR_DOWNWARD);
      break;
    default:
      break;
    }
  }

}

void  GridDataset::writeout_connectivity(mscomplex_t *msgraph)
{

  addCriticalPointsToMSComplex
      (msgraph,m_critical_cells.begin(),m_critical_cells.end());

  msgraph->m_cp_fns.resize(m_critical_cells.size());


  for (cellid_list_t::iterator it = m_critical_cells.begin() ;
  it != m_critical_cells.end();++it)
  {
    cellid_t c = *it;

    uint cp_idx = msgraph->m_id_cp_map[c];

    msgraph->m_cp_fns[cp_idx] = get_cell_fn(c);

    if(getCellDim(c) == 1)
    {
      cellid_t f[4],cf[4];

      uint f_ct = getCellFacets(c,f);
      uint cf_ct = getCellCofacets(c,cf);

      for(uint i = 0 ; i < f_ct;++i)
      {
        cellid_t f_own_cp = (*m_cell_own)(f[i]);

        if(f_own_cp != cellid_t(-1,-1))
          connectCps(msgraph,c,f_own_cp);
      }

      for(uint i = 0 ; i < cf_ct;++i)
      {
        cellid_t cf_own_cp = (*m_cell_own)(cf[i]);

        if(cf_own_cp != cellid_t(-1,-1))
          connectCps(msgraph,c,cf_own_cp);
      }
    }

    if(!isCellPaired(c))  continue;

    msgraph->m_cps[cp_idx]->isBoundryCancelable = true;

    msgraph->m_cps[cp_idx]->pair_idx =
        msgraph->m_id_cp_map[getCellPairId(c)];
  }
}

void GridDataset::getCellCoord (cellid_t c,double &x,double &y,double &z)
{
  x = c[0];
  y = 0;
  z = c[1];

  cellid_t pts[20];

  if(m_ext_rect.contains(c))
  {
    y= get_cell_fn(c);

  }
}

void GridDataset::log_flags()
{
  for (cell_coord_t y = m_ext_rect.bottom(); y <= m_ext_rect.top();y += 1)
  {
    for (cell_coord_t x = m_ext_rect.left(); x <= m_ext_rect.right();x += 1)
    {
      cellid_t c(x,y);

      int val = (*m_cell_flags)(c);

      std::cout<<val<<" ";
    }
    std::cout<<std::endl;
  }
}

void GridDataset::log_pairs()
{
  for (cell_coord_t y = m_ext_rect.bottom(); y <= m_ext_rect.top();y += 1)
  {
    for (cell_coord_t x = m_ext_rect.left(); x <= m_ext_rect.right();x += 1)
    {
      cellid_t c(x,y);
      std::cout<<(*m_cell_pairs)(c)<<" ";
    }
    std::cout<<std::endl;
  }
}
