#include <cmath>
#include <queue>
#include <iostream>

#include <grid_mscomplex.h>
#include <grid_dataset.h>
#include <fstream>
#include <limits>

using namespace std;

namespace grid
{
  inline std::string idx_to_string(mscomplex_t *msc,uint i)
  {
    return to_string(msc->m_cps[i]->cellid);
  }

  inline std::string edge_to_string(mscomplex_t *msc,uint_pair_t e)
  {
    std::stringstream ss;

    ss<<idx_to_string(msc,e[0])<<"----"<<idx_to_string(msc,e[1]);

    return ss.str();
  }

  critpt_t::critpt_t(cellid_t c,uchar i,cell_fn_t f, cellid_t v):
      cellid(c),index(i),fn(f),vert_cell(v)
  {
    is_cancelled          = false;
    is_paired             = false;
    pair_idx              = -1;
  }

  critpt_t::critpt_t(const critpt_t &c)
  {
    cellid       = c.cellid;
    index        = c.index;
    vert_cell    = c.vert_cell;
    fn           = c.fn;
    vert_cell    = c.vert_cell;

    is_cancelled = false;
    is_paired    = false;
    pair_idx     = -1;
  }

  mscomplex_t::mscomplex_t(rect_t r,rect_t e)
    :m_rect(r),m_ext_rect(e)
  {
  }

  mscomplex_t::~mscomplex_t()
  {
    clear();
  }

  int mscomplex_t::add_critpt(cellid_t c,uchar i,cell_fn_t f,cellid_t v)
  {
    ASSERT(m_id_cp_map.count(c) == 0);

    critpt_t * cp  = new critpt_t(c,i,f,v);
    m_id_cp_map.insert(std::make_pair(c,m_cps.size()));
    m_cps.push_back(cp);

    return (m_cps.size()-1);
  }

  int mscomplex_t::add_critpt(const critpt_t &c)
  {
    ASSERT(m_id_cp_map.count(c.cellid) == 0);

    critpt_t * cp  = new critpt_t(c);
    m_id_cp_map.insert(std::make_pair(cp->cellid,m_cps.size()));
    m_cps.push_back(cp);

    return (m_cps.size()-1);
  }

  void mscomplex_t::connect_cps(cellid_t c0,cellid_t c1)
  {
    try
    {
      ASSERT(m_id_cp_map.count(c0) == 1);
      ASSERT(m_id_cp_map.count(c1) == 1);

      connect_cps(uint_pair_t(m_id_cp_map[c0],m_id_cp_map[c1]));
    }
    catch(assertion_error e)
    {
      e.push(_FFL);
      e.push(SVAR2(c0,m_id_cp_map.count(c0)));
      e.push(SVAR2(c1,m_id_cp_map.count(c1)));
      throw;
    }
  }

  void mscomplex_t::connect_cps(uint_pair_t e)
  {
    order_pr_by_cp_index(this,e);

    ASSERT(m_cps[e[1]]->index+1 == m_cps[e[0]]->index);

    for(uint dir = 0 ; dir < GRADDIR_COUNT;++dir)
      if(m_cps[e[dir]]->conn[dir].count(e[dir^1]) == 2)
        return;

    for(uint dir = 0 ; dir < GRADDIR_COUNT;++dir)
      m_cps[e[dir]]->conn[dir].insert(e[dir^1]);

  }

  void mscomplex_t::dir_connect_cps(cellid_t c1,cellid_t c2)
  {
    ASSERT(m_id_cp_map.count(c1) == 1);
    ASSERT(m_id_cp_map.count(c2) == 1);

    dir_connect_cps(uint_pair_t(m_id_cp_map[c1],m_id_cp_map[c2]));
  }

  void mscomplex_t::dir_connect_cps(uint_pair_t pr)
  {
    critpt_t * c = m_cps[pr[0]];
    critpt_t * p = m_cps[pr[1]];

    try
    {
      ASSERT(c->is_paired != p->is_paired);
      ASSERT(abs(c->index-p->index) == 1);

      if(p->is_paired)
      {
        std::swap(c,p);
        std::swap(pr[0],pr[1]);
      }

      if(c->index > p->index)
        c->conn[0].insert(pr[1]);
      else
        c->conn[1].insert(pr[1]);
    }
    catch(assertion_error e)
    {
      e.push(_FFL);
      e.push(SVAR2(c->cellid,p->cellid));
      e.push(SVAR2(c->is_paired,p->is_paired));
      e.push(SVAR2((int)c->index,(int)p->index));
      throw;
    }
  }

  void mscomplex_t::pair_cps(cellid_t c1,cellid_t c2)
  {
    ASSERT(m_id_cp_map.count(c1) == 1);
    ASSERT(m_id_cp_map.count(c2) == 1);

    pair_cps(uint_pair_t(m_id_cp_map[c1],m_id_cp_map[c2]));
  }

  void mscomplex_t::pair_cps(uint_pair_t p)
  {
    m_cps[p[0]]->pair_idx = p[1];
    m_cps[p[1]]->pair_idx = p[0];

    m_cps[p[0]]->is_paired = true;
    m_cps[p[1]]->is_paired = true;
  }

  void cancelPairs ( mscomplex_t *msc,uint_pair_t e,
                     uint_pair_list_t * new_edges = NULL)
  {

    order_pr_by_cp_index(msc,e);

    critpt_t * cp[] = {msc->m_cps[e[0]],msc->m_cps[e[1]]};

    ASSERT(cp[0]->index == cp[1]->index+1);
    ASSERT(cp[0]->conn[0].count(e[1]) == 1);
    ASSERT(cp[0]->conn[0].count(e[1]) == cp[1]->conn[1].count(e[0]));

    conn_iter_t it[GRADDIR_COUNT];

    for(uint dir = 0 ; dir < 2;++dir)
      cp[dir]->conn[dir].erase(e[dir^1]);

    // cps in lower of u except l
    for(it[0] = cp[0]->conn[0].begin();it[0] != cp[0]->conn[0].end();++it[0])
      for(it[1] = cp[1]->conn[1].begin();it[1] != cp[1]->conn[1].end();++it[1])
      {
        ASSERT(msc->m_cps[*it[0]]->is_cancelled == false);
        ASSERT(msc->m_cps[*it[1]]->is_cancelled == false);

        msc->connect_cps(uint_pair_t(*it[0],*it[1]));
      }

    if(new_edges)
      for(it[0] = cp[0]->conn[0].begin();it[0] != cp[0]->conn[0].end();++it[0])
        for(it[1] = cp[1]->conn[1].begin();it[1] != cp[1]->conn[1].end();++it[1])
            new_edges->push_back(uint_pair_t(*it[1],*it[0]));

    for(uint dir = 0 ; dir<2;++dir)
     for(conn_iter_t it = cp[dir]->conn[dir].begin();it != cp[dir]->conn[dir].end();++it)
      {
        msc->m_cps[*it]->conn[dir^1].erase(e[dir]);
      }

    for(uint dir = 0 ; dir<2;++dir)
      for(conn_iter_t it = cp[dir]->conn[dir^1].begin();it != cp[dir]->conn[dir^1].end();++it)
      {
        msc->m_cps[*it]->conn[dir].erase(e[dir]);
      }

    for(uint dir = 0 ; dir < 2;++dir)
      cp[dir]->conn[dir].insert(e[dir^1]);

    for(uint dir = 0 ; dir < 2;++dir)
      cp[dir]->is_cancelled = true;

    for(uint dir = 0 ; dir < 2;++dir)
      cp[dir]->conn[dir^1].clear();
  }

  void uncancel_pair( mscomplex_t  *msc,uint_pair_t e)
  {
    order_pr_by_cp_index(msc,e);

    critpt_t * cp[] = {msc->m_cps[e[0]],msc->m_cps[e[1]]};

    ASSERT(cp[0]->is_cancelled == true && cp[1]->is_cancelled ==true);
    ASSERT(cp[0]->index == cp[1]->index + 1);
    ASSERT(cp[0]->conn[0].count(e[1]) != 0);
    ASSERT(cp[0]->conn[0].count(e[1]) == cp[1]->conn[1].count(e[0]));

    conn_iter_t i_it,j_it;

    for(uint dir = 0 ; dir <2 ; ++dir)
    {
      conn_set_t new_conn;

      for(i_it = cp[dir]->conn[dir].begin();i_it != cp[dir]->conn[dir].end() ; ++i_it )
      {
        if(*i_it == e[dir^1])
          continue;

        critpt_t *ccp = msc->m_cps[*i_it];

        if(ccp->is_paired == false)
        {
          new_conn.insert(*i_it);
          continue;
        }

        ASSERT(ccp->is_cancelled == false);

        critpt_t *ccp_pr = msc->m_cps[ccp->pair_idx];

        try
        {
          for(j_it = ccp_pr->conn[dir].begin(); j_it != ccp_pr->conn[dir].end() ; ++j_it )
          {
            ASSERT(msc->m_cps[*j_it]->is_paired == false);
            ASSERT(abs<int>(msc->m_cps[*j_it]->index-cp[dir]->index) == 1);

            if(new_conn.count(*j_it) == 0)
              new_conn.insert(*j_it);
          }
        }
        catch(assertion_error ex)
        {
          ex.push(_FFL);

          ex.push(SVAR(msc->m_cps[*j_it]->cellid));
          ex.push(SVAR(msc->m_cps[*j_it]->is_paired));

          if(msc->m_cps[*j_it]->is_paired)
            ex.push(SVAR(msc->m_cps[msc->m_cps[*j_it]->pair_idx]->cellid));

          ex.push("");
          ex.push(SVAR(ccp->cellid));
          ex.push(SVAR(ccp_pr->cellid));

          ex.push("");
          ex.push(SVAR(dir));
          ex.push(SVAR(cp[dir]->cellid));
          ex.push(SVAR(cp[dir^1]->cellid));

          throw;
        }
      }
      msc->m_cps[e[dir]]->is_cancelled = false;
      cp[dir]->conn[dir].clear();
      cp[dir]->conn[dir].insert(new_conn.begin(),new_conn.end());
    }
  }

  void mscomplex_t::merge_up
      (const mscomplex_t& msc1,
       const mscomplex_t& msc2,
       const rect_t& bnd)
  {
    typedef conn_set_t::const_iterator cconn_it_t;

    try
    {
      ASSERT(bnd.eff_dim() == gc_grid_dim-1);
      ASSERT(m_ext_rect.intersection(bnd) == bnd);

      const mscomplex_t *msc[] = {&msc1,&msc2};

      for (int i = 0 ; i < 2 ; ++i)
      {
        for(int j = 0; j < msc[i]->m_cps.size();++j)
        {
          const critpt_t *cp = msc[i]->m_cps[j];

          try
          {
            if(cp->is_cancelled)
              continue;

            ASSERT((bnd.contains(cp->cellid) == false) || (msc[i^1]->m_id_cp_map.count(cp->cellid) == 1));
            ASSERT((bnd.contains(cp->cellid) == true)  || (msc[i^1]->m_id_cp_map.count(cp->cellid) == 0));

            if((i == 1) && (bnd.contains(cp->cellid)))
              continue;

            add_critpt(*cp);
          }
          catch(assertion_error e)
          {
            e.push(_FFL);
            e.push(SVAR(cp->cellid));
            e.push(SVAR(j));
            e.push(SVAR(i));
            throw;
          }
        }
      }

      for (int i = 0 ; i < 2 ; ++i)
      {
        for(int j = 0; j < msc[i]->m_cps.size();++j)
        {
          const critpt_t *cp = msc[i]->m_cps[j];

          if(cp->is_cancelled)
            continue;

          bool is_cp_in_bnd = bnd.contains(cp->cellid);

          for(cconn_it_t it = cp->conn[0].begin();it != cp->conn[0].end();++it)
          {
            const critpt_t *ccp = msc[i]->m_cps[*it];

            if(ccp->is_cancelled)
              continue;

            if((i == 1) && is_cp_in_bnd && bnd.contains(ccp->cellid))
              continue;

            connect_cps(cp->cellid,ccp->cellid);
          }

          if(!cp->is_paired)
            continue;

          pair_cps(cp->cellid,msc[i]->m_cps[cp->pair_idx]->cellid);
        }
      }

      rect_t ixn = m_rect.intersection(bnd);

      for(cellid_t c = ixn.lower_corner() ; c[2] <= ixn[2][1]; ++c[2])
      {
        for(c[1] = ixn[1][0] ; c[1] <= ixn[1][1]; ++c[1])
        {
          for(c[0] = ixn[0][0] ; c[0] <= ixn[0][1]; ++c[0])
          {
            if (m_id_cp_map.count(c) == 0)
              continue;

            int cp_idx = m_id_cp_map[c];

            critpt_t *cp = m_cps[cp_idx];

            if(!cp->is_paired)
              continue;

            if(bnd.contains(m_cps[cp->pair_idx]->cellid))
              continue;

            cancelPairs(this,uint_pair_t(cp_idx,cp->pair_idx));
          }
        }
      }
    }
    catch(assertion_error e)
    {
      e.push(_FFL);
      e.push(SVAR(bnd));
      e.push(SVAR2(m_rect,m_ext_rect));
      e.push(SVAR2(msc1.m_rect,msc1.m_ext_rect));
      e.push(SVAR2(msc2.m_rect,msc2.m_ext_rect));
      throw;
    }
  }

  void mscomplex_t::merge_down
      (mscomplex_t& msc1,
       mscomplex_t& msc2,
       const rect_t& bnd)
  {
    mscomplex_t *msc_arr[] = {&msc1,&msc2};

    try
    {
      ASSERT(bnd.eff_dim() == gc_grid_dim-1);
      ASSERT(m_ext_rect.intersection(bnd) == bnd);

      rect_t ixn = m_rect.intersection(bnd);

      for(cellid_t c = ixn.upper_corner() ; c[2] >= ixn[2][0]; --c[2])
      {
        for(c[1] = ixn[1][1] ; c[1] >= ixn[1][0]; --c[1])
        {
          for(c[0] = ixn[0][1] ; c[0] >= ixn[0][0]; --c[0])
          {
            if (m_id_cp_map.count(c) == 0)
              continue;

            int cp_idx = m_id_cp_map[c];

            critpt_t * cp= m_cps[cp_idx];

            if(!cp->is_paired)
              continue;

            critpt_t * cp_pr =  m_cps[cp->pair_idx];

            if(bnd.contains(cp_pr->cellid))
              continue;

            uncancel_pair(this,uint_pair_t(cp_idx,cp->pair_idx));
          }
        }
      }
    }
    catch(assertion_error e)
    {
      e.push(_FFL);
      e.push(SVAR(bnd));
      e.push(SVAR2(m_rect,m_ext_rect));
      e.push(SVAR2(msc1.m_rect,msc1.m_ext_rect));
      e.push(SVAR2(msc2.m_rect,msc2.m_ext_rect));
      throw;
    }

    vector<int> cpi_to_mcpi(m_cps.size());

    for(int m = 0 ; m < 2; ++m)
    {
      mscomplex_t &msc = *msc_arr[m];

      for(int cpi = 0 ; cpi < m_cps.size(); ++cpi)
      {
        critpt_t * cp = m_cps[cpi];

        if(msc.m_id_cp_map.count(cp->cellid) == 0 )
          cpi_to_mcpi[cpi] = -1;
        else
          cpi_to_mcpi[cpi] = msc.m_id_cp_map[cp->cellid];
      }

      for(int cpi = 0 ; cpi < m_cps.size(); ++cpi)
      {
        critpt_t * cp = m_cps[cpi];

        if(cpi_to_mcpi[cpi] == -1 )
          continue;

        int mcpi       = cpi_to_mcpi[cpi];
        critpt_t * mcp = msc.m_cps[mcpi];

        mcp->conn[0].clear();
        mcp->conn[1].clear();

        if(cp->is_paired == false)
          continue;

        critpt_t * cp_pr =  m_cps[cp->pair_idx];

        if(cpi_to_mcpi[cp->pair_idx] != -1)
          continue;

        cpi_to_mcpi[cp->pair_idx] = msc.add_critpt(*cp_pr);
        msc.pair_cps(cp->cellid,cp_pr->cellid);
      }

      for(int cpi = 0 ; cpi < m_cps.size(); ++cpi)
      {
        critpt_t * cp = m_cps[cpi];

        if(cpi_to_mcpi[cpi] == -1 )
          continue;

        if(cp->is_paired == false)
          continue;

        critpt_t * cp_pr =  m_cps[cp->pair_idx];

        int d = (cp_pr->index - cp->index + 1)/2;

        if(d == 0 && (msc.m_rect.contains(cp->cellid)))
          continue;

        if(d == 1 && (msc.m_rect.contains(cp_pr->cellid) == false))
          continue;

        for(conn_iter_t it = cp->conn[d].begin();it != cp->conn[d].end();++it)
        {
          int ccpi = *it;

          critpt_t * ccp  = m_cps[ccpi];

          try
          {
            ASSERT(ccp->is_paired == false);

            if(cpi_to_mcpi[ccpi] == -1)
              cpi_to_mcpi[ccpi] = msc.add_critpt(*ccp);

            uint_pair_t mpr(cpi_to_mcpi[cpi],cpi_to_mcpi[ccpi]);

            msc.dir_connect_cps(mpr);
          }
          catch(assertion_error e)
          {
            e.push(_FFL);
            e.push(SVAR2(cp->cellid,ccp->cellid));
            e.push(SVAR2(cp->cellid,cp_pr->cellid));
            e.push(SVAR2(cp->is_paired,ccp->is_paired));
            e.push(SVAR2(msc.m_rect,msc.m_ext_rect));
            e.push(SVAR2(m_rect,m_ext_rect));
            throw;
          }
        }
      }
    }
  }

  void mscomplex_t::clear()
  {
    std::for_each(m_cps.begin(),m_cps.end(),&delete_ftor<critpt_t>);
    m_cps.clear();
    m_id_cp_map.clear();
  }

  struct persistence_comparator_t
  {
    mscomplex_t *m_msc;

    persistence_comparator_t(mscomplex_t *m):m_msc(m){}

    bool operator()(const uint_pair_t & p0, const uint_pair_t &p1)
    {
      return cmp_lt(p1,p0);
    }

    bool cmp_lt(uint_pair_t p0, uint_pair_t p1)
    {
      order_pr_by_cp_index(m_msc,p0);
      order_pr_by_cp_index(m_msc,p1);

      cellid_t v00 = m_msc->m_cps[p0[0]]->vert_cell;
      cellid_t v01 = m_msc->m_cps[p0[1]]->vert_cell;
      cellid_t v10 = m_msc->m_cps[p1[0]]->vert_cell;
      cellid_t v11 = m_msc->m_cps[p1[1]]->vert_cell;

      cellid_t c00 = m_msc->m_cps[p0[0]]->cellid;
      cellid_t c01 = m_msc->m_cps[p0[1]]->cellid;
      cellid_t c10 = m_msc->m_cps[p1[0]]->cellid;
      cellid_t c11 = m_msc->m_cps[p1[1]]->cellid;

      if( (v00 == v01 ) != (v10 == v11))
        return (v00 == v01 );

      if( (v00 == v01 ) &&(v10 == v11))
      {
        if(v00 == v10)
        {
          if(c00 != c10)
            return c00 < c10;
          else
            return c01 < c11;
        }
        else
        {
          return (v00 < v10);
        }
      }

      cell_fn_t f00 = m_msc->m_cps[p0[0]]->fn;
      cell_fn_t f01 = m_msc->m_cps[p0[1]]->fn;
      cell_fn_t f10 = m_msc->m_cps[p1[0]]->fn;
      cell_fn_t f11 = m_msc->m_cps[p1[1]]->fn;

      cell_fn_t d1 = std::abs(f01-f00);
      cell_fn_t d2 = std::abs(f11-f10);

      if(d1 != d2)
        return d1 < d2;

      if(c00 != c10)
        return c00 < c10;

      return c01 < c11;
    }
  };

  bool is_valid_canc_edge(mscomplex_t *msc,uint_pair_t e )
  {
    order_pr_by_cp_index(msc,e);

    critpt_t * cp[] = {msc->m_cps[e[0]],msc->m_cps[e[1]]};

    for(uint dir = 0 ; dir < 2; ++dir)
      if(cp[dir]->is_cancelled)
        return false;

    if(msc->m_rect.isOnBoundry(cp[0]->cellid) !=
       msc->m_rect.isOnBoundry(cp[1]->cellid))
      return false;

    for(uint dir = 0 ; dir < 2; ++dir)
      if(cp[dir]->conn[dir].count(e[dir^1]) != 1)
        return false;

    return true;
  }

  bool is_epsilon_persistent(mscomplex_t *msc,uint_pair_t e )
  {
    return (msc->m_cps[e[0]]->vert_cell == msc->m_cps[e[1]]->vert_cell);
  }

  void mscomplex_t::simplify(uint_pair_list_t & canc_pairs_list,
                               double simplification_treshold)
  {
    typedef std::priority_queue
        <uint_pair_t,uint_pair_list_t,persistence_comparator_t>
        canc_pair_priq_t;

    persistence_comparator_t comp(this);

    canc_pair_priq_t  canc_pair_priq(comp);



    cell_fn_t max_val = std::numeric_limits<cell_fn_t>::min();
    cell_fn_t min_val = std::numeric_limits<cell_fn_t>::max();

    for(uint i = 0 ;i < m_cps.size();++i)
    {
      critpt_t *cp = m_cps[i];

      max_val = std::max(max_val,m_cps[i]->fn);

      min_val = std::min(min_val,m_cps[i]->fn);

      for(const_conn_iter_t it = cp->conn[0].begin();it != cp->conn[0].end() ;++it)
      {
        if(is_valid_canc_edge(this,uint_pair_t(i,*it)))
          canc_pair_priq.push(uint_pair_t(i,*it));
      }
    }

    double max_persistence = max_val - min_val;

    uint num_cancellations = 0;

    uint num_cancellations_eps = 0;

    while (canc_pair_priq.size() !=0)
    {
      uint_pair_t pr = canc_pair_priq.top();

      canc_pair_priq.pop();

      double persistence = std::abs(m_cps[pr[0]]->fn-m_cps[pr[1]]->fn)/max_persistence;

      if(is_valid_canc_edge(this,pr) == false)
        continue;

      if(is_epsilon_persistent(this,pr) == false)
      {
        if(persistence >= simplification_treshold)
          break;
      }
      else
      {
        num_cancellations_eps++;
      }

//      std::cout
//             <<   "no = "<<num_cancellations<<" "
//             << "pers = "<<persistence<<" "
//             <<"index = ("<<(int)m_cps[pr[0]]->index<<","<<(int)m_cps[pr[1]]->index<<") "
//             << "edge = "<<edge_to_string(this,pr)<<" "
//             <<std::endl;

      uint_pair_list_t new_edges;

      cancelPairs ( this,pr,&new_edges);
      num_cancellations++;

      pair_cps(pr);

      canc_pairs_list.push_back(pr);

      for(uint i = 0 ; i < new_edges.size(); i++)
      {
        canc_pair_priq.push(new_edges[i]);
      }
    }
    std::cout<<"num_cancellations::"<<(num_cancellations)<<std::endl;
    std::cout<<"num_cancellations_eps::"<<(num_cancellations_eps)<<std::endl;
  }

  void mscomplex_t::un_simplify(const uint_pair_list_t &canc_pairs_list)
  {
    for(uint_pair_list_t::const_reverse_iterator it = canc_pairs_list.rbegin();
    it != canc_pairs_list.rend() ; ++it)
      uncancel_pair(this,*it);
  }

  void mscomplex_t::simplify_un_simplify(double simplification_treshold)
  {
    uint_pair_list_t canc_pairs_list;

    simplify(canc_pairs_list,simplification_treshold);

    un_simplify(canc_pairs_list);
  }

  void mscomplex_t::add_disc_tracking_seed_cps()
  {
    for(uint i = 0 ; i < m_cps.size(); ++i)
    {
      if(!m_cps[i]->is_paired)
      {
        for(uint dir = 0 ; dir < 2 ;++dir)
        {
          m_cps[i]->disc[dir].push_back(m_cps[i]->cellid);
          m_cps[i]->contrib[dir].push_back(i);
        }
        continue;
      }

      uint_pair_t e(i,m_cps[i]->pair_idx);

      order_pr_by_cp_index(this,e);

      if(e[0] != i) continue;

      critpt_t * cp[] = {m_cps[e[0]],m_cps[e[1]]};

      ASSERT(cp[1]->index+1 == cp[0]->index);
      ASSERT(cp[0]->is_paired == true && cp[1]->is_paired == true);
      ASSERT(cp[0]->pair_idx  == e[1] && cp[1]->pair_idx  == e[0]);

      for(uint dir = 0 ; dir < 2 ;++dir)
      {
        bool need_disc = false;

        for(conn_iter_t it  = cp[dir]->conn[dir].begin();
                        it != cp[dir]->conn[dir].end(); ++it)
        {
          ASSERT(m_cps[*it]->is_paired == false);

          m_cps[*it]->contrib[dir^1].push_back(e[dir^1]);

          need_disc = true;
        }

        if(need_disc)
          cp[dir^1]->disc[dir^1].push_back(cp[dir^1]->cellid);

      }
    }
  }

  void mscomplex_t::write_manifolds(std::ostream &os)
  {
    const char *mfold_txt[] = {"des","asc"};

    for(uint i = 0 ; i < m_cps.size();++i)
    {
      critpt_t * cp = m_cps[i];

      if(cp->is_paired == true)
        continue;

      for(uint dir = 0 ; dir <2 ;++dir)
      {
        std::set<cellid_t> cset;

        for(uint j = 0 ; j < cp->contrib[dir].size();++j)
        {
          critpt_t *cp_contrib = m_cps[cp->contrib[dir][j]];

          for(uint k = 0; k < cp_contrib->disc[dir].size(); ++k)
          {
            cellid_t c = cp_contrib->disc[dir][k];

            if(cset.count(c) == 0)
              cset.insert(c);
          }
        }

        os<<mfold_txt[dir]<<" "<<cp->cellid<<" "<<cset.size()<<"\n";

        std::copy(cset.begin(),cset.end(),
                  std::ostream_iterator<cellid_t>(os,"\n"));

        os<<std::endl;
      }
    }
  }

  void mscomplex_t::write_manifolds(const std::string &fn)
  {
    std::fstream os(fn.c_str(),std::ios::out);

    ensure(os.is_open(),"failed to open file");

    write_manifolds((std::ostream&)os);

    os.close();
  }

  cellid_t get_cp_cellid(mscomplex_t *msc,uint idx)
  {
    return msc->m_cps[idx]->cellid;
  }

  void mscomplex_t::write_graph(std::ostream & os) const
  {
    using namespace std;

    os<<"# Num Cps"<<endl;
    os<<m_cps.size()<<endl;

    os<<"# grid "<<endl;
    os<<m_rect<<endl;

    os<<"# SL.No idx isPaired pairIdx cpCell vertCell fn"<<endl;

    for(int i = 0 ; i < m_cps.size();++i)
    {
      critpt_t * cp = m_cps[i];

      os<<i<<" ";
      os<<(int)cp->index<<" ";
      os<<(bool)cp->is_paired<<" ";
      os<<(int)cp->pair_idx<<" ";
      os<<cp->cellid<<" ";
      os<<cp->vert_cell<<" ";
      os<<cp->fn<<" ";
      os<<endl;
    }

    os<<"# SL.No  numDes numAsc connList"<<std::endl;

    for(uint i = 0 ; i < m_cps.size();++i)
    {
      critpt_t * cp = m_cps[i];

      os<<(int)i<<" ";
      os<<(int)cp->conn[0].size()<<" ";
      os<<(int)cp->conn[1].size()<<" ";

      for(uint dir = 0 ; dir <2 ;++dir)
      {
        conn_set_t &conn = cp->conn[dir];

        for(conn_iter_t it = conn.begin(); it != conn.end(); ++it)
        {
          os<<*it<<" ";
        }
      }
      os<<endl;
    }

  }

  void mscomplex_t::write_graph(const std::string &fn) const
  {
    std::fstream os(fn.c_str(),std::ios::out);

    ensure(os.is_open(),"failed to open file");

    write_graph(os);

    os.close();
  }

}


//#include <boost/serialization/vector.hpp>
//#include <boost/serialization/array.hpp>
//#include <boost/serialization/map.hpp>
//#include <boost/serialization/set.hpp>
//#include <boost/serialization/base_object.hpp>
//#include <boost/serialization/binary_object.hpp>
//
//#include <boost/archive/binary_iarchive.hpp>
//#include <boost/archive/binary_oarchive.hpp>
//
//using namespace grid;
//
//namespace boost
//{
//  namespace serialization
//  {
//    template<class Archive>
//    void serialize(Archive & ar, rect_range_t & r, const unsigned int )
//    {
//      typedef boost::array<rect_range_t::value_type,rect_range_t::static_size>
//          rect_range_base_t;
//
//      ar & boost::serialization::base_object<rect_range_base_t>(r);
//    }
//
//    template<class Archive>
//    void serialize(Archive & ar, rect_point_t & p, const unsigned int )
//    {
//      typedef boost::array<rect_point_t::value_type,rect_point_t::static_size>
//          rect_point_base_t;
//
//      ar & boost::serialization::base_object<rect_point_base_t>(p);
//    }
//
//    template<class Archive>
//    void serialize(Archive & ar, rect_t & r, const unsigned int )
//    {
//      typedef boost::array<rect_t::value_type,rect_t::static_size>
//          rect_base_t;
//
//      ar & boost::serialization::base_object<rect_base_t>(r);
//    }
//
//    template<class Archive>
//    void serialize(Archive & ar, critpt_t & c, const unsigned int )
//    {
//      ar & c.cellid;
//      for(uint dir = 0 ;dir <0 ;++dir)
//        ar& c.conn[dir];
//      ar & c.is_paired;
//      ar & c.isCancelled;
//      ar & c.pair_idx;
//      ar & c.fn;
//    }
//
//
//    template<class Archive>
//    void serialize(Archive & ar, mscomplex_t & g, const unsigned int )
//    {
//      ar & g.m_rect;
//      ar & g.m_ext_rect;
//      ar & g.m_id_cp_map;
//      ar & g.m_cps;
//    }
//
//    //    template<class Archive>
//    //    void serialize(Archive & ar, GridDataset & ds, const unsigned int )
//    //    {
//    //       ar & ds.m_rect;
//    //       ar & ds.m_ext_rect;
//    //
//    //       GridDataset::rect_size_t ext_sz = ds.m_ext_rect.size();
//    //       uint num_data_items = (ext_sz[0]+1)*(ext_sz[1]+1);
//    //
//    //       if(Archive::is_loading::value)
//    //         ds.init(NULL);
//    //
//    //       ar & make_binary_object(ds.(*m_cell_flags).data(),num_data_items*sizeof(GridDataset::cell_flag_t));
//    //       ar & make_binary_object(ds.m_cell_pairs.data(),num_data_items*sizeof(GridDataset::cellid_t));
//    //    }
//  }
//}
//
////// without the explicit instantiations below, the program will
////// fail to link for lack of instantiantiation of the above function
////// The impls are visible only in this file to save compilation time..
////
////template void boost::serialization::serialize<boost::archive::text_iarchive>(
////    boost::archive::text_iarchive & ar,
////    GridDataset & g,
////    const unsigned int file_version
////);
////
////template void boost::serialization::serialize<boost::archive::text_oarchive>(
////    boost::archive::text_oarchive & ar,
////    GridDataset & g,
////    const unsigned int file_version
////);
//
//
//// without the explicit instantiations below, the program will
//// fail to link for lack of instantiantiation of the above function
//// The impls are visible only in this file to save compilation time..
//
//template void boost::serialization::serialize<boost::archive::binary_iarchive>(
//    boost::archive::binary_iarchive & ar,
//    mscomplex_t & g,
//    const unsigned int file_version
//    );
//template void boost::serialization::serialize<boost::archive::binary_oarchive>(
//    boost::archive::binary_oarchive & ar,
//    mscomplex_t & g,
//    const unsigned int file_version
//    );
//
//
