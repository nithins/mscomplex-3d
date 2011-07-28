#include <cmath>
#include <queue>
#include <iostream>

#include <grid_mscomplex.h>
#include <grid_mscomplex_ensure.h>
#include <grid_dataset.h>
#include <fstream>
#include <limits>

namespace grid
{
  void mark_cancel_pair(mscomplex_t *msc,uint_pair_t e)
  {
    critpt_t * cp1 = msc->m_cps[e[0]];
    critpt_t * cp2 = msc->m_cps[e[1]];

    cp1->is_paired = true;
    cp2->is_paired = true;

    cp1->pair_idx  = e[1];
    cp2->pair_idx  = e[0];

  }

  void mscomplex_t::add_critpt(cellid_t c,uchar i,cell_fn_t f,cellid_t vert_cell)
  {
    critpt_t * cp  = new critpt_t;
    cp->cellid     = c;
    cp->index      = i;
    cp->fn         = f;
    cp->vert_cell  = vert_cell;
    m_id_cp_map.insert(std::make_pair(c,m_cps.size()));
    m_cps.push_back(cp);
  }

  void mscomplex_t::connect_cps(cellid_t c0,cellid_t c1)
  {
    ensure_cellid_critical(this,c0);
    ensure_cellid_critical(this,c1);

    connect_cps(uint_pair_t(m_id_cp_map[c0],m_id_cp_map[c1]));
  }

  bool is_saddle_extremum_pair(mscomplex_t * msc,uint_pair_t e)
  {
    order_pr_by_cp_index(msc,e);

    ensure_index_one_separation(msc,e);

    if(msc->m_cps[e[0]]->index == gc_grid_dim ||
       msc->m_cps[e[1]]->index == 0)
      return true;

    return false;

  }

  void mscomplex_t::connect_cps(uint_pair_t e)
  {
    order_pr_by_cp_index(this,e);

    ensure_ordered_index_one_separation(this,e);

    for(uint dir = 0 ; dir < GRADDIR_COUNT;++dir)
      if(m_cps[e[dir]]->conn[dir].count(e[dir^1]) == 2)
        return;

    for(uint dir = 0 ; dir < GRADDIR_COUNT;++dir)
      m_cps[e[dir]]->conn[dir].insert(e[dir^1]);

  }

  void cancelPairs ( mscomplex_t *msc,uint_pair_t e,
                     uint_pair_list_t * new_edges = NULL)
  {

    order_pr_by_cp_index(msc,e);

    ensure_single_connectivity(msc,e);

    critpt_t * cp[] = {msc->m_cps[e[0]],msc->m_cps[e[1]]};

    conn_iter_t it[GRADDIR_COUNT];

    for(uint dir = 0 ; dir < 2;++dir)
      cp[dir]->conn[dir].erase(e[dir^1]);

    // cps in lower of u except l
    for(it[0] = cp[0]->conn[0].begin();it[0] != cp[0]->conn[0].end();++it[0])
      for(it[1] = cp[1]->conn[1].begin();it[1] != cp[1]->conn[1].end();++it[1])
      {
        ensure_not_cancelled(msc,*it[0]);
        ensure_not_cancelled(msc,*it[1]);

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
      cp[dir]->isCancelled = true;

    for(uint dir = 0 ; dir < 2;++dir)
      cp[dir]->conn[dir^1].clear();
  }

  void uncancel_pair( mscomplex_t  *msc,uint_pair_t e)
  {
    order_pr_by_cp_index(msc,e);

    ensure_ordered_connectivity(msc,e);

    ensure_cp_is_cancelled(msc,e[0]);
    ensure_cp_is_cancelled(msc,e[1]);

    critpt_t * cp[] = {msc->m_cps[e[0]],msc->m_cps[e[1]]};

    conn_iter_t i_it,j_it;

    for(uint dir = 0 ; dir <2 ; ++dir)
    {
      conn_set_t new_conn;

      for(i_it = cp[dir]->conn[dir].begin();i_it != cp[dir]->conn[dir].end() ; ++i_it )
      {
        if(*i_it == e[dir^1]) continue;

        critpt_t *conn_cp = msc->m_cps[*i_it];

        if(!conn_cp->is_paired)
        {
          new_conn.insert(*i_it);
          continue;
        }

        ensure_cp_is_not_cancelled(msc,*i_it);

        ensure_cp_is_paired(msc,*i_it);

        critpt_t *cp_pr = msc->m_cps[conn_cp->pair_idx];

        for(j_it = cp_pr->conn[dir].begin(); j_it != cp_pr->conn[dir].end() ; ++j_it )
        {
          ensure_cp_is_not_paired(msc,*j_it);

          ensure_index_one_separation(msc,uint_pair_t(*j_it,e[dir]));

          if(new_conn.count(*j_it) == 0)
            new_conn.insert(*j_it);
        }
      }

      msc->m_cps[e[dir]]->isCancelled = false;
      cp[dir]->conn[dir].clear();
      cp[dir]->conn[dir].insert(new_conn.begin(),new_conn.end());
    }
  }

  mscomplex_t::~mscomplex_t()
  {
    clear();
  }

//  void shallow_replicate_cp(mscomplex_t &msc,const critpt_t &cp)
//  {
//    if(msc.m_id_cp_map.count(cp.cellid) != 0)
//      throw std::logic_error("this cp is present in msc");
//
//    critpt_t * dest_cp = new critpt_t;
//
//    dest_cp->isCancelled              = cp.isCancelled;
//    dest_cp->is_paired                = cp.is_paired;
//    dest_cp->cellid                   = cp.cellid;
//    dest_cp->fn                       = cp.fn;
//    dest_cp->index                    = cp.index;
//
//    msc.m_id_cp_map[dest_cp->cellid]  = msc.m_cps.size();
//    msc.m_cps.push_back(dest_cp);
//  }
//
//  mscomplex_t * mscomplex_t::merge_up
//      (const mscomplex_t& msc1,
//       const mscomplex_t& msc2)
//  {
//
//    // form the intersection rect
//    rect_t ixn;
//
//    if (!msc2.m_rect.intersection (msc1.m_rect,ixn))
//      throw std::logic_error ("rects should intersect for merge");
//
//    if ( msc2.m_rect.eff_dim() != gc_grid_dim-1)
//      throw std::logic_error ("rects must merge along a d-1 manifold");
//
//    // TODO: ensure that the union of  rects is not including anything extra
//
//    rect_t r = msc1.m_rect.bounding_box(msc1.m_rect);
//
//    rect_t e = msc1.m_ext_rect.bounding_box(msc1.m_ext_rect);
//
//    mscomplex_t * out_msc = new mscomplex_t(r,e);
//
//    const mscomplex_t* msc_arr[] = {&msc1,&msc2};
//
//    // make a union of the critical points in this
//    for (uint i = 0 ; i <2;++i)
//    {
//      const mscomplex_t * msc = msc_arr[i];
//
//      for (uint j = 0 ; j <msc->m_cps.size();++j)
//      {
//        const critpt_t *src_cp = msc->m_cps[j];
//
//        // if it is contained or not
//        if (i == 1 && (out_msc->m_id_cp_map.count(src_cp->cellid) == 1))
//          continue;
//
//        if(src_cp->isCancelled)
//          continue;
//
//        shallow_replicate_cp(*out_msc,*src_cp);
//
//      }
//    }
//
//    for (uint i = 0 ; i <2;++i)
//    {
//      const mscomplex_t * msc = msc_arr[i];
//
//      // copy over connectivity information
//      for (uint j = 0 ; j <msc->m_cps.size();++j)
//      {
//        const critpt_t *src_cp = msc->m_cps[j];
//
//        if(src_cp->isCancelled)
//          continue;
//
//        critpt_t *dest_cp = out_msc->m_cps[out_msc->m_id_cp_map[src_cp->cellid]];
//
//        if(src_cp->is_paired)
//        {
//          critpt_t *src_pair_cp = msc->m_cps[src_cp->pair_idx];
//
//          dest_cp->pair_idx = out_msc->m_id_cp_map[src_pair_cp->cellid];
//        }
//
//        const conn_t *acdc_src[]  = {&src_cp->asc, &src_cp->des};
//
//        conn_t *acdc_dest[] = {&dest_cp->asc,&dest_cp->des};
//
//        bool is_src_cmn_bndry = (ixn.contains(src_cp->cellid) && i == 1);
//
//        for (uint j = 0 ; j < 2; ++j)
//        {
//          for (const_conn_iter_t it = acdc_src[j]->begin();
//          it != acdc_src[j]->end();++it)
//          {
//            const critpt_t *conn_cp = msc->m_cps[*it];
//
//            // common boundry connections would have been found along the boundry
//            if( is_src_cmn_bndry && ixn.contains(conn_cp->cellid))
//              continue;
//
//            if (conn_cp->isCancelled)
//              continue;
//
//            acdc_dest[j]->insert (out_msc->m_id_cp_map[conn_cp->cellid]);
//          }
//        }
//      }
//    }
//
//    static_assert(gc_grid_dim == 3&&"defined for 3-manifolds only");
//
//    cellid_t c;
//
//    for(c[2] = ixn[2][0] ; c[2] <= ixn[2][1]; ++c[2])
//    {
//      for(c[1] = ixn[1][0] ; c[1] <= ixn[1][1]; ++c[1])
//      {
//        for(c[0] = ixn[0][0] ; c[0] <= ixn[0][1]; ++c[0])
//        {
//
//          if(out_msc->m_id_cp_map.count(c) != 1)
//            throw std::logic_error("missing common bndry cp");
//
//          u_int src_idx = out_msc->m_id_cp_map[c];
//
//          critpt_t *src_cp = out_msc->m_cps[src_idx];
//
//          if(src_cp->is_paired || !src_cp->is_paired)
//            continue;
//
//          u_int pair_idx = src_cp->pair_idx;
//
//          cellid_t p = out_msc->m_cps[pair_idx]->cellid;
//
//          if(!out_msc->m_rect.isInInterior(c)&& !out_msc->m_ext_rect.isOnBoundry(c))
//            continue;
//
//          if(!out_msc->m_rect.isInInterior(p)&& !out_msc->m_ext_rect.isOnBoundry(p))
//            continue;
//
//          cancelPairs(out_msc,uint_pair_t(src_idx,pair_idx));
//        }
//      }
//    }
//
//    return out_msc;
//  }
//
//  void mscomplex_t::merge_down(mscomplex_t& msc1,mscomplex_t& msc2)
//  {
//    // form the intersection rect
//    rect_t ixn;
//
//    if (!msc2.m_rect.intersection (msc1.m_rect,ixn))
//      throw std::logic_error ("rects should intersect for merge");
//
//    if ( msc2.m_rect.eff_dim() != gc_grid_dim-1)
//      throw std::logic_error ("rects must merge along a d-1 manifold");
//
//    static_assert(gc_grid_dim == 3&&"defined for 3-manifolds only");
//
//    cellid_t c;
//
//    for(c[2] = ixn[2][1] ; c[2] >= ixn[2][0]; --c[2])
//    {
//      for(c[1] = ixn[1][1] ; c[1] >= ixn[1][0]; --c[1])
//      {
//        for(c[0] = ixn[0][1] ; c[0] >= ixn[0][0]; --c[0])
//        {
//
//          if(this->m_id_cp_map.count(c) != 1)
//            throw std::logic_error("missing common bndry cp");
//
//          u_int src_idx = this->m_id_cp_map[c];
//
//          critpt_t *src_cp = this->m_cps[src_idx];
//
//          if(!src_cp->isCancelled )
//            continue;
//
//          u_int pair_idx = src_cp->pair_idx;
//
//          cellid_t p = this->m_cps[pair_idx]->cellid;
//
//          if(!this->m_rect.isInInterior(c)&& !this->m_ext_rect.isOnBoundry(c))
//            continue;
//
//          if(!this->m_rect.isInInterior(p)&& !this->m_ext_rect.isOnBoundry(p))
//            continue;
//
//          uncancel_pairs(this,uint_pair_t(src_idx,pair_idx));
//        }
//      }
//    }
//
//    // identify and copy the results to msc1 and msc2
//
//    mscomplex_t* msc_arr[] = {&msc1,&msc2};
//
//    for (uint i = 0 ; i <2;++i)
//    {
//      mscomplex_t * msc = msc_arr[i];
//
//      // adjust connections for uncancelled cps in msc
//      for(uint j = 0 ; j < m_cps.size();++j)
//      {
//        critpt_t * src_cp = m_cps[j];
//
//        if(src_cp->isCancelled)
//          throw std::logic_error("all cps ought to be uncancelled by now");
//
//        if(!src_cp->is_paired)
//          continue;
//
//        critpt_t * src_pair_cp = m_cps[src_cp->pair_idx];
//
//        bool src_in_msc      = (msc->m_id_cp_map.count(src_cp->cellid) != 0);
//        bool src_pair_in_msc = (msc->m_id_cp_map.count(src_pair_cp->cellid) != 0);
//
//        if(!src_in_msc && !src_pair_in_msc)
//          continue;
//
//        if(!src_in_msc)
//        {
//          shallow_replicate_cp(*msc,*src_cp);
//        }
//
//        if(!src_pair_in_msc)
//        {
//          shallow_replicate_cp(*msc,*src_pair_cp);
//        }
//
//        uint dest_cp_idx = msc->m_id_cp_map[src_cp->cellid];
//
//        critpt_t *dest_cp = msc->m_cps[dest_cp_idx];
//
//        if(!src_in_msc || !src_pair_in_msc || !dest_cp->is_paired)
//        {
//          uint dest_pair_cp_idx  = msc->m_id_cp_map[src_pair_cp->cellid];
//          critpt_t *dest_pair_cp = msc->m_cps[dest_pair_cp_idx];
//
//          dest_cp->is_paired      = true;
//          dest_pair_cp->is_paired = true;
//
//          dest_cp->pair_idx      = dest_pair_cp_idx;
//          dest_pair_cp->pair_idx = dest_cp_idx;
//        }
//
//        conn_t *src_acdc[] = {&src_cp->asc,&src_cp->des};
//        conn_t *dest_acdc[] = {&dest_cp->asc,&dest_cp->des};
//
//        for(uint k = 0 ; k < 2;++k)
//        {
//          dest_acdc[k]->clear();
//
//          for(conn_iter_t it = src_acdc[k]->begin(); it!=src_acdc[k]->end();++it)
//          {
//            critpt_t *src_conn_cp = m_cps[*it];
//
//            if(src_conn_cp->is_paired == true)
//              throw std::logic_error("only non cancellable cps must be remaining");
//
//            if(msc->m_id_cp_map.count(src_conn_cp->cellid) == 0)
//            {
//              shallow_replicate_cp(*msc,*src_conn_cp);
//            }
//
//            dest_acdc[k]->insert(msc->m_id_cp_map[src_conn_cp->cellid]);
//          }// end it
//        }// end k
//      }// end j
//
//      // adjust connections for non uncancelled cps in msc
//      for(uint j = 0 ; j < m_cps.size();++j)
//      {
//        critpt_t * src_cp = m_cps[j];
//
//        if(src_cp->is_paired)
//          continue;
//
//        if(msc->m_id_cp_map.count(src_cp->cellid) != 1)
//          continue;
//
//        critpt_t *dest_cp = msc->m_cps[msc->m_id_cp_map[src_cp->cellid]];
//
//        conn_t *src_acdc[] = {&src_cp->asc,&src_cp->des};
//        conn_t *dest_acdc[] = {&dest_cp->asc,&dest_cp->des};
//
//        for(uint k = 0 ; k < 2;++k)
//        {
//          dest_acdc[k]->clear();
//
//          for(conn_iter_t it = src_acdc[k]->begin(); it!=src_acdc[k]->end();++it)
//          {
//            critpt_t *src_conn_cp = m_cps[*it];
//
//            if(src_conn_cp->is_paired == true)
//              throw std::logic_error("only non cancellable cps must be remaining 1");
//
//
//            if(msc->m_id_cp_map.count(src_conn_cp->cellid) == 0)
//              continue;
//
//            dest_acdc[k]->insert(msc->m_id_cp_map[src_conn_cp->cellid]);
//          }// end it
//        }// end k
//      }// end j
//    }//end i
//  }

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
      if(cp[dir]->isCancelled)
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

      std::cout
             <<   "no = "<<num_cancellations<<" "
             << "pers = "<<persistence<<" "
             <<"index = ("<<(int)m_cps[pr[0]]->index<<","<<(int)m_cps[pr[1]]->index<<") "
             << "edge = "<<edge_to_string(this,pr)<<" "
             <<std::endl;

      uint_pair_list_t new_edges;

      cancelPairs ( this,pr,&new_edges);
      num_cancellations++;

      mark_cancel_pair(this,pr);

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

      ensure_ordered_index_one_separation(this,e);

      ensure_pairing(this,e);

      for(uint dir = 0 ; dir < 2 ;++dir)
      {
        bool need_disc = false;

        for(conn_iter_t it  = cp[dir]->conn[dir].begin();
                        it != cp[dir]->conn[dir].end(); ++it)
        {
          ensure_cp_is_not_paired(this,*it);

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

    if(os.is_open() == false)
      throw std::runtime_error(std::string("failed to open file :")+fn);

    write_manifolds((std::ostream&)os);

    os.close();
  }

  cellid_t get_cp_cellid(mscomplex_t *msc,uint idx)
  {
    return msc->m_cps[idx]->cellid;
  }

  void mscomplex_t::write_graph(std::ostream & os)
  {
    const char *conn_txt[] = {"des","asc"};

    for(uint i = 0 ; i < m_cps.size();++i)
    {
      critpt_t * cp = m_cps[i];

      if(cp->is_paired)
        continue;

      if(cp->is_paired == true)
        continue;

      for(uint dir = 0 ; dir <2 ;++dir)
      {
        os<<conn_txt[dir]<<" "<<cp->cellid<<" "<<cp->conn[dir].size()<<"\n";

        std::transform(cp->conn[dir].begin(),cp->conn[dir].end(),
                       std::ostream_iterator<cellid_t>(os,"\n"),
                       boost::bind1st(get_cp_cellid,this));

        os<<std::endl;
      }
      os<<std::endl;
    }
  }

  void mscomplex_t::write_graph(const std::string &fn)
  {
    std::fstream os(fn.c_str(),std::ios::out);

    if(os.is_open() == false)
      throw std::runtime_error(std::string("failed to open file :")+fn);

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
