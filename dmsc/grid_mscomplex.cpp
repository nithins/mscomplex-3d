#include <cmath>
#include <queue>

#include <grid_mscomplex.h>
#include <grid_dataset.h>
#include <fstream>
#include <limits>

namespace grid
{

  inline void order_pr_by_cp_index(mscomplex_t *msc,uint_pair_t &e)
  {
    if(msc->m_cps[e[0]]->index > msc->m_cps[e[1]]->index)
      std::swap(e[0],e[1]);
  }

  inline void ensure_index_one_separation(mscomplex_t *msc,uint_pair_t e)
  {
    if(msc->m_cps[e[0]]->index +1 != msc->m_cps[e[1]]->index)
      throw std::logic_error("index one separation violated");
  }

  inline void ensure_connectivity(mscomplex_t *msc,uint_pair_t e)
  {
    ensure_index_one_separation(msc,e);

    if(msc->m_cps[e[0]]->asc.count(e[1]) == 0 ||
       msc->m_cps[e[1]]->des.count(e[0]) == 0)
      throw std::logic_error("connectivity violated");
  }

  inline void ensure_single_connectivity(mscomplex_t *msc,uint_pair_t e)
  {
    ensure_index_one_separation(msc,e);

    if(msc->m_cps[e[0]]->asc.count(e[1]) != 1 ||
       msc->m_cps[e[1]]->des.count(e[0]) != 1)
      throw std::logic_error("single connectivity violated");
  }

  void ensure_cellid_critical(mscomplex_t * msc,cellid_t c)
  {
    if(msc->m_id_cp_map.count(c) == 0)
      throw std::logic_error("cellid not entered as critical in msc");
  }

  inline void ensure_cp_is_cancelled(mscomplex_t *msc,uint i)
  {
    if(!msc->m_cps[i]->isCancelled)
      throw std::logic_error("failed to ensure cp is not canceled");
  }

  inline void ensure_cp_is_not_cancelled(mscomplex_t *msc,uint i)
  {
    if(msc->m_cps[i]->isCancelled)
      throw std::logic_error("failed to ensure cp is canceled");
  }

  inline void ensure_pairing(mscomplex_t *msc,uint_pair_t e)
  {
    if(!msc->m_cps[e[0]]->is_paired||
       !msc->m_cps[e[1]]->is_paired||
       msc->m_cps[e[1]]->pair_idx != e[0]||
       msc->m_cps[e[0]]->pair_idx != e[1])
      throw std::logic_error("failed to ensure that edge forms a sane pairing ");
  }

  inline void ensure_cp_is_paired(mscomplex_t *msc,uint c)
  {
    if(!msc->m_cps[c]->is_paired)
      throw std::logic_error("failed to ensure cell is paired ");

    ensure_pairing(msc,uint_pair_t(c,msc->m_cps[c]->pair_idx));
  }

  inline void ensure_cp_is_not_paired(mscomplex_t *msc,uint c)
  {
    if(msc->m_cps[c]->is_paired)
      throw std::logic_error("failed to ensure cell is not paired ");
  }

  void mark_cancel_pair(mscomplex_t *msc,uint_pair_t e)
  {
    critpt_t * cp1 = msc->m_cps[e[0]];
    critpt_t * cp2 = msc->m_cps[e[1]];

    cp1->is_paired = true;
    cp2->is_paired = true;

    cp1->pair_idx  = e[1];
    cp2->pair_idx  = e[0];

  }

  void mscomplex_t::add_critpt(cellid_t c,uchar i,cell_fn_t f)
  {
    critpt_t * cp = new critpt_t;
    cp->cellid    = c;
    cp->index     = i;
    cp->fn        = f;
    m_id_cp_map.insert(std::make_pair(c,m_cps.size()));
    m_cps.push_back(cp);
  }

  void mscomplex_t::connect_cps(cellid_t c0,cellid_t c1)
  {
    ensure_cellid_critical(this,c0);
    ensure_cellid_critical(this,c1);

    uint_pair_t e(m_id_cp_map[c0],m_id_cp_map[c1]);

    order_pr_by_cp_index(this,e);

    ensure_index_one_separation(this,e);

    m_cps[e[0]]->asc.insert(e[1]);
    m_cps[e[1]]->des.insert(e[0]);
  }

  void cancelPairs ( mscomplex_t *msc,uint_pair_t e,
                     uint_pair_list_t * new_edges = NULL)
  {

    order_pr_by_cp_index(msc,e);

    ensure_single_connectivity(msc,e);

    critpt_t * cp0 = msc->m_cps[e[0]];
    critpt_t * cp1 = msc->m_cps[e[1]];

    for ( conn_iter_t asc0_it = cp0->asc.begin();asc0_it != cp0->asc.end(); ++asc0_it )
    {
      if(*asc0_it == e[1])
        continue;

      for ( conn_iter_t des1_it = cp1->des.begin();des1_it != cp1->des.end(); ++des1_it )
      {
        if(*des1_it == e[0])
          continue;

        msc->m_cps[*asc0_it]->des.insert(*des1_it);
        msc->m_cps[*des1_it]->asc.insert(*asc0_it);

        if(new_edges != NULL)
          new_edges->push_back(uint_pair_t(*des1_it,*asc0_it ));
      }
    }

    // merge here the des disc of cp1_ind to the des disc of asc0_it;
    for ( conn_iter_t asc0_it  = cp0->asc.begin(); asc0_it != cp0->asc.end(); ++asc0_it )
      if ( *asc0_it != e[1] )
        msc->m_cps[*asc0_it]->des.erase ( e[0] );

    for ( conn_iter_t  des0_it = cp0->des.begin();des0_it != cp0->des.end(); ++des0_it )
      msc->m_cps[*des0_it]->asc.erase ( e[0] );

    for ( conn_iter_t asc1_it  = cp1->asc.begin();asc1_it != cp1->asc.end(); ++asc1_it )
      msc->m_cps[*asc1_it]->des.erase ( e[1] );

    // merge here the asc disc of cp0_ind to the asc disc of des1_it;
    for ( conn_iter_t des1_it = cp1->des.begin(); des1_it != cp1->des.end(); ++des1_it )
      if ( *des1_it!= e[0] )
        msc->m_cps[*des1_it]->asc.erase ( e[1] );


    cp0->isCancelled = true;
    cp0->des.clear();

    cp1->isCancelled = true;
    cp1->asc.clear();
  }

  void uncancel_pairs( mscomplex_t  *msc,uint_pair_t e)
  {
    order_pr_by_cp_index(msc,e);

    ensure_connectivity(msc,e);

    conn_t * ad_conns [] = {&msc->m_cps[e[0]]->asc,&msc->m_cps[e[1]]->des};

    for(uint ad = 0 ; ad <2 ; ++ad)
    {
      conn_t new_ad;

      for(conn_iter_t i_it = ad_conns[ad]->begin(); i_it != ad_conns[ad]->end() ; ++i_it )
      {
        if(*i_it == e[(ad+1)%2]) continue;

        critpt_t *cp = msc->m_cps[*i_it];

        if(!cp->is_paired)
        {
          new_ad.insert(*i_it);
          continue;
        }

        ensure_cp_is_not_cancelled(msc,*i_it);

        ensure_cp_is_paired(msc,*i_it);

        critpt_t *cp_pr = msc->m_cps[cp->pair_idx];

        conn_t *pr_ad_conns = (ad == 0)?(&cp_pr->asc):(&cp_pr->des);

        for(conn_iter_t j_it = pr_ad_conns->begin(); j_it != pr_ad_conns->end() ; ++j_it )
        {
          ensure_cp_is_not_paired(msc,*j_it);

          new_ad.insert(*j_it);
        }
      }

      msc->m_cps[e[ad]]->isCancelled = false;
      ad_conns[ad]->clear();
      ad_conns[ad]->insert(new_ad.begin(),new_ad.end());
    }
  }

  void print_cp_connections(std::ostream & os,const mscomplex_t &msc,
                            const conn_t &conn)
  {

    os<<"{ ";
    for(conn_iter_t it = conn.begin(); it != conn.end(); ++it)
    {
      if(msc.m_cps[*it]->is_paired)
        os<<"*";
      os<<msc.m_cps[*it]->cellid;
      os<<", ";
    }
    os<<"}";
  }


  void mscomplex_t::print_connections(std::ostream & os)
  {
    for(uint i = 0 ; i < m_cps.size();++i)
    {
      os<<"des(";
      if(m_cps[i]->is_paired)
        os<<"*";
      os<<m_cps[i]->cellid<<") = ";
      print_cp_connections(os,*this,m_cps[i]->des);
      os<<std::endl;

      os<<"asc(";
      if(m_cps[i]->is_paired)
        os<<"*";
      os<<m_cps[i]->cellid<<") = ";
      print_cp_connections(os,*this,m_cps[i]->asc);
      os<<std::endl;
      os<<std::endl;
    }
  }

  mscomplex_t::~mscomplex_t()
  {
    clear();
  }

  void shallow_replicate_cp(mscomplex_t &msc,const critpt_t &cp)
  {
    if(msc.m_id_cp_map.count(cp.cellid) != 0)
      throw std::logic_error("this cp is present in msc");

    critpt_t * dest_cp = new critpt_t;

    dest_cp->isCancelled              = cp.isCancelled;
    dest_cp->is_paired                = cp.is_paired;
    dest_cp->cellid                   = cp.cellid;
    dest_cp->fn                       = cp.fn;
    dest_cp->index                    = cp.index;

    msc.m_id_cp_map[dest_cp->cellid]  = msc.m_cps.size();
    msc.m_cps.push_back(dest_cp);
  }

  mscomplex_t * mscomplex_t::merge_up
      (const mscomplex_t& msc1,
       const mscomplex_t& msc2)
  {

    // form the intersection rect
    rect_t ixn;

    if (!msc2.m_rect.intersection (msc1.m_rect,ixn))
      throw std::logic_error ("rects should intersect for merge");

    if ( msc2.m_rect.eff_dim() != gc_grid_dim-1)
      throw std::logic_error ("rects must merge along a d-1 manifold");

    // TODO: ensure that the union of  rects is not including anything extra

    rect_t r = msc1.m_rect.bounding_box(msc1.m_rect);

    rect_t e = msc1.m_ext_rect.bounding_box(msc1.m_ext_rect);

    mscomplex_t * out_msc = new mscomplex_t(r,e);

    const mscomplex_t* msc_arr[] = {&msc1,&msc2};

    // make a union of the critical points in this
    for (uint i = 0 ; i <2;++i)
    {
      const mscomplex_t * msc = msc_arr[i];

      for (uint j = 0 ; j <msc->m_cps.size();++j)
      {
        const critpt_t *src_cp = msc->m_cps[j];

        // if it is contained or not
        if (i == 1 && (out_msc->m_id_cp_map.count(src_cp->cellid) == 1))
          continue;

        if(src_cp->isCancelled)
          continue;

        shallow_replicate_cp(*out_msc,*src_cp);

      }
    }

    for (uint i = 0 ; i <2;++i)
    {
      const mscomplex_t * msc = msc_arr[i];

      // copy over connectivity information
      for (uint j = 0 ; j <msc->m_cps.size();++j)
      {
        const critpt_t *src_cp = msc->m_cps[j];

        if(src_cp->isCancelled)
          continue;

        critpt_t *dest_cp = out_msc->m_cps[out_msc->m_id_cp_map[src_cp->cellid]];

        if(src_cp->is_paired)
        {
          critpt_t *src_pair_cp = msc->m_cps[src_cp->pair_idx];

          dest_cp->pair_idx = out_msc->m_id_cp_map[src_pair_cp->cellid];
        }

        const conn_t *acdc_src[]  = {&src_cp->asc, &src_cp->des};

        conn_t *acdc_dest[] = {&dest_cp->asc,&dest_cp->des};

        bool is_src_cmn_bndry = (ixn.contains(src_cp->cellid) && i == 1);

        for (uint j = 0 ; j < 2; ++j)
        {
          for (const_conn_iter_t it = acdc_src[j]->begin();
          it != acdc_src[j]->end();++it)
          {
            const critpt_t *conn_cp = msc->m_cps[*it];

            // common boundry connections would have been found along the boundry
            if( is_src_cmn_bndry && ixn.contains(conn_cp->cellid))
              continue;

            if (conn_cp->isCancelled)
              continue;

            acdc_dest[j]->insert (out_msc->m_id_cp_map[conn_cp->cellid]);
          }
        }
      }
    }

    static_assert(gc_grid_dim == 3&&"defined for 3-manifolds only");

    cellid_t c;

    for(c[2] = ixn[2][0] ; c[2] <= ixn[2][1]; ++c[2])
    {
      for(c[1] = ixn[1][0] ; c[1] <= ixn[1][1]; ++c[1])
      {
        for(c[0] = ixn[0][0] ; c[0] <= ixn[0][1]; ++c[0])
        {

          if(out_msc->m_id_cp_map.count(c) != 1)
            throw std::logic_error("missing common bndry cp");

          u_int src_idx = out_msc->m_id_cp_map[c];

          critpt_t *src_cp = out_msc->m_cps[src_idx];

          if(src_cp->is_paired || !src_cp->is_paired)
            continue;

          u_int pair_idx = src_cp->pair_idx;

          cellid_t p = out_msc->m_cps[pair_idx]->cellid;

          if(!out_msc->m_rect.isInInterior(c)&& !out_msc->m_ext_rect.isOnBoundry(c))
            continue;

          if(!out_msc->m_rect.isInInterior(p)&& !out_msc->m_ext_rect.isOnBoundry(p))
            continue;

          cancelPairs(out_msc,uint_pair_t(src_idx,pair_idx));
        }
      }
    }

    return out_msc;
  }

  void mscomplex_t::merge_down(mscomplex_t& msc1,mscomplex_t& msc2)
  {
    // form the intersection rect
    rect_t ixn;

    if (!msc2.m_rect.intersection (msc1.m_rect,ixn))
      throw std::logic_error ("rects should intersect for merge");

    if ( msc2.m_rect.eff_dim() != gc_grid_dim-1)
      throw std::logic_error ("rects must merge along a d-1 manifold");

    static_assert(gc_grid_dim == 3&&"defined for 3-manifolds only");

    cellid_t c;

    for(c[2] = ixn[2][1] ; c[2] >= ixn[2][0]; --c[2])
    {
      for(c[1] = ixn[1][1] ; c[1] >= ixn[1][0]; --c[1])
      {
        for(c[0] = ixn[0][1] ; c[0] >= ixn[0][0]; --c[0])
        {

          if(this->m_id_cp_map.count(c) != 1)
            throw std::logic_error("missing common bndry cp");

          u_int src_idx = this->m_id_cp_map[c];

          critpt_t *src_cp = this->m_cps[src_idx];

          if(!src_cp->isCancelled )
            continue;

          u_int pair_idx = src_cp->pair_idx;

          cellid_t p = this->m_cps[pair_idx]->cellid;

          if(!this->m_rect.isInInterior(c)&& !this->m_ext_rect.isOnBoundry(c))
            continue;

          if(!this->m_rect.isInInterior(p)&& !this->m_ext_rect.isOnBoundry(p))
            continue;

          uncancel_pairs(this,uint_pair_t(src_idx,pair_idx));
        }
      }
    }

    // identify and copy the results to msc1 and msc2

    mscomplex_t* msc_arr[] = {&msc1,&msc2};

    for (uint i = 0 ; i <2;++i)
    {
      mscomplex_t * msc = msc_arr[i];

      // adjust connections for uncancelled cps in msc
      for(uint j = 0 ; j < m_cps.size();++j)
      {
        critpt_t * src_cp = m_cps[j];

        if(src_cp->isCancelled)
          throw std::logic_error("all cps ought to be uncancelled by now");

        if(!src_cp->is_paired)
          continue;

        critpt_t * src_pair_cp = m_cps[src_cp->pair_idx];

        bool src_in_msc      = (msc->m_id_cp_map.count(src_cp->cellid) != 0);
        bool src_pair_in_msc = (msc->m_id_cp_map.count(src_pair_cp->cellid) != 0);

        if(!src_in_msc && !src_pair_in_msc)
          continue;

        if(!src_in_msc)
        {
          shallow_replicate_cp(*msc,*src_cp);
        }

        if(!src_pair_in_msc)
        {
          shallow_replicate_cp(*msc,*src_pair_cp);
        }

        uint dest_cp_idx = msc->m_id_cp_map[src_cp->cellid];

        critpt_t *dest_cp = msc->m_cps[dest_cp_idx];

        if(!src_in_msc || !src_pair_in_msc || !dest_cp->is_paired)
        {
          uint dest_pair_cp_idx  = msc->m_id_cp_map[src_pair_cp->cellid];
          critpt_t *dest_pair_cp = msc->m_cps[dest_pair_cp_idx];

          dest_cp->is_paired      = true;
          dest_pair_cp->is_paired = true;

          dest_cp->pair_idx      = dest_pair_cp_idx;
          dest_pair_cp->pair_idx = dest_cp_idx;
        }

        conn_t *src_acdc[] = {&src_cp->asc,&src_cp->des};
        conn_t *dest_acdc[] = {&dest_cp->asc,&dest_cp->des};

        for(uint k = 0 ; k < 2;++k)
        {
          dest_acdc[k]->clear();

          for(conn_iter_t it = src_acdc[k]->begin(); it!=src_acdc[k]->end();++it)
          {
            critpt_t *src_conn_cp = m_cps[*it];

            if(src_conn_cp->is_paired == true)
              throw std::logic_error("only non cancellable cps must be remaining");

            if(msc->m_id_cp_map.count(src_conn_cp->cellid) == 0)
            {
              shallow_replicate_cp(*msc,*src_conn_cp);
            }

            dest_acdc[k]->insert(msc->m_id_cp_map[src_conn_cp->cellid]);
          }// end it
        }// end k
      }// end j

      // adjust connections for non uncancelled cps in msc
      for(uint j = 0 ; j < m_cps.size();++j)
      {
        critpt_t * src_cp = m_cps[j];

        if(src_cp->is_paired)
          continue;

        if(msc->m_id_cp_map.count(src_cp->cellid) != 1)
          continue;

        critpt_t *dest_cp = msc->m_cps[msc->m_id_cp_map[src_cp->cellid]];

        conn_t *src_acdc[] = {&src_cp->asc,&src_cp->des};
        conn_t *dest_acdc[] = {&dest_cp->asc,&dest_cp->des};

        for(uint k = 0 ; k < 2;++k)
        {
          dest_acdc[k]->clear();

          for(conn_iter_t it = src_acdc[k]->begin(); it!=src_acdc[k]->end();++it)
          {
            critpt_t *src_conn_cp = m_cps[*it];

            if(src_conn_cp->is_paired == true)
              throw std::logic_error("only non cancellable cps must be remaining 1");


            if(msc->m_id_cp_map.count(src_conn_cp->cellid) == 0)
              continue;

            dest_acdc[k]->insert(msc->m_id_cp_map[src_conn_cp->cellid]);
          }// end it
        }// end k
      }// end j
    }//end i
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

    bool operator()(const uint_pair_t & p1, const uint_pair_t &p2)
    {
      cell_fn_t f1 = m_msc->m_cps[p1[0]]->fn;
      cell_fn_t f2 = m_msc->m_cps[p1[1]]->fn;
      cell_fn_t f3 = m_msc->m_cps[p2[0]]->fn;
      cell_fn_t f4 = m_msc->m_cps[p2[1]]->fn;

      cell_fn_t d1 = std::abs(f2-f1);
      cell_fn_t d2 = std::abs(f4-f3);

      if(d1 != d2)
        return d1>d2;

      cellid_t c1 = m_msc->m_cps[p1[0]]->cellid;
      cellid_t c2 = m_msc->m_cps[p1[1]]->cellid;

      cellid_t c3 = m_msc->m_cps[p1[0]]->cellid;
      cellid_t c4 = m_msc->m_cps[p1[1]]->cellid;

      d1 = (c1-c2)*(c1-c2);
      d2 = (c3-c4)*(c3-c4);

      if(d1 != d2)
        return d1>d2;

      if(c1 > c2)
        std::swap(c1,c2);

      if(c3 > c4)
        std::swap(c3,c4);

      if(c1 != c3)
        return c1 > c3;

      return c2 > c4;
    }

  };

  bool is_valid_canc_edge(mscomplex_t *msc,uint_pair_t e )
  {
    order_pr_by_cp_index(msc,e);

    critpt_t *cp0 = msc->m_cps[e[0]];
    critpt_t *cp1 = msc->m_cps[e[1]];

    if(cp0->isCancelled ||
       cp1->isCancelled )
      return false;

    if((msc->m_rect.isOnBoundry(cp1->cellid) && !msc->m_rect.isOnBoundry(cp0->cellid))||
       (msc->m_rect.isOnBoundry(cp0->cellid) && !msc->m_rect.isOnBoundry(cp1->cellid)))
      return false;

    if(cp0->asc.count(e[1]) != 1 ||
       cp1->des.count(e[0]) != 1)
      return false;

    return true;
  }

  void mscomplex_t::simplify(uint_pair_list_t & canc_pairs_list,
                             double simplification_treshold)
  {
    typedef std::priority_queue
        <uint_pair_t,uint_pair_list_t,persistence_comparator_t>
        canc_pair_priq_t;

    persistence_comparator_t comp(this);

    canc_pair_priq_t  canc_pair_priq(comp);

    // add every edge in the descending manifold of the critical point

    cell_fn_t max_persistence = 0.0;

    cell_fn_t max_val = std::numeric_limits<cell_fn_t>::min();
    cell_fn_t min_val = std::numeric_limits<cell_fn_t>::max();

    for(uint i = 0 ;i < m_cps.size();++i)
    {
      critpt_t *cp = m_cps[i];

      max_val = std::max(max_val,m_cps[i]->fn);

      min_val = std::min(min_val,m_cps[i]->fn);

      for(const_conn_iter_t it = cp->des.begin();it != cp->des.end() ;++it)
        if(is_valid_canc_edge(this,uint_pair_t(i,*it)))
          canc_pair_priq.push(uint_pair_t(i,*it));
    }

    max_persistence = max_val - min_val;

    while (canc_pair_priq.size() !=0)
    {
      uint_pair_t pr = canc_pair_priq.top();

      canc_pair_priq.pop();

      if(is_valid_canc_edge(this,pr) == false)
        continue;

      cell_fn_t persistence = std::abs(m_cps[pr[0]]->fn-m_cps[pr[1]]->fn);

      if((double)persistence/(double)max_persistence > simplification_treshold)
        break;

      std::vector<uint_pair_t> new_edges;

      cancelPairs ( this,pr,&new_edges);

      mark_cancel_pair(this,pr);

      canc_pairs_list.push_back(pr);

      for(uint i = 0 ; i < new_edges.size(); i++)
        canc_pair_priq.push(new_edges[i]);
    }
  }

  void mscomplex_t::un_simplify(const uint_pair_list_t &canc_pairs_list)
  {
    for(uint_pair_list_t::const_reverse_iterator it = canc_pairs_list.rbegin();
    it != canc_pairs_list.rend() ; ++it)
      uncancel_pairs(this,*it);
  }

  void mscomplex_t::simplify_un_simplify(double simplification_treshold)
  {
    uint_pair_list_t canc_pairs_list;

    simplify(canc_pairs_list,simplification_treshold);

    un_simplify(canc_pairs_list);
  }

  void write_disc(const critpt_disc_t *disc,
                  const std::string &prefix,
                  const cellid_t &c )
  {

    std::stringstream ss;
    ss<<prefix;
    ((std::ostream&)ss)<<c;

    std::ofstream os;

    os.open(ss.str().c_str(),std::ios::out|std::ios::ate|std::ios::app|std::ios::binary);


    if(os.is_open() == false)
      throw std::runtime_error("asc/des disc file not writeable");
    os.write((const char*)(const void*)disc->data(),disc->size()*sizeof(cellid_t));

    os.close();
  }

  void mscomplex_t::write_discs(const std::string &fn_prefix)
  {
    critpt_t * cp ;

    for(uint i = 0 ; i < m_cps.size();++i)
    {
      cp = m_cps[i];

      if(cp->is_paired == true)
        continue;

      if(cp->asc_disc.size() != 0 )
        write_disc(&cp->asc_disc,fn_prefix+"asc_",cp->cellid);

      if(cp->des_disc.size() != 0 )
        write_disc(&cp->des_disc,fn_prefix+"des_",cp->cellid);
    }
  }
}


#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/binary_object.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

using namespace grid;

namespace boost
{
  namespace serialization
  {
    template<class Archive>
    void serialize(Archive & ar, rect_range_t & r, const unsigned int )
    {
      typedef boost::array<rect_range_t::value_type,rect_range_t::static_size>
          rect_range_base_t;

      ar & boost::serialization::base_object<rect_range_base_t>(r);
    }

    template<class Archive>
    void serialize(Archive & ar, rect_point_t & p, const unsigned int )
    {
      typedef boost::array<rect_point_t::value_type,rect_point_t::static_size>
          rect_point_base_t;

      ar & boost::serialization::base_object<rect_point_base_t>(p);
    }

    template<class Archive>
    void serialize(Archive & ar, rect_t & r, const unsigned int )
    {
      typedef boost::array<rect_t::value_type,rect_t::static_size>
          rect_base_t;

      ar & boost::serialization::base_object<rect_base_t>(r);
    }

    template<class Archive>
    void serialize(Archive & ar, critpt_t & c, const unsigned int )
    {
      ar & c.cellid;
      ar & c.asc;
      ar & c.des;
      ar & c.is_paired;
      ar & c.isCancelled;
      ar & c.pair_idx;
      ar & c.fn;
    }


    template<class Archive>
    void serialize(Archive & ar, mscomplex_t & g, const unsigned int )
    {
      ar & g.m_rect;
      ar & g.m_ext_rect;
      ar & g.m_id_cp_map;
      ar & g.m_cps;
    }

    //    template<class Archive>
    //    void serialize(Archive & ar, GridDataset & ds, const unsigned int )
    //    {
    //       ar & ds.m_rect;
    //       ar & ds.m_ext_rect;
    //
    //       GridDataset::rect_size_t ext_sz = ds.m_ext_rect.size();
    //       uint num_data_items = (ext_sz[0]+1)*(ext_sz[1]+1);
    //
    //       if(Archive::is_loading::value)
    //         ds.init(NULL);
    //
    //       ar & make_binary_object(ds.(*m_cell_flags).data(),num_data_items*sizeof(GridDataset::cell_flag_t));
    //       ar & make_binary_object(ds.m_cell_pairs.data(),num_data_items*sizeof(GridDataset::cellid_t));
    //    }
  }
}

//// without the explicit instantiations below, the program will
//// fail to link for lack of instantiantiation of the above function
//// The impls are visible only in this file to save compilation time..
//
//template void boost::serialization::serialize<boost::archive::text_iarchive>(
//    boost::archive::text_iarchive & ar,
//    GridDataset & g,
//    const unsigned int file_version
//);
//
//template void boost::serialization::serialize<boost::archive::text_oarchive>(
//    boost::archive::text_oarchive & ar,
//    GridDataset & g,
//    const unsigned int file_version
//);


// without the explicit instantiations below, the program will
// fail to link for lack of instantiantiation of the above function
// The impls are visible only in this file to save compilation time..

template void boost::serialization::serialize<boost::archive::binary_iarchive>(
    boost::archive::binary_iarchive & ar,
    mscomplex_t & g,
    const unsigned int file_version
    );
template void boost::serialization::serialize<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive & ar,
    mscomplex_t & g,
    const unsigned int file_version
    );


