#ifndef GRID_MSCOMPLEX_ENSURE_H_INCLUDED
#define GRID_MSCOMPLEX_ENSURE_H_INCLUDED

#include <grid_mscomplex.h>

// bunch of predicates that throw when I suspect something could be logically
// wrong with the state of the MS complex.. to be disabled in release builds

namespace grid
{

  inline void ensure_index_one_separation(mscomplex_t *msc,uint_pair_t e)
  {
    if((msc->m_cps[e[0]]->index +1 != msc->m_cps[e[1]]->index)&&
       (msc->m_cps[e[1]]->index +1 != msc->m_cps[e[0]]->index))
      throw std::logic_error("index one separation violated");
  }

  inline void ensure_ordered_index_one_separation(mscomplex_t *msc,uint_pair_t e)
  {
    if(msc->m_cps[e[0]]->index +1 != msc->m_cps[e[1]]->index)
      throw std::logic_error("ordered index one separation violated");
  }

  inline void ensure_connectivity(mscomplex_t *msc,uint_pair_t e)
  {
    ensure_ordered_index_one_separation(msc,e);

    if(msc->m_cps[e[0]]->asc.count(e[1]) == 0 ||
       msc->m_cps[e[1]]->des.count(e[0]) == 0)
      throw std::logic_error("connectivity violated");
  }

  inline void ensure_single_connectivity(mscomplex_t *msc,uint_pair_t e)
  {
    ensure_ordered_index_one_separation(msc,e);

    if(msc->m_cps[e[0]]->asc.count(e[1]) != 1 ||
       msc->m_cps[e[1]]->des.count(e[0]) != 1)
      throw std::logic_error("single connectivity violated");
  }

  inline void ensure_cellid_critical(mscomplex_t * msc,cellid_t c)
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
}
#endif
