inline void init
  ( dataset_t ds,
    rect_t ex_rect,
    __read_only image3d_t flag_img,
    __global int * ex_own_buf
  )
{
  int N = num_cells2(ex_rect);

  for( int i = get_global_id(0) ; i < N; i += get_global_size(0))
  {
    cell_t c  = i_to_c2(ex_rect,i);
    flag_t fg = read_imageui(flag_img, flag_sampler, to_int4(c-ds.e.lc)).x;

    cell_t o = c;

    if(!is_critical(fg))
    {
      cell_t p = flag_to_pair(c,fg);
      o  = p + p-c;
    }

    ex_own_buf[c_to_i2(ex_rect,c)] = c_to_i2(ex_rect,o);
  }
}
inline void propagate
  ( dataset_t ds,
    rect_t ex_rect,
    __global int * ex_own_buf1,
    __global int * ex_own_buf2,
    __global int * is_updated
  )
{
  int N = num_cells2(ex_rect);

  int updated = 0;

  for( int i = get_global_id(0) ; i < N; i += get_global_size(0))
  {
    int  o = ex_own_buf1[i];
    int oo = ex_own_buf1[o];
    ex_own_buf2[i] = oo;

    if(o != oo)
      updated = 1;
  }

  if(updated == 1)
    *is_updated = 1;
}
inline void finalize
  ( dataset_t ds,
    rect_t ex_rect,
    __read_only image3d_t  flag_img,
    __global int * ex_own_buf1,
    __global int * ex_own_buf2)
{
  int N = num_cells2(ex_rect);

  for( int i = get_global_id(0) ; i < N; i += get_global_size(0))
  {
    cell_t c  = i_to_c2(ex_rect,i);
    flag_t fg = read_imageui(flag_img, flag_sampler, to_int4(c-ds.e.lc)).x;

    if(!is_critical(fg))
    {
      int  o = ex_own_buf1[i];
      int oo = ex_own_buf1[o];
      ex_own_buf2[i] = oo;
    }
  }
}

__kernel void init_minima
( cell_pair_t rct,
  cell_pair_t ext,
  cell_pair_t dom,
  __read_only image3d_t  flag_img,
  __global int    * ex_own_buf
)
{
  dataset_t ds = make_dataset2(rct,ext,dom);
  init(ds,ds.r,flag_img,ex_own_buf);
}
__kernel void propagate_minima
( cell_pair_t rct,
  cell_pair_t ext,
  cell_pair_t dom,
  __global int * own_buf1,
  __global int * own_buf2,
  __global int * is_updated

)
{
  dataset_t ds = make_dataset2(rct,ext,dom);
  propagate(ds,ds.r,own_buf1,own_buf2,is_updated);
}
__kernel void finalize_minima
( cell_pair_t rct,
  cell_pair_t ext,
  cell_pair_t dom,
  __read_only image3d_t  flag_img,
  __global int * own_buf1,
  __global int * own_buf2
)
{
  dataset_t ds = make_dataset2(rct,ext,dom);
  finalize(ds,ds.r,flag_img,own_buf1,own_buf2);
}




__kernel void init_maxima
( cell_pair_t rct,
  cell_pair_t ext,
  cell_pair_t dom,
  __read_only image3d_t  flag_img,
  __global int    * ex_own_buf
)
{
  dataset_t ds = make_dataset2(rct,ext,dom);
  rect_t ex_rect = make_rect(ds.r.lc+to_cell(1,1,1),ds.r.uc-to_cell(1,1,1));
  init(ds,ex_rect,flag_img,ex_own_buf);
}
__kernel void propagate_maxima
( cell_pair_t rct,
  cell_pair_t ext,
  cell_pair_t dom,
  __global int * own_buf1,
  __global int * own_buf2,
  __global int * is_updated
)
{
  dataset_t ds = make_dataset2(rct,ext,dom);
  rect_t ex_rect = make_rect(ds.r.lc+to_cell(1,1,1),ds.r.uc-to_cell(1,1,1));
  propagate(ds,ex_rect,own_buf1,own_buf2,is_updated);
}
__kernel void finalize_maxima
( cell_pair_t rct,
  cell_pair_t ext,
  cell_pair_t dom,
  __read_only image3d_t  flag_img,
  __global int * own_buf1,
  __global int * own_buf2
)
{
  dataset_t ds = make_dataset2(rct,ext,dom);
  rect_t ex_rect = make_rect(ds.r.lc+to_cell(1,1,1),ds.r.uc-to_cell(1,1,1));
  finalize(ds,ex_rect,flag_img,own_buf1,own_buf2);
}
