const sampler_t func_sampler  = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;
const sampler_t flag_sampler  = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;

inline bool compare_verts
( __read_only image3d_t func_img,
  __read_only image3d_t flag_img,
  cell_t v1, cell_t v2)
{
  int4 p1 = {(v1.x/2),(v1.y/2),(v1.z/2),(0)};
  int4 p2 = {(v2.x/2),(v2.y/2),(v2.z/2),(0)};

  func_t f1 = read_imagef(func_img, func_sampler, p1).x;
  func_t f2 = read_imagef(func_img, func_sampler, p2).x;

  if( f1 != f2)
    return f1 < f2;

  return __compare_cells(v1,v2);
}

//inline bool compare_edges
//( __read_only image3d_t func_img,
//  __read_only image3d_t flag_img,
//  cell_t e1, cell_t e2)
//{
//  flag_t flg1 = read_imageui(flag_img, flag_sampler, cell_to_int4(e1)).x;
//  flag_t flg2 = read_imageui(flag_img, flag_sampler, cell_to_int4(e2)).x;

//  cell_t v1 = flag_to_mxfct(e1,flg1);
//  cell_t v2 = flag_to_mxfct(e2,flg2);

//  if( __cell_equal(v1,v2))
//  {
//    v1 = second_max_facet(e1,v1);
//    v2 = second_max_facet(e2,v2);
//  }

//  return compare_verts(func_img,flag_img,v1,v2);
//}

//inline bool compare_faces(const dataset_t ds,cell_t f1, cell_t f2)
//{
//  int4 if1 = {(f1.x),(f1.y),(f1.z),(0)};
//  int4 if2 = {(f2.x),(f2.y),(f2.z),(0)};

//  flag_t flg1 = read_imageui(ds.m_func, flag_sampler, if1).x;
//  flag_t flg2 = read_imageui(ds.m_func, flag_sampler, if2).x;

//  cell_t e1 = flag_to_mxfct(f1,flg1);
//  cell_t e2 = flag_to_mxfct(f2,flg2);

//  if( __cell_equal(e1,e2))
//  {
//    e1 = second_max_facet(f1,e1);
//    e2 = second_max_facet(f2,e2);
//  }

//  return compare_edges(ds,e1,e2);
//}

//inline bool compare_cubes(const dataset_t ds,cell_t c1, cell_t c2)
//{
//  int4 ic1 = {(c1.x),(c1.y),(c1.z),(0)};
//  int4 ic2 = {(c2.x),(c2.y),(c2.z),(0)};

//  flag_t flg1 = read_imageui(ds.m_func, flag_sampler, ic1).x;
//  flag_t flg2 = read_imageui(ds.m_func, flag_sampler, ic2).x;

//  cell_t f1 = flag_to_mxfct(c1,flg1);
//  cell_t f2 = flag_to_mxfct(c2,flg2);

//  if( __cell_equal(f1,f2))
//  {
//    f1 = second_max_facet(c1,f1);
//    f2 = second_max_facet(c2,f2);
//  }

//  return compare_faces(ds,f1,f2);
//}

__kernel void assign_max_facet_edge
(
  __read_only  image3d_t  g_func,
  __read_only image3d_t   g_flag,
  cell_t ext_rect_lc,
  cell_t ext_rect_uc,
  cell_t edir,
  __global flag_t   * g_flag_out
)
{
  int tid      = get_global_id(0);
  int num_thds = get_global_size(0);

  rect_t  ext_rect = make_rect(ext_rect_lc,ext_rect_uc);
  rect_t cell_rect = make_rect(ext_rect_lc+edir,ext_rect_uc-edir);

  int N = num_cells2(cell_rect);

  for( int i = tid ; i < N; i += num_thds)
  {
    cell_t e = i_to_c2(cell_rect,i);

    cell_t v1 = e - edir;
    cell_t v2 = e + edir;

    cell_t v = v1;

    if( compare_verts(g_func,g_flag,v,v2))
      v = v2;

    g_flag_out[c_to_i(ext_rect,e)] = mxfct_to_flag(e,v);
  }
}

__kernel void assign_zero
(
  cell_t ext_rect_lc,
  cell_t ext_rect_uc,
  __global flag_t   * g_flag_out
)
{
  int tid      = get_global_id(0);
  int num_thds = get_global_size(0);

  rect_t  ext_rect = make_rect(ext_rect_lc,ext_rect_uc);
  int N = num_cells(ext_rect);

  for( int i = tid ; i < N; i += num_thds)
  {
    g_flag_out[i] = 0;
  }
}
