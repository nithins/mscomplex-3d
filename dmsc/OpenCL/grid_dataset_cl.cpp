#define __CL_ENABLE_EXCEPTIONS
#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.hpp>
#else
#include <cl.hpp>
#endif
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include <stdexcept>
#include <cpputils.h>

#include<grid_dataset_cl.h>

using namespace grid;
using namespace std;

typedef cl_short4 cl_cell_t;

inline cl_cell_t to_cell(const cellid_t & b)
{
  cl_cell_t a;

  a.x = b[0];
  a.y = b[1];
  a.z = b[2];

  return a;
}

inline cl::size_t<3> to_size(const cellid_t & b)
{
  cl::size_t<3> a;

  a.assign(b.begin(),b.end());

  return a;
}

inline cl::size_t<3> to_size(int x,int y,int z)
{
  cl::size_t<3> a;

  a[0] = x;
  a[1] = y;
  a[2] = z;

  return a;
}


cl::Context      s_context;
cl::CommandQueue s_queue;

cl::KernelFunctor s_assign_max_facet_edge;
cl::KernelFunctor s_assign_max_facet_face;
cl::KernelFunctor s_assign_max_facet_cube;
cl::KernelFunctor s_assign_pairs;
cl::KernelFunctor s_mark_cps;
cl::KernelFunctor s_save_cps;

cl::KernelFunctor s_scan_local_sums;
cl::KernelFunctor s_scan_group_sums;
cl::KernelFunctor s_scan_update_sums;

const char * s_header_file = "/home/nithin/projects/mscomplex-3d/dmsc/OpenCL/grid_dataset.clh";
const char * s_source_file = "/home/nithin/projects/mscomplex-3d/dmsc/OpenCL/grid_dataset.cl";

const int WI_SIZE = 256;
const int WG_NUM  = 32;
const int WG_SIZE = WG_NUM*WI_SIZE;

void init_opencl(void)
{
  cl::Program             program;
  std::vector<cl::Device> devices;

  try
  {
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);

    if (platforms.size() == 0)
        throw std::runtime_error("cl Platform size 0\n");

    cl_context_properties properties[] =
       { CL_CONTEXT_PLATFORM, (cl_context_properties)(platforms[0])(), 0};
    s_context = cl::Context(CL_DEVICE_TYPE_GPU, properties);

    devices = s_context.getInfo<CL_CONTEXT_DEVICES>();

    std::ifstream sourceFile(s_source_file);
    ensure(sourceFile.is_open(),"unable to open file");

    std::ifstream headerFile(s_header_file);
    ensure(headerFile.is_open(),"unable to open file");

    std::string headerCode(std::istreambuf_iterator<char>(headerFile),(std::istreambuf_iterator<char>()));
    std::string sourceCode(std::istreambuf_iterator<char>(sourceFile),(std::istreambuf_iterator<char>()));

    cl::Program::Sources sources;
    sources.push_back(std::make_pair(headerCode.c_str(),headerCode.size()));
    sources.push_back(std::make_pair(sourceCode.c_str(),sourceCode.size()));

    program = cl::Program(s_context, sources);
    program.build(devices);

    s_queue = cl::CommandQueue(s_context, devices[0]);

    s_assign_max_facet_edge =  cl::Kernel(program, "assign_max_facet_edge").
        bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

    s_assign_max_facet_face =  cl::Kernel(program, "assign_max_facet_face").
        bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

    s_assign_max_facet_cube =  cl::Kernel(program, "assign_max_facet_cube").
        bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

    s_assign_pairs = cl::Kernel(program, "assign_pairs").
        bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

    s_mark_cps = cl::Kernel(program, "mark_cps").
        bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

    s_save_cps = cl::Kernel(program, "save_cps").
        bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

    s_scan_local_sums = cl::Kernel(program, "scan_local_sums").
        bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

    s_scan_group_sums= cl::Kernel(program, "scan_group_sums").
        bind(s_queue,cl::NullRange,cl::NDRange(WI_SIZE),cl::NDRange(WI_SIZE));

    s_scan_update_sums = cl::Kernel(program, "scan_update_sums").
        bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

  }
  catch (cl::Error err)
  {
   std::cerr
      << "ERROR: "
      << err.what()
      << "("
      << err.err()
      << ")"
      << std::endl;

   std::cerr
      <<program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0])
      <<std::endl;

   throw;
  }
}

template<typename T>
void log_buffer(cl::Buffer buf,int n,std::ostream &os=cout)
{
  std::vector<T> buf_cpu(n);

  s_queue.enqueueReadBuffer(buf,true,0,sizeof(T)*n,buf_cpu.data());

  for( int i = 0 ; i < buf_cpu.size(); ++i)
    os << buf_cpu[i]<< " ";

  os<<endl;
}

void assign_gradient_opencl(rect_t rct,rect_t ext,rect_t dom,
                            cell_fn_t *func, cell_flag_t *flag)
{
  cl_cell_t rct_lc = to_cell(rct.lc());
  cl_cell_t rct_uc = to_cell(rct.uc());

  cl_cell_t ext_lc = to_cell(ext.lc());
  cl_cell_t ext_uc = to_cell(ext.uc());

  cl_cell_t dom_lc = to_cell(dom.lc());
  cl_cell_t dom_uc = to_cell(dom.uc());


  cl::size_t<3> func_size = to_size(ext.span()/2+ 1);
  cl::size_t<3> flag_size = to_size(ext.span()  + 1);

  int flag_size_bytes = flag_size[0]*flag_size[1]*flag_size[2]*sizeof(cell_flag_t);

  cl::Buffer  flag_buf;

  try
  {
    cl::Image3D  func_img(s_context,CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,
                          cl::ImageFormat(CL_R,CL_FLOAT),
                          func_size[0],func_size[1],func_size[2],0,0,func);

    cl::Image3D  flag_img(s_context,CL_MEM_READ_ONLY,
                          cl::ImageFormat(CL_R,CL_UNSIGNED_INT8),
                          flag_size[0],flag_size[1],flag_size[2],0,0,flag);

    flag_buf = cl::Buffer(s_context,CL_MEM_READ_WRITE,flag_size_bytes);

    s_assign_max_facet_edge
        (func_img,flag_img,rct_lc,rct_uc,ext_lc,ext_uc,dom_lc,dom_uc,flag_buf);
    s_queue.finish();

    s_queue.enqueueCopyBufferToImage
        (flag_buf,flag_img,0,to_size(0,0,0),flag_size);
    s_queue.finish();

    s_assign_max_facet_face
        (func_img,flag_img,rct_lc,rct_uc,ext_lc,ext_uc,dom_lc,dom_uc,flag_buf);
    s_queue.finish();

    s_queue.enqueueCopyBufferToImage
        (flag_buf,flag_img,0,to_size(0,0,0),flag_size);
    s_queue.finish();

    s_assign_max_facet_cube
        (func_img,flag_img,rct_lc,rct_uc,ext_lc,ext_uc,dom_lc,dom_uc,flag_buf);
    s_queue.finish();

    s_queue.enqueueCopyBufferToImage
        (flag_buf,flag_img,0,to_size(0,0,0),flag_size);
    s_queue.finish();

    s_assign_pairs
        (func_img,flag_img,rct_lc,rct_uc,ext_lc,ext_uc,dom_lc,dom_uc,flag_buf);
    s_queue.finish();
  }
  catch(cl::Error err)
  {
    std::cerr<<_FFL<<std::endl;
    std::cerr<< "ERROR: "<< err.what()<< "("<< err.err()<< ")"<< std::endl;
    throw;
  }

  try
  {
    cl::Buffer local_sums_buf(s_context,CL_MEM_READ_WRITE,sizeof(int)*WG_SIZE);
    cl::Buffer group_sums_buf(s_context,CL_MEM_READ_WRITE,sizeof(int)*(WG_NUM));

    s_mark_cps(rct_lc,rct_uc,ext_lc,ext_uc,dom_lc,dom_uc,flag_buf,local_sums_buf);
    s_scan_local_sums(local_sums_buf,group_sums_buf);
    s_scan_group_sums(group_sums_buf);
    s_scan_update_sums(local_sums_buf,group_sums_buf);
    s_queue.finish();

    int num_cps;

    s_queue.enqueueReadBuffer(group_sums_buf,true,sizeof(int)*(WG_NUM-1),sizeof(int),&num_cps);

    cl::Buffer cp_cellid_buf(s_context,CL_MEM_READ_WRITE,sizeof(cellid_t)*num_cps);
    cl::Buffer cp_pair_idx_buf(s_context,CL_MEM_READ_WRITE,sizeof(int)*num_cps);


    s_save_cps(rct_lc,rct_uc,ext_lc,ext_uc,dom_lc,dom_uc,flag_buf,local_sums_buf,
               cp_cellid_buf,cp_pair_idx_buf);


    s_queue.enqueueReadBuffer(flag_buf,true,0,flag_size_bytes,flag);
    s_queue.finish();

  }
  catch(cl::Error err)
  {
    std::cerr<<_FFL<<std::endl;
    std::cerr<< "ERROR: "<< err.what()<< "("<< err.err()<< ")"<< std::endl;
    throw;
  }

}

#include <grid_dataset.h>

namespace utls
{
  template <> inline std::string to_string (const unsigned char & t)
  {
    char t_str[9];

    t_str[0] = '0'+char((t>>7)&1);
    t_str[1] = '0'+char((t>>6)&1);
    t_str[2] = '0'+char((t>>5)&1);
    t_str[3] = '0'+char((t>>4)&1);
    t_str[4] = '0'+char((t>>3)&1);
    t_str[5] = '0'+char((t>>2)&1);
    t_str[6] = '0'+char((t>>1)&1);
    t_str[7] = '0'+char((t>>0)&1);
    t_str[8] =  0;
    return t_str;
  }
}


void check_assign_gradient_opencl
  (dataset_ptr_t ds,int dim,rect_t check_rect,cell_flag_t mask)
{
  rect_size_t   span = ds->m_ext_rect.span() + 1;
  rect_point_t   bl  = ds->m_ext_rect.lower_corner();

  dataset_t::cellflag_array_t flag(span,boost::fortran_storage_order());

  flag.reindex(bl);

  assign_gradient_opencl(ds->m_rect,ds->m_ext_rect,ds->m_domain_rect,
                          ds->m_vert_fns.data(),flag.data());

  for(int d = 0 ; d <= dim; ++d)
  {
    cellid_t c,s(0,0,0),stride(2,2,2);

    for(int i = 0 ;  i < d; ++i)
      s[i] = 1;

    while(true)
    {
      rect_t rect  = rect_t(check_rect.lc()+s,check_rect.uc()-s);

      int n = c_to_i(rect.uc(),rect,stride) + 1;

      for( int i = 0; i < n; i ++)
      {
        c = i_to_c(i,rect,stride);

        if((flag(c)&mask) != (ds->m_cell_flags(c)&mask))
        {
          cell_flag_t f_gpu =flag(c)&mask;
          cell_flag_t f_cpu =(ds->m_cell_flags(c)&mask) &mask;

          cout<<SVAR(c);
          if(ds->isCellPaired(c))
          {
            cout<<"     "<<SVAR(ds->getCellPairId(c));
          }
          cout<<endl;

          cout<<hex<<SVAR(f_cpu)<<" "<<SVAR(f_gpu);
          cout<<endl;

        }
      }

      if(!next_permutation(s.rbegin(),s.rend()))
        break;
    }
  }
}


