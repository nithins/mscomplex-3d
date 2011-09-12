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

const char * s_header_file = "/home/nithin/projects/mscomplex-3d/dmsc/OpenCL/grid_dataset.clh";
const char * s_source_file = "/home/nithin/projects/mscomplex-3d/dmsc/OpenCL/grid_dataset.cl";

const int WI_SIZE = 256;
const int WG_SIZE = 64*WI_SIZE;

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

   std::cerr<<program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0])<<std::endl;
   throw;
  }
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

  try
  {
    cl::Image3D func_img(s_context,CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,
                         cl::ImageFormat(CL_R,CL_FLOAT),
                         func_size[0],func_size[1],func_size[2],0,0,func);

    cl::Image3D flag_img(s_context,CL_MEM_READ_ONLY,
                         cl::ImageFormat(CL_R,CL_UNSIGNED_INT8),
                         flag_size[0],flag_size[1],flag_size[2],0,0,flag);

    cl::Buffer  flag_buf(s_context,CL_MEM_READ_WRITE,flag_size_bytes);

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

    s_queue.enqueueReadBuffer(flag_buf,true,0,flag_size_bytes,flag);
    s_queue.finish();
  }
  catch(cl::Error err)
  {

    std::cerr
       << "ERROR: "
       << err.what()
       << "("
       << err.err()
       << ")"
       << std::endl;

    throw;
  }
}
