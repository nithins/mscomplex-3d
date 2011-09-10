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

cl_cell_t to_cell(const cellid_t & b)
{
  cl_cell_t a;

  a.x = b[0];
  a.y = b[1];
  a.z = b[2];

  return a;
}


cl::Context      s_context;
cl::CommandQueue s_queue;

cl::KernelFunctor s_assign_max_facet_edge;
cl::KernelFunctor s_assign_zero;

const char * s_header_file = "/home/nithin/projects/mscomplex-3d/dmsc/OpenCL/grid_dataset.clh";
const char * s_source_file = "/home/nithin/projects/mscomplex-3d/dmsc/OpenCL/grid_dataset.cl";

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
        bind(s_queue,cl::NullRange,cl::NDRange(128),cl::NDRange(32));

    s_assign_zero          =  cl::Kernel(program, "assign_zero").
        bind(s_queue,cl::NullRange,cl::NDRange(128),cl::NDRange(32));


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

void assign_max_facet_opencl(rect_t ext_rect,cell_fn_t *func, cell_flag_t *flag)
{
  cellid_t func_size = ext_rect.span()/2+ 1;
  cellid_t flag_size = ext_rect.span()  + 1;

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

    s_assign_zero(to_cell(ext_rect.lc()),to_cell(ext_rect.uc()),flag_buf);

    s_queue.finish();

    s_assign_max_facet_edge
        (func_img,flag_img,to_cell(ext_rect.lc()),to_cell(ext_rect.uc()),
         to_cell(cellid_t(1,0,0)),flag_buf);

    s_assign_max_facet_edge
        (func_img,flag_img,to_cell(ext_rect.lc()),to_cell(ext_rect.uc()),
         to_cell(cellid_t(0,1,0)),flag_buf);

    s_assign_max_facet_edge
        (func_img,flag_img,to_cell(ext_rect.lc()),to_cell(ext_rect.uc()),
         to_cell(cellid_t(0,0,1)),flag_buf);

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






