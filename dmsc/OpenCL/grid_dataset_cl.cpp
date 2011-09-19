#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include <stdexcept>
#include <cpputils.h>

#include<grid_dataset_cl.h>

#include<grid_dataset.h>
#include<grid_mscomplex.h>

const int WI_SIZE = 256;
const int WG_NUM  = 32;
const int WG_SIZE = WG_NUM*WI_SIZE;

cl::Context       s_context;
cl::CommandQueue  s_queue;

cl::KernelFunctor s_assign_max_facet_edge;
cl::KernelFunctor s_assign_max_facet_face;
cl::KernelFunctor s_assign_max_facet_cube;
cl::KernelFunctor s_assign_pairs;

cl::KernelFunctor s_mark_cps;
cl::KernelFunctor s_mark_boundry_cps;
cl::KernelFunctor s_count_cps;
cl::KernelFunctor s_count_boundry_cps;
cl::KernelFunctor s_save_boundry_cps;
cl::KernelFunctor s_save_cps;

cl::KernelFunctor s_scan_local_sums;
cl::KernelFunctor s_scan_group_sums;
cl::KernelFunctor s_scan_update_sums;

cl::KernelFunctor s_init_propagate;
cl::KernelFunctor s_propagate;
cl::KernelFunctor s_init_update_to_cp_no;
cl::KernelFunctor s_update_to_cp_no;
cl::KernelFunctor s_update_to_surviving_cp_no;


const char * s_header_file =
    "/home/nithin/projects/mscomplex-3d/dmsc/OpenCL/grid_dataset.clh";
const char * s_source1_file =
    "/home/nithin/projects/mscomplex-3d/dmsc/OpenCL/grid_dataset_assigngradient.cl";
const char * s_source2_file =
    "/home/nithin/projects/mscomplex-3d/dmsc/OpenCL/grid_dataset_markandcollect.cl";
const char * s_source3_file =
    "/home/nithin/projects/mscomplex-3d/dmsc/OpenCL/grid_dataset_ownerextrema.cl";


using namespace std;

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

namespace grid
{
  namespace opencl
  {
    typedef cl_short4 cell_t;
    typedef cl_short8   cell_pair_t;

    inline cell_t to_cell(const cellid_t & b)
    {
      cell_t a;

      a.x = b[0];
      a.y = b[1];
      a.z = b[2];

      return a;
    }

    inline cellid_t from_cell(const cell_t & b)
    {
      cellid_t a;

      a[0] = b.x;
      a[1] = b.y;
      a[2] = b.z;

      return a;
    }


    inline cell_pair_t to_cell_pair(const rect_t & b)
    {
      cell_pair_t a;

      a.lo = to_cell(b.lc());
      a.hi = to_cell(b.uc());

      return a;
    }

    inline rect_t from_cell_pair(const cell_pair_t & r)
    {
      return rect_t(from_cell(r.lo),from_cell(r.hi));
    }

    inline cell_t to_cell(int x , int y , int z)
    {
      cell_t a;

      a.x = x;
      a.y = y;
      a.z = z;

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

    inline cl::size_t<3> get_size(const cell_pair_t &p)
    {
      cl::size_t<3> s;

      s[0] = p.hi.x -p.lo.x +1;
      s[1] = p.hi.y -p.lo.y +1;
      s[2] = p.hi.z -p.lo.z +1;

      return s;
    }

    inline int num_cells(const cell_pair_t &p)
    {
      return grid::num_cells(from_cell_pair(p));
    }

    template<typename T>
    void log_buffer(cl::Buffer buf,int n,int nlrepeat = -1,std::ostream &os=cout)
    {
      std::vector<T> buf_cpu(n);

      s_queue.finish();

      s_queue.enqueueReadBuffer(buf,true,0,sizeof(T)*n,buf_cpu.data());

      for( int i = 0 ; i < buf_cpu.size(); ++i)
      {
        if( (nlrepeat > 0)  && (i%nlrepeat == 0))
          cout<<endl;

        os << buf_cpu[i]<< " ";
      }

      os<<endl;
    }

    void init(void)
    {
      std::vector<cl::Device> devices;
      cl::Program             program1;
      cl::Program             program2;
      cl::Program             program3;

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

        s_queue = cl::CommandQueue(s_context, devices[0]);

      }
      catch (cl::Error err)
      {
       cerr<< "SETUP QUEUE ERROR: "<< err.what()<< "("<< err.err()<< ")"<< endl;
       throw;
      }

      try
      {
        std::ifstream sourceFile(s_source1_file);
        ensure(sourceFile.is_open(),"unable to open file");

        std::ifstream headerFile(s_header_file);
        ensure(headerFile.is_open(),"unable to open file");

        std::string headerCode(std::istreambuf_iterator<char>(headerFile),(std::istreambuf_iterator<char>()));
        std::string sourceCode(std::istreambuf_iterator<char>(sourceFile),(std::istreambuf_iterator<char>()));

        cl::Program::Sources sources;
        sources.push_back(std::make_pair(headerCode.c_str(),headerCode.size()));
        sources.push_back(std::make_pair(sourceCode.c_str(),sourceCode.size()));

        program1 = cl::Program(s_context, sources);
        program1.build(devices);

        s_assign_max_facet_edge =  cl::Kernel(program1, "assign_max_facet_edge").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_assign_max_facet_face =  cl::Kernel(program1, "assign_max_facet_face").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_assign_max_facet_cube =  cl::Kernel(program1, "assign_max_facet_cube").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_assign_pairs = cl::Kernel(program1, "assign_pairs").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));
      }
      catch (cl::Error err)
      {
       cerr<< "PROGRAM1 ERROR: "<< err.what()<< "("<< err.err()<< ")"<< endl;
       cerr<<program1.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0])<<endl;

       throw;
      }

      try
      {
        std::ifstream sourceFile(s_source2_file);
        ensure(sourceFile.is_open(),"unable to open file");

        std::ifstream headerFile(s_header_file);
        ensure(headerFile.is_open(),"unable to open file");

        std::string headerCode(std::istreambuf_iterator<char>(headerFile),(std::istreambuf_iterator<char>()));
        std::string sourceCode(std::istreambuf_iterator<char>(sourceFile),(std::istreambuf_iterator<char>()));

        cl::Program::Sources sources;
        sources.push_back(std::make_pair(headerCode.c_str(),headerCode.size()));
        sources.push_back(std::make_pair(sourceCode.c_str(),sourceCode.size()));

        program2 = cl::Program(s_context, sources);
        program2.build(devices);

        s_mark_cps = cl::Kernel(program2, "mark_cps").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_count_cps = cl::Kernel(program2, "count_cps").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_save_cps = cl::Kernel(program2, "save_cps").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_mark_boundry_cps = cl::Kernel(program2, "mark_boundry_cps").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_count_boundry_cps = cl::Kernel(program2, "count_boundry_cps").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_save_boundry_cps = cl::Kernel(program2, "save_boundry_cps").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_scan_local_sums = cl::Kernel(program2, "scan_local_sums").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_scan_group_sums= cl::Kernel(program2, "scan_group_sums").
            bind(s_queue,cl::NullRange,cl::NDRange(WI_SIZE),cl::NDRange(WI_SIZE));

        s_scan_update_sums = cl::Kernel(program2, "scan_update_sums").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

      }
      catch (cl::Error err)
      {
       cerr<< "PROGRAM2 ERROR: "<< err.what()<< "("<< err.err()<< ")"<< endl;
       cerr<<program2.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0])<<endl;

       throw;
      }

      try
      {
        std::ifstream sourceFile(s_source3_file);
        ensure(sourceFile.is_open(),"unable to open file");

        std::ifstream headerFile(s_header_file);
        ensure(headerFile.is_open(),"unable to open file");

        std::string headerCode(std::istreambuf_iterator<char>(headerFile),(std::istreambuf_iterator<char>()));
        std::string sourceCode(std::istreambuf_iterator<char>(sourceFile),(std::istreambuf_iterator<char>()));

        cl::Program::Sources sources;
        sources.push_back(std::make_pair(headerCode.c_str(),headerCode.size()));
        sources.push_back(std::make_pair(sourceCode.c_str(),sourceCode.size()));

        program3 = cl::Program(s_context, sources);
        program3.build(devices);

        s_init_propagate =  cl::Kernel(program3, "init_propagate").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_propagate =  cl::Kernel(program3, "propagate").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_init_update_to_cp_no =  cl::Kernel(program3, "init_update_to_cp_no").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_update_to_cp_no = cl::Kernel(program3, "update_to_cp_no").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

//        s_update_to_surviving_cp_no =  cl::Kernel(program3, "update_to_surviving_cp_no").
//            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

      }
      catch (cl::Error err)
      {
       cerr<< "PROGRAM3 ERROR: "<< err.what()<< "("<< err.err()<< ")"<< endl;
       cerr<<program3.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0])<<endl;

       throw;
      }
    }

    void __assign_gradient
      ( cell_pair_t rct,
        cell_pair_t ext,
        cell_pair_t dom,
        cl::Image3D &func_img,
        cl::Image3D &flag_img,
        cell_fn_t   *h_func,
        cell_flag_t *h_flag)
    {
      cl::size_t<3> func_size = to_size(from_cell_pair(ext).span()/2+ 1);
      cl::size_t<3> flag_size = get_size(ext);

      int cell_ct = num_cells(ext);

      try
      {
        func_img = cl::Image3D(s_context,CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,
                               cl::ImageFormat(CL_R,CL_FLOAT),
                               func_size[0],func_size[1],func_size[2],0,0,h_func);

        flag_img = cl::Image3D(s_context,CL_MEM_READ_ONLY,
                               cl::ImageFormat(CL_R,CL_UNSIGNED_INT8),
                               flag_size[0],flag_size[1],flag_size[2],0,0);

        cl::Buffer flag_buf(s_context,CL_MEM_READ_WRITE,cell_ct*sizeof(cell_flag_t));

        s_assign_max_facet_edge
            (func_img,flag_img,rct.lo,rct.hi,ext.lo,ext.hi,dom.lo,dom.hi,flag_buf);
        s_queue.finish();

        s_queue.enqueueCopyBufferToImage
            (flag_buf,flag_img,0,to_size(0,0,0),flag_size);
        s_queue.finish();

        s_assign_max_facet_face
            (func_img,flag_img,rct.lo,rct.hi,ext.lo,ext.hi,dom.lo,dom.hi,flag_buf);
        s_queue.finish();

        s_queue.enqueueCopyBufferToImage
            (flag_buf,flag_img,0,to_size(0,0,0),flag_size);
        s_queue.finish();

        s_assign_max_facet_cube
            (func_img,flag_img,rct.lo,rct.hi,ext.lo,ext.hi,dom.lo,dom.hi,flag_buf);
        s_queue.finish();

        s_queue.enqueueCopyBufferToImage
            (flag_buf,flag_img,0,to_size(0,0,0),flag_size);
        s_queue.finish();

        s_assign_pairs
            (func_img,flag_img,rct.lo,rct.hi,ext.lo,ext.hi,dom.lo,dom.hi,flag_buf);
        s_queue.finish();

        s_mark_cps(rct,ext,dom,flag_buf);

        rect_list_t bnds;

        get_boundry_rects(from_cell_pair(rct),from_cell_pair(ext),bnds);

        for( int i = 0 ; i < bnds.size(); ++i)
        {
          cell_pair_t bnd = to_cell_pair(bnds[i]);
          cell_t  bnd_dir = to_cell(bnds[i].get_normal());

          s_mark_boundry_cps(rct,ext,dom,bnd,bnd_dir,flag_buf);
        }

        s_queue.enqueueCopyBufferToImage(flag_buf,flag_img,0,to_size(0,0,0),flag_size);
        s_queue.enqueueReadBuffer(flag_buf,false,0,cell_ct*sizeof(cell_flag_t),h_flag);
      }
      catch(cl::Error err)
      {
        std::cerr<<_FFL<<std::endl;
        std::cerr<< "ERROR: "<< err.what()<< "("<< err.err()<< ")"<< std::endl;
        throw;
      }
    }

    void __count_and_scan_cps
    ( cell_pair_t rct,
      cell_pair_t ext,
      cell_pair_t dom,
      cl::Image3D &flag_img,
      cl::Buffer  &cp_count_buf,
      int &num_cps)
    {
      try
      {
        cp_count_buf = cl::Buffer(s_context,CL_MEM_READ_WRITE,sizeof(int)*WG_SIZE);
        cl::Buffer group_sums_buf(s_context,CL_MEM_READ_WRITE,sizeof(int)*(WG_NUM));

        s_count_cps(rct,ext,dom,flag_img,cp_count_buf);

        rect_list_t bnds;

        get_boundry_rects(from_cell_pair(rct),from_cell_pair(ext),bnds);

        for( int i = 0 ; i < bnds.size(); ++i)
        {
          cell_pair_t bnd = to_cell_pair(bnds[i]);
          cell_t  bnd_dir = to_cell(bnds[i].get_normal());

          s_count_boundry_cps(rct,ext,dom,bnd,bnd_dir,flag_img,cp_count_buf);
        }

        s_scan_local_sums(cp_count_buf,group_sums_buf);
        s_scan_group_sums(group_sums_buf);
        s_scan_update_sums(cp_count_buf,group_sums_buf);
        s_queue.finish();
        s_queue.enqueueReadBuffer(group_sums_buf,false,sizeof(int)*(WG_NUM-1),sizeof(int),&num_cps);
        s_queue.finish();

      }
      catch(cl::Error err)
      {
        std::cerr<<_FFL<<std::endl;
        std::cerr<< "ERROR: "<< err.what()<< "("<< err.err()<< ")"<< std::endl;
        throw;
      }
    }

    void __save_cps
      ( cell_pair_t rct,
        cell_pair_t ext,
        cell_pair_t dom,
        cl::Image3D &func_img,
        cl::Image3D &flag_img,
        cl::Buffer  &cp_offset_buf,
        int num_cps,
        cellid_t   *h_cellid,
        cellid_t   *h_vertid,
        int        *h_pair_idx,
        char       *h_index,
        cell_fn_t  *h_func)
    {
      try
      {
        cl::Buffer cp_cellid_buf(s_context,CL_MEM_READ_WRITE,sizeof(cellid_t)*num_cps);
        cl::Buffer cp_vertid_buf(s_context,CL_MEM_READ_WRITE,sizeof(cellid_t)*num_cps);
        cl::Buffer cp_pair_idx_buf(s_context,CL_MEM_READ_WRITE,sizeof(int)*num_cps);
        cl::Buffer cp_index_buf(s_context,CL_MEM_READ_WRITE,sizeof(char)*num_cps);
        cl::Buffer cp_func_buf(s_context,CL_MEM_READ_WRITE,sizeof(cell_fn_t)*num_cps);

        s_save_cps(rct,ext,dom,func_img,flag_img,cp_offset_buf,cp_cellid_buf,
                   cp_index_buf,cp_pair_idx_buf,cp_vertid_buf,cp_func_buf);

        rect_list_t bnds;

        get_boundry_rects(from_cell_pair(rct),from_cell_pair(ext),bnds);

        for( int i = 0 ; i < bnds.size(); ++i)
        {
          cell_pair_t bnd = to_cell_pair(bnds[i]);
          cell_t  bnd_dir = to_cell(bnds[i].get_normal());


          s_save_boundry_cps(rct,ext,dom,bnd,bnd_dir,func_img,flag_img,
                             cp_offset_buf,cp_cellid_buf,cp_index_buf,
                             cp_pair_idx_buf,cp_vertid_buf,cp_func_buf);
        }

        s_queue.finish();

        s_queue.enqueueReadBuffer(cp_cellid_buf,false,0,sizeof(cellid_t)*num_cps,h_cellid);
        s_queue.enqueueReadBuffer(cp_vertid_buf,false,0,sizeof(cellid_t)*num_cps,h_vertid);
        s_queue.enqueueReadBuffer(cp_pair_idx_buf,false,0,sizeof(int)*num_cps,h_pair_idx);
        s_queue.enqueueReadBuffer(cp_index_buf,false,0,sizeof(char)*num_cps,h_index);
        s_queue.enqueueReadBuffer(cp_func_buf,false,0,sizeof(cell_fn_t)*num_cps,h_func);

        s_queue.finish();
      }
      catch(cl::Error err)
      {
        std::cerr<<_FFL<<std::endl;
        std::cerr<< "ERROR: "<< err.what()<< "("<< err.err()<< ")"<< std::endl;
        throw;
      }
    }

    void worker::assign_gradient(dataset_ptr_t ds)
    {
      cl::Image3D  func_img;
      cell_pair_t rct = to_cell_pair(ds->m_rect);
      cell_pair_t ext = to_cell_pair(ds->m_ext_rect);
      cell_pair_t dom = to_cell_pair(ds->m_domain_rect);

      __assign_gradient(rct,ext,dom,func_img,flag_img,
                        ds->m_vert_fns.data(),ds->m_cell_flags.data());
    }


    void worker::assign_gradient(dataset_ptr_t ds, mscomplex_ptr_t msc)
    {
      cl::Image3D  func_img;
      cell_pair_t rct = to_cell_pair(ds->m_rect);
      cell_pair_t ext = to_cell_pair(ds->m_ext_rect);
      cell_pair_t dom = to_cell_pair(ds->m_domain_rect);

      __assign_gradient(rct,ext,dom,func_img,flag_img,
                        ds->m_vert_fns.data(),ds->m_cell_flags.data());


      cl::Buffer cp_offset_buf;
      int          num_cps;

      __count_and_scan_cps(rct,ext,dom,flag_img,cp_offset_buf,num_cps);

      msc->resize(num_cps);

      __save_cps(rct,ext,dom,func_img,flag_img,cp_offset_buf,num_cps,
                 msc->m_cp_cellid.data(),msc->m_cp_vertid.data(),
                 msc->m_cp_pair_idx.data(),msc->m_cp_index.data(),
                 msc->m_cp_fn.data());
    }

    void __owner_extrema
    ( cell_pair_t rct,
      cell_pair_t ext,
      cell_pair_t dom,
      cell_pair_t ex_rect,
      cl::Image3D &flag_img,
      cellid_t *h_cp_cellid,
      cl_int    num_cps,
      int      *h_result)
    {
      int num_cells = num_cells2(from_cell_pair(ex_rect));

      try
      {
        cl::Buffer own_buf1(s_context,CL_MEM_READ_WRITE,num_cells*sizeof(int));
        cl::Buffer own_buf2(s_context,CL_MEM_READ_WRITE,num_cells*sizeof(int));
        cl::Buffer is_updated_buf(s_context,CL_MEM_READ_WRITE,sizeof(int));
        cl::Buffer cp_cellid_buf(s_context,CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,
                                 num_cps*sizeof(cellid_t),h_cp_cellid);

        s_init_propagate(rct,ext,dom,ex_rect,flag_img,own_buf1);

        int is_updated;

        do
        {
          is_updated = 0;

          s_queue.enqueueWriteBuffer(is_updated_buf,true,0,sizeof(int),&is_updated);
          s_propagate(rct,ext,dom,ex_rect,own_buf1,own_buf2,is_updated_buf);
          s_queue.finish();

          s_queue.enqueueReadBuffer(is_updated_buf,true,0,sizeof(int),&is_updated);
          s_queue.finish();

          std::swap(own_buf1,own_buf2);
        }
        while(is_updated == 1);

        s_init_update_to_cp_no(rct,ext,dom,ex_rect,cp_cellid_buf,num_cps,own_buf1);
        s_update_to_cp_no(rct,ext,dom,ex_rect,flag_img,own_buf1,own_buf2);
        s_queue.finish();

        s_queue.enqueueReadBuffer(own_buf2,false,0,num_cells*sizeof(int),h_result);
        s_queue.finish();

      }
      catch(cl::Error err)
      {
        std::cerr<<_FFL<<std::endl;
        std::cerr<< "ERROR: "<< err.what()<< "("<< err.err()<< ")"<< std::endl;
        throw;
      }
    }

    void worker::owner_extrema(dataset_ptr_t ds,mscomplex_ptr_t msc)
    {
      cell_pair_t rct = to_cell_pair(ds->m_rect);
      cell_pair_t ext = to_cell_pair(ds->m_ext_rect);
      cell_pair_t dom = to_cell_pair(ds->m_domain_rect);

      cell_pair_t max_rect = to_cell_pair(rect_t(ds->m_rect.lc()+1,ds->m_rect.uc()-1));
      cell_pair_t min_rect = to_cell_pair(ds->m_rect);

      __owner_extrema(rct,ext,dom,max_rect,flag_img,
                      msc->m_cp_cellid.data(),msc->get_num_critpts(),
                      ds->m_owner_maxima.data());

      __owner_extrema(rct,ext,dom,min_rect,flag_img,
                      msc->m_cp_cellid.data(),msc->get_num_critpts(),
                      ds->m_owner_minima.data());

      s_queue.finish();
    }
  }
}

//    void check_assign_gradient_opencl
//      (dataset_ptr_t ds,int dim,rect_t check_rect,cell_flag_t mask)
//    {
//      rect_size_t   span = ds->m_ext_rect.span() + 1;
//      rect_point_t   bl  = ds->m_ext_rect.lower_corner();

//      dataset_t::cellflag_array_t flag(span,boost::fortran_storage_order());

//      flag.reindex(bl);

////      __assign_gradient(ds->m_rect,ds->m_ext_rect,ds->m_domain_rect,
////                              ds->m_vert_fns.data(),flag.data());

//      for(int d = 0 ; d <= dim; ++d)
//      {
//        cellid_t c,s(0,0,0),stride(2,2,2);

//        for(int i = 0 ;  i < d; ++i)
//          s[i] = 1;

//        while(true)
//        {
//          rect_t rect  = rect_t(check_rect.lc()+s,check_rect.uc()-s);

//          int n = c_to_i(rect.uc(),rect,stride) + 1;

//          for( int i = 0; i < n; i ++)
//          {
//            c = i_to_c(i,rect,stride);

//            if((flag(c)&mask) != (ds->m_cell_flags(c)&mask))
//            {
//              cell_flag_t f_gpu =flag(c)&mask;
//              cell_flag_t f_cpu =(ds->m_cell_flags(c)&mask) &mask;

//              cout<<SVAR(c);
//              if(ds->isCellPaired(c))
//              {
//                cout<<"     "<<SVAR(ds->getCellPairId(c));
//              }
//              cout<<endl;

//              cout<<hex<<SVAR(f_cpu)<<" "<<SVAR(f_gpu);
//              cout<<endl;

//            }
//          }

//          if(!next_permutation(s.rbegin(),s.rend()))
//            break;
//        }
//      }
//    }




