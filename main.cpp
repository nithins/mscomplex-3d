#include <exception>
#include <string>

#include <grid_viewer_mainwindow.h>
#include <grid_datamanager.h>

#include <cpputils.h>
#include <iostream>

#include <boost/program_options.hpp>

#include <stdexcept>

using namespace std;

namespace bpo = boost::program_options ;

int main(int ac , char **av)
{
  string filename;

  grid::cellid_t size;

  bool   use_ocl = false;

  double   simp_tresh= 0.0;

  bool   gui = false;

  grid::rect_t roi ;

  bpo::options_description desc("Allowed options");
  desc.add_options()
      ("help,h", "produce help message")
      ("file,f",bpo::value<std::string >(), "grid file name")
      ("size,d", bpo::value<grid::cellid_t >(), "size of grid entered as [x,y,z]")
      ("cl","use OpenCL ")
      ("simp-tresh,t",bpo::value<double>(),"simplification treshold")
      ("gui,g","show gui")
      ("roi,r", bpo::value<grid::rect_t >(),
       "visualize a regoin of interset \nentered as [[x0,x1],[y0,y1],[z0,z1]]")
      ;


  bpo::variables_map vm;
  bpo::store(bpo::parse_command_line(ac, av, desc), vm);
  bpo::notify(vm);

  if (vm.count("help"))
  {
    std::cout << desc << "\n";
    return 1;
  }

  if (vm.count("size"))
    size = vm["size"].as<grid::cellid_t >();
  else
    throw std::invalid_argument("no dim specified");

  if (vm.count("file"))
    filename = vm["file"].as<std::string>();
  else
    throw std::invalid_argument("no filename specified");

  if (vm.count("cl"))
    use_ocl = true;

  if (vm.count("simp-tresh"))
    simp_tresh = vm["simp-tresh"].as<double>();

  if (vm.count("gui"))
    gui = true;

  if (vm.count("roi"))
    roi = vm["roi"].as<grid::rect_t >();

  grid::data_manager_t * gdm = new grid::data_manager_t
                          (filename,size,
                           use_ocl,
                           simp_tresh);

  gdm->work();

  if(gui)
  {
    QApplication application(ac,av);

    grid::viewer_mainwindow gvmw(gdm,roi);

    gvmw.setWindowTitle("ms complex vis");

    gvmw.show();

    application.exec();
  }
  else
  {
    delete gdm;
  }
}
