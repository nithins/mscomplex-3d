#include <exception>
#include <string>

#include <config.h>

#ifdef BUILD_EXEC_GUI
#include <grid_viewer_mainwindow.h>
#endif

#include <grid_datamanager.h>
#include <cpputils.h>

#include <boost/program_options.hpp>

#include <stdexcept>
#include <iostream>

using namespace std;

namespace bpo = boost::program_options ;

int main(int ac , char **av)
{
  string filename;

  grid::cellid_t dim;

  bool   use_ocl = false;

  double   simp_tresh= 0.0;

#ifdef BUILD_EXEC_GUI
  bool   gui = false;
#endif

  bpo::options_description desc("Allowed options");
  desc.add_options()
      ("help,h", "produce help message")
      ("file,f",bpo::value<string >(), "grid file name")
      ("dim,d", bpo::value<grid::cellid_t >(), "dim of grid entered as (x,y,z)")
      ("cl","use OpenCL ")
      ("simp-tresh,t",bpo::value<double>(),"simplification treshold")
#ifdef BUILD_EXEC_GUI
      ("gui,g","show gui")
#endif
      ;

  bpo::variables_map vm;
  bpo::store(bpo::parse_command_line(ac, av, desc), vm);
  bpo::notify(vm);

  if (vm.count("help"))
  {
    cout << desc << "\n";
    return 1;
  }

  if (vm.count("dim"))
    dim = vm["dim"].as<grid::cellid_t >();
  else
    throw invalid_argument("no dim specified");

  if (vm.count("file"))
    filename = vm["file"].as<string>();
  else
    throw invalid_argument("no filename specified");

  if (vm.count("cl"))
    use_ocl = true;

  if (vm.count("simp-tresh"))
    simp_tresh = vm["simp-tresh"].as<double>();

#ifdef BUILD_EXEC_GUI
  if (vm.count("gui"))
    gui = true;
#endif

  grid::data_manager_t * gdm = new grid::data_manager_t
      (filename,dim,use_ocl,simp_tresh);

  gdm->work();

#ifdef BUILD_EXEC_GUI
  if(gui)
  {
    QApplication application(ac,av);

    grid::viewer_mainwindow gvmw(gdm);

    gvmw.setWindowTitle("ms complex vis");

    gvmw.show();

    application.exec();
  }
  else
  {
    delete gdm;
  }
#else
  delete gdm;
#endif
}
