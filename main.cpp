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
  string         filename;
  grid::cellid_t size;
  bool           use_ocl;
  double         simp_tresh;
  int            levels;

#ifdef BUILD_EXEC_GUI
  bool           gui;
#endif

  bpo::options_description desc("Allowed options");
  desc.add_options()
      ("help,h", "produce help message")
      ("file,f",bpo::value<string >(&filename)->required(), "grid file name")
      ("dim,d", bpo::value<grid::cellid_t>(&size)->required(), "dim of grid entered as (x,y,z)")
      ("levels,l",bpo::value<int>(&levels)->default_value(0),"number of subdivision levels")
      ("simp-tresh,t",bpo::value<double>(&simp_tresh)->default_value(0.0),"simplification treshold")
#ifdef BUILD_EXEC_GUI
      ("gui,g",bpo::value<bool>(&gui)->default_value(false),"show gui")
#endif
      ;

  bpo::variables_map vm;
  bpo::store(bpo::parse_command_line(ac, av, desc), vm);

  if (vm.count("help"))
  {
    cout << desc << endl;
    return 0;
  }
  try
  {
    bpo::notify(vm);
  }
  catch(bpo::required_option e)
  {
    cout<<e.what()<<endl;
    cout<<desc<<endl;
    return 1;
  }

  grid::data_manager_t * gdm = new grid::data_manager_t
      (filename,size,levels,simp_tresh);

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
