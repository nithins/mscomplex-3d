#include <iostream>

#include <boost/python.hpp>
#include <boost/python/object.hpp>
#include <boost/python/stl_iterator.hpp>


#include <grid.h>
#include <grid_datamanager.h>

namespace bp = boost::python;

namespace grid
{
  cellid_t tup_to_cellid(bp::object pobj)
  {
    bp::tuple dim_tup = bp::extract<bp::tuple>(pobj);

    cellid_t c;

    for(int i = 0 ; i < c.size() ; ++i)
      c[i] = bp::extract<cell_coord_t>(dim_tup[i]);

    return c;
  }

  boost::shared_ptr<data_manager_t> make_data_manager
      (std::string f,bp::object dim,int ml,double t)
  {
    return boost::shared_ptr<data_manager_t>
        (new data_manager_t(f,tup_to_cellid(dim),ml,t));
  }
}

using namespace grid;

BOOST_PYTHON_MODULE(pymscomplex3d)
{

  bp::class_<data_manager_t, boost::shared_ptr<data_manager_t> >
      ("data_manager_t",bp::no_init)
      .def("__init__", bp::make_constructor(make_data_manager))
      .def("work", &data_manager_t::work);
}
