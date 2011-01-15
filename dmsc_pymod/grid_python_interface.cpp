#include <boost/python.hpp>
#include <boost/python/object.hpp>
#include <boost/python/stl_iterator.hpp>


#include <grid.h>
#include <grid_datamanager.h>

using namespace boost::python;
using namespace grid;

void cellid_assign(cellid_t& c, object o)
{
  stl_input_iterator<cell_coord_t> b(o), e;

  cellid_t::iterator o_b = c.begin(),o_e = c.end();

  while ( b != e && o_b != o_e)
    *o_b++ = *b++;
}

BOOST_PYTHON_MODULE(pymscomplex3d)
{
  class_<cellid_t>("cellid_t", init<>())
      .def(init<cell_coord_t,cell_coord_t,cell_coord_t>())
      .def("__iter__", iterator<cellid_t >())
      .def("assign", &cellid_assign)
      .def(self + other<cellid_t>())
      .def(self - other<cellid_t>())
      .def(self * other<cellid_t>())
      .def(self / other<cellid_t>())
      .def(self < other<cellid_t>())
      .def(self + other<cell_coord_t>())
      .def(self - other<cell_coord_t>())
      .def(self * other<cell_coord_t>())
      .def(self / other<cell_coord_t>())
      .def("__str__",&boost::lexical_cast<std::string,cellid_t>)
      ;

  class_<data_manager_t>("data_manager_t",init<std::string,cellid_t,bool,double>())
      .def("work", &data_manager_t::work)
      ;
}
