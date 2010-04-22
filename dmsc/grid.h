#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

#include <cpputils.h>

#include <boost/static_assert.hpp>
#define static_assert BOOST_STATIC_ASSERT

namespace grid
{
  template <typename coord_type,
  coord_type invalid_value = -1> class   rectangle_complex{

  public:

    struct point_def:public three_tuple_t<coord_type>
    {
      typedef three_tuple_t<coord_type> base_t;

      point_def ( const coord_type &x,
                  const coord_type &y,
                  const coord_type &z) :
      three_tuple_t<coord_type>(x,y,z){}

      point_def () :
          three_tuple_t<coord_type>(invalid_value,invalid_value,invalid_value){}

      inline point_def operator/(const coord_type s) const
      {
        point_def ret;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          ret[i] = (*this)[i]/s;

        return ret;
      }

      inline point_def operator-(const point_def & o) const
      {
        point_def ret;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          ret[i] = (*this)[i]-o[i];

        return ret;
      }

      inline point_def operator+(const point_def & o) const
      {
        point_def ret;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          ret[i] = (*this)[i]-o[i];

        return ret;
      }

      inline coord_type operator*(const point_def & o) const
      {
        coord_type ret = 0;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          ret += (*this)[i] * o[i];

        return ret;

      }

      inline point_def operator*(const coord_type & o) const
      {
        point_def ret;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          ret [i]= (*this)[i] * o;

        return ret;

      }

      static point_def zero;

      static point_def one;
    };

    struct range_def:public two_tuple_t<coord_type>
    {
      range_def ( const coord_type &l,const coord_type &u)
      {

        (*this)[0] = std::min ( l,u);
        (*this)[1] = std::max ( l,u);
      }

      range_def ():two_tuple_t<coord_type>(invalid_value,invalid_value){}

      inline bool isInOpen(const coord_type &c) const
      {
        return (( (*this)[0] < c ) && (  c < (*this)[1] ));
      }

      inline bool isInClosed(const coord_type &c) const
      {
        return (( (*this)[0] <= c ) && (  c <= (*this)[1] ));
      }

      inline bool isOnBndry(const coord_type &c) const
      {
        return (( (*this)[0] == c ) || (  c == (*this)[1] ));
      }

      inline bool contains(const range_def & r) const
      {
        return isInOpen(r[0]) && isInOpen(r[1]);
      }

      inline bool intersects(const range_def & r) const
      {
        return !((r[0] > (*this)[1]) || ((*this)[0] > r[1]));
      }

      inline bool intersection(const range_def & r,range_def & i) const
      {
        i = range_def(std::max(r[0],(*this)[0]),std::min(r[1],(*this)[1]));

        return intersects(r);
      }

      inline range_def range_union(const range_def & r) const
      {
        return range_def(std::min(r[0],(*this)[0]),std::max(r[1],(*this)[1]));
      }
    };

    struct rectangle_def:public three_tuple_t<range_def>
    {
      typedef three_tuple_t<range_def> base_t;

      rectangle_def
          (
              const coord_type & start_x,
              const coord_type & end_x,
              const coord_type & start_y,
              const coord_type & end_y,
              const coord_type & start_z,
              const coord_type & end_z
              )

      {
        (*this)[0] = range_def(start_x,end_x);
        (*this)[1] = range_def(start_y,end_y);
        (*this)[2] = range_def(start_z,end_z);
      }

      rectangle_def
          (
              const range_def &r1,
              const range_def &r2,
              const range_def &r3
              )
      {
        (*this)[0] = r1;
        (*this)[1] = r2;
        (*this)[2] = r3;
      }

      rectangle_def
          (
              const point_def &p1,
              const point_def &p2
              )
      {
        (*this)[0] = range_def(p1[0],p2[0]);
        (*this)[1] = range_def(p1[1],p2[1]);
        (*this)[2] = range_def(p1[2],p2[2]);

      }

      rectangle_def(){}

      inline point_def size() const
      {
        point_def ret;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          ret[i] = (*this)[i][1]-(*this)[i][0];

        return ret;
      }


      bool isInInterior ( const point_def & p ) const
      {
        bool ret = true;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          ret &= (*this)[i].isInOpen(p[i]);

        return ret;
      }

      bool contains ( const point_def & p ) const
      {
        bool ret = true;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          ret &= (*this)[i].isInClosed(p[i]);

        return ret;
      }

      bool isOnBoundry ( const point_def & p ) const
      {
        return ( contains ( p ) && !isInInterior ( p ) );
      }

      bool contains ( const rectangle_def &r ) const
      {
        bool ret = true;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          ret &= (*this)[i].contains(r[i]);

        return ret;
      }

      bool intersects ( const rectangle_def &r ) const
      {
        bool ret = true;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          ret &= (*this)[i].intersects(r[i]);

        return ret;
      }

      bool intersection(const rectangle_def & r,rectangle_def &ixn) const
      {
        bool ret = true;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          ret &= (*this)[i].intersection(r[i],ixn[i]);

        return ret;
      }

      rectangle_def bounding_box(const rectangle_def & r) const
      {
        rectangle_def ret;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          ret[i] = (*this)[i].range_union(r[i]);

        return ret;
      }

      point_def lower_corner() const
      {
        point_def c;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          c[i]= (*this)[i][0];

        return c;
      }

      point_def upper_corner() const
      {
        point_def c;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          c[i]= (*this)[i][1];

        return c;
      }

//      friend std::ostream& operator<< ( std::ostream& o, const rectangle_def& r )
//      {
//
//        o<<r[0];
//
//        for(size_t i = 1 ; i < rectangle_def::static_size;++i )
//          o<<"x"<<r[i];
//
//        return o;
//      }

      coord_type eff_dim() const
      {
        coord_type d;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          d += ((*this)[i][1] != (*this)[i][0]) ?(1):(0);

        return d;
      }
    };
  };

  template <typename coord_type,coord_type invalid_value>
      typename rectangle_complex<coord_type,invalid_value>::point_def
      rectangle_complex<coord_type,invalid_value>::point_def::zero =
      rectangle_complex<coord_type,invalid_value>::point_def(0,0,0);

  template <typename coord_type,coord_type invalid_value>
      typename rectangle_complex<coord_type,invalid_value>::point_def
      rectangle_complex<coord_type,invalid_value>::point_def::one =
      rectangle_complex<coord_type,invalid_value>::point_def(1,1,1);

  typedef int16_t                              cell_coord_t;
  typedef float                                cell_fn_t;
  typedef rectangle_complex<cell_coord_t>      rect_cmplx_t;
  typedef rect_cmplx_t::rectangle_def          rect_t;
  typedef rect_cmplx_t::point_def              cellid_t;
  typedef rect_cmplx_t::point_def              rect_point_t;
  typedef rect_cmplx_t::point_def              rect_size_t;
  typedef rect_cmplx_t::range_def              rect_range_t;
  typedef std::vector<cellid_t>                cellid_list_t;

  const uint gc_grid_dim = rect_t::base_t::static_size;

  enum eGradDirection
  {
    GRADIENT_DIR_DOWNWARD,
    GRADIENT_DIR_UPWARD,
    GRADIENT_DIR_COUNT,
  };

}

#endif
