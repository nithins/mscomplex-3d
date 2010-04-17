/***************************************************************************
 *   Copyright (C) 2009 by nithin,,,   *
 *   nithin@gauss   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef __RECTANGLE_COMPLEX_H__
#define __RECTANGLE_COMPLEX_H__

#include <vector>
#include <algorithm>
#include <iostream>
#include <set>
#include <map>

#include <cpputils.h>
#include <logutil.h>


template <typename coord_type>
    class rectangle_complex
{

public:

  typedef three_tuple_t<coord_type> size_def;

  struct point_def:public three_tuple_t<coord_type>
  {
    point_def ( const coord_type &x,
                const coord_type &y,
                const coord_type &z,) :
    three_tuple_t<coord_type>(x,y,z){}

    point_def () :three_tuple_t<coord_type>(-1,-1,-1){}

    inline point_def operator/(const coord_type s) const
    {
      point_def ret;

      for(size_t i = 0 ; i < static_size;++i )
        ret[i] = (*this)[i]/s;

      return ret;
    }

    inline size_def operator-(const point_def & o) const
    {
      point_def ret;

      for(size_t i = 0 ; i < static_size;++i )
        ret[i] = (*this)[i]-s[i];

      return ret;
    }
  };

  struct range_def:public two_tuple_t
  {
    point_def ( const coord_type &x,
                const coord_type &y)
    {

      (*this)[0] = std::min ( x,y);
      (*this)[1] = std::max ( x,y);
    }

    inline bool isInOpen(const coord_type &c)
    {
      return (( (*this)[0] < c ) && (  c < (*this)[1] ));
    }

    inline bool isInClosed(const coord_type &c)
    {
      return (( (*this)[0] <= c ) && (  c <= (*this)[1] ));
    }

    inline bool isOnBndry(const coord_type &c)
    {
      return (( (*this)[0] == c ) || (  c == (*this)[1] ));
    }

    inline bool contains(const range_def & r)
    {
      return isInOpen(r[0]) && isInOpen(r[1]);
    }

    inline bool intersects(const range_def & r)
    {
      return !((r[0] > (*this)[1]) || ((*this)[0] > r[1]));
    }

    inline bool intersection(const range_def & r,range_def & i)
    {
      i = range_def(std::max(r[0],(*this)[0]),std::min(r[1],(*this)[1]));

      return intersects(r);
    }
  };

  struct rectangle_def:public three_tuple_t<range_def>
  {
    rectangle_def
        (
            const coord_type & start_x,
            const coord_type & end_x,
            const coord_type & start_y,
            const coord_type & end_y,
            const coord_type & start_z,
            const coord_type & end_z,
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

    inline size_def size() const
    {
      size_def ret;

      for(size_t i = 0 ; i < static_size;++i )
        ret[i] = (*this)[i][1]-(*this)[i][0];

      return ret;
    }


    bool isInInterior ( const point_def & p ) const
    {
      bool ret = true;

      for(size_t i = 0 ; i < static_size;++i )
        ret &= (*this)[i].isInOpen(p[i]);

      return ret;
    }

    bool contains ( const point_def & p ) const
    {
      bool ret = true;

      for(size_t i = 0 ; i < static_size;++i )
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

      for(size_t i = 0 ; i < static_size;++i )
        ret &= (*this)[i].contains(r[i]);

      return ret;
    }

    bool intersects ( const rectangle_def &r ) const
    {
      bool ret = true;

      for(size_t i = 0 ; i < static_size;++i )
        ret &= (*this)[i].intersects(r[i]);

      return ret;
    }

    bool intersection(const rectangle_def & r,rectangle_def &ixn) const
    {
      bool ret = true;

      for(size_t i = 0 ; i < static_size;++i )
        ret &= (*this)[i].intersection(r[i],ixn[i]);

      return ret;
    }

    //    coord_type left() const
    //    {
    //      return bl[0];
    //    }
    //
    //    coord_type right() const
    //    {
    //      return tr[0];
    //    }
    //
    //    coord_type top() const
    //    {
    //      return tr[1];
    //    }
    //
    //    coord_type bottom() const
    //    {
    //      return bl[1];
    //    }
    //
    //    point_def bottom_left() const
    //    {
    //      return point_def ( left(), bottom() );
    //    }
    //
    //    point_def top_left() const
    //    {
    //      return point_def ( left(),top() );
    //    }
    //
    //    point_def top_right() const
    //    {
    //      return point_def ( right(),top() );
    //    }
    //    point_def bottom_right() const
    //    {
    //      return point_def ( right(),bottom() );
    //    }

    friend std::ostream& operator<< ( std::ostream& o, const rectangle_def& r )
    {

      o<<r[0];

      for(size_t i = 1 ; i < rectangle_def::static_size-1;++i )
        o<<r[i]<<"x";

      return o;
    }
  };
};

#endif
