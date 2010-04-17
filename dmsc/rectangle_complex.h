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

  typedef two_tuple_t<coord_type> size_def;

  struct point_def:public two_tuple_t<coord_type>
  {
    point_def ( const coord_type &x,const coord_type &y ) :
        two_tuple_t<coord_type>(x,y){}

    point_def () :two_tuple_t<coord_type>(-1,-1){}

    inline point_def operator/(const coord_type s) const
    {
      return point_def((*this)[0]/s,(*this)[1]/s);
    }

    inline size_def operator-(const point_def & o) const
    {
      return size_def((*this)[0]-o[0],(*this)[1]-o[1]);
    }
  };

  struct rectangle_def
  {
    typedef rectangle_def _SELF_NAME;

    point_def bl;
    point_def tr;


    inline size_def size() const
    {
      return tr - bl;
    }


    rectangle_def
        (
            const coord_type & start_x,
            const coord_type & start_y,
            const coord_type & end_x,
            const coord_type & end_y
            )

    {

      bl[0] = std::min ( start_x,end_x);
      bl[1] = std::min ( start_y,end_y);

      tr[0] = std::max ( start_x,end_x);
      tr[1] = std::max ( start_y,end_y);

    }

    rectangle_def
        (
            const point_def &p1,
            const point_def &p2
            )
    {
      bl[0] = std::min ( p1[0],p2[0] );
      bl[1] = std::min ( p1[1],p2[1] );

      tr[0] = std::max ( p1[0],p2[0] );
      tr[1] = std::max ( p1[1],p2[1] );

    }

    rectangle_def(){}


    bool isInInterior ( const point_def & p ) const
    {
      return
          (
              ( p[0] > bl[0] ) &&
              ( p[1] > bl[1] ) &&
              ( p[0] < tr[0] ) &&
              ( p[1] < tr[1] )
              );

    }

    bool contains ( const point_def & p ) const
    {
      return
          (
              ( p[0] >= bl[0] ) &&
              ( p[1] >= bl[1] ) &&
              ( p[0] <= tr[0] ) &&
              ( p[1] <= tr[1] )
              );
    }

    bool isOnBoundry ( const point_def & p ) const
    {
      return ( contains ( p ) && !isInInterior ( p ) );
    }

    bool contains ( const _SELF_NAME &rec ) const
    {
      return
          (
              ( bl[0] < rec.bl[0] ) &&
              ( bl[1] < rec.bl[1] ) &&
              ( tr[0] > rec.tr[0] ) &&
              ( tr[1] > rec.tr[1] )
              );
    }

    bool intersects ( const _SELF_NAME &rec ) const
    {
      return
          ! (
              ( ( tr[0] ) < ( rec.bl[0] ) ) ||
              ( ( rec.tr[0] ) < ( bl[0] ) ) ||
              ( ( tr[1] ) < ( rec.bl[1] ) ) ||
              ( ( rec.tr[1] ) < ( bl[1] ) )
              );
    }

    bool intersection(const rectangle_def & rect,rectangle_def &i) const
    {
      if(!intersects(rect))
        return false;

      coord_type l = std::max(rect.left(),left());
      coord_type r = std::min(rect.right(),right());
      coord_type b = std::max(rect.bottom(),bottom());
      coord_type t = std::min(rect.top(),top());

      i = rectangle_def(l,b,r,t);
      return true;
    }

    coord_type left() const
    {
      return bl[0];
    }

    coord_type right() const
    {
      return tr[0];
    }

    coord_type top() const
    {
      return tr[1];
    }

    coord_type bottom() const
    {
      return bl[1];
    }

    point_def bottom_left() const
    {
      return point_def ( left(), bottom() );
    }

    point_def top_left() const
    {
      return point_def ( left(),top() );
    }

    point_def top_right() const
    {
      return point_def ( right(),top() );
    }
    point_def bottom_right() const
    {
      return point_def ( right(),bottom() );
    }

    friend std::ostream& operator<< ( std::ostream& o, const rectangle_def& r )
    {
      return o<<"[ bl="<<r.bl<<" tr= "<<r.tr<<"]";
    }
  };
};

#endif
