/* 
* 
*
*  Copyright (c) 2015, Greg Landrum
*  All rights reserved.
*
*  This file is part of the RDKit.
*  The contents are covered by the terms of the BSD license
*  which is included in the file license.txt, found at the root
*  of the RDKit source tree.
*
* taken from this stackexchange answer: http://stackoverflow.com/a/12768318
*/
%{
#include <boost/tuple/tuple.hpp>
%}

namespace boost {
  template <typename T1=void, typename T2=void, typename T3=void> 
  struct tuple;

  template <>
  struct tuple<void,void,void> {
  };

  template <typename T1>
  struct tuple<T1, void, void> {
    tuple(T1);
    %extend {
      T1 first() const {
        return boost::get<0>(*$self);
      }
    }
  };

  template <typename T1, typename T2>
  struct tuple <T1, T2, void> {
    tuple(T1,T2);
    %extend {
      T1 first() const {
        return boost::get<0>(*$self);
      }
      T2 second() const { 
        return boost::get<1>(*$self);
      }
    }
  };

  template <typename T1, typename T2, typename T3> 
  struct tuple <T1,T2,T3> {
    tuple(T1,T2,T3);
    %extend {
      T1 first() const {
        return boost::get<0>(*$self);
      }
      T2 second() const {
        return boost::get<1>(*$self);
      }
      T3 third() const {
        return boost::get<2>(*$self);
      }
    }
  };
}
