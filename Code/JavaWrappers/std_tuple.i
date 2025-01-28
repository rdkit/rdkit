/* 
* 
*
*  Copyright (c) 2015-2022, Greg Landrum
*  All rights reserved.
*
*  This file is part of the RDKit.
*  The contents are covered by the terms of the BSD license
*  which is included in the file license.txt, found at the root
*  of the RDKit source tree.
*
* adapted from this stackexchange answer: https://stackoverflow.com/a/72849351
*/
// [1]
%{
#include <tuple>
#include <utility>
%}

// [2]    
#define make_getter(pos, type) const type& get##pos() const { return std::get<pos>(*$self); }
#define make_setter(pos, type) void set##pos(const type& val) { std::get<pos>(*$self) = val; }
#define make_ctorargN(pos, type) , type v##pos
#define make_ctorarg(first, ...) const first& v0 FOR_EACH(make_ctorargN, __VA_ARGS__)

// [3]
#define FE_0(...)
#define FE_1(action,a1) action(0,a1)
#define FE_2(action,a1,a2) action(0,a1) action(1,a2)
#define FE_3(action,a1,a2,a3) action(0,a1) action(1,a2) action(2,a3)
#define FE_4(action,a1,a2,a3,a4) action(0,a1) action(1,a2) action(2,a3) action(3,a4)
#define FE_5(action,a1,a2,a3,a4,a5) action(0,a1) action(1,a2) action(2,a3) action(3,a4) action(4,a5)

#define GET_MACRO(_1,_2,_3,_4,_5,NAME,...) NAME
%define FOR_EACH(action,...)
  GET_MACRO(__VA_ARGS__, FE_5, FE_4, FE_3, FE_2, FE_1, FE_0)(action,__VA_ARGS__)
%enddef

// [4]
%define %std_tuple(Name, ...)
%rename(Name) std::tuple<__VA_ARGS__>;
namespace std {
    struct tuple<__VA_ARGS__> {
        // [5]
        tuple(make_ctorarg(__VA_ARGS__));
        %extend {
            // [6]
            FOR_EACH(make_getter, __VA_ARGS__)
            FOR_EACH(make_setter, __VA_ARGS__)
            size_t __len__() const { return std::tuple_size<std::decay_t<decltype(*$self)>>{}; }
            %pythoncode %{
                # [7]
                def __getitem__(self, n):
                    if n >= len(self): raise IndexError()
                    return getattr(self, 'get%d' % n)()
                def __setitem__(self, n, val):
                    if n >= len(self): raise IndexError()
                    getattr(self, 'set%d' % n)(val)
            %}
        }
    };
}
%enddef