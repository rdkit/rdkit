//
//  Copyright (C) 2003-2006 greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _SEQS_HPP_
#define _SEQS_HPP_

#include <GraphMol/RDKitBase.h>
#include <iostream>
#include <boost/python.hpp>
namespace python = boost::python;

namespace RDKit{
  // Note: T1 should be some iterator type,
  //       T2 is the value when we dereference T
  template<class T1, class T2> class ReadOnlySeq{
  private:
    T1 _start,_end,_pos;
    int _size;
  public:
    //ReadOnlySeq() {};
    ~ReadOnlySeq() {
      //std::cerr<<" GC: "<<this<<std::endl;
    }
    ReadOnlySeq(T1 start,T1 end) : _start(start), _end(end), _pos(start), _size(-1) {};
    ReadOnlySeq(const ReadOnlySeq<T1,T2> &other) :  _start(other._start), _end(other._end),
                                                    _pos(other._pos), _size(other._size) {};
    void reset() {
      //std::cerr << "**** reset ****" << this<< std::endl;
      _pos = _start;
    }
    ReadOnlySeq<T1,T2> *__iter__() {
      //std::cerr << "**** ITER ****" << this << std::endl;
      reset();
      //std::cerr << "  finish ****" << this << std::endl;
      return this;
    };
    T2 next() {
      //std::cerr << "\tnext: " << _pos._pos << " " << _end._pos << std::endl;
      if(_pos == _end){
	PyErr_SetString(PyExc_StopIteration,"End of sequence hit");
	throw python::error_already_set();
      }      
      T2 res = *_pos;
      ++_pos;
      return res;
    }
    T2 get_item(int which) {
      //std::cerr << "get_item: " <<which<< std::endl;
      if(which >= len() ){
	PyErr_SetString(PyExc_IndexError,"End of sequence hit");
	throw python::error_already_set();
      }
      T1 it=_start;
      for(int i=0;i<which;i++){
	++it;
      }
      return *it;
    }
    int len() {
      //std::cerr << "len " << std::endl;
      if(_size<0){
	_size = 0;
	for(T1 tmp=_start;tmp!=_end;tmp++){
	  //std::cerr << "\tincr: "  <<  std::endl;
	  _size++;
	}
      }
      //std::cerr << "\tret" << std::endl;
      return _size;
    }
  };

  typedef ReadOnlySeq<ROMol::AtomIterator,Atom*> AtomIterSeq;
  typedef ReadOnlySeq<ROMol::QueryAtomIterator,Atom*> QueryAtomIterSeq;
  typedef ReadOnlySeq<ROMol::BondIterator,Bond*> BondIterSeq;
}
#endif
