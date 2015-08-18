//
// Copyright (c) 2003-2005 greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_EXCEPTIONS_H
#define _RD_EXCEPTIONS_H
#include <stdexcept>
#include <string>

//! \brief Class to allow us to throw an \c IndexError from C++ and have
//!         it make it back to Python
//!
class IndexErrorException : public std::runtime_error
{
public:
  IndexErrorException(int i) : std::runtime_error("IndexErrorException"),_idx(i) {};
  int index () const { return _idx; };
  ~IndexErrorException () throw () {};
private:
  int _idx;
};

//! \brief Class to allow us to throw a \c ValueError from C++ and have
//!         it make it back to Python
//!
class ValueErrorException : public std::runtime_error
{
public:
  ValueErrorException(const std::string i) : std::runtime_error("ValueErrorException"), _value(i) {};
  ValueErrorException(const char *msg) : std::runtime_error("ValueErrorException"), _value(msg) {};
  std::string message () const { return _value; };
  ~ValueErrorException () throw () {};
private:
  std::string _value;
};


//! \brief Class to allow us to throw a \c KeyError from C++ and have
//!         it make it back to Python
//!
class KeyErrorException : public std::runtime_error
{
public:
  KeyErrorException(std::string key) : std::runtime_error("KeyErrorException"), _key(key) {};
  std::string key() const { return _key; };
  ~KeyErrorException () throw () {};
private:
  std::string _key;
};

#endif
