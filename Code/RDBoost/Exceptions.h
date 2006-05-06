//
// Copyright (c) 2003-2005 greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _RD_EXCEPTIONS_H
#define _RD_EXCEPTIONS_H
#include <exception>
#include <string>

//! \brief Class to allow us to throw an \c IndexError from C++ and have
//!         it make it back to Python
//!
class IndexErrorException : public std::exception
{
public:
  IndexErrorException(int i) : _idx(i) {};
  int index () const { return _idx; };
  ~IndexErrorException () throw () {};
private:
  int _idx;
};

//! \brief Class to allow us to throw a \c ValueError from C++ and have
//!         it make it back to Python
//!
class ValueErrorException : public std::exception
{
public:
  ValueErrorException(const std::string i) : _value(i) {};
  ValueErrorException(const char *msg) : _value(msg) {};
  std::string message () const { return _value; };
  ~ValueErrorException () throw () {};
private:
  std::string _value;
};


//! \brief Class to allow us to throw a \c KeyError from C++ and have
//!         it make it back to Python
//!
class KeyErrorException : public std::exception
{
public:
  KeyErrorException(std::string key) : _key(key) {};
  std::string key() const { return _key; };
  ~KeyErrorException () throw () {};
private:
  std::string _key;
};

#endif
