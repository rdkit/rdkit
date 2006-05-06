//
//  Copyright (C) 2005-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//

#ifndef _DATASTRUCTS_EXCEPTION_H_20050126
#define _DATASTRUCTS_EXCEPTION_H_20050126

class DatastructsException : public std::exception {
 public:
  //! construct with an error message
  DatastructsException(const char *msg) : _msg(msg) {};
  //! construct with an error message
  DatastructsException(const std::string msg) : _msg(msg) {};
  //! get the error message
  const char *message () const { return _msg.c_str(); };
  ~DatastructsException () throw () {};
  private:
  std::string _msg;
};

#endif
