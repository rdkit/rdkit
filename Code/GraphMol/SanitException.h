//
//  Copyright (C) 2002-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#ifndef _RD_SANITEXCEPTION_H
#define _RD_SANITEXCEPTION_H

#include <RDGeneral/types.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>

#include <string>
#include <vector>
#include <exception>

namespace RDKit {
  
  //! class for flagging sanitization errors
  class MolSanitizeException : public std::exception {
    public :
      MolSanitizeException(const char *msg) : _msg(msg) {};
      MolSanitizeException(const std::string msg) : _msg(msg) {};
      const char *message () const { return _msg.c_str(); };
      ~MolSanitizeException () throw () {};
    
    private :
      std::string _msg;
  };
}

#endif
