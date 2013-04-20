//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_GASTEIGERPARAMS_H
#define _RD_GASTEIGERPARAMS_H

#include <RDGeneral/types.h>
#include <string>
#include <map>

namespace RDKit {
  extern std::string paramData;
  extern std::string additionalParamData;

  // this is a constant used during the iteration procedure for the hydrogen atoms
  // for the remaining atoms it is computed on the fly 
  const double IONXH = 20.02; 

  const double DAMP_SCALE = 0.5;
  const double DAMP = 0.5;
    
  class GasteigerParams {
    /* \brief Container for all the partial charge paramters
     * 
     * It is filled by the paramData string defined in GasteigerParams.cpp
     * The main data member is a STL map that take a pair<std::string, std::string>
     * of element name and mode (hybridization or bonding mode) and return a vector
     * of three parameters, used int eh ierative partial charges euqlization procedure
     */

  public:

    static const GasteigerParams *getParams(const std::string &paramData="");
    
    ~GasteigerParams() {
      d_paramMap.clear();
    }

    DOUBLE_VECT getParams(std::string elem, std::string mode,bool throwOnFailure=false) const {
      std::pair<std::string, std::string> query(elem, mode);
      std::map<std::pair<std::string, std::string>, DOUBLE_VECT>::const_iterator iter;
      iter=d_paramMap.find(query);
      if (iter != d_paramMap.end()) {
        return iter->second;
      }
      else {
	if(throwOnFailure){
	  std::string message = "ERROR: No Gasteiger Partial Charge parameters for Element: ";
	  message += elem;
	  message += " Mode: ";
	  message += mode;
	  throw message.c_str();
	} else {
          iter=d_paramMap.find(std::make_pair(std::string("X"),std::string("*")));
          if (iter != d_paramMap.end()) {
            return iter->second;
          } else {
            std::string message = "ERROR: Default Gasteiger Partial Charge parameters are missing";
            throw message.c_str();
          }
	}
      }
    }

    GasteigerParams(std::string paramData="");
  private:
    std::map<std::pair<std::string, std::string>, DOUBLE_VECT> d_paramMap;

    static class GasteigerParams *ds_instance;
  };
};

#endif
