//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _RD_GASTEIGERPARAMS_H
#define _RD_GASTEIGERPARAMS_H

#include <RDGeneral/types.h>
#include <string>
#include <map>

namespace RDKit {
  extern std::string paramData;

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

    static GasteigerParams *getParams();
    
    ~GasteigerParams() {
      d_paramMap.clear();
    }

    DOUBLE_VECT getParams(std::string elem, std::string mode,bool throwOnFailure=false) {
      std::pair<std::string, std::string> query(elem, mode);
      if (d_paramMap.find(query) != d_paramMap.end()) {
        return d_paramMap[query];
      }
      else {
	if(throwOnFailure){
	  std::string message = "ERROR: No Gasteiger Partial Charge parameters for Element: ";
	  message += elem;
	  message += " Mode: ";
	  message += mode;
	  throw message.c_str();
	} else {
	  return d_paramMap[std::make_pair<std::string,std::string>("X","*")];
	}
      }
    }

  private:
    GasteigerParams();
    std::map<std::pair<std::string, std::string>, DOUBLE_VECT> d_paramMap;

    static class GasteigerParams *ds_instance;
  };
};

#endif
