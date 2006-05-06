//  $Id$
//
//  Copyright (C) 2003 Rational Discovery LLC
//   All Rights Reserved
//
#ifndef _TEMPLENUM_H_
#define _TEMPLENUM_H_
#include <exception>
#include <vector>
#include <GraphMol/RDKitBase.h>

namespace TemplateEnum {
  using namespace RDKit;

  class EnumException : public std::exception {
  public :
    EnumException(const char *msg) : _msg(msg) {};
    EnumException(const std::string msg) : _msg(msg) {};
    const char *message () const { return _msg.c_str(); };
    ~EnumException () throw () {};
  private :
    std::string _msg;
  };

  void orientSidechain(RWMol *mol,RWMol *sidechain,
		       int molAttachIdx,int sidechainAttachIdx);

  typedef std::vector< RWMOL_SPTR_VECT > VECT_RWMOL_SPTR_VECT;
  void markAttachmentPoints(RWMOL_SPTR *mol,char frontMarker='X');
  void markAttachmentPoints(RWMol *mol,char frontMarker='X');
  void prepareSidechains(RWMOL_SPTR_VECT *sidechains,char frontMarker='Y');
  RWMOL_SPTR_VECT enumerateLibrary(RWMol *mol,VECT_RWMOL_SPTR_VECT &sidechains,
				   bool orientSidechains=true);
  RWMOL_SPTR_VECT enumFromFiles(const char *templateName,
				std::vector<const char *> &sidechainName);

} // end of TemplateEnum namespace
#endif

