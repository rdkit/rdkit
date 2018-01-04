//
//  Copyright (C) 2018 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifdef _MSC_VER
#pragma warning(disable:4503)
#endif


#include <RDGeneral/Invariant.h>
#include <GraphMol/MolInterchange/MolInterchange.h>
#include <RDGeneral/FileParseException.h>

#include <sstream>
#include <exception>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>

namespace property_tree = boost::property_tree;


namespace RDKit {

namespace MolInterchange {

namespace {


} // end of anonymous namespace

std::vector<boost::shared_ptr<RWMol> > JSONDataStreamToMols(std::istream *inStream){
  PRECONDITION(inStream,"no stream");

#if 0
  property_tree::ptree pt;
  property_tree::read_json(*inStream, pt);

  if(pt.get<int>("moljson-header.version") != 10)
    throw FileParseException("bad version");
  property_tree::ptree mols=pt.get_child("molecules");
#else
rapidjson::IStreamWrapper isw(*inStream);
rapidjson::Document doc;
doc.ParseStream(isw);

#endif

  std::vector<boost::shared_ptr<RWMol> > res;
  return res;
}
std::vector<boost::shared_ptr<RWMol> >  JSONDataToMols(const std::string &jsonBlock){
  std::istringstream inStream(jsonBlock);
  return JSONDataStreamToMols(&inStream);
}
} // end of namespace MolInterchange
} // end of namespace RDKit
