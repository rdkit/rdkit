// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/Wrap.h>

#include <boost/python.hpp>

#include <GraphMol/FragCatalog/FragCatGenerator.h>
#include <GraphMol/FragCatalog/FragCatParams.h>
#include <GraphMol/FragCatalog/FragCatalogEntry.h>


namespace python = boost::python;
namespace RDKit{

  struct fragcatalog_pickle_suite : python::pickle_suite
  {
    static python::tuple
    getinitargs(const FragCatalog& self)
    {
      std::string res;
      res = self.Serialize();
      return python::make_tuple(python::object(python::handle<>(PyBytes_FromStringAndSize(res.c_str(),res.length()))));
    };
  };
  unsigned int GetBitEntryId(const FragCatalog *self,unsigned int idx){
    if(idx > self->getFPLength())
      throw_index_error(idx);
    return self->getIdOfEntryWithBitId(idx);
  }

  unsigned int GetEntryBitId(const FragCatalog *self,unsigned int idx){
    if(idx > self->getNumEntries())
      throw_index_error(idx);
    return self->getEntryWithIdx(idx)->getBitId();
  }
  std::string GetEntryDescription(const FragCatalog *self,unsigned int idx){
    if(idx > self->getNumEntries())
      throw_index_error(idx);
    return self->getEntryWithIdx(idx)->getDescription();
  }
  std::string GetBitDescription(const FragCatalog *self,unsigned int idx){
    if(idx > self->getFPLength())
      throw_index_error(idx);
    return self->getEntryWithBitId(idx)->getDescription();
  }
  unsigned int GetEntryOrder(const FragCatalog *self,unsigned int idx){
    if(idx > self->getNumEntries())
      throw_index_error(idx);
    return self->getEntryWithIdx(idx)->getOrder();
  }
  unsigned int GetBitOrder(const FragCatalog *self,unsigned int idx){
    if(idx > self->getFPLength())
      throw_index_error(idx);
    return self->getEntryWithBitId(idx)->getOrder();
  }
  INT_VECT GetEntryFuncGroupIds(const FragCatalog *self,unsigned int idx){
    if(idx > self->getNumEntries())
      throw_index_error(idx);
    INT_VECT res;
    INT_INT_VECT_MAP gps = self->getEntryWithIdx(idx)->getFuncGroupMap();
    for(INT_INT_VECT_MAP::const_iterator i=gps.begin();i!=gps.end();i++){
      for(INT_VECT_CI ivci=i->second.begin();ivci!=i->second.end();ivci++){
	res.push_back(*ivci);
      }
    }
    return res;
  }
  INT_VECT GetBitFuncGroupIds(const FragCatalog *self,unsigned int idx){
    if(idx > self->getFPLength())
      throw_index_error(idx);
    INT_VECT res;
    INT_INT_VECT_MAP gps = self->getEntryWithBitId(idx)->getFuncGroupMap();
    for(INT_INT_VECT_MAP::const_iterator i=gps.begin();i!=gps.end();i++){
      for(INT_VECT_CI ivci=i->second.begin();ivci!=i->second.end();ivci++){
	res.push_back(*ivci);
      }
    }
    return res;
  }
  INT_VECT GetEntryDownIds(const FragCatalog *self,unsigned int idx){
    if(idx > self->getNumEntries())
      throw_index_error(idx);
    return self->getDownEntryList(idx);
  }


  DOUBLE_VECT GetBitDiscrims(const FragCatalog *self,unsigned int idx){
    if(idx > self->getFPLength())
      throw_index_error(idx);
    DOUBLE_VECT res;
    const FragCatalogEntry *entry=self->getEntryWithBitId(idx);
    Subgraphs::DiscrimTuple tmp=entry->getDiscrims();
    res.push_back(tmp.get<0>());
    res.push_back(tmp.get<1>());
    res.push_back(tmp.get<2>());
    return res;
  }


  
  struct fragcat_wrapper {
    static void wrap() {

      // FIX: none of the functions giving access to the entries in the catalog
      // are  being exposed to python
      // right now, adding entries for example should happen through the
      // FragCatGenerator
      
      python::class_<FragCatalog>("FragCatalog", python::init<FragCatParams *>())
	.def(python::init<const std::string &>())
	.def("GetNumEntries", &FragCatalog::getNumEntries)
	.def("GetFPLength", &FragCatalog::getFPLength)
	.def("GetCatalogParams", (FragCatParams* (FragCatalog::*)())&FragCatalog::getCatalogParams,
	     python::return_value_policy<python::reference_existing_object>())
	.def("Serialize", &FragCatalog::Serialize)

	.def("GetBitDescription", &GetBitDescription)
	.def("GetBitOrder", &GetBitOrder)
	.def("GetBitFuncGroupIds", &GetBitFuncGroupIds)
	.def("GetBitEntryId", &GetBitEntryId)

	.def("GetEntryBitId", &GetEntryBitId)
	.def("GetEntryDescription", &GetEntryDescription)
	.def("GetEntryOrder", &GetEntryOrder)
	.def("GetEntryFuncGroupIds", &GetEntryFuncGroupIds)
	.def("GetEntryDownIds", &GetEntryDownIds)


	.def("GetBitDiscrims", &GetBitDiscrims)

	// enable pickle support
	.def_pickle(fragcatalog_pickle_suite())

	;
    };
  };

} //end of namespace

void wrap_fragcat() {
  RDKit::fragcat_wrapper::wrap();
}



