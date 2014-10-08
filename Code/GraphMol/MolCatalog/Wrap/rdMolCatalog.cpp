// $Id$
//
//  Copyright (C) 2006 Greg Landrum
//
#include "rdMolCatalog.h"
#include <boost/python.hpp>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolCatalog/MolCatalog.h>
#include <GraphMol/MolCatalog/MolCatalogEntry.h>
#include <GraphMol/MolCatalog/MolCatalogParams.h>

namespace python = boost::python;
using namespace RDKit;
namespace {
  struct molcatalog_pickle_suite : python::pickle_suite
  {
    static python::tuple
    getinitargs(const MolCatalog& self)
    {
      std::string res;
      res = self.Serialize();
      return python::make_tuple(python::object(python::handle<>(PyBytes_FromStringAndSize(res.c_str(),res.length()))));
    };
  };

  struct molcatalogentry_pickle_suite : python::pickle_suite
  {
    static python::tuple
    getinitargs(const MolCatalogEntry& self)
    {
      std::string res;
      res = self.Serialize();
      return python::make_tuple(python::object(python::handle<>(PyBytes_FromStringAndSize(res.c_str(),res.length()))));
    };
  };


  unsigned int GetBitEntryId(const MolCatalog *self,unsigned int idx){
    if(idx > self->getFPLength())
      throw_index_error(idx);
    return self->getIdOfEntryWithBitId(idx);
  }

  unsigned int GetEntryBitId(const MolCatalog *self,unsigned int idx){
    if(idx > self->getNumEntries())
      throw_index_error(idx);
    return self->getEntryWithIdx(idx)->getBitId();
  }
  std::string GetEntryDescription(const MolCatalog *self,unsigned int idx){
    if(idx > self->getNumEntries())
      throw_index_error(idx);
    return self->getEntryWithIdx(idx)->getDescription();
  }
  std::string GetBitDescription(const MolCatalog *self,unsigned int idx){
    if(idx > self->getFPLength())
      throw_index_error(idx);
    return self->getEntryWithBitId(idx)->getDescription();
  }
  INT_VECT GetEntryDownIds(const MolCatalog *self,unsigned int idx){
    if(idx > self->getNumEntries())
      throw_index_error(idx);
    return self->getDownEntryList(idx);
  }

  unsigned int AddEntry(MolCatalog *self,MolCatalogEntry *entry){
    MolCatalogEntry *cpy=new MolCatalogEntry(*entry);
    return self->addEntry(cpy);
    //return self->addEntry(entry);
  }

  void catalogEntrySetMol(MolCatalogEntry *self,const ROMol *mol){
    ROMol *cpy = new ROMol(*mol);
    self->setMol(cpy);
  }
  
  const ROMol &catalogEntryGetMol(MolCatalogEntry &self){
    return *self.getMol();
  }

  MolCatalog *createMolCatalog(){
    return new MolCatalog(new MolCatalogParams());
  }
  struct MolCatalog_wrapper {
    static void wrap() {

      python::class_<MolCatalog>("MolCatalog",python::init<const std::string &>())
	.def("GetNumEntries", &MolCatalog::getNumEntries)
	.def("GetFPLength", &MolCatalog::getFPLength)
	.def("Serialize", &MolCatalog::Serialize)

	.def("GetBitDescription", GetBitDescription)
	.def("GetBitEntryId", GetBitEntryId)

	.def("GetEntryBitId", GetEntryBitId)
	.def("GetEntryDescription", GetEntryDescription)
	.def("GetEntryDownIds", GetEntryDownIds)

	.def("AddEntry",AddEntry)
	.def("AddEdge",&MolCatalog::addEdge)
      
	// enable pickle support
	.def_pickle(molcatalog_pickle_suite())
	;
      python::def("CreateMolCatalog",createMolCatalog,
		  python::return_value_policy<python::manage_new_object>());
    };
  };
  struct MolCatalogEntry_wrapper {
    static void wrap() {

      python::class_<MolCatalogEntry>("MolCatalogEntry", python::init<>())
	.def(python::init<const std::string &>())
	.def("GetDescription", &MolCatalogEntry::getDescription)
	.def("SetDescription", &MolCatalogEntry::setDescription)
	.def("GetMol", catalogEntryGetMol,
	     python::return_internal_reference<1>())
	.def("SetMol", catalogEntrySetMol)
	.def("GetOrder", &MolCatalogEntry::getOrder)
	.def("SetOrder", &MolCatalogEntry::setOrder)


	// enable pickle support
	.def_pickle(molcatalogentry_pickle_suite())

	;
    };
  };
}


BOOST_PYTHON_MODULE(rdMolCatalog)
{
  python::register_exception_translator<IndexErrorException>(&translate_index_error);
  python::register_exception_translator<ValueErrorException>(&translate_value_error);
  MolCatalog_wrapper::wrap();
  MolCatalogEntry_wrapper::wrap();
}
