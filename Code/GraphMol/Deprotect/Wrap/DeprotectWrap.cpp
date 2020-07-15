//
//  Copyright (C) 2020 Brian P Kelley
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/python.h>
#include <RDBoost/Wrap.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Deprotect/Deprotect.h>

namespace python = boost::python;
namespace RDKit {
  // boost::python doesn't support unique_ptr so we convert to shared for
  //  python.
  template< typename T >
  inline
  std::vector< T > to_std_vector( const python::object& iterable )
  {
    return std::vector< T >( python::stl_input_iterator< T >( iterable ),
                             python::stl_input_iterator< T >( ) );
  }

  boost::shared_ptr<ROMol> DeprotectVectWrap(
		  const ROMol &mol,
		  const python::object &iterable) {
    //		  const std::vector<DeprotectData> &deprotections) {
    auto deprotections = to_std_vector<DeprotectData>(iterable);
    auto res = deprotect(mol, deprotections);
    auto m = boost::shared_ptr<ROMol>(res.get());
    res.release();
    return m;
  }

  boost::shared_ptr<ROMol> DeprotectWrap(const ROMol &mol) {
    auto res = deprotect(mol, getDeprotections());
    auto m = boost::shared_ptr<ROMol>(res.get());
    res.release();
    return m;
  }

  //! Make a copy so we don't try and change a const vector
  std::vector<DeprotectData> GetDeprotectionsWrap() {
    return getDeprotections();
  }
}

struct deprotect_wrap {
  static void wrap() {
    const char * constructor_doc = 
      "Construct a new DeprotectData instance.\n"
      "  >>> reaction_class = \"amine\"\n"
      "  >>> reaction_smarts = \"[C;R0][C;R0]([C;R0])([O;R0][C;R0](=[O;R0])[NX3;H0,H1:1])C>>[N:1]\"\n"
      "  >>> abbreviation = \"Boc\"\n"
      "  >>> full_name = \"tert-butyloxycarbonyl\"\n"
      "  >>> data = DeprotectData(reaction_class, reaction_smarts, abbreviation, full_name)\n"
      "  >>> assert data.isValid()\n"
      "\n";

    const char * deprotect_doc_string = 
      "DeprotectData class, contains a single deprotection reaction and information\n"
      "\n"
      " deprotectdata.deprotection_class - functional group being protected\n"
      " deprotectdata.reaction_smarts - reaction smarts used for deprotection\n"
      " deprotectdata.abbreviation - common abbreviation for the protecting group\n"
      " deprotectdata.full_name - full IUPAC name for the protecting group\n"
      "\n"
      "\n"
      ;
      
    python::class_<std::vector<RDKit::DeprotectData> >("DeprotectDataVect").
      def(python::vector_indexing_suite<std::vector<RDKit::DeprotectData>>());
							       
    python::class_<RDKit::DeprotectData>(
	    "DeprotectData",
	    deprotect_doc_string,
	    python::init<std::string, std::string,
	    std::string, std::string>(constructor_doc))
      .def_readonly("deprotection_class", &RDKit::DeprotectData::deprotection_class)
      .def_readonly("full_name", &RDKit::DeprotectData::full_name)
      .def_readonly("abbreviation", &RDKit::DeprotectData::abbreviation)
      .def_readonly("reaction_smarts", &RDKit::DeprotectData::reaction_smarts)
      .def("isValid", &RDKit::DeprotectData::isValid,
	   "Returns True if the DeprotectData has a valid reaction");

	  python::def("getDeprotections", &RDKit::GetDeprotectionsWrap,
		      "Return the default list of deprotections");
    
    python::def("Deprotect", &RDKit::DeprotectWrap,
		python::arg("mol"),
		"Return the deprotected version of the molecule.");
    
    python::def("Deprotect", &RDKit::DeprotectVectWrap,
		python::arg("mol"), python::arg("deprotections"),
		"Given a list of deprotections, return the deprotected version of the molecule.");
	  }
	  
};

void wrap_deprotect() { deprotect_wrap::wrap(); }
