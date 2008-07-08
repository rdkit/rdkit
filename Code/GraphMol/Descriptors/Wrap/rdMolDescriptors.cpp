// $Id$
//
//  Copyright (C) 2007 Greg Landrum
//
//   @@ All Rights Reserved  @@
//
#include <RDBoost/Wrap.h>
#include <GraphMol/Atom.h>
#include <GraphMol/GraphMol.h>

#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/AtomPairs.h>

#include <vector>

namespace python = boost::python;

namespace {
  std::vector<unsigned int> atomPairTypes(RDKit::Descriptors::AtomPairs::atomNumberTypes,
					  RDKit::Descriptors::AtomPairs::atomNumberTypes+sizeof(RDKit::Descriptors::AtomPairs::atomNumberTypes)/sizeof(unsigned int));
  python::tuple computeASAContribs(const RDKit::ROMol &mol,bool includeHs=true,
				   bool force=false){
    std::vector<double> contribs(mol.getNumAtoms());
    double hContrib=0.0;
    RDKit::Descriptors::getLabuteAtomContribs(mol,contribs,hContrib,includeHs,force);
    python::tuple pycontribs(contribs);
    return python::make_tuple(contribs,hContrib);
  }
  python::list computeCrippenContribs(const RDKit::ROMol &mol,
				       bool force=false){
    std::vector<double> logpContribs(mol.getNumAtoms());
    std::vector<double> mrContribs(mol.getNumAtoms());

    RDKit::Descriptors::getCrippenAtomContribs(mol,logpContribs,mrContribs,force);
    python::list pycontribs;
    for(unsigned int i=0;i<mol.getNumAtoms();++i){
      pycontribs.append(python::make_tuple(logpContribs[i],mrContribs[i]));
    }
    return pycontribs;
  }
  python::tuple calcCrippenDescriptors(const RDKit::ROMol &mol,bool includeHs=true,
				       bool force=false){
    double logp,mr;
    RDKit::Descriptors::CalcCrippenDescriptors(mol,logp,mr,includeHs,force);
    return python::make_tuple(logp,mr);
  }
}

BOOST_PYTHON_MODULE(rdMolDescriptors) {
  python::scope().attr("__doc__") =
    "Module containing functions to compute molecular descriptors"
    ;

  python::register_exception_translator<IndexErrorException>(&translate_index_error);
  python::register_exception_translator<ValueErrorException>(&translate_value_error);

  std::string docString = "";

  python::class_<python::object>("AtomPairsParameters")
    .def_readonly("version",
		  RDKit::Descriptors::AtomPairs::atomPairsVersion)
    .def_readonly("numTypeBits",
		  RDKit::Descriptors::AtomPairs::numTypeBits)
    .def_readonly("numPiBits",
		  RDKit::Descriptors::AtomPairs::numPiBits)
    .def_readonly("numBranchBits",
		  RDKit::Descriptors::AtomPairs::numBranchBits)
    .def_readonly("codeSize",
		  RDKit::Descriptors::AtomPairs::codeSize)
    .def_readonly("atomTypes",
		  atomPairTypes)
    .def_readonly("numPathBits",
		  RDKit::Descriptors::AtomPairs::numPathBits)
    .def_readonly("numAtomPairFingerprintBits",
		  RDKit::Descriptors::AtomPairs::numAtomPairFingerprintBits)
    ;
  docString="Returns the atom code (hash) for an atom";
  python::def("GetAtomPairAtomCode", RDKit::Descriptors::AtomPairs::getAtomCode,
              (python::arg("atom"), python::arg("branchSubtract")=0),
              docString.c_str());
  docString="Returns the atom-pair fingerprint for a molecule as an IntSparseIntVect";
  python::def("GetAtomPairFingerprint",
	      RDKit::Descriptors::AtomPairs::getAtomPairFingerprint,
	      python::arg("mol"),
              docString.c_str(),
	      python::return_value_policy<python::manage_new_object>());
  docString="Returns the topological-torsion fingerprint for a molecule as a LongIntSparseIntVect";
  python::def("GetTopologicalTorsionFingerprint",
	      RDKit::Descriptors::AtomPairs::getTopologicalTorsionFingerprint,
	      (python::arg("mol"),python::arg("targetSize")=4),
              docString.c_str(),
	      python::return_value_policy<python::manage_new_object>());


  docString="returns (as a list of 2-tuples) the contributions of each atom to\n"
    "the Wildman-Cripppen logp and mr value";
  python::def("_CalcCrippenContribs",
	      computeCrippenContribs,
	      (python::arg("mol"),
	       python::arg("force")=false),
              docString.c_str());
  docString="returns a 2-tuple with the Wildman-Crippen logp,mr values";
  python::def("CalcCrippenDescriptors",
	      calcCrippenDescriptors,
	      (python::arg("mol"),python::arg("includeHs")=true,
	       python::arg("force")=false),
              docString.c_str());
  python::scope().attr("__CalcCrippenDescriptors_version__")=
    RDKit::Descriptors::crippenVersion;
  
  docString="returns the Labute ASA value for a molecule";
  python::def("CalcLabuteASA",
	      RDKit::Descriptors::calcLabuteASA,
	      (python::arg("mol"),python::arg("includeHs")=true,
	       python::arg("force")=false),
              docString.c_str());
  python::scope().attr("__CalcLabuteASA_version__")=RDKit::Descriptors::labuteASAVersion;

  docString="returns a list of atomic contributions to the Labute ASA";
  python::def("_CalcLabuteASAContribs",
	      computeASAContribs,
	      (python::arg("mol"),python::arg("includeHs")=true,
	       python::arg("force")=false),
              docString.c_str());

  docString="returns the molecule's molecular weight";
  python::def("_CalcMolWt",
	      RDKit::Descriptors::CalcAMW,
	      (python::arg("mol"),python::arg("onlyHeavy")=false),
              docString.c_str());
  python::scope().attr("__CalcMolWt_version__")="1.0.0";

}
