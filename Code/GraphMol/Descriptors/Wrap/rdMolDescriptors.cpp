// $Id$
//
//  Copyright (C) 2007-2008 Greg Landrum
//
//   @@ All Rights Reserved  @@
//
#include <RDBoost/Wrap.h>
#include <GraphMol/Atom.h>
#include <GraphMol/GraphMol.h>

#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/AtomPairs.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/Fingerprints.h>

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

  RDKit::SparseIntVect<boost::int32_t> *GetAtomPairFingerprint(const RDKit::ROMol &mol,
                                                               unsigned int minLength,
                                                               unsigned int maxLength,
                                                               python::object fromAtoms){
    std::vector<boost::uint32_t> *vect=0;
    if(fromAtoms){
      vect = new std::vector<boost::uint32_t>;
      unsigned int nFrom=python::extract<unsigned int>(fromAtoms.attr("__len__")());
      for(unsigned int i=0;i<nFrom;++i){
        boost::uint32_t v=python::extract<boost::uint32_t>(fromAtoms[i]);
        if(v>=mol.getNumAtoms()){
          throw_value_error("atom index specified that is larger than the number of atoms");
        }
        vect->push_back(v);
      }
    }

    RDKit::SparseIntVect<boost::int32_t> *res;
    res = RDKit::Descriptors::AtomPairs::getAtomPairFingerprint(mol,minLength,maxLength,
                                                                vect);
    if(vect) delete vect;
    return res;
  }

  RDKit::SparseIntVect<boost::int64_t> *GetTopologicalTorsionFingerprint(const RDKit::ROMol &mol,
                                                                         unsigned int targetSize,
                                                                         python::object fromAtoms){
    std::vector<boost::uint32_t> *vect=0;
    if(fromAtoms){
      vect = new std::vector<boost::uint32_t>;
      unsigned int nFrom=python::extract<unsigned int>(fromAtoms.attr("__len__")());
      for(unsigned int i=0;i<nFrom;++i){
        boost::uint32_t v=python::extract<boost::uint32_t>(fromAtoms[i]);
        if(v>=mol.getNumAtoms()){
          throw_value_error("atom index specified that is larger than the number of atoms");
        }
        vect->push_back(v);
      }
    }

    if(targetSize*RDKit::Descriptors::AtomPairs::codeSize>64){
      std::ostringstream errout;
      errout << "Maximum supported topological torsion path length is " << 64/RDKit::Descriptors::AtomPairs::codeSize<<std::endl;
      throw_value_error(errout.str());
    }
    
    RDKit::SparseIntVect<boost::int64_t> *res;
    res = RDKit::Descriptors::AtomPairs::getTopologicalTorsionFingerprint(mol,targetSize,vect);
    if(vect) delete vect;
    return res;
  }

  RDKit::SparseIntVect<boost::int64_t> *GetHashedTopologicalTorsionFingerprint(const RDKit::ROMol &mol,
                                                                               unsigned int nBits,
                                                                         unsigned int targetSize,
                                                                         python::object fromAtoms){
    std::vector<boost::uint32_t> *vect=0;
    if(fromAtoms){
      vect = new std::vector<boost::uint32_t>;
      unsigned int nFrom=python::extract<unsigned int>(fromAtoms.attr("__len__")());
      for(unsigned int i=0;i<nFrom;++i){
        boost::uint32_t v=python::extract<boost::uint32_t>(fromAtoms[i]);
        if(v>=mol.getNumAtoms()){
          throw_value_error("atom index specified that is larger than the number of atoms");
        }
        vect->push_back(v);
      }
    }
    RDKit::SparseIntVect<boost::int64_t> *res;
    res = RDKit::Descriptors::AtomPairs::getHashedTopologicalTorsionFingerprint(mol,nBits,targetSize,vect);
    if(vect) delete vect;
    return res;
  }

  RDKit::SparseIntVect<boost::uint32_t> *GetMorganFingerprint(const RDKit::ROMol &mol,
                                                              int radius,
                                                              python::object invariants,
                                                              python::object fromAtoms){
    std::vector<boost::uint32_t> *invars=0;
    if(invariants){
      unsigned int nInvar=python::extract<unsigned int>(invariants.attr("__len__")());
      if(nInvar){
        if(nInvar!=mol.getNumAtoms()){
          throw_value_error("length of invariant vector != number of atoms");
        }
        invars = new std::vector<boost::uint32_t>(mol.getNumAtoms());
        for(unsigned int i=0;i<mol.getNumAtoms();++i){
          (*invars)[i] = python::extract<boost::uint32_t>(invariants[i]);
        }
      }
    }
    std::vector<boost::uint32_t> *froms=0;
    if(fromAtoms){
      unsigned int nFrom=python::extract<unsigned int>(fromAtoms.attr("__len__")());
      if(nFrom){
        froms = new std::vector<boost::uint32_t>();
        for(unsigned int i=0;i<nFrom;++i){
          froms->push_back(python::extract<boost::uint32_t>(fromAtoms[i]));
        }
      }
    }
    RDKit::SparseIntVect<boost::uint32_t> *res;
    res = RDKit::MorganFingerprints::getFingerprint(mol,
                                                    static_cast<unsigned int>(radius),
                                                    invars,froms);
    if(invars) delete invars;
    if(froms) delete froms;
    return res;
  }
  ExplicitBitVect *wrapLayeredFingerprint(const RDKit::ROMol &mol,unsigned int layerFlags,
                                          unsigned int minPath,unsigned int maxPath,
                                          unsigned int fpSize,
                                          double tgtDensity,
                                          unsigned int minSize,
                                          python::list atomCounts,
                                          ExplicitBitVect *includeOnlyBits){
    std::vector<unsigned int> *atomCountsV=0;
    if(atomCounts){
      atomCountsV = new std::vector<unsigned int>;
      unsigned int nAts=python::extract<unsigned int>(atomCounts.attr("__len__")());
      if(nAts<mol.getNumAtoms()){
        throw_value_error("atomCounts shorter than the number of atoms");
      }
      atomCountsV->resize(nAts);
      for(unsigned int i=0;i<nAts;++i){
        (*atomCountsV)[i] = python::extract<unsigned int>(atomCounts[i]);
      }
    }

    ExplicitBitVect *res;
    res = RDKit::LayeredFingerprintMol(mol,layerFlags,minPath,maxPath,fpSize,tgtDensity,minSize,atomCountsV,includeOnlyBits);

    if(atomCountsV){
      for(unsigned int i=0;i<atomCountsV->size();++i){
        atomCounts[i] = (*atomCountsV)[i];
      }
      delete atomCountsV;
    }
    
    return res;
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
  python::def("GetAtomPairFingerprint", GetAtomPairFingerprint,
              (python::arg("mol"),
               python::arg("minLength")=1,
               python::arg("maxLength")=RDKit::Descriptors::AtomPairs::maxPathLen-1,
               python::arg("fromAtoms")=python::list()),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  python::def("GetHashedAtomPairFingerprint",
	      (RDKit::SparseIntVect<boost::int32_t> *(*)(const RDKit::ROMol&,unsigned int,unsigned int,unsigned int))RDKit::Descriptors::AtomPairs::getHashedAtomPairFingerprint,
	      (python::arg("mol"),
               python::arg("nBits")=2048,
               python::arg("minLength")=1,
               python::arg("maxLength")=RDKit::Descriptors::AtomPairs::maxPathLen-1),
              docString.c_str(),
	      python::return_value_policy<python::manage_new_object>());
  docString="Returns the topological-torsion fingerprint for a molecule as a LongIntSparseIntVect";
  python::def("GetTopologicalTorsionFingerprint",
	      GetTopologicalTorsionFingerprint,
	      (python::arg("mol"),python::arg("targetSize")=4,
               python::arg("fromAtoms")=0),
              docString.c_str(),
	      python::return_value_policy<python::manage_new_object>());
  python::def("GetHashedTopologicalTorsionFingerprint",
	      GetHashedTopologicalTorsionFingerprint,
	      (python::arg("mol"),
               python::arg("nBits")=2048,
               python::arg("targetSize")=4,
               python::arg("fromAtoms")=0),
              docString.c_str(),
	      python::return_value_policy<python::manage_new_object>());

  docString="Returns a Morgan fingerprint for a molecule";
  python::def("GetMorganFingerprint", GetMorganFingerprint,
              (python::arg("mol"),python::arg("radius"),
               python::arg("invariants")=python::list(),
               python::arg("fromAtoms")=python::list()),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  python::scope().attr("__MorganFingerprint_version__")=
    RDKit::MorganFingerprints::morganFingerprintVersion;

      // ------------------------------------------------------------------------
      docString="Returns an RDKit topological fingerprint for a molecule\n\
\n\
  Explanation of the algorithm below.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - minPath: (optional) minimum number of bonds to include in the subgraphs\n\
      Defaults to 1.\n\
\n\
    - maxPath: (optional) maximum number of bonds to include in the subgraphs\n\
      Defaults to 7.\n\
\n\
    - fpSize: (optional) number of bits in the fingerprint\n\
      Defaults to 2048.\n\
\n\
    - nBitsPerPath: (optional) number of bits to set per path\n\
      Defaults to 4.\n\
\n\
    - useHs: (optional) include information about number of Hs on each\n\
      atom when calculating path hashes.\n\
      Defaults to 1.\n\
\n\
    - tgtDensity: (optional) fold the fingerprint until this minimum density has\n\
      been reached\n\
      Defaults to 0.\n\
\n\
    - minSize: (optional) the minimum size the fingerprint will be folded to when\n\
      trying to reach tgtDensity\n\
      Defaults to 128.\n\
\n\
  RETURNS: a DataStructs.ExplicitBitVect with _fpSize_ bits\n\
\n\
  ALGORITHM:\n\
\n\
   This algorithm functions by find all paths between minPath and maxPath in\n \
   length.  For each path:\n\
\n\
     1) A hash is calculated.\n\
\n\
     2) The hash is used to seed a random-number generator\n\
\n\
     3) _nBitsPerPath_ random numbers are generated and used to set the corresponding\n\
        bits in the fingerprint\n\
\n\
\n";
      python::def("RDKFingerprint", RDKit::RDKFingerprintMol,
                  (python::arg("mol"),python::arg("minPath")=1,
                   python::arg("maxPath")=7,python::arg("fpSize")=2048,
                   python::arg("nBitsPerHash")=4,python::arg("useHs")=true,
                   python::arg("tgtDensity")=0.0,python::arg("minSize")=128),
                  docString.c_str(),python::return_value_policy<python::manage_new_object>());
  python::scope().attr("__RDKFingerprint_version__")=
    RDKit::RDKFingerprintMolVersion;

      // ------------------------------------------------------------------------
      docString="Returns a layered fingerprint for a molecule\n\
\n\
  NOTE: This function is experimental. The API or results may change from\n\
    release to release.\n\
\n\
  Explanation of the algorithm below.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - layerFlags: (optional) which layers to include in the fingerprint\n\
      See below for definitions. Defaults to all.\n\
\n\
    - minPath: (optional) minimum number of bonds to include in the subgraphs\n\
      Defaults to 1.\n\
\n\
    - maxPath: (optional) maximum number of bonds to include in the subgraphs\n\
      Defaults to 7.\n\
\n\
    - fpSize: (optional) number of bits in the fingerprint\n\
      Defaults to 2048.\n\
\n\
    - tgtDensity: (optional) fold the fingerprint until this minimum density has\n\
      been reached\n\
      Defaults to 0.\n\
\n\
    - minSize: (optional) the minimum size the fingerprint will be folded to when\n\
      trying to reach tgtDensity\n\
      Defaults to 128.\n\
\n\
  RETURNS: a DataStructs.ExplicitBitVect with _fpSize_ bits\n\
\n\
  Layer definitions:\n\
     - 0x01: pure topology\n\
     - 0x02: bond order\n\
     - 0x04: atom types\n\
     - 0x08: presence of rings\n\
     - 0x10: ring sizes\n\
\n\
\n";
      python::def("LayeredFingerprint", wrapLayeredFingerprint,
                  (python::arg("mol"),
                   python::arg("layerFlags")=0xFFFFFFFF,
                   python::arg("minPath")=1,
                   python::arg("maxPath")=7,python::arg("fpSize")=2048,
                   python::arg("tgtDensity")=0.0,python::arg("minSize")=128,
                   python::arg("atomCounts")=python::list(),
                   python::arg("setOnlyBits")=(ExplicitBitVect *)0),
                  docString.c_str(),python::return_value_policy<python::manage_new_object>());
  python::scope().attr("__LayeredFingerprint_version__")=
    RDKit::LayeredFingerprintMolVersion;

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
