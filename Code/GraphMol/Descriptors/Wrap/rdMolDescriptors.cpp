// $Id$
//
//  Copyright (C) 2007-2010 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/Wrap.h>
#include <GraphMol/Atom.h>
#include <GraphMol/GraphMol.h>
#include <numpy/arrayobject.h>
#include <boost/foreach.hpp>

#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <DataStructs/BitVects.h>

#include <vector>

namespace python = boost::python;

namespace {
  std::vector<unsigned int> atomPairTypes(RDKit::AtomPairs::atomNumberTypes,
					  RDKit::AtomPairs::atomNumberTypes+sizeof(RDKit::AtomPairs::atomNumberTypes)/sizeof(unsigned int));
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
    res = RDKit::AtomPairs::getAtomPairFingerprint(mol,minLength,maxLength,
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

    if(targetSize*RDKit::AtomPairs::codeSize>64){
      std::ostringstream errout;
      errout << "Maximum supported topological torsion path length is " << 64/RDKit::AtomPairs::codeSize<<std::endl;
      throw_value_error(errout.str());
    }
    
    RDKit::SparseIntVect<boost::int64_t> *res;
    res = RDKit::AtomPairs::getTopologicalTorsionFingerprint(mol,targetSize,vect);
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
    res = RDKit::AtomPairs::getHashedTopologicalTorsionFingerprint(mol,nBits,targetSize,vect);
    if(vect) delete vect;
    return res;
  }

  ExplicitBitVect *GetHashedTopologicalTorsionFingerprintAsBitVect(const RDKit::ROMol &mol,
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
    ExplicitBitVect *res;
    res = RDKit::AtomPairs::getHashedTopologicalTorsionFingerprintAsBitVect(mol,nBits,targetSize,vect);
    if(vect) delete vect;
    return res;
  }

  RDKit::SparseIntVect<boost::uint32_t> *GetMorganFingerprint(const RDKit::ROMol &mol,
                                                              int radius,
                                                              python::object invariants,
                                                              python::object fromAtoms,
                                                              bool useChirality,
                                                              bool useBondTypes,
                                                              bool useFeatures){
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
    } else if(useFeatures){
      invars = new std::vector<boost::uint32_t>(mol.getNumAtoms());
      RDKit::MorganFingerprints::getFeatureInvariants(mol,*invars);
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
                                                    invars,froms,useChirality,useBondTypes);
    if(invars) delete invars;
    if(froms) delete froms;
    return res;
  }

  ExplicitBitVect *GetMorganFingerprintBV(const RDKit::ROMol &mol,
                                          int radius,
                                          unsigned int nBits,
                                          python::object invariants,
                                          python::object fromAtoms,
                                          bool useChirality,
                                          bool useBondTypes,
                                          bool useFeatures){
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
    } else if(useFeatures){
      invars = new std::vector<boost::uint32_t>(mol.getNumAtoms());
      RDKit::MorganFingerprints::getFeatureInvariants(mol,*invars);
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
    ExplicitBitVect *res;
    res = RDKit::MorganFingerprints::getFingerprintAsBitVect(mol,
                                                             static_cast<unsigned int>(radius),
                                                             nBits,
                                                             invars,froms,useChirality,
                                                             useBondTypes);
    if(invars) delete invars;
    if(froms) delete froms;
    return res;
  }

  python::list GetConnectivityInvariants(const RDKit::ROMol &mol,bool includeRingMembership){
    std::vector<boost::uint32_t> invars(mol.getNumAtoms());
    RDKit::MorganFingerprints::getConnectivityInvariants(mol,invars,includeRingMembership);
    python::list res;
    for(std::vector<boost::uint32_t>::const_iterator iv=invars.begin();iv!=invars.end();++iv){
      res.append(python::long_(*iv));
    }
    return res;
  }
  python::list GetFeatureInvariants(const RDKit::ROMol &mol){
    std::vector<boost::uint32_t> invars(mol.getNumAtoms());
    RDKit::MorganFingerprints::getFeatureInvariants(mol,invars);
    python::list res;
    for(std::vector<boost::uint32_t>::const_iterator iv=invars.begin();iv!=invars.end();++iv){
      res.append(python::long_(*iv));
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
    .setattr("version",
		  RDKit::AtomPairs::atomPairsVersion)
    .setattr("numTypeBits",
		  RDKit::AtomPairs::numTypeBits)
    .setattr("numPiBits",
		  RDKit::AtomPairs::numPiBits)
    .setattr("numBranchBits",
		  RDKit::AtomPairs::numBranchBits)
    .setattr("codeSize",
		  RDKit::AtomPairs::codeSize)
    .setattr("atomTypes",
		  atomPairTypes)
    .setattr("numPathBits",
		  RDKit::AtomPairs::numPathBits)
    .setattr("numAtomPairFingerprintBits",
		  RDKit::AtomPairs::numAtomPairFingerprintBits)
    ;
  docString="Returns the atom code (hash) for an atom";
  python::def("GetAtomPairAtomCode", RDKit::AtomPairs::getAtomCode,
              (python::arg("atom"), python::arg("branchSubtract")=0),
              docString.c_str());
  docString="Returns the atom-pair fingerprint for a molecule as an IntSparseIntVect";
  python::def("GetAtomPairFingerprint", GetAtomPairFingerprint,
              (python::arg("mol"),
               python::arg("minLength")=1,
               python::arg("maxLength")=RDKit::AtomPairs::maxPathLen-1,
               python::arg("fromAtoms")=python::list()),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  python::def("GetHashedAtomPairFingerprint",
	      (RDKit::SparseIntVect<boost::int32_t> *(*)(const RDKit::ROMol&,unsigned int,unsigned int,unsigned int))RDKit::AtomPairs::getHashedAtomPairFingerprint,
	      (python::arg("mol"),
               python::arg("nBits")=2048,
               python::arg("minLength")=1,
               python::arg("maxLength")=RDKit::AtomPairs::maxPathLen-1),
              docString.c_str(),
	      python::return_value_policy<python::manage_new_object>());

  docString="Returns the atom-pair fingerprint for a molecule as an ExplicitBitVect";
  python::def("GetHashedAtomPairFingerprintAsBitVect",
	      RDKit::AtomPairs::getHashedAtomPairFingerprintAsBitVect,
	      (python::arg("mol"),
               python::arg("nBits")=2048,
               python::arg("minLength")=1,
               python::arg("maxLength")=RDKit::AtomPairs::maxPathLen-1),
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
  docString="Returns the topological-torsion fingerprint for a molecule as an ExplicitBitVect";
  python::def("GetHashedTopologicalTorsionFingerprintAsBitVect",
	      GetHashedTopologicalTorsionFingerprintAsBitVect,
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
               python::arg("fromAtoms")=python::list(),
               python::arg("useChirality")=false,
               python::arg("useBondTypes")=true,
               python::arg("useFeatures")=false),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  docString="Returns a Morgan fingerprint for a molecule as a bit vector";
  python::def("GetMorganFingerprintAsBitVect", GetMorganFingerprintBV,
              (python::arg("mol"),python::arg("radius"),python::arg("nBits")=2048,
               python::arg("invariants")=python::list(),
               python::arg("fromAtoms")=python::list(),
               python::arg("useChirality")=false,
               python::arg("useBondTypes")=true,
               python::arg("useFeatures")=false),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  python::scope().attr("_MorganFingerprint_version")=
    RDKit::MorganFingerprints::morganFingerprintVersion;
  docString="Returns connectivity invariants (ECFP-like) for a molecule.";
  python::def("GetConnectivityInvariants", GetConnectivityInvariants,
              (python::arg("mol"),python::arg("includeRingMembership")=true),
              docString.c_str());
  python::scope().attr("_ConnectivityInvariants_version")=
    RDKit::MorganFingerprints::morganConnectivityInvariantVersion;

  docString="Returns feature invariants (FCFP-like) for a molecule.";
  python::def("GetFeatureInvariants", GetFeatureInvariants,
              (python::arg("mol")),
              docString.c_str());
  python::scope().attr("_FeatureInvariants_version")=
    RDKit::MorganFingerprints::morganFeatureInvariantVersion;

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
  python::scope().attr("_CalcCrippenDescriptors_version")=
    RDKit::Descriptors::crippenVersion;
  
  docString="returns the Labute ASA value for a molecule";
  python::def("CalcLabuteASA",
	      RDKit::Descriptors::calcLabuteASA,
	      (python::arg("mol"),python::arg("includeHs")=true,
	       python::arg("force")=false),
              docString.c_str());
  python::scope().attr("_CalcLabuteASA_version")=RDKit::Descriptors::labuteASAVersion;

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
  python::scope().attr("_CalcMolWt_version")="1.0.0";

}
