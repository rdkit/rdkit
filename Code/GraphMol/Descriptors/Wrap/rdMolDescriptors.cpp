// $Id$
//
//  Copyright (C) 2007-2011 Greg Landrum
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
#include <GraphMol/Fingerprints/MACCS.h>
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
  python::tuple computeTPSAContribs(const RDKit::ROMol &mol,
				   bool force=false){
    std::vector<double> contribs(mol.getNumAtoms());
    RDKit::Descriptors::getTPSAAtomContribs(mol,contribs,force);
    python::tuple pycontribs(contribs);
    return pycontribs;
  }

  python::list computeCrippenContribs(const RDKit::ROMol &mol,
                                      bool force=false,
                                      python::list atomTypes=python::list(),
                                      python::list atomTypeLabels=python::list()){
    std::vector<unsigned int> *tAtomTypes=0;
    std::vector<std::string> *tAtomTypeLabels=0;
    if(python::extract<unsigned int>(atomTypes.attr("__len__")())!=0){
      if(python::extract<unsigned int>(atomTypes.attr("__len__")())!=mol.getNumAtoms()){
        throw_value_error("if atomTypes vector is provided, it must be as long as the number of atoms");
      } else {
        tAtomTypes = new std::vector<unsigned int>(mol.getNumAtoms(),0);
      }
    }
    if(python::extract<unsigned int>(atomTypeLabels.attr("__len__")())!=0){
      if(python::extract<unsigned int>(atomTypeLabels.attr("__len__")())!=mol.getNumAtoms()){
        throw_value_error("if atomTypeLabels vector is provided, it must be as long as the number of atoms");
      } else {
        tAtomTypeLabels = new std::vector<std::string>(mol.getNumAtoms(),"");
      }
    }

    std::vector<double> logpContribs(mol.getNumAtoms());
    std::vector<double> mrContribs(mol.getNumAtoms());
    
    RDKit::Descriptors::getCrippenAtomContribs(mol,logpContribs,mrContribs,force,tAtomTypes,tAtomTypeLabels);
    python::list pycontribs;
    for(unsigned int i=0;i<mol.getNumAtoms();++i){
      pycontribs.append(python::make_tuple(logpContribs[i],mrContribs[i]));
    }
    if(tAtomTypes){
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        atomTypes[i] = (*tAtomTypes)[i];
      }
      delete tAtomTypes;
    }
    if(tAtomTypeLabels){
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        atomTypeLabels[i] = (*tAtomTypeLabels)[i];
      }
      delete tAtomTypeLabels;
    }
    return pycontribs;
  }
  python::tuple calcCrippenDescriptors(const RDKit::ROMol &mol,bool includeHs=true,
				       bool force=false){
    double logp,mr;
    RDKit::Descriptors::calcCrippenDescriptors(mol,logp,mr,includeHs,force);
    return python::make_tuple(logp,mr);
  }

  RDKit::SparseIntVect<boost::int32_t> *GetAtomPairFingerprint(const RDKit::ROMol &mol,
                                                               unsigned int minLength,
                                                               unsigned int maxLength,
                                                               python::object fromAtoms,
                                                               python::object ignoreAtoms,
                                                               python::object atomInvariants,
                                                               bool includeChirality,
                                                               bool use2D){
    std::vector<boost::uint32_t> *fvect=pythonObjectToVect(fromAtoms,mol.getNumAtoms());
    std::vector<boost::uint32_t> *ivect=pythonObjectToVect(ignoreAtoms,mol.getNumAtoms());
    std::vector<boost::uint32_t> *invvect=pythonObjectToVect(atomInvariants,static_cast<unsigned int>(1<<RDKit::AtomPairs::codeSize));
    RDKit::SparseIntVect<boost::int32_t> *res;
    res = RDKit::AtomPairs::getAtomPairFingerprint(mol,minLength,maxLength,
                                                   fvect,ivect,invvect,includeChirality,use2D);
    if(fvect) delete fvect;
    if(ivect) delete ivect;
    if(invvect) delete invvect;
    return res;
  }
  RDKit::SparseIntVect<boost::int32_t> *GetHashedAtomPairFingerprint(const RDKit::ROMol &mol,
                                                                     unsigned int nBits,
                                                                     unsigned int minLength,
                                                                     unsigned int maxLength,
                                                                     python::object fromAtoms,
                                                                     python::object ignoreAtoms,
                                                                     python::object atomInvariants,
                                                                     bool includeChirality,
                                                                     bool use2D){
    std::vector<boost::uint32_t> *fvect=pythonObjectToVect(fromAtoms,mol.getNumAtoms());
    std::vector<boost::uint32_t> *ivect=pythonObjectToVect(ignoreAtoms,mol.getNumAtoms());
    std::vector<boost::uint32_t> *invvect=pythonObjectToVect(atomInvariants,static_cast<unsigned int>(1<<RDKit::AtomPairs::codeSize));
    RDKit::SparseIntVect<boost::int32_t> *res;
    res = RDKit::AtomPairs::getHashedAtomPairFingerprint(mol,nBits,minLength,maxLength,
                                                         fvect,ivect,invvect,includeChirality,
                                                         use2D);
    if(fvect) delete fvect;
    if(ivect) delete ivect;
    if(invvect) delete invvect;
    return res;
  }

  RDKit::SparseIntVect<boost::int64_t> *GetTopologicalTorsionFingerprint(const RDKit::ROMol &mol,
                                                                         unsigned int targetSize,
                                                                         python::object fromAtoms,
                                                                         python::object ignoreAtoms,
                                                                         python::object atomInvariants,
                                                                         bool includeChirality){
    std::vector<boost::uint32_t> *fvect=pythonObjectToVect(fromAtoms,mol.getNumAtoms());
    std::vector<boost::uint32_t> *ivect=pythonObjectToVect(ignoreAtoms,mol.getNumAtoms());
    std::vector<boost::uint32_t> *invvect=pythonObjectToVect(atomInvariants,static_cast<unsigned int>(1<<RDKit::AtomPairs::codeSize));
    if(targetSize*RDKit::AtomPairs::codeSize>64){
      std::ostringstream errout;
      errout << "Maximum supported topological torsion path length is " << 64/RDKit::AtomPairs::codeSize<<std::endl;
      throw_value_error(errout.str());
    }
    
    RDKit::SparseIntVect<boost::int64_t> *res;
    res = RDKit::AtomPairs::getTopologicalTorsionFingerprint(mol,targetSize,fvect,ivect,invvect,includeChirality);
    if(fvect) delete fvect;
    if(ivect) delete ivect;
    if(invvect) delete invvect;
    return res;
  }

  RDKit::SparseIntVect<boost::int64_t> *GetHashedTopologicalTorsionFingerprint(const RDKit::ROMol &mol,
                                                                               unsigned int nBits,
                                                                               unsigned int targetSize,
                                                                               python::object fromAtoms,
                                                                               python::object ignoreAtoms,
                                                                               python::object atomInvariants,
                                                                               bool includeChirality){
    std::vector<boost::uint32_t> *fvect=pythonObjectToVect(fromAtoms,mol.getNumAtoms());
    std::vector<boost::uint32_t> *ivect=pythonObjectToVect(ignoreAtoms,mol.getNumAtoms());
    std::vector<boost::uint32_t> *invvect=pythonObjectToVect(atomInvariants,static_cast<unsigned int>(1<<RDKit::AtomPairs::codeSize));
    RDKit::SparseIntVect<boost::int64_t> *res;
    res = RDKit::AtomPairs::getHashedTopologicalTorsionFingerprint(mol,nBits,targetSize,fvect,ivect,
                                                                   invvect,includeChirality);
    if(fvect) delete fvect;
    if(ivect) delete ivect;
    if(invvect) delete invvect;
    return res;
  }

  ExplicitBitVect *GetHashedTopologicalTorsionFingerprintAsBitVect(const RDKit::ROMol &mol,
                                                                   unsigned int nBits,
                                                                   unsigned int targetSize,
                                                                   python::object fromAtoms,
                                                                   python::object ignoreAtoms,
                                                                   python::object atomInvariants,
                                                                   unsigned int nBitsPerEntry,
                                                                   bool includeChirality
                                                                   ){
    std::vector<boost::uint32_t> *fvect=pythonObjectToVect(fromAtoms,mol.getNumAtoms());
    std::vector<boost::uint32_t> *ivect=pythonObjectToVect(ignoreAtoms,mol.getNumAtoms());
    std::vector<boost::uint32_t> *invvect=pythonObjectToVect(atomInvariants,static_cast<unsigned int>(1<<RDKit::AtomPairs::codeSize));
    ExplicitBitVect *res;
    res = RDKit::AtomPairs::getHashedTopologicalTorsionFingerprintAsBitVect(mol,nBits,targetSize,fvect,ivect,
                                                                            invvect,
                                                                            nBitsPerEntry,
                                                                            includeChirality);
    if(fvect) delete fvect;
    if(ivect) delete ivect;
    if(invvect) delete invvect;
    return res;
  }

  ExplicitBitVect *GetHashedAtomPairFingerprintAsBitVect(const RDKit::ROMol &mol,
                                                         unsigned int nBits,
                                                         unsigned int minLength,
                                                         unsigned int maxLength,
                                                         python::object fromAtoms,
                                                         python::object ignoreAtoms,
                                                         python::object atomInvariants,
                                                         unsigned int nBitsPerEntry,
                                                         bool includeChirality,
                                                         bool use2D
                                                         ){
    std::vector<boost::uint32_t> *fvect=pythonObjectToVect(fromAtoms,mol.getNumAtoms());
    std::vector<boost::uint32_t> *ivect=pythonObjectToVect(ignoreAtoms,mol.getNumAtoms());
    std::vector<boost::uint32_t> *invvect=pythonObjectToVect(atomInvariants,static_cast<unsigned int>(1<<RDKit::AtomPairs::codeSize));
    ExplicitBitVect *res;
    res = RDKit::AtomPairs::getHashedAtomPairFingerprintAsBitVect(mol,nBits,
                                                                  minLength,maxLength,
                                                                  fvect,ivect,invvect,
                                                                  nBitsPerEntry,includeChirality,
                                                                  use2D);
    if(fvect) delete fvect;
    if(ivect) delete ivect;
    if(invvect) delete invvect;
    return res;
  }

  namespace {
    double kappaHelper(double (* fn)(const RDKit::ROMol &,std::vector<double>*),const RDKit::ROMol &mol,python::object atomContribs){
      std::vector<double> *lContribs=0;
      if(atomContribs!=python::object()){
        // make sure the optional argument actually was a list
        python::list typecheck=python::extract<python::list>(atomContribs);

        if(python::extract<unsigned int>(typecheck.attr("__len__")())!=mol.getNumAtoms()){
          throw_value_error("length of atomContribs list != number of atoms");
        }          

        lContribs = new std::vector<double>(mol.getNumAtoms());
      }
      double res=fn(mol,lContribs);
      if(lContribs){
        python::list acl=python::extract<python::list>(atomContribs);
        for(unsigned int i=0;i<mol.getNumAtoms();++i){
          acl[i] = (*lContribs)[i];
        }
        delete lContribs;
      }
      return res;
    }
    double hkAlphaHelper(const RDKit::ROMol &mol,python::object atomContribs){
      return kappaHelper(RDKit::Descriptors::calcHallKierAlpha,mol,atomContribs);
    }

    RDKit::SparseIntVect<boost::uint32_t> *MorganFingerprintHelper(const RDKit::ROMol &mol,
                                 int radius,
                                 int nBits,
                                 python::object invariants,
                                 python::object fromAtoms,
                                 bool useChirality,
                                 bool useBondTypes,
                                 bool useFeatures,
                                 bool useCounts,
                                 python::object bitInfo){
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
      RDKit::MorganFingerprints::BitInfoMap *bitInfoMap=0;
      if(bitInfo!=python::object()){
        // make sure the optional argument actually was a dictionary
        python::dict typecheck=python::extract<python::dict>(bitInfo);
        bitInfoMap=new RDKit::MorganFingerprints::BitInfoMap();
      }
      RDKit::SparseIntVect<boost::uint32_t> *res;
      if(nBits<0){
        res = RDKit::MorganFingerprints::getFingerprint(mol,
                                                        static_cast<unsigned int>(radius),
                                                        invars,froms,useChirality,
                                                        useBondTypes,useCounts,false,bitInfoMap);
      } else {
        res = RDKit::MorganFingerprints::getHashedFingerprint(mol,
                                                        static_cast<unsigned int>(radius),
                                                        static_cast<unsigned int>(nBits),                                                                                                     invars,froms,useChirality,useBondTypes,
                                                        false,bitInfoMap);
      }
      if(bitInfoMap){
        bitInfo.attr("clear")();
        for(RDKit::MorganFingerprints::BitInfoMap::const_iterator iter=bitInfoMap->begin();
            iter!=bitInfoMap->end();++iter){
          const std::vector<std::pair<boost::uint32_t,boost::uint32_t> > &v=iter->second;
          python::list localL;
          for(std::vector<std::pair<boost::uint32_t,boost::uint32_t> >::const_iterator vIt=v.begin();
              vIt!=v.end();++vIt){
            localL.append(python::make_tuple(vIt->first,vIt->second));
          }
          bitInfo[iter->first]=python::tuple(localL);
        }
        delete bitInfoMap;
      }
      if(invars) delete invars;
      if(froms) delete froms;
      return res;
    }
  }
  RDKit::SparseIntVect<boost::uint32_t> *GetMorganFingerprint(const RDKit::ROMol &mol,
                                                              int radius,
                                                              python::object invariants,
                                                              python::object fromAtoms,
                                                              bool useChirality,
                                                              bool useBondTypes,
                                                              bool useFeatures,
                                                              bool useCounts,
                                                              python::object bitInfo){
    return MorganFingerprintHelper(mol,radius,-1,invariants,fromAtoms,useChirality,useBondTypes,
                                   useFeatures,useCounts,bitInfo);
  }
  RDKit::SparseIntVect<boost::uint32_t> *GetHashedMorganFingerprint(const RDKit::ROMol &mol,
                                                                    int radius,
                                                                    int nBits,
                                                              python::object invariants,
                                                              python::object fromAtoms,
                                                              bool useChirality,
                                                              bool useBondTypes,
                                                              bool useFeatures,
                                                              python::object bitInfo){
    return MorganFingerprintHelper(mol,radius,nBits,invariants,fromAtoms,useChirality,useBondTypes,
                                   useFeatures,true,bitInfo);
  }

  ExplicitBitVect *GetMorganFingerprintBV(const RDKit::ROMol &mol,
                                          int radius,
                                          unsigned int nBits,
                                          python::object invariants,
                                          python::object fromAtoms,
                                          bool useChirality,
                                          bool useBondTypes,
                                          bool useFeatures,
                                          python::object bitInfo){
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
    RDKit::MorganFingerprints::BitInfoMap *bitInfoMap=0;
    if(bitInfo!=python::object()){
      // make sure the optional argument actually was a dictionary
      python::dict typecheck=python::extract<python::dict>(bitInfo);
      bitInfoMap=new RDKit::MorganFingerprints::BitInfoMap();
    }
    ExplicitBitVect *res;
    res = RDKit::MorganFingerprints::getFingerprintAsBitVect(mol,
                                                             static_cast<unsigned int>(radius),
                                                             nBits,
                                                             invars,froms,useChirality,
                                                             useBondTypes,
                                                             false,bitInfoMap);
    if(bitInfoMap){
      bitInfo.attr("clear")();
      for(RDKit::MorganFingerprints::BitInfoMap::const_iterator iter=bitInfoMap->begin();
          iter!=bitInfoMap->end();++iter){
        const std::vector<std::pair<boost::uint32_t,boost::uint32_t> > &v=iter->second;
        python::list localL;
        for(std::vector<std::pair<boost::uint32_t,boost::uint32_t> >::const_iterator vIt=v.begin();
            vIt!=v.end();++vIt){
          localL.append(python::make_tuple(vIt->first,vIt->second));
        }
        bitInfo[iter->first]=python::tuple(localL);
      }
      delete bitInfoMap;
    }
    if(invars) delete invars;
    if(froms) delete froms;
    return res;
  }

  python::list GetConnectivityInvariants(const RDKit::ROMol &mol,bool includeRingMembership){
    std::vector<boost::uint32_t> invars(mol.getNumAtoms());
    RDKit::MorganFingerprints::getConnectivityInvariants(mol,invars,includeRingMembership);
    python::list res;
    BOOST_FOREACH(boost::uint32_t iv,invars){
      res.append(python::long_(iv));
    }
    return res;
  }
  python::list GetFeatureInvariants(const RDKit::ROMol &mol){
    std::vector<boost::uint32_t> invars(mol.getNumAtoms());
    RDKit::MorganFingerprints::getFeatureInvariants(mol,invars);
    python::list res;
    BOOST_FOREACH(boost::uint32_t iv,invars){
      res.append(python::long_(iv));
    }
    return res;
  }

  python::list CalcSlogPVSA(const RDKit::ROMol &mol,
                            python::object bins,
                            bool force){
    std::vector<double> *lbins=0;
    if(bins){
      unsigned int nBins=python::extract<unsigned int>(bins.attr("__len__")());
      if(nBins){
        lbins = new std::vector<double>(nBins,0.0);
        for(unsigned int i=0;i<nBins;++i){
          (*lbins)[i] = python::extract<double>(bins[i]);
        }
      }
    }
    std::vector<double> res;
    res=RDKit::Descriptors::calcSlogP_VSA(mol,lbins,force);

    python::list pyres;
    BOOST_FOREACH(double dv,res){
      pyres.append(dv);
    }
    return pyres;
  }
  python::list CalcSMRVSA(const RDKit::ROMol &mol,
                            python::object bins,
                            bool force){
    std::vector<double> *lbins=0;
    if(bins){
      unsigned int nBins=python::extract<unsigned int>(bins.attr("__len__")());
      if(nBins){
        lbins = new std::vector<double>(nBins,0.0);
        for(unsigned int i=0;i<nBins;++i){
          (*lbins)[i] = python::extract<double>(bins[i]);
        }
      }
    }
    std::vector<double> res;
    res=RDKit::Descriptors::calcSMR_VSA(mol,lbins,force);

    python::list pyres;
    BOOST_FOREACH(double dv,res){
      pyres.append(dv);
    }
    return pyres;
  }
  python::list CalcPEOEVSA(const RDKit::ROMol &mol,
                            python::object bins,
                            bool force){
    std::vector<double> *lbins=0;
    if(bins){
      unsigned int nBins=python::extract<unsigned int>(bins.attr("__len__")());
      if(nBins){
        lbins = new std::vector<double>(nBins,0.0);
        for(unsigned int i=0;i<nBins;++i){
          (*lbins)[i] = python::extract<double>(bins[i]);
        }
      }
    }
    std::vector<double> res;
    res=RDKit::Descriptors::calcPEOE_VSA(mol,lbins,force);

    python::list pyres;
    BOOST_FOREACH(double dv,res){
      pyres.append(dv);
    }
    return pyres;
  }
  python::list CalcMQNs(const RDKit::ROMol &mol,
                       bool force){
    std::vector<unsigned int> res;
    res=RDKit::Descriptors::calcMQNs(mol,force);

    python::list pyres;
    BOOST_FOREACH(unsigned int iv,res){
      pyres.append(iv);
    }
    return pyres;
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
              (python::arg("atom"), python::arg("branchSubtract")=0,
               python::arg("includeChirality")=false),
              docString.c_str());
  docString="Returns the atom-pair fingerprint for a molecule as an IntSparseIntVect";
  python::def("GetAtomPairFingerprint", GetAtomPairFingerprint,
              (python::arg("mol"),
               python::arg("minLength")=1,
               python::arg("maxLength")=RDKit::AtomPairs::maxPathLen-1,
               python::arg("fromAtoms")=0,
               python::arg("ignoreAtoms")=0,
               python::arg("atomInvariants")=0,
               python::arg("includeChirality")=false,
               python::arg("use2D")=true),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString="Returns the hashed atom-pair fingerprint for a molecule as an IntSparseIntVect";
  python::def("GetHashedAtomPairFingerprint",
	      GetHashedAtomPairFingerprint,
	      (python::arg("mol"),
               python::arg("nBits")=2048,
               python::arg("minLength")=1,
               python::arg("maxLength")=RDKit::AtomPairs::maxPathLen-1,
               python::arg("fromAtoms")=0,
               python::arg("ignoreAtoms")=0,
               python::arg("atomInvariants")=0,
               python::arg("includeChirality")=false,
               python::arg("use2D")=true),
              docString.c_str(),
	      python::return_value_policy<python::manage_new_object>());

  docString="Returns the atom-pair fingerprint for a molecule as an ExplicitBitVect";
  python::def("GetHashedAtomPairFingerprintAsBitVect",
	      GetHashedAtomPairFingerprintAsBitVect,
	      (python::arg("mol"),
               python::arg("nBits")=2048,
               python::arg("minLength")=1,
               python::arg("maxLength")=RDKit::AtomPairs::maxPathLen-1,
               python::arg("fromAtoms")=0,
               python::arg("ignoreAtoms")=0,
               python::arg("atomInvariants")=0,
               python::arg("nBitsPerEntry")=4,
               python::arg("includeChirality")=false,
               python::arg("use2D")=true),
              docString.c_str(),
	      python::return_value_policy<python::manage_new_object>());

  docString="Returns the topological-torsion fingerprint for a molecule as a LongIntSparseIntVect";
  python::def("GetTopologicalTorsionFingerprint",
	      GetTopologicalTorsionFingerprint,
	      (python::arg("mol"),python::arg("targetSize")=4,
               python::arg("fromAtoms")=0,
               python::arg("ignoreAtoms")=0,
               python::arg("atomInvariants")=0,
               python::arg("includeChirality")=false),
              docString.c_str(),
	      python::return_value_policy<python::manage_new_object>());
  docString="Returns the hashed topological-torsion fingerprint for a molecule as a LongIntSparseIntVect";
  python::def("GetHashedTopologicalTorsionFingerprint",
	      GetHashedTopologicalTorsionFingerprint,
	      (python::arg("mol"),
               python::arg("nBits")=2048,
               python::arg("targetSize")=4,
               python::arg("fromAtoms")=0,
               python::arg("ignoreAtoms")=0,
               python::arg("atomInvariants")=0,
               python::arg("includeChirality")=false),
              docString.c_str(),
	      python::return_value_policy<python::manage_new_object>());
  docString="Returns the topological-torsion fingerprint for a molecule as an ExplicitBitVect";
  python::def("GetHashedTopologicalTorsionFingerprintAsBitVect",
	      GetHashedTopologicalTorsionFingerprintAsBitVect,
	      (python::arg("mol"),
               python::arg("nBits")=2048,
               python::arg("targetSize")=4,
               python::arg("fromAtoms")=0,
               python::arg("ignoreAtoms")=0,
               python::arg("atomInvariants")=0,
               python::arg("nBitsPerEntry")=4,
               python::arg("includeChirality")=false),
              docString.c_str(),
	      python::return_value_policy<python::manage_new_object>());

  docString="Returns a Morgan fingerprint for a molecule";
  python::def("GetMorganFingerprint", GetMorganFingerprint,
              (python::arg("mol"),python::arg("radius"),
               python::arg("invariants")=python::list(),
               python::arg("fromAtoms")=python::list(),
               python::arg("useChirality")=false,
               python::arg("useBondTypes")=true,
               python::arg("useFeatures")=false,
               python::arg("useCounts")=true,
               python::arg("bitInfo")=python::object()),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  docString="Returns a hashed Morgan fingerprint for a molecule";
  python::def("GetHashedMorganFingerprint", GetHashedMorganFingerprint,
              (python::arg("mol"),python::arg("radius"),
               python::arg("nBits")=2048,
               python::arg("invariants")=python::list(),
               python::arg("fromAtoms")=python::list(),
               python::arg("useChirality")=false,
               python::arg("useBondTypes")=true,
               python::arg("useFeatures")=false,
               python::arg("bitInfo")=python::object()),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  docString="Returns a Morgan fingerprint for a molecule as a bit vector";
  python::def("GetMorganFingerprintAsBitVect", GetMorganFingerprintBV,
              (python::arg("mol"),python::arg("radius"),python::arg("nBits")=2048,
               python::arg("invariants")=python::list(),
               python::arg("fromAtoms")=python::list(),
               python::arg("useChirality")=false,
               python::arg("useBondTypes")=true,
               python::arg("useFeatures")=false,
               python::arg("bitInfo")=python::object()),
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
	       python::arg("force")=false,
               python::arg("atomTypes")=python::list(),
               python::arg("atomTypeLabels")=python::list()),
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

  docString="returns the TPSA value for a molecule";
  python::def("CalcTPSA",
	      RDKit::Descriptors::calcTPSA,
	      (python::arg("mol"),
	       python::arg("force")=false),
              docString.c_str());
  python::scope().attr("_CalcTPSA_version")=RDKit::Descriptors::tpsaVersion;

  docString="returns a list of atomic contributions to the TPSA";
  python::def("_CalcTPSAContribs",
	      computeTPSAContribs,
	      (python::arg("mol"),
	       python::arg("force")=false),
              docString.c_str());

  docString="returns the molecule's molecular weight";
  python::def("_CalcMolWt",
	      RDKit::Descriptors::calcAMW,
	      (python::arg("mol"),python::arg("onlyHeavy")=false),
              docString.c_str());
  python::scope().attr("_CalcMolWt_version")="1.0.0";

  docString="returns the molecule's exact molecular weight";
  python::def("CalcExactMolWt",
	      RDKit::Descriptors::calcExactMW,
	      (python::arg("mol"),python::arg("onlyHeavy")=false),
              docString.c_str());
  python::scope().attr("_CalcExactMolWt_version")="1.0.0";

  docString="returns the molecule's formula";
  python::def("CalcMolFormula",
	      RDKit::Descriptors::calcMolFormula,
	      (python::arg("mol"),
               python::arg("separateIsotopes")=false,
               python::arg("abbreviateHIsotopes")=true),
              docString.c_str());
  python::scope().attr("_CalcMolFormula_version")="1.3.0";

  docString="returns the number of Lipinski H-bond donors for a molecule";
  python::def("CalcNumLipinskiHBD",
	      RDKit::Descriptors::calcLipinskiHBD,
	      (python::arg("mol")),
              docString.c_str());
  python::scope().attr("_CalcNumLipinskiHBD_version")=RDKit::Descriptors::lipinskiHBDVersion;
  docString="returns the number of Lipinski H-bond acceptors for a molecule";
  python::def("CalcNumLipinskiHBA",
	      RDKit::Descriptors::calcLipinskiHBA,
	      (python::arg("mol")),
              docString.c_str());
  python::scope().attr("_CalcNumLipinskiHBA_version")=RDKit::Descriptors::lipinskiHBAVersion;
  docString="returns the number of H-bond donors for a molecule";
  python::def("CalcNumHBD",
	      RDKit::Descriptors::calcNumHBD,
	      (python::arg("mol")),
              docString.c_str());
  python::scope().attr("_CalcNumHBD_version")=RDKit::Descriptors::NumHBDVersion;
  docString="returns the number of H-bond acceptors for a molecule";
  python::def("CalcNumHBA",
	      RDKit::Descriptors::calcNumHBA,
	      (python::arg("mol")),
              docString.c_str());
  python::scope().attr("_CalcNumHBA_version")=RDKit::Descriptors::NumHBAVersion;
  docString="returns the number of rotatable bonds for a molecule. The strict option (default) does not count things like amide or ester bonds";
  python::def("CalcNumRotatableBonds",
	      RDKit::Descriptors::calcNumRotatableBonds,
	      (python::arg("mol"),
               python::arg("strict")=true),
              docString.c_str());
  python::scope().attr("_CalcNumRotatableBonds_version")=RDKit::Descriptors::NumRotatableBondsVersion;

  docString="returns the number of rings for a molecule";
  python::def("CalcNumRings",
	      RDKit::Descriptors::calcNumRings,
	      (python::arg("mol")),
              docString.c_str());
  python::scope().attr("_CalcNumRings_version")=RDKit::Descriptors::NumRingsVersion;

  docString="returns the number of aromatic rings for a molecule";
  python::def("CalcNumAromaticRings",
	      RDKit::Descriptors::calcNumAromaticRings,
	      (python::arg("mol")),
              docString.c_str());
  python::scope().attr("_CalcNumAromaticRings_version")=RDKit::Descriptors::NumAromaticRingsVersion;

  docString="returns the number of saturated rings for a molecule";
  python::def("CalcNumSaturatedRings",
	      RDKit::Descriptors::calcNumSaturatedRings,
	      (python::arg("mol")),
              docString.c_str());
  python::scope().attr("_CalcNumSaturatedRings_version")=RDKit::Descriptors::NumSaturatedRingsVersion;

  docString="returns the number of heterocycles for a molecule";
  python::def("CalcNumHeterocycles",
	      RDKit::Descriptors::calcNumHeterocycles,
	      (python::arg("mol")),
              docString.c_str());
  python::scope().attr("_CalcNumHeterocycles_version")=RDKit::Descriptors::NumHeterocyclesVersion;

  docString="returns the number of aromatic heterocycles for a molecule";
  python::def("CalcNumAromaticHeterocycles",
	      RDKit::Descriptors::calcNumAromaticHeterocycles,
	      (python::arg("mol")),
              docString.c_str());
  python::scope().attr("_CalcNumAromaticHeterocycles_version")=RDKit::Descriptors::NumAromaticHeterocyclesVersion;

  docString="returns the number of aromatic carbocycles for a molecule";
  python::def("CalcNumAromaticCarbocycles",
	      RDKit::Descriptors::calcNumAromaticCarbocycles,
	      (python::arg("mol")),
              docString.c_str());
  python::scope().attr("_CalcNumAromaticCarbocycles_version")=RDKit::Descriptors::NumAromaticCarbocyclesVersion;

  docString="returns the number of saturated heterocycles for a molecule";
  python::def("CalcNumSaturatedHeterocycles",
	      RDKit::Descriptors::calcNumSaturatedHeterocycles,
	      (python::arg("mol")),
              docString.c_str());
  python::scope().attr("_CalcNumSaturatedHeterocycles_version")=RDKit::Descriptors::NumSaturatedHeterocyclesVersion;

  docString="returns the number of saturated carbocycles for a molecule";
  python::def("CalcNumSaturatedCarbocycles",
	      RDKit::Descriptors::calcNumSaturatedCarbocycles,
	      (python::arg("mol")),
              docString.c_str());
  python::scope().attr("_CalcNumSaturatedCarbocycles_version")=RDKit::Descriptors::NumSaturatedCarbocyclesVersion;

  docString="returns the number of aliphatic (containing at least one non-aromatic bond) rings for a molecule";
  python::def("CalcNumAliphaticRings",
	      RDKit::Descriptors::calcNumAliphaticRings,
	      (python::arg("mol")),
              docString.c_str());
  python::scope().attr("_CalcNumAliphaticRings_version")=RDKit::Descriptors::NumAliphaticRingsVersion;

  docString="returns the number of aliphatic (containing at least one non-aromatic bond) heterocycles for a molecule";
  python::def("CalcNumAliphaticHeterocycles",
	      RDKit::Descriptors::calcNumAliphaticHeterocycles,
	      (python::arg("mol")),
              docString.c_str());
  python::scope().attr("_CalcNumAliphaticHeterocycles_version")=RDKit::Descriptors::NumAliphaticHeterocyclesVersion;

  docString="returns the number of aliphatic (containing at least one non-aromatic bond) carbocycles for a molecule";
  python::def("CalcNumAliphaticCarbocycles",
	      RDKit::Descriptors::calcNumAliphaticCarbocycles,
	      (python::arg("mol")),
              docString.c_str());
  python::scope().attr("_CalcNumAliphaticCarbocycles_version")=RDKit::Descriptors::NumAliphaticCarbocyclesVersion;

  docString="returns the number of heteroatoms for a molecule";
  python::def("CalcNumHeteroatoms",
	      RDKit::Descriptors::calcNumHeteroatoms,
	      (python::arg("mol")),
              docString.c_str());
  python::scope().attr("_CalcNumHeteroatoms_version")=RDKit::Descriptors::NumHeteroatomsVersion;
  docString="returns the number of amide bonds in a molecule";
  python::def("CalcNumAmideBonds",
	      RDKit::Descriptors::calcNumAmideBonds,
	      (python::arg("mol")),
              docString.c_str());
  python::scope().attr("_CalcNumAmideBonds_version")=RDKit::Descriptors::NumAmideBondsVersion;

  docString="returns the fraction of C atoms that are SP3 hybridized";
  python::def("CalcFractionCSP3",
	      RDKit::Descriptors::calcFractionCSP3,
	      (python::arg("mol")),
              docString.c_str());
  python::scope().attr("_CalcFractionCSP3_version")=RDKit::Descriptors::FractionCSP3Version;

  docString="returns the SlogP VSA contributions for a molecule";
  python::def("SlogP_VSA_",
              CalcSlogPVSA,
              (python::arg("mol"),
               python::arg("bins")=python::list(),
               python::arg("force")=false));
  docString="returns the SMR VSA contributions for a molecule";
  python::def("SMR_VSA_",
              CalcSMRVSA,
              (python::arg("mol"),
               python::arg("bins")=python::list(),
               python::arg("force")=false));
  docString="returns the PEOE VSA contributions for a molecule";
  python::def("PEOE_VSA_",
              CalcPEOEVSA,
              (python::arg("mol"),
               python::arg("bins")=python::list(),
               python::arg("force")=false));
  docString="returns the MQN descriptors for a molecule";
  python::def("MQNs_",
              CalcMQNs,
              (python::arg("mol"),
               python::arg("force")=false));

  docString="From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991)";
  python::def("CalcChiNv",
              RDKit::Descriptors::calcChiNv,
              (python::arg("mol"),
               python::arg("n"),
               python::arg("force")=false));
  python::scope().attr("_CalcChiNv_version")=RDKit::Descriptors::chiNvVersion;
  docString="From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991)";
  python::def("CalcChi0v",
              RDKit::Descriptors::calcChi0v,
              (python::arg("mol"),
               python::arg("force")=false));
  python::scope().attr("_CalcChi0v_version")=RDKit::Descriptors::chi0vVersion;
  docString="From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991)";
  python::def("CalcChi1v",
              RDKit::Descriptors::calcChi1v,
              (python::arg("mol"),
               python::arg("force")=false));
  python::scope().attr("_CalcChi1v_version")=RDKit::Descriptors::chi1vVersion;
  docString="From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991)";
  python::def("CalcChi2v",
              RDKit::Descriptors::calcChi2v,
              (python::arg("mol"),
               python::arg("force")=false));
  python::scope().attr("_CalcChi2v_version")=RDKit::Descriptors::chi2vVersion;
  docString="From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991)";
  python::def("CalcChi3v",
              RDKit::Descriptors::calcChi3v,
              (python::arg("mol"),
               python::arg("force")=false));
  python::scope().attr("_CalcChi3v_version")=RDKit::Descriptors::chi3vVersion;
  docString="From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991)";
  python::def("CalcChi4v",
              RDKit::Descriptors::calcChi4v,
              (python::arg("mol"),
               python::arg("force")=false));
  python::scope().attr("_CalcChi4v_version")=RDKit::Descriptors::chi4vVersion;
  docString="Similar to ChiXv, but uses uses nVal instead of valence. This makes a big difference after we get out of the first row.";
  python::def("CalcChiNn",
              RDKit::Descriptors::calcChiNn,
              (python::arg("mol"),
               python::arg("n"),
               python::arg("force")=false));
  python::scope().attr("_CalcChiNn_version")=RDKit::Descriptors::chiNnVersion;

  docString="Similar to ChiXv, but uses uses nVal instead of valence. This makes a big difference after we get out of the first row.";
  python::def("CalcChi0n",
              RDKit::Descriptors::calcChi0n,
              (python::arg("mol"),
               python::arg("force")=false));
  python::scope().attr("_CalcChi0n_version")=RDKit::Descriptors::chi0nVersion;

  docString="Similar to ChiXv, but uses uses nVal instead of valence. This makes a big difference after we get out of the first row.";
  python::def("CalcChi1n",
              RDKit::Descriptors::calcChi1n,
              (python::arg("mol"),
               python::arg("force")=false));
  python::scope().attr("_CalcChi1n_version")=RDKit::Descriptors::chi1nVersion;

  docString="Similar to ChiXv, but uses uses nVal instead of valence. This makes a big difference after we get out of the first row.";
  python::def("CalcChi2n",
              RDKit::Descriptors::calcChi2n,
              (python::arg("mol"),
               python::arg("force")=false));
  python::scope().attr("_CalcChi2n_version")=RDKit::Descriptors::chi2nVersion;

  docString="Similar to ChiXv, but uses uses nVal instead of valence. This makes a big difference after we get out of the first row.";
  python::def("CalcChi3n",
              RDKit::Descriptors::calcChi3n,
              (python::arg("mol"),
               python::arg("force")=false));
  python::scope().attr("_CalcChi3n_version")=RDKit::Descriptors::chi3nVersion;

  docString="Similar to ChiXv, but uses uses nVal instead of valence. This makes a big difference after we get out of the first row.";
  python::def("CalcChi4n",
              RDKit::Descriptors::calcChi4n,
              (python::arg("mol"),
               python::arg("force")=false));
  python::scope().attr("_CalcChi4n_version")=RDKit::Descriptors::chi4nVersion;

  docString="From equation (58) of Rev. Comp. Chem. vol 2, 367-422, (1991)";
  python::def("CalcHallKierAlpha",
              hkAlphaHelper,
              (python::arg("mol"),
               python::arg("atomContribs")=python::object()));
  python::scope().attr("_CalcHallKierAlpha_version")=RDKit::Descriptors::hallKierAlphaVersion;

  docString="From equations (58) and (59) of Rev. Comp. Chem. vol 2, 367-422, (1991)";
  python::def("CalcKappa1",
              RDKit::Descriptors::calcKappa1,
              (python::arg("mol")));
  python::scope().attr("_CalcKappa1_version")=RDKit::Descriptors::kappa1Version;
  docString="From equations (58) and (60) of Rev. Comp. Chem. vol 2, 367-422, (1991)";
  python::def("CalcKappa2",
              RDKit::Descriptors::calcKappa2,
              (python::arg("mol")));
  python::scope().attr("_CalcKappa2_version")=RDKit::Descriptors::kappa2Version;
  docString="From equations (58), (61) and (62) of Rev. Comp. Chem. vol 2, 367-422, (1991)";
  python::def("CalcKappa3",
              RDKit::Descriptors::calcKappa3,
              (python::arg("mol")));
  python::scope().attr("_CalcKappa3_version")=RDKit::Descriptors::kappa3Version;

  docString="Returns the MACCS keys for a molecule as an ExplicitBitVect";
  python::def("GetMACCSKeysFingerprint",
	      RDKit::MACCSFingerprints::getFingerprintAsBitVect,
	      (python::arg("mol")),
              docString.c_str(),
	      python::return_value_policy<python::manage_new_object>());
  
}
