// $Id$
//
//  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#define NO_IMPORT_ARRAY
#include "rdmolops.h"
#include <boost/python.hpp>
#include <numpy/arrayobject.h>
#include <string>
#include <math.h>

#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <RDBoost/Wrap.h>

namespace python = boost::python;

namespace RDKit{
  ROMol *addHs(const ROMol &orig,bool explicitOnly=false,bool addCoords=false){
    return MolOps::addHs(orig,explicitOnly,addCoords);
  }
  ROMol *removeHs(const ROMol &orig,bool implicitOnly=false){
    return MolOps::removeHs(orig,implicitOnly);
  }
  int getSSSR(ROMol &mol) {
    VECT_INT_VECT rings;
    int nr = MolOps::findSSSR(mol, rings);
    return nr;
  }

  void addRecursiveQuery(ROMol &mol,
                         const ROMol &query,
                         unsigned int atomIdx,
                         bool preserveExistingQuery){
    if(atomIdx>=mol.getNumAtoms()){
      throw_value_error("atom index exceeds mol.GetNumAtoms()");
    }
    RecursiveStructureQuery *q = new RecursiveStructureQuery(new ROMol(query));

    Atom *oAt=mol.getAtomWithIdx(atomIdx);
    if(!oAt->hasQuery()){
      QueryAtom qAt(*oAt);
      static_cast<RWMol &>(mol).replaceAtom(atomIdx,&qAt);
      oAt = mol.getAtomWithIdx(atomIdx);
    }

    
    if(!preserveExistingQuery){
      // FIX: this leaks a query atom on oAt
      oAt->setQuery(q);
    } else {
      oAt->expandQuery(q,Queries::COMPOSITE_AND);
    }
    
  }



  void sanitizeMol(ROMol &mol) {
    RWMol &wmol = static_cast<RWMol &>(mol);
    MolOps::sanitizeMol(wmol);
  }

  RWMol *getEditable(const ROMol &mol) {
    RWMol *res = static_cast<RWMol *>(new ROMol(mol,false));
    return res;
  }

  ROMol *getNormal(const RWMol &mol) {
    ROMol *res = static_cast<ROMol *>(new RWMol(mol));
    return res;
  }

  void kekulizeMol(ROMol &mol,bool clearAromaticFlags=false) {
    RWMol &wmol = static_cast<RWMol &>(mol);
    MolOps::Kekulize(wmol,clearAromaticFlags);
  }

  void setAromaticityMol(ROMol &mol){
    RWMol &wmol = static_cast<RWMol &>(mol);
    MolOps::setAromaticity(wmol);
  }

  VECT_INT_VECT getSymmSSSR(ROMol &mol) {
    VECT_INT_VECT rings;
    MolOps::symmetrizeSSSR(mol, rings);
    return rings;
  }
  
  PyObject *getDistanceMatrix(ROMol &mol, bool useBO=false,
                              bool useAtomWts=false,bool force=false,
                              const char *prefix=0) {
    int nats = mol.getNumAtoms();
    npy_intp dims[2];
    dims[0] = nats;
    dims[1] = nats;
    double *distMat;
    
    distMat = MolOps::getDistanceMat(mol, useBO, useAtomWts,force,prefix);
    
    PyArrayObject *res = (PyArrayObject *)PyArray_SimpleNew(2,dims,NPY_DOUBLE);
    
    memcpy(static_cast<void *>(res->data),
         static_cast<void *>(distMat),nats*nats*sizeof(double));
    
    return PyArray_Return(res);
  }

  PyObject *getAdjacencyMatrix(ROMol &mol, bool useBO=false,
                               int emptyVal=0,bool force=false,
                               const char *prefix=0) {
    int nats = mol.getNumAtoms();
    npy_intp  dims[2];
    dims[0] = nats;
    dims[1] = nats;

    double *tmpMat = MolOps::getAdjacencyMatrix(mol, useBO, emptyVal,force,prefix);
    
    PyArrayObject *res;
    if(useBO){
      // if we're using valence, the results matrix is made up of doubles
      res = (PyArrayObject *)PyArray_SimpleNew(2,dims,
                                              NPY_DOUBLE);
      memcpy(static_cast<void *>(res->data),
             static_cast<void *>(tmpMat),nats*nats*sizeof(double));
    } else {
      res = (PyArrayObject *)PyArray_SimpleNew(2,dims,
                                              NPY_INT);
      int *data = (int *)res->data;
      for(int i=0;i<nats;i++){
        for(int j=0;j<nats;j++){
          data[i*nats+j] = (int)round(tmpMat[i*nats+j]);
        }
      }
    }
    return PyArray_Return(res);
  }

  python::tuple GetMolFrags(const ROMol &mol,bool asMols){
    python::list res;

    if(!asMols){
      VECT_INT_VECT frags;
      MolOps::getMolFrags(mol,frags);

      for(unsigned int i=0;i<frags.size();++i){
        python::list tpl;
        for(unsigned int j=0;j<frags[i].size();++j){
          tpl.append(frags[i][j]);
        }
        res.append(python::tuple(tpl));
      }
    } else {
      std::vector<boost::shared_ptr<ROMol> > frags;
      frags=MolOps::getMolFrags(mol);
      for(unsigned int i=0;i<frags.size();++i){
        res.append(frags[i]);
      }
    }
    return python::tuple(res);
  }


  struct molops_wrapper {
    static void wrap() {
      std::string docString;

      // ------------------------------------------------------------------------
      docString="Kekulize, check valencies, set aromaticity, conjugation and hybridization\n\
\n\
    - The molecule is modified in place.\n\
\n\
    - If sanitization fails, an exception will be thrown\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
  NOTES:\n\
\n";
      python::def("SanitizeMol", sanitizeMol,
                  docString.c_str());

      // ------------------------------------------------------------------------
      docString="Get the smallest set of simple rings for a molecule.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use.\n\
\n\
  RETURNS: the number of rings found\n\
         This will be equal to NumBonds-NumAtoms+1 for single-fragment molecules.\n\
\n";
      python::def("GetSSSR", getSSSR, 
                  docString.c_str());
      
      // ------------------------------------------------------------------------
      docString="Get a symmetrized SSSR for a molecule.\n\
\n\
  The symmetrized SSSR is at least as large as the SSSR for a molecule.\n\
  In certain highly-symmetric cases (e.g. cubane), the symmetrized SSSR can be\n\
  a bit larger (i.e. the number of symmetrized rings is >= NumBonds-NumAtoms+1).\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use.\n\
\n\
  RETURNS: the number of rings found\n\
\n";
      python::def("GetSymmSSSR", getSymmSSSR,
                  docString.c_str());

      // ------------------------------------------------------------------------
      docString="Adds hydrogens to the graph of a molecule.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - explicitOnly: (optional) if this toggle is set, only explicit Hs will\n\
      be added to the molecule.  Default value is 0 (add implicit and explicit Hs).\n\
\n\
    - addCoords: (optional) if this toggle is set, The Hs will have 3D coordinates\n\
      set.  Default value is 0 (no 3D coords).\n\
\n\
  RETURNS: a new molecule with added Hs\n\
\n\
  NOTES:\n\
\n\
    - The original molecule is *not* modified.\n\
\n\
    - Much of the code assumes that Hs are not included in the molecular\n\
      topology, so be *very* careful with the molecule that comes back from\n\
      this function.\n\
\n";
      python::def("AddHs", addHs,
                  (python::arg("mol"),python::arg("explicitOnly")=false,
                   python::arg("addCoords")=false),
                  docString.c_str(),
                  python::return_value_policy<python::manage_new_object>());

      // ------------------------------------------------------------------------
      docString="Removes any hydrogens from the graph of a molecule.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - implicitOnly: (optional) if this toggle is set, only implicit Hs will\n\
      be removed from the graph.  Default value is 0 (remove implicit and explicit Hs).\n\
\n\
  RETURNS: a new molecule with the Hs removed\n\
\n\
  NOTES:\n\
\n\
    - The original molecule is *not* modified.\n\
\n";
      python::def("RemoveHs", removeHs,
                  (python::arg("mol"),python::arg("implicitOnly")=false),
                  docString.c_str(),
                  python::return_value_policy<python::manage_new_object>());

      python::def("MergeQueryHs", MolOps::mergeQueryHs,
                  (python::arg("mol")),
                  "merges hydrogens into their neighboring atoms as queries",
                  python::return_value_policy<python::manage_new_object>());

      // ------------------------------------------------------------------------
      docString="Returns the molecule's topological distance matrix.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - useBO: (optional) toggles use of bond orders in calculating the distance matrix.\n\
      Default value is 0.\n\
\n\
    - useAtomWts: (optional) toggles using atom weights for the diagonal elements of the\n\
      matrix (to return a \"Balaban\" distance matrix).\n\
      Default value is 0.\n\
\n\
    - force: (optional) forces the calculation to proceed, even if there is a cached value.\n\
      Default value is 0.\n\
\n\
    - prefix: (optional, internal use) sets the prefix used in the property cache\n\
      Default value is "".\n\
\n\
  RETURNS: a Numeric array of floats with the distance matrix\n\
\n";
      python::def("GetDistanceMatrix", getDistanceMatrix,
                  (python::arg("mol"),python::arg("useBO")=false,
                   python::arg("useAtomWts")=false,
                   python::arg("force")=false,
                   python::arg("prefix")=""),
                  docString.c_str());
      // ------------------------------------------------------------------------
      docString="Returns the molecule's adjacency matrix.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - useBO: (optional) toggles use of bond orders in calculating the matrix.\n\
      Default value is 0.\n\
\n\
    - emptyVal: (optional) sets the elements of the matrix between non-adjacent atoms\n\
      Default value is 0.\n\
\n\
    - force: (optional) forces the calculation to proceed, even if there is a cached value.\n\
      Default value is 0.\n\
\n\
    - prefix: (optional, internal use) sets the prefix used in the property cache\n\
      Default value is "".\n\
\n\
  RETURNS: a Numeric array of floats containing the adjacency matrix\n\
\n";
      python::def("GetAdjacencyMatrix", getAdjacencyMatrix, 
                  (python::arg("mol"), python::arg("useBO")=false,
                   python::arg("emptyVal")=0,
                   python::arg("force")=false,
                   python::arg("prefix")=""),
                  docString.c_str());

      
      // ------------------------------------------------------------------------
      docString="Kekulizes the molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - clearAromaticFlags: (optional) if this toggle is set, all atoms and bonds in the \n\
      molecule will be marked non-aromatic following the kekulization.\n\
      Default value is 0.\n\
\n\
  NOTES:\n\
\n\
    - The molecule is modified in place.\n\
\n";
      python::def("Kekulize", kekulizeMol,
                  (python::arg("mol"),python::arg("clearAromaticFlags")=false),
                  docString.c_str());
      

      // ------------------------------------------------------------------------
      docString="does aromaticity perception\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
  NOTES:\n\
\n\
    - The molecule is modified in place.\n\
\n";
      python::def("SetAromaticity", setAromaticityMol,
                  (python::arg("mol")),
                  docString.c_str());
      
   
      // ------------------------------------------------------------------------
      docString="Finds the disconnected fragments from a molecule.\n\
\n\
  For example, for the molecule 'CC(=O)[O-].[NH3+]C' GetMolFrags() returns\n\
  ((0, 1, 2, 3), (4, 5))\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
    - asMols: (optional) if this is provided and true, the fragments\n\
      will be returned as molecules instead of atom ids.\n\
\n\
  RETURNS: a tuple of tuples with IDs for the atoms in each fragment\n\
           or a tuple of molecules.\n\
\n";
      python::def("GetMolFrags", &GetMolFrags,
                  (python::arg("mol"),python::arg("asMols")=false),
                  docString.c_str());

      // ------------------------------------------------------------------------
      docString="Returns the formal charge for the molecule.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n";
      python::def("GetFormalCharge", &MolOps::getFormalCharge,docString.c_str());


      // ------------------------------------------------------------------------
      docString="Does the CIP stereochemistry assignment \n\
  for the molecule's atoms (R/S) and double bond (Z/E).\n\
  Chiral atoms will have a property '_CIPCode' indicating\n\
  their chiral code.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
    - cleanIt: (optional) if provided, atoms with a chiral specifier that aren't\n\
      actually chiral (e.g. atoms with duplicate substituents or only 2 substituents,\n\
      etc.) will have their chiral code set to CHI_UNSPECIFIED\n\
    - force: (optional) causes the calculation to be repeated, even if it has already\n\
      been done\n\
\n";
      python::def("AssignStereochemistry", MolOps::assignStereochemistry,
                  (python::arg("mol"),python::arg("cleanIt")=false,python::arg("force")=false),
                  docString.c_str());

      // ------------------------------------------------------------------------
      docString="Removes all stereochemistry info from the molecule.\n\
\n";
      python::def("RemoveStereochemistry", MolOps::removeStereochemistry,
                  (python::arg("mol")),
                  docString.c_str());

      // ------------------------------------------------------------------------
      docString="Sets the chiral tags on a molecule's atoms based on \n\
  a 3D conformation.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
    - confId: the conformer id to use, -1 for the default \n\
    - replaceExistingTags: if True, existing stereochemistry information will be cleared \n\
                           before running the calculation. \n\
\n";
      python::def("AssignAtomChiralTagsFromStructure", MolOps::assignChiralTypesFrom3D,
                  (python::arg("mol"),python::arg("confId")=-1,python::arg("replaceExistingTags")=true),
                  docString.c_str());

      docString="Set the wedging on single bonds in a molecule.\n \
   The wedging scheme used is that from Mol files.\n \
\n\
  ARGUMENTS:\n\
\n\
    - molecule: the molecule to update\n \
\n\
\n";
      python::def("WedgeMolBonds", WedgeMolBonds,
                  docString.c_str());

      // ------------------------------------------------------------------------
      docString="Adds a recursive query to an atom\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - query: the molecule to be used as the recursive query (this will be copied)\n\
\n\
    - atomIdx: the atom to modify\n\
\n\
    - preserveExistingQuery: (optional) if this is set, existing query information on the atom will be preserved\n\
\n\
  RETURNS: None\n\
\n";
      python::def("AddRecursiveQuery", addRecursiveQuery,
                  (python::arg("mol"),python::arg("query"),
                   python::arg("atomIdx"),python::arg("preserveExistingQuery")=true),
                  docString.c_str());



    };
  };
}

void wrap_molops() {
  RDKit::molops_wrapper::wrap();
}

