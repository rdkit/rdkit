// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#define NO_IMPORT_ARRAY
#include "rdmolops.h"
#include <boost/python.hpp>
#include <Numeric/arrayobject.h>
#include <string>
#include <math.h>

#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>

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

  PyObject* replaceSubstructures(const ROMol &orig,
                                 const ROMol &query,
                                 const ROMol &replacement,
                                 bool replaceAll=false) {
    std::vector<ROMOL_SPTR> v=replaceSubstructs(orig, query,
                                                replacement, replaceAll);
    PyObject *res=PyTuple_New(v.size());
    for(unsigned int i=0;i<v.size();++i){
      PyTuple_SetItem(res,i,
                      python::converter::shared_ptr_to_python(v[i]));
    }
    return res;
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

  VECT_INT_VECT getSymmSSSR(ROMol &mol) {
    VECT_INT_VECT rings;
    MolOps::symmetrizeSSSR(mol, rings);
    return rings;
  }
  
  PyObject *getDistanceMatrix(ROMol &mol, bool useBO=false,
                              bool useAtomWts=false,bool force=false,
                              const char *prefix=0) {
    int dims[2];
    int nats = mol.getNumAtoms();
    dims[0] = nats;
    dims[1] = nats;
    double *distMat;
    
    distMat = MolOps::getDistanceMat(mol, useBO, useAtomWts,force,prefix);
    
    PyArrayObject *res = (PyArrayObject *)PyArray_FromDims(2,dims,PyArray_DOUBLE);
    
    memcpy(static_cast<void *>(res->data),
         static_cast<void *>(distMat),nats*nats*sizeof(double));
    
    return (PyObject *) res;
  }

  PyObject *getAdjacencyMatrix(ROMol &mol, bool useBO=false,
                               int emptyVal=0,bool force=false,
                               const char *prefix=0) {
    int dims[2];
    int nats = mol.getNumAtoms();
    dims[0] = nats;
    dims[1] = nats;

    double *tmpMat = MolOps::getAdjacencyMatrix(mol, useBO, emptyVal,force,prefix);
    
    PyArrayObject *res;
    if(useBO){
      // if we're using valence, the results matrix is made up of doubles
      res = (PyArrayObject *)PyArray_FromDims(2,dims,
                                              PyArray_DOUBLE);
      memcpy(static_cast<void *>(res->data),
             static_cast<void *>(tmpMat),nats*nats*sizeof(double));
    } else {
      res = (PyArrayObject *)PyArray_FromDims(2,dims,
                                              PyArray_INT);
      int *data = (int *)res->data;
      for(int i=0;i<nats;i++){
        for(int j=0;j<nats;j++){
          data[i*nats+j] = (int)round(tmpMat[i*nats+j]);
        }
      }
    }
    return (PyObject *) res;
  }

  VECT_INT_VECT GetMolFrags(const ROMol &mol){
    VECT_INT_VECT frags;
    MolOps::getMolFrags(mol,frags);
    return frags;
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

      // ------------------------------------------------------------------------
      docString="Removes atoms matching a substructure query from a molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - query: the molecule to be used as a substructure query\n\
\n\
    - onlyFrags: (optional) if this toggle is set, atoms will only be removed if\n\
      the entire fragment in which they are found is matched by the query.\n\
      See below for examples.\n\
      Default value is 0 (remove the atoms whether or not the entire fragment matches)\n\
\n\
  RETURNS: a new molecule with the substructure removed\n\
\n\
  NOTES:\n\
\n\
    - The original molecule is *not* modified.\n\
\n\
  EXAMPLES:\n\
\n\
   The following examples substitute SMILES/SMARTS strings for molecules, you'd have\n\
   to actually use molecules:\n\
\n\
    - DeleteSubstructs('CCOC','OC') -> 'CC'\n\
\n\
    - DeleteSubstructs('CCOC','OC',1) -> 'CCOC'\n\
\n\
    - DeleteSubstructs('CCOCCl.Cl','Cl',1) -> 'CCOCCl'\n\
\n\
    - DeleteSubstructs('CCOCCl.Cl','Cl') -> 'CCOC'\n\
\n";
      python::def("DeleteSubstructs", deleteSubstructs,
                  (python::arg("mol"),python::arg("query"),
                   python::arg("onlyFrags")=false),
                  docString.c_str(),
                  python::return_value_policy<python::manage_new_object>());

      // ------------------------------------------------------------------------
      docString="Replaces atoms matching a substructure query in a molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - query: the molecule to be used as a substructure query\n\
\n\
    - replacement: the molecule to be used as the replacement\n\
\n\
    - replaceAll: (optional) if this toggle is set, all substructures matching\n\
      the query will be replaced in a single result, otherwise each result will\n\
      contain a separate replacement.\n\
      Default value is False (return multiple replacements)\n\
\n\
  RETURNS: a tuple of new molecules with the substructures replaced removed\n\
\n\
  NOTES:\n\
\n\
    - The original molecule is *not* modified.\n\
\n\
  EXAMPLES:\n\
\n\
   The following examples substitute SMILES/SMARTS strings for molecules, you'd have\n\
   to actually use molecules:\n\
\n\
    - ReplaceSubstructs('CCOC','OC','NC') -> ('CCNC',)\n\
\n\
    - ReplaceSubstructs('COCCOC','OC','NC') -> ('COCCNC','CNCCOC')\n\
\n\
    - ReplaceSubstructs('COCCOC','OC','NC',True) -> ('CNCCNC',)\n\
\n";
      python::def("ReplaceSubstructs", replaceSubstructures,
                  (python::arg("mol"),python::arg("query"),
                   python::arg("replacement"),
                   python::arg("replaceAll")=false),
                  docString.c_str());

      // ------------------------------------------------------------------------
      docString="Returns the molecule's distance matrix.\n\
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
      docString="Finds all subgraphs of a particular length in a molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - length: an integer with the target number of bonds for the subgraphs.\n\
\n\
    - useHs: (optional) toggles whether or not bonds to Hs that are part of the graph\n\
      should be included in the results.\n\
      Defaults to 0.\n\
\n\
    - verbose: (optional, internal use) toggles verbosity in the search algorithm.\n\
      Defaults to 0.\n\
\n\
  RETURNS: a tuple of 2-tuples with bond IDs\n\
\n\
  NOTES: \n\
\n\
   - Difference between _subgraphs_ and _paths_ :: \n\
\n\
       Subgraphs are potentially branched, whereas paths (in our \n\
       terminology at least) cannot be.  So, the following graph: \n\
\n\
            C--0--C--1--C--3--C\n\
                  |\n\
                  2\n\
                  |\n\
                  C\n\
  has 3 _subgraphs_ of length 3: (0,1,2),(0,1,3),(2,1,3)\n\
  but only 2 _paths_ of length 3: (0,1,3),(2,1,3)\n\
\n";
      python::def("FindAllSubgraphsOfLengthN", &findAllSubgraphsOfLengthN,
                  (python::arg("mol"),python::arg("length"),
                   python::arg("useHs")=false),
                  docString.c_str());
      // ------------------------------------------------------------------------
      docString="Finds unique subgraphs of a particular length in a molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - length: an integer with the target number of bonds for the subgraphs.\n\
\n\
    - useHs: (optional) toggles whether or not bonds to Hs that are part of the graph\n\
      should be included in the results.\n\
      Defaults to 0.\n\
\n\
    - useBO: (optional) Toggles use of bond orders in distinguishing one subgraph from\n\
      another.\n\
      Defaults to 1.\n\
\n\
  RETURNS: a tuple of tuples with bond IDs\n\
\n\
\n";
      python::def("FindUniqueSubgraphsOfLengthN", &findUniqueSubgraphsOfLengthN, 
                  (python::arg("mol"),python::arg("length"),
                   python::arg("useHs")=false,python::arg("useBO")=true),
                  docString.c_str());
                  
      // ------------------------------------------------------------------------
      docString="Finds all paths of a particular length in a molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - length: an integer with the target length for the paths.\n\
\n\
    - useBonds: (optional) toggles the use of bond indices in the paths.\n\
      Otherwise atom indices are used.  *Note* this behavior is different\n\
      from that for subgraphs.\n\
      Defaults to 1.\n\
\n\
  RETURNS: a tuple of tuples with IDs for the bonds.\n\
\n\
  NOTES: \n\
\n\
   - Difference between _subgraphs_ and _paths_ :: \n\
\n\
       Subgraphs are potentially branched, whereas paths (in our \n\
       terminology at least) cannot be.  So, the following graph: \n\
\n\
            C--0--C--1--C--3--C\n\
                  |\n\
                  2\n\
                  |\n\
                  C\n\
\n\
       has 3 _subgraphs_ of length 3: (0,1,2),(0,1,3),(2,1,3)\n\
       but only 2 _paths_ of length 3: (0,1,3),(2,1,3)\n\
\n";
      python::def("FindAllPathsOfLengthN", &findAllPathsOfLengthN, 
                  (python::arg("mol"),python::arg("length"),
                   python::arg("useBonds")=true,python::arg("useHs")=false),
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
\n\
  RETURNS: a tuple of tuples with IDs for the atoms in each fragment.\n\
\n";
      python::def("GetMolFrags", &GetMolFrags,docString.c_str());

      // ------------------------------------------------------------------------
      docString="Returns the formal charge for the molecule.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n";
      python::def("GetFormalCharge", &MolOps::getFormalCharge,docString.c_str());



      // ------------------------------------------------------------------------
      docString="Does the CIP chirality assignment (R/S) \n\
  for the molecule's atoms.\n\
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
      python::def("AssignAtomChiralCodes", MolOps::assignAtomChiralCodes,
                  (python::arg("mol"),python::arg("cleanIt")=false,python::arg("force")=false),
                  docString.c_str());


      // ------------------------------------------------------------------------
      docString="Does the CIP stereochemistry assignment (Z/E)\n\
   for the molecule's bonds .\n\
  Qualifying bonds will have a property '_CIPCode' indicating\n\
  their stereochemistry.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
    - cleanIt: (optional) ignored\n\
    - force: (optional) causes the calculation to be repeated, even if it has already\n\
      been done\n\
\n";
      python::def("AssignBondStereoCodes", MolOps::assignBondStereoCodes,
                  (python::arg("mol"),python::arg("cleanIt")=false,python::arg("force")=false),
                  docString.c_str());


      // ------------------------------------------------------------------------
      docString="Returns a \"Daylight\"-type fingerprint for a molecule\n\
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
  RETURNS: a DataStructs.ExplicitBitVect with _fpSize_ bits\n\
\n\
  ALGORITHM:\n\
\n\
   This algorithm functions by find all paths between minPath and maxPath in\n \
   length.  For each path:\n\
\n\
     1) The Balaban J value is calculated.\n\
\n\
     2) The 32 bit Balaban J value is used to seed a random-number generator\n\
\n\
     3) _nBitsPerPath_ random numbers are generated and used to set the corresponding\n\
        bits in the fingerprint\n\
\n\
\n";
      python::def("DaylightFingerprint", DaylightFingerprintMol,
                  (python::arg("mol"),python::arg("minPath")=1,
                   python::arg("maxPath")=7,python::arg("fpSize")=2048,
                   python::arg("nBitsPerHash")=4,python::arg("useHs")=true,
                   python::arg("tgtDensity")=0.0,python::arg("minSize")=128),
                  docString.c_str(),python::return_value_policy<python::manage_new_object>());



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
      docString="Replaces sidechains in a molecule with dummy atoms for their attachment points.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - coreQuery: the molecule to be used as a substructure query for recognizing the core\n\
\n\
  RETURNS: a new molecule with the sidechains removed\n\
\n\
  NOTES:\n\
\n\
    - The original molecule is *not* modified.\n\
\n\
  EXAMPLES:\n\
\n\
   The following examples substitute SMILES/SMARTS strings for molecules, you'd have\n\
   to actually use molecules:\n\
\n\
    - ReplaceSidechains('CCC1CCC1','C1CCC1') -> '[Xa]C1CCC1'\n\
\n\
    - ReplaceSidechains('CCC1CC1','C1CCC1') -> ''\n\
\n\
    - ReplaceSidechains('C1CC2C1CCC2','C1CCC1') -> '[Xa]C1CCC1[Xb]'\n\
\n";
      python::def("ReplaceSidechains", replaceSidechains,
                  (python::arg("mol"),python::arg("coreQuery")),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());

      // ------------------------------------------------------------------------
      docString="Removes the core of a molecule and labels the sidechains with dummy atoms.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - coreQuery: the molecule to be used as a substructure query for recognizing the core\n\
\n\
    - replaceDummies: toggles replacement of atoms that match dummies in the query\n\
\n\
  RETURNS: a new molecule with the core removed\n\
\n\
  NOTES:\n\
\n\
    - The original molecule is *not* modified.\n\
\n\
  EXAMPLES:\n\
\n\
   The following examples substitute SMILES/SMARTS strings for molecules, you'd have\n\
   to actually use molecules:\n\
\n\
    - ReplaceCore('CCC1CCC1','C1CCC1') -> 'CC[Xa]'\n\
\n\
    - ReplaceCore('CCC1CC1','C1CCC1') -> ''\n\
\n\
    - ReplaceCore('C1CC2C1CCC2','C1CCC1') -> '[Xa]C1CCC1[Xb]'\n\
\n\
    - ReplaceCore('C1CNCC1','N') -> '[Xa]CCCC[Xb]'\n\
\n\
    - ReplaceCore('C1CCC1CN','C1CCC1[*]',False) -> '[Xa]CN'\n\
\n";
      python::def("ReplaceCore", replaceCore,
                  (python::arg("mol"),python::arg("coreQuery"),
                   python::arg("replaceDummies")=true),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());

    };
  };
}

void wrap_molops() {
  RDKit::molops_wrapper::wrap();
}

