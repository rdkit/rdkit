%module RDKFuncs
%include "std_string.i"

%{
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
  using namespace RDKit;
#include <GraphMol/ROMol.h>
#include <GraphMol/Atom.h>
#include "RDKFuncs.h"

%}

//%ignore RDKitMol::dp_mol;
//%include "RDKFuncs.h"
%ignore getAllAtomsWithBookmark;
%ignore getAtomBookmarks;
%ignore getAllBondsWithBookmark;
%ignore getBondBookmarks;
%ignore getAtomNeighbors;
%ignore getAtomBonds;
%ignore getAtomPMap;
%ignore getBondPMap;
%ignore getVertices;
%ignore getEdges;
%ignore getTopology;
%ignore debugMol;
%ignore beginAtoms;
%ignore endAtoms;
%ignore beginAromaticAtoms;
%ignore endAromaticAtoms;
%ignore beginHeteros;
%ignore endHeteros;
%ignore beginQueryAtoms;
%ignore endQueryAtoms;
%ignore beginBonds;
%ignore endBonds;
%ignore getPropList;

%ignore getRingInfo;

%ignore getConformer;
%ignore addConformer;
%ignore beginConformers;
%ignore endConformers;
%include <GraphMol/ROMol.h>


%ignore copy;
%ignore getOwningMol;
%ignore setIdx;
%ignore setFormalCharge;
%ignore setNoImplicit;
%ignore setNumExplicitHs;
%ignore setIsAromatic;
%ignore setMass;
%ignore setDativeFlag;
%ignore clearDativeFlag;
%ignore setChiralTag;
%ignore setHybridizationType;
%ignore hasQuery;
%ignore setQuery;
%ignore getQuery;
%ignore expandQuery;
%ignore Match;
%ignore getPropList;
%ignore getPerturbationOrder;
%include <GraphMol/Atom.h>

%ignore setIsConjugated;
%ignore setBeginAtom;
%ignore setEndAtom;
%ignore setBondDir;
%ignore setStereo;
%ignore getStereoAtoms;
%include <GraphMol/Bond.h>


%include "RDKFuncs.h"


