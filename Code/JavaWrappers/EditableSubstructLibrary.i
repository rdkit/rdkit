%{
//#include <GraphMol/SubstructLibrary/ChunkedHitlist.h>
//#include <GraphMol/SubstructLibrary/SubstructLibrary.h>
#include <GraphMol/SubstructLibrary/EditableSubstructLibrary.h>
%}


//%template(ChunkedHitlist) RDKit::ChunkedHitlist<std::string>;
//%template(EditableSubstructLibraryTrustedSmilesWithPattern) RDKit::EditableSubstructLibrary<std::string, RDKit::CachedTrustedSmilesMolHolder, RDKit::PatternHolder>;

%include <GraphMol/SubstructLibrary/EditableSubstructLibrary.h>
//%include <GraphMol/SubstructLibrary/ChunkedHistlist.h>