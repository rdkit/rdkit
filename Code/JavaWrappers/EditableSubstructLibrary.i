
%{
#include <RDGeneral/export.h>
#include <GraphMol/SubstructLibrary/ChunkedHitlist.h>
#include <GraphMol/SubstructLibrary/SubstructLibrary.h>
#include <GraphMol/SubstructLibrary/EditableSubstructLibrary.h>
#include <GraphMol/SubstructLibrary/EditableSubstructLibraryTrustedSmilesWithPattern.h>
%}

%include <RDGeneral/export.h>

%include <GraphMol/SubstructLibrary/SubstructLibrary.h>
%include <GraphMol/SubstructLibrary/EditableSubstructLibrary.h>
%include <GraphMol/SubstructLibrary/ChunkedHitlist.h>

%template(EditableSubstructLibraryTrustedSmilesWithPatternBase) RDKit::EditableSubstructLibrary<std::string, RDKit::CachedTrustedSmilesMolHolder, RDKit::PatternHolder>;
%template(StringChunkedHitlist) RDKit::ChunkedHitlist<std::string>;
%template(StringChunkedHitlistPtr) shared_ptr<RDKit::ChunkedHitlist<std::string>>;                                                                                                                                                                                                                                   

%include <GraphMol/SubstructLibrary/EditableSubstructLibraryTrustedSmilesWithPattern.h>
