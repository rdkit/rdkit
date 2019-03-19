%import "SubstructLibrary.i"
%include "std_vector.i"

%{
#include <RDGeneral/export.h>
#include <GraphMol/SubstructLibrary/ChunkedHitlist.h>
#include <GraphMol/SubstructLibrary/SubstructLibrary.h>
#include <GraphMol/SubstructLibrary/EditableSubstructLibrary.h>
#include <GraphMol/SubstructLibrary/EditableSubstructLibraryTrustedSmilesWithPattern.h>
%}

%shared_ptr(RDKit::ChunkedHitlist<std::string>)
%include <RDGeneral/export.h>

%include <GraphMol/SubstructLibrary/SubstructLibrary.h>
%include <GraphMol/SubstructLibrary/EditableSubstructLibrary.h>
%include <GraphMol/SubstructLibrary/ChunkedHitlist.h>

%template(EditableSubstructLibraryTrustedSmilesWithPatternBase) RDKit::EditableSubstructLibrary<std::string, RDKit::CachedTrustedSmilesMolHolder, RDKit::PatternHolder>;
%template(StringChunkedHitlist) RDKit::ChunkedHitlist<std::string>;
%template(StringChunkedHitlistPtr) shared_ptr<RDKit::ChunkedHitlist<std::string>>;                                                                                                                                                                                                                                   

#ifdef SWIGJAVA
%typemap(jni) std::string RDKit::EditableSubstructLibraryTrustedSmilesWithPattern::makeStringFingerprint "jbyteArray"
%typemap(jtype) std::string RDKit::EditableSubstructLibraryTrustedSmilesWithPattern::makeStringFingerprint "byte[]"
%typemap(jstype) std::string RDKit::EditableSubstructLibraryTrustedSmilesWithPattern::makeStringFingerprint "byte[]"
%typemap(javaout) std::string RDKit::EditableSubstructLibraryTrustedSmilesWithPattern::makeStringFingerprint {
  return $jnicall;
}
%typemap(out) std::string RDKit::EditableSubstructLibraryTrustedSmilesWithPattern::makeStringFingerprint {
  $result = JCALL1(NewByteArray, jenv, $1.size());
  JCALL4(SetByteArrayRegion, jenv, $result, 0, $1.size(), (const jbyte*)$1.c_str());
}

%apply(char *STRING, size_t LENGTH) { (char *str, size_t len) };

%extend RDKit::EditableSubstructLibraryTrustedSmilesWithPattern {
    int addSmiles(const std::string smiles, char *str, size_t len) {
         std::string fp(str, len);
         return $self->addSmiles(smiles, fp);
    }
};
#endif

%ignore RDKit::EditableSubstructLibraryTrustedSmilesWithPattern::addSmiles;

%include <GraphMol/SubstructLibrary/EditableSubstructLibraryTrustedSmilesWithPattern.h>