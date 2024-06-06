/*
*
*  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
*  All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are
* met:
*
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above
*       copyright notice, this list of conditions and the following
*       disclaimer in the documentation and/or other materials provided
*       with the distribution.
*     * Neither the name of Novartis Institutes for BioMedical Research Inc.
*       nor the names of its contributors may be used to endorse or promote
*       products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
* OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
%module RDKFuncs

/* Suppress the unimportant warnings */
#pragma SWIG nowarn=503,516

%include <boost_shared_ptr.i>
%{
    #include <boost/shared_ptr.hpp>
    #include <boost/shared_array.hpp>
%}
// The actual definition isn't in the top level hpp file!
// The next two lines are to work around a problem caused by the fact that older versions of
// SWIG don't work with newer versions of boost.
#define BOOST_SP_NOEXCEPT
#define BOOST_SP_NOEXCEPT_WITH_ASSERT
#define BOOST_NOEXCEPT
#define BOOST_NO_CXX11_RVALUE_REFERENCES
#define BOOST_NO_CXX11_NULLPTR
%include <boost/smart_ptr/shared_array.hpp>

/* undefine RDKIT_<LIBNAME>_EXPORT macros */
%include <RDGeneral/RDExportMacros.h>
%include <RDGeneral/export.h>
/* Include the base types before anything that will utilize them */
#ifdef SWIGWIN
%include "../msvc_stdint.i"
#else
%include "../stdint.i"
#endif
%include "std_string.i"
%include "std_list.i"
%include "extend_std_vector.i"
%include "std_map.i"
%include "std_pair.i"
%include "carrays.i"

/*
 * Custom handler for longs.  The problem is described in swig-Bugs-2965875
 * and most of this solution is taken from the proposed patch in that bug report.
 * -------------------------------------------------------------------------
 * Define typemaps for `long`
 *
 * This is complicated by the fact `long` is 32-bits on some platforms
 * but is 64-bits on other platforms. We're just going to override the
 * important ones here.
 */
#if defined(SWIGWORDSIZE64)
typedef long long int		int64_t;
typedef unsigned long long int	uint64_t;
typedef long long int		int_least64_t;
typedef unsigned long long int	uint_least64_t;
typedef long long int		int_fast64_t;
typedef unsigned long long int	uint_fast64_t;
typedef long long int		intmax_t;
typedef unsigned long long int	uintmax_t;
%apply long long { long };
%apply const long long & { const long & };
%apply unsigned long long { unsigned long };
%apply const unsigned long long & { const unsigned long & };
/*
#elif defined(SWIGWORDSIZE32)
%apply int { long };
%apply const int & { const long & };
%apply unsigned int { unsigned long };
%apply const unsigned int & { const unsigned long & };
#else
#error "Neither SWIGWORDSIZE64 nor SWIGWORDSIZE32 is defined"
*/
#endif

%shared_ptr(std::exception)
%shared_ptr(RDKit::RDProps)
%shared_ptr(RDKit::Conformer)
%shared_ptr(RDKit::ROMol)
%shared_ptr(RDKit::RWMol)
%shared_ptr(RDKit::Atom)
%shared_ptr(RDKit::Bond)
%shared_ptr(RDKit::PeriodicTable)
%shared_ptr(Canon::MolStackElem)
%shared_ptr(RDKit::QueryAtom)
%shared_ptr(RDKit::QueryBond)
%shared_ptr(RDKit::QueryOps)
%shared_ptr(RDKit::MolBundle)
%shared_ptr(RDKit::FixedMolSizeMolBundle)
%shared_ptr(RDKit::MolSanitizeException)
%shared_ptr(RDKit::AtomSanitizeException)
%shared_ptr(RDKit::AtomValenceException)
%shared_ptr(RDKit::AtomKekulizeException)
%shared_ptr(RDKit::KekulizeException)
%shared_ptr(RDKit::SmilesParseException)
%shared_ptr(RDKit::MolPicklerException)
%shared_ptr(RDKit::RingInfo)
%shared_ptr(RDKit::ChemicalReaction)
%shared_ptr(ForceFields::ForceFieldContrib);
%shared_ptr(ForceFields::UFF::AngleBendContrib);
%shared_ptr(ForceFields::UFF::BondStretchContrib);
%shared_ptr(ForceFields::UFF::DistanceConstraintContrib);
%shared_ptr(ForceFields::UFF::vdWContrib);
%shared_ptr(ForceFields::UFF::TorsionAngleContrib);
%shared_ptr(ForceFields::UFF::InversionContrib);
%shared_ptr(RDKit::FilterCatalogEntry);

/* Some utility classes for passing arrays in and out */
%array_class(double, Double_Array);

/* Since documentation management is deprecated in SWIG 1.3, we're using the suggested workarounds.  Apply them
   here so that can be removed easily later */
// Documentation
%include "../Atom_doc.i"
%include "../Bond_doc.i"
%include "../BondIterators_doc.i"
%include "../ChemReactions_doc.i"
%include "../Conformer_doc.i"
%include "../ExplicitBitVect_doc.i"
%include "../PeriodicTable_doc.i"
%include "../Point3D_doc.i"
%include "../QueryAtom_doc.i"
%include "../QueryBond_doc.i"
%include "../RDKFuncs_doc.i"
%include "../RingInfo_doc.i"
%include "../ROMol_doc.i"
%include "../RWMol_doc.i"
%include "../SDMolSupplier_doc.i"
%include "../SmilesMolSupplier_doc.i"
%include "../SmilesWriter_doc.i"
%include "../TDTMolSupplier_doc.i"
%include "../TDTWriter_doc.i"
%include "../Transform2D_doc.i"
%include "../Transform3D_doc.i"
%include "../FilterCatalog_doc.i"
%include "../FilterCatalogParams_doc.i"
%include "../FilterCatalogs_doc.i"

// DO THIS BEFORE ANY OF THE OTHER INCLUDES
%include "../RDKitExceptions.i"

%include "../point.i"
// Need the types wrapper or we get undefined errors for STR_VECT
%include "../types.i"
// Conformer seems to need to come before ROMol
%include "../Conformer.i"
%include "../Dict.i"
%include "../RDLogger.i"
%include "../RDProps.i"
%include "../StereoGroup.i"
%include "../ROMol.i"
%include "../RWMol.i"
%include "../Bond.i"
%include "../BondIterators.i"
%include "../Atom.i"
%include "../AtomIterators.i"
%include "../AtomPairs.i"
%include "../Canon.i"
%include "../Conformer.i"
%include "../QueryAtom.i"
%include "../QueryBond.i"
%include "../QueryOps.i"
%include "../MolBundle.i"
%include "../MonomerInfo.i"
%include "../PeriodicTable.i"
%include "../SanitException.i"
%include "../SmilesParse.i"
%include "../SmilesWrite.i"
%include "../SmartsWrite.i"
%include "../MolOps.i"
%include "../MolSupplier.i"
%include "../MolWriters.i"
%include "../RingInfo.i"
%include "../ChemReactions.i"
%include "../BitOps.i"
%include "../ExplicitBitVect.i"
%include "../Fingerprints.i"
%include "../MorganFingerprints.i"
%include "../ReactionFingerprints.i"
%include "../Rings.i"
%include "../transforms.i"
%include "../DistGeom.i"
%include "../ForceField.i"
%include "../ChemTransforms.i"
%include "../Subgraphs.i"
%include "../MolTransforms.i"
%include "../FMCS.i"
%include "../MolDraw2D.i"
%include "../FilterCatalog.i"
%include "../Trajectory.i"
%include "../MolStandardize.i"
%include "../SubstructLibrary.i"
%include "../RGroupDecomposition.i"
%include "../ScaffoldNetwork.i"
%include "../TautomerQuery.i"
%include "../SubstanceGroup.i"
%include "../MolEnumerator.i"
%include "../MolHash.i"
%include "../Abbreviations.i"
%include "../Streams.i"
%include "../GeneralizedSubstruct.i"
%include "../RascalMCES.i"
%include "../Queries.i"

// Create a class to throw various sorts of errors for testing.  Required for unit tests in ErrorHandlingTests.java
#ifdef INCLUDE_ERROR_GENERATOR
%include "../ErrorGenerator.i"
#endif

/* Done explicitly in BitOps.i
%include <DataStructs/BitOps.h>
%template(TanimotoSimilarityEBV) TanimotoSimilarity<ExplicitBitVect,ExplicitBitVect>;
%template(DiceSimilarityEBV) DiceSimilarity<ExplicitBitVect,ExplicitBitVect>;
*/
%template(DiceSimilarity) RDKit::DiceSimilarity<boost::uint32_t>;

/* vector */
%template(Int_Vect) std::vector<int>;
%template(Byte_Vect) std::vector<signed char>;
%template(Double_Vect) std::vector<double>;
%template(UInt_Vect) std::vector<boost::uint32_t>;
%template(Str_Vect)  std::vector<std::string> ;
%template(Point_Vect) std::vector<RDGeom::Point *>;
%template(Point2D_Vect) std::vector<RDGeom::Point2D *>;
%template(Point3D_Vect) std::vector<RDGeom::Point3D *>;
%template(Atomic_Params_Vect) std::vector<const ForceFields::UFF::AtomicParams *>;

/* pair */
%template(Int_Pair) std::pair<boost::int32_t, int >;
%template(Double_Pair) std::pair<double,double>;
%template(UInt_Pair) std::pair<boost::uint32_t, int >;
%template(Long_Pair) std::pair<boost::int64_t,int>;

/* map */
%template(String_String_Map) std::map<std::string,std::string>;
%template(Int_Int_Map) std::map<int,int>;
%template(Int_Point2D_Map) std::map<int, RDGeom::Point2D>;
%template(Int_Point3D_Map) std::map<int, RDGeom::Point3D>;
%template(Int_Int_Vect_List_Map) std::map<int,std::list<std::vector<int> > >;

/* vector pair */
%template(UInt_Pair_Vect) std::vector<std::pair<boost::uint32_t,int> >;
%template(Match_Vect) std::vector<std::pair<boost::int32_t,int> >;
%template(Long_Pair_Vect) std::vector<std::pair<boost::int64_t,int> >;

/* vector vector */
%template(Int_Vect_Vect) std::vector<std::vector<int> >;

/* list */
%template(Int_Vect_List) std::list<std::vector<int> >;
%template(Int_List) std::list<int>;
%template(UInt_List) std::list<unsigned int>;


/* other */
%template(Match_Vect_Vect) std::vector<std::vector<std::pair<int,int> > >;
%template(Flagged_Atomic_Params_Vect) std::pair<std::vector<const ForceFields::UFF::AtomicParams *>,bool>;
%template(Shared_Int_Array) boost::shared_array<int>;
%template(Shared_Double_Array) boost::shared_array<double>;

// Methods to get at elements of shared arrays
%extend boost::shared_array<double> {
  double getElement(int i) {
    return (*($self))[i];
  }
  void setElement(int i, double value) {
    (*($self))[i] = value;
  }
}
%extend boost::shared_array<int> {
  int getElement(int i) {
    return (*($self))[i];
  }
  void setElement(int i, int value) {
    (*($self))[i] = value;
  }
}

%include "../Descriptors.i"

#ifdef RDK_BUILD_AVALON_SUPPORT
%include "../AvalonLib.i"
#endif
#ifdef RDK_BUILD_INCHI_SUPPORT
%include "../Inchi.i"
#endif

%include "../DiversityPick.i"

%{
#include <RDGeneral/versions.h>
%}

%immutable RDKit::rdkitVersion;
%immutable RDKit::boostVersion;
%immutable RDKit::rdkitBuild;

%include <RDGeneral/versions.h>
