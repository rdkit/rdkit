find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(CMAKE_SWIG_FLAGS "")

if(WIN32)
	set(RDKit_Wrapper_Libs
        FileParsers SmilesParse
	Depictor SubstructMatch ChemReactions Fingerprints ChemTransforms
        Descriptors MolTransforms  PartialCharges 
	Subgraphs GraphMol DataStructs
	DistGeomHelpers DistGeometry ForceFieldHelpers ForceField
	EigenSolvers
	Optimizer
	MolAlign
	Alignment
	RDGeometryLib RDGeneral )
        if(RDK_BUILDING_INCHI_SUPPORT)
          set(RDKit_Wrapper_Libs RDInchiLib Inchi ${RDKit_Wrapper_Libs})
        endif(RDK_BUILDING_INCHI_SUPPORT)
else()
	set(RDKit_Wrapper_Libs 
        FileParsers_static SmilesParse_static
	Depictor_static SubstructMatch_static ChemReactions_static Fingerprints_static ChemTransforms_static
	Subgraphs_static GraphMol_static DataStructs_static Descriptors_static PartialCharges_static MolTransforms_static
	DistGeomHelpers_static DistGeometry_static ForceFieldHelpers_static ForceField_static
	EigenSolvers_static
	Optimizer_static
	MolAlign_static
	Alignment_static
	RDGeometryLib_static RDGeneral_static )
        if(RDK_BUILDING_INCHI_SUPPORT)
          set(RDKit_Wrapper_Libs RDInchiLib_static Inchi_static ${RDKit_Wrapper_Libs})
        endif(RDK_BUILDING_INCHI_SUPPORT)


endif()

add_subdirectory(gmwrapper)
