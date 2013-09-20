
if(RDK_TEST_MMFF_COMPLIANCE)
 set(T_LIBS ${Boost_LIBRARIES})
 find_package(Boost 1.39.0 COMPONENTS iostreams REQUIRED)
 set(Boost_LIBRARIES ${T_LIBS} ${Boost_LIBRARIES})
 if(MSVC)
   ADD_DEFINITIONS("/D YY_NO_UNISTD_H")
   # as of cmake 2.8.0 and boost 1.43 we need to do this to avoid
   # forcing the user to build both static and dynamic versions
   # of the regex library:
   ADD_DEFINITIONS("/D BOOST_ALL_DYN_LINK") 
 endif(MSVC)

 rdkit_test(testMMFFForceField testMMFFForceField.cpp 
	    LINK_LIBRARIES 
	    DistGeomHelpers DistGeometry ForceFieldHelpers ForceField 
	    FileParsers MolTransforms SmilesParse SubstructMatch MolAlign Optimizer EigenSolvers Alignment GraphMol RDGeometryLib RDGeneral
 ${Boost_IOSTREAMS_LIBRARY})
endif(RDK_TEST_MMFF_COMPLIANCE)
#add_subdirectory(Wrap)
