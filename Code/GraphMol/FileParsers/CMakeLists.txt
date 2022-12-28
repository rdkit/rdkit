

if(RDK_USE_BOOST_IOSTREAMS)
if(WIN32)
  find_package(Boost ${RDK_BOOST_VERSION} COMPONENTS system iostreams REQUIRED)
  set (link_iostreams ${Boost_LIBRARIES})
  if (NOT Boost_USE_STATIC_LIBS)
     add_definitions("-DBOOST_IOSTREAMS_DYN_LINK")
  endif()
  find_package(ZLIB)
  include_directories(${ZLIB_INCLUDE_DIRS})
  find_package(Boost ${RDK_BOOST_VERSION} COMPONENTS zlib)
  if(Boost_zlib_FOUND)
    set(zlib_lib Boost::zlib)
  else()
    set(zlib_lib ${ZLIB_LIBRARIES})
  endif()


else()
  find_package(Boost ${RDK_BOOST_VERSION} COMPONENTS system iostreams REQUIRED)
  set (link_iostreams ${Boost_LIBRARIES})
  if (NOT Boost_USE_STATIC_LIBS)
     add_definitions("-DBOOST_IOSTREAMS_DYN_LINK")
  endif()
  find_package(ZLIB)
  include_directories(${ZLIB_INCLUDE_DIRS})
  set(zlib_lib ${ZLIB_LIBRARIES})
endif()
ADD_DEFINITIONS("-DRDK_USE_BOOST_IOSTREAMS")
endif(RDK_USE_BOOST_IOSTREAMS)

if(RDK_BUILD_MAEPARSER_SUPPORT)
    include_directories(${maeparser_INCLUDE_DIRS})
    set (maesupplier MaeMolSupplier.cpp)
endif()

if(RDK_USE_BOOST_REGEX)
  add_definitions(-DRDKIT_USE_BOOST_REGEX)
  find_package(Boost ${RDK_BOOST_VERSION} COMPONENTS regex REQUIRED)
  set(regex_lib Boost::regex)
endif()

rdkit_library(FileParsers
              CDXMLParser.cpp
              Mol2FileParser.cpp MolFileParser.cpp
              MolSGroupParsing.cpp MolSGroupWriting.cpp
              MolFileStereochem.cpp MolFileWriter.cpp
              ForwardSDMolSupplier.cpp SDMolSupplier.cpp SDWriter.cpp
              SmilesMolSupplier.cpp SmilesWriter.cpp
              TDTMolSupplier.cpp TDTWriter.cpp
              TplFileParser.cpp TplFileWriter.cpp
              PDBParser.cpp PDBWriter.cpp PDBSupplier.cpp
              CMLWriter.cpp XYZFileWriter.cpp XYZFileParser.cpp
              ${maesupplier}
              ProximityBonds.cpp
              SequenceParsers.cpp SequenceWriters.cpp
              SVGParser.cpp PNGParser.cpp
							MultithreadedMolSupplier.cpp
							MultithreadedSmilesMolSupplier.cpp
							MultithreadedSDMolSupplier.cpp
              LINK_LIBRARIES GenericGroups Depictor SmilesParse ChemTransforms GraphMol ${RDK_MAEPARSER_LIBS} ${regex_lib} ${link_iostreams} ${zlib_lib})
target_compile_definitions(FileParsers PRIVATE RDKIT_FILEPARSERS_BUILD)

rdkit_headers(CDXMLParser.h
              FileParsers.h
              FileParserUtils.h
              MolFileStereochem.h
              MolSupplier.h
              MolWriters.h
              SequenceParsers.h SequenceWriters.h
              GeneralFileReader.h
							MultithreadedMolSupplier.h
							MultithreadedSmilesMolSupplier.h
							MultithreadedSDMolSupplier.h
              PNGParser.h
              DEST GraphMol/FileParsers)

rdkit_test(fileParsersTest1 test1.cpp
           LINK_LIBRARIES FileParsers SubstructMatch )



rdkit_test(testMolSupplier testMolSupplier.cpp
           LINK_LIBRARIES  FileParsers RDStreams)

rdkit_test(testGeneralFileReader testGeneralFileReader.cpp
					LINK_LIBRARIES FileParsers RDStreams)

if(RDK_TEST_MULTITHREADED)
rdkit_test(testMultithreadedMolSupplier testMultithreadedMolSupplier.cpp
					LINK_LIBRARIES FileParsers Fingerprints RDStreams)
endif(RDK_TEST_MULTITHREADED)

rdkit_test(testMolWriter testMolWriter.cpp LINK_LIBRARIES FileParsers )

rdkit_test(testTplParser testTpls.cpp LINK_LIBRARIES FileParsers )

rdkit_test(testMol2ToMol testMol2ToMol.cpp LINK_LIBRARIES FileParsers )

rdkit_test(testSequence testSequence.cpp LINK_LIBRARIES FileParsers )

rdkit_test(testExtendedStereoParsing testExtendedStereoParsing.cpp
           LINK_LIBRARIES FileParsers )

rdkit_catch_test(fileParsersCatchTest file_parsers_catch.cpp 
           LINK_LIBRARIES FileParsers )

rdkit_catch_test(testPropertyLists testPropertyLists.cpp
           LINK_LIBRARIES FileParsers )

rdkit_catch_test(XYZFileParserCatchTest XYZFileParserCatchTest.cpp
           LINK_LIBRARIES FileParsers )

rdkit_catch_test(cdxmlParserCatchTest cdxml_parser_catch.cpp
           LINK_LIBRARIES FileParsers )

rdkit_catch_test(molfileStereoCatchTest molfile_stereo_catch.cpp 
           LINK_LIBRARIES FileParsers Subgraphs)

