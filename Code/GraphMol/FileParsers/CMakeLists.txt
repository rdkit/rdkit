rdkit_library(FileParsers
              Mol2FileParser.cpp  
              MolFileParser.cpp MolFileStereochem.cpp MolFileWriter.cpp 
              ForwardSDMolSupplier.cpp SDMolSupplier.cpp SDWriter.cpp
              SmilesMolSupplier.cpp SmilesWriter.cpp
              TDTMolSupplier.cpp TDTWriter.cpp
              TplFileParser.cpp TplFileWriter.cpp
              PDBParser.cpp PDBWriter.cpp PDBSupplier.cpp
              ProximityBonds.cpp
              SequenceParsers.cpp SequenceWriters.cpp
              LINK_LIBRARIES SmilesParse GraphMol)
              
rdkit_headers(FileParsers.h
              FileParserUtils.h
              MolFileStereochem.h
              MolSupplier.h
              MolWriters.h
              SequenceParsers.h SequenceWriters.h
              DEST GraphMol/FileParsers)

rdkit_test(fileParsersTest1 test1.cpp 
           LINK_LIBRARIES FileParsers SmilesParse Depictor SubstructMatch GraphMol RDGeneral RDGeometryLib )

if(RDK_BUILD_TEST_GZIP)
    find_package(Boost 1.39.0 COMPONENTS iostreams REQUIRED)
    set (link_iostreams ${Boost_LIBRARIES})
    add_definitions("-DTEST_GZIP_SD")
    if (NOT Boost_USE_STATIC_LIBS)
      add_definitions("-DBOOST_IOSTREAMS_DYN_LINK")
    endif()
else()
    set (link_iostreams "")
endif()

rdkit_test(testMolSupplier testMolSupplier.cpp 
           LINK_LIBRARIES FileParsers SmilesParse Depictor SubstructMatch GraphMol RDGeneral RDGeometryLib ${link_iostreams})

rdkit_test(testMolWriter testMolWriter.cpp LINK_LIBRARIES FileParsers SmilesParse GraphMol RDGeneral RDGeometryLib )

rdkit_test(testTplParser testTpls.cpp LINK_LIBRARIES FileParsers SmilesParse GraphMol RDGeneral RDGeometryLib )

rdkit_test(testMol2ToMol testMol2ToMol.cpp LINK_LIBRARIES FileParsers SmilesParse GraphMol RDGeneral RDGeometryLib )

rdkit_test(testSequence testSequence.cpp LINK_LIBRARIES FileParsers SmilesParse GraphMol RDGeneral RDGeometryLib )

