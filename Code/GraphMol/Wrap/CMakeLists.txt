rdkit_python_extension(rdchem 
                       rdchem.cpp Table.cpp Atom.cpp Bond.cpp Mol.cpp
                       Conformer.cpp RingInfo.cpp EditableMol.cpp MonomerInfo.cpp
                       DEST Chem 
                       LINK_LIBRARIES 
                       SmilesParse ChemTransforms SubstructMatch GraphMol RDGeometryLib RDGeneral RDBoost)

rdkit_python_extension(rdmolops
                       rdmolops.cpp MolOps.cpp 
                       DEST Chem 
                       LINK_LIBRARIES 
                       FileParsers SubstructMatch Fingerprints ChemTransforms
                       Catalogs Subgraphs SmilesParse GraphMol RDGeometryLib DataStructs RDGeneral RDBoost  )

rdkit_python_extension(rdqueries
                       rdqueries.cpp Queries.cpp
                       DEST Chem 
                       LINK_LIBRARIES 
                       GraphMol RDGeometryLib DataStructs RDGeneral RDBoost  )

if(RDK_BUILD_COMPRESSED_SUPPLIERS)
find_package(Boost 1.39.0 COMPONENTS iostreams REQUIRED)
set(rdmolfiles_sources rdmolfiles.cpp 
                       ForwardSDMolSupplier.cpp
                       SDMolSupplier.cpp TDTMolSupplier.cpp
                       SmilesMolSupplier.cpp SmilesWriter.cpp SDWriter.cpp
                       TDTWriter.cpp CompressedSDMolSupplier.cpp
                       PDBWriter.cpp )
else(RDK_BUILD_COMPRESSED_SUPPLIERS)
set(rdmolfiles_sources rdmolfiles.cpp
                       ForwardSDMolSupplier.cpp
                       SDMolSupplier.cpp TDTMolSupplier.cpp
                       SmilesMolSupplier.cpp SmilesWriter.cpp SDWriter.cpp
                       TDTWriter.cpp
                       PDBWriter.cpp )
endif(RDK_BUILD_COMPRESSED_SUPPLIERS)

rdkit_python_extension(rdmolfiles
                       ${rdmolfiles_sources}
                       DEST Chem 
                       LINK_LIBRARIES SmilesParse FileParsers GraphMol
                                      RDGeometryLib RDGeneral RDBoost)
if(RDK_BUILD_COMPRESSED_SUPPLIERS)
  set_target_properties(rdmolfiles PROPERTIES DEFINE_SYMBOL SUPPORT_COMPRESSED_SUPPLIERS )
endif(RDK_BUILD_COMPRESSED_SUPPLIERS)

add_pytest(pyGraphMolWrap
         ${CMAKE_CURRENT_SOURCE_DIR}/rough_test.py)

add_pytest(pyTestConformerWrap
         ${CMAKE_CURRENT_SOURCE_DIR}/testConformer.py)
