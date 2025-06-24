if(NOT RDK_BUILD_FREESASA_SUPPORT)
  return()
endif(NOT RDK_BUILD_FREESASA_SUPPORT)

if(NOT DEFINED FREESASA_DIR)
  set(FREESASA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/freesasa")
  set(needDownload "TRUE")
  if(EXISTS "${FREESASA_DIR}/src/freesasa.h")
    set(needDownload "FALSE")
  endif()
endif()

if(needDownload)
  # don't actually use the md5 here
  set(FREESASA_VERSION "2.0.3")
  set(FREESASA_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/freesasa-${FREESASA_VERSION}")
  if(NOT EXISTS "${FREESASA_SRC_DIR}/src")
    downloadAndCheckMD5("https://github.com/mittinatten/freesasa/releases/download/${FREESASA_VERSION}/freesasa-${FREESASA_VERSION}.tar.gz"
          "${CMAKE_CURRENT_SOURCE_DIR}/master.tar.gz"
          "")
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar zxf
      ${CMAKE_CURRENT_SOURCE_DIR}/master.tar.gz
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    message("Copying freesasa2.c to ${FREESASA_SRC_DIR}/src/freesasa2.c")
    file(COPY freesasa2.c DESTINATION ${FREESASA_SRC_DIR}/src/)
  else()
    message("FreeSASA found in ${FREESASA_SRC_DIR}, no download required")
  endif()
endif()

set(CMAKE_C_STANDARD 99)

set (freesasa_clib_srcs
    ${FREESASA_SRC_DIR}/src/classifier.c
    ${FREESASA_SRC_DIR}/src/classifier_naccess.c
    ${FREESASA_SRC_DIR}/src/classifier_oons.c
    ${FREESASA_SRC_DIR}/src/classifier_protor.c
    ${FREESASA_SRC_DIR}/src/coord.c
    ${FREESASA_SRC_DIR}/src/freesasa2.c
    ${FREESASA_SRC_DIR}/src/lexer.c
    ${FREESASA_SRC_DIR}/src/log.c
    ${FREESASA_SRC_DIR}/src/nb.c
    ${FREESASA_SRC_DIR}/src/node.c
    ${FREESASA_SRC_DIR}/src/parser.c
    ${FREESASA_SRC_DIR}/src/pdb.c
    ${FREESASA_SRC_DIR}/src/rsa.c
    ${FREESASA_SRC_DIR}/src/sasa_lr.c
    ${FREESASA_SRC_DIR}/src/sasa_sr.c
    ${FREESASA_SRC_DIR}/src/selection.c
    ${FREESASA_SRC_DIR}/src/structure.c
    ${FREESASA_SRC_DIR}/src/util.c
    )
if((MSVC AND RDK_INSTALL_DLLS_MSVC) OR ((NOT MSVC) AND WIN32))
  set(freesasa_additional_exports ${CMAKE_CURRENT_SOURCE_DIR}/additional_exports.def)
  file(WRITE ${freesasa_additional_exports}
    "EXPORTS\n"
    "freesasa_calc_structure\n"
    "freesasa_calc_coord\n"
    "freesasa_calc_tree\n"
    "freesasa_result_classes\n"
    "freesasa_result_free\n"
    "freesasa_classifier_from_file\n"
    "freesasa_classifier_free\n"
    "freesasa_classifier_radius\n"
    "freesasa_classifier_class\n"
    "freesasa_classifier_class2str\n"
    "freesasa_classifier_name\n"
    "freesasa_selection_new\n"
    "freesasa_selection_free\n"
    "freesasa_selection_name\n"
    "freesasa_selection_command\n"
    "freesasa_selection_area\n"
    "freesasa_set_verbosity\n"
    "freesasa_get_verbosity\n"
    "freesasa_set_err_out\n"
    "freesasa_get_err_out\n"
    "freesasa_structure_new\n"
    "freesasa_structure_free\n"
    "freesasa_structure_from_pdb\n"
    "freesasa_structure_array\n"
    "freesasa_structure_add_atom\n"
    "freesasa_structure_add_atom_wopt\n"
    "freesasa_structure_get_chains\n"
    "freesasa_structure_chain_labels\n"
    "freesasa_structure_n\n"
    "freesasa_structure_n_residues\n"
    "freesasa_structure_n_chains\n"
    "freesasa_structure_radius\n"
    "freesasa_structure_set_radius\n"
    "freesasa_structure_atom_name\n"
    "freesasa_structure_atom_res_name\n"
    "freesasa_structure_atom_res_number\n"
    "freesasa_structure_atom_chain\n"
    "freesasa_structure_atom_symbol\n"
    "freesasa_structure_atom_radius\n"
    "freesasa_structure_atom_set_radius\n"
    "freesasa_structure_residue_name\n"
    "freesasa_structure_residue_number\n"
    "freesasa_structure_residue_chain\n"
    "freesasa_structure_model\n"
    "freesasa_structure_coord_array\n"
    "freesasa_structure_residue_atoms\n"
    "freesasa_structure_chain_atoms\n"
    "freesasa_structure_chain_residues\n"
    "freesasa_structure_classifier_name\n"
    "freesasa_tree_new\n"
    "freesasa_tree_init\n"
    "freesasa_tree_add_result\n"
    "freesasa_tree_join\n"
    "freesasa_tree_export\n"
    "freesasa_node_free\n"
    "freesasa_node_area\n"
    "freesasa_node_children\n"
    "freesasa_node_next\n"
    "freesasa_node_parent\n"
    "freesasa_node_type\n"
    "freesasa_node_name\n"
    "freesasa_node_classified_by\n"
    "freesasa_node_atom_is_polar\n"
    "freesasa_node_atom_is_mainchain\n"
    "freesasa_node_atom_radius\n"
    "freesasa_node_atom_pdb_line\n"
    "freesasa_node_residue_number\n"
    "freesasa_node_residue_n_atoms\n"
    "freesasa_node_residue_reference\n"
    "freesasa_node_chain_n_residues\n"
    "freesasa_node_structure_n_chains\n"
    "freesasa_node_structure_n_atoms\n"
    "freesasa_node_structure_chain_labels\n"
    "freesasa_node_structure_model\n"
    "freesasa_node_structure_result\n"
    "freesasa_node_structure_selections\n"
    "freesasa_node_structure_add_selection\n"
    "freesasa_node_result_parameters\n"
    "freesasa_select_area\n"
    "freesasa_default_parameters DATA\n"
    "freesasa_protor_classifier DATA\n"
    "freesasa_naccess_classifier DATA\n"
    "freesasa_oons_classifier DATA\n"
  )
  set (freesasa_clib_srcs
    ${freesasa_clib_srcs}
    ${freesasa_additional_exports}
  )
endif()

set (freesasa_h ${FREESASA_SRC_DIR}/src/freesasa.h)
file(READ ${freesasa_h} freesasa_h_data)
if (NOT freesasa_h_data MATCHES "RDGeneral/export")
  string(REGEX REPLACE "(#include <stdio.h>)" "\\1\n#include <RDGeneral/export.h>" freesasa_h_data "${freesasa_h_data}")
endif()
if (MSVC)
  set (RDKIT_FREESASA_CLIB_EXPORT "RDKIT_FREESASA_CLIB_EXPORT ")
else()
  set (RDKIT_FREESASA_CLIB_EXPORT "")
endif()
string(REGEX REPLACE "(RDKIT_FREESASA_CLIB_EXPORT )?(extern const)" "${RDKIT_FREESASA_CLIB_EXPORT}\\2" freesasa_h_data "${freesasa_h_data}")
file(WRITE ${freesasa_h} "${freesasa_h_data}")

add_definitions(-DUSE_THREADS=0)
add_definitions(-DUSE_JSON=0)
add_definitions(-DUSE_XML=0)
rdkit_library(freesasa_clib ${freesasa_clib_srcs})
target_compile_definitions(freesasa_clib PRIVATE RDKIT_FREESASALIB_BUILD)

if((MSVC AND RDK_INSTALL_DLLS_MSVC) OR ((NOT MSVC) AND WIN32))
  set_target_properties(freesasa_clib PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
endif()
include_directories("${FREESASA_SRC_DIR}/src")

rdkit_library(FreeSASALib RDFreeSASA.cpp SHARED 
    LINK_LIBRARIES freesasa_clib GraphMol )
target_compile_definitions(FreeSASALib PRIVATE RDKIT_FREESASALIB_BUILD)

rdkit_headers(RDFreeSASA.h DEST GraphMol)

rdkit_test(testFreeSASALib testFreeSASA.cpp
           LINK_LIBRARIES FreeSASALib FileParsers SmilesParse )

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
