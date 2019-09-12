/*
 * This file is part of the RingDecomposerLib, licensed
 * under BSD New license (see LICENSE in the root directory).
 * Copyright (c) 2016
 * University of Hamburg, ZBH - Center for Bioinformatics
 * Niek Andresen, Florian Flachsenberg, Matthias Rarey
 * 
 * Please cite:
 * 
 * Kolodzik, A.; Urbaczek, S.; Rarey, M.
 * Unique Ring Families: A Chemically Meaningful Description
 * of Molecular Ring Topologies.
 * J. Chem. Inf. Model., 2012, 52 (8), pp 2013-2021
 * 
 * Flachsenberg, F.; Andresen, N.; Rarey, M.
 * RingDecomposerLib: An Open-Source Implementation of
 * Unique Ring Families and Other Cycle Bases.
 * J. Chem. Inf. Model., 2017, 57 (2), pp 122-126
 */

/**
 * @file
 *
 * @brief Runtime bechmark tool for the UniqueRingFamily library and RDKit
 *
 * Start the tool with
 *
 *     RDKitURF -i <filename> -r <rumber of repeats> [--rc]
 *
 * The input file must be in SMILES, SDF or DIMCAS format.
 * Number of repeats specifies the number of times the benchmark
 * is executed. The last parameter is an optional flag.
 * It enables benchmarking of RC calculation (as long as the
 * number is less than 4e6.
 *
 * It outputs some molecule characteristics (number of
 * bonds, atoms, relevant cycles and URFs) as well
 * as the mean runtime for a number of experiments.
 */

#include <algorithm>
#include <chrono>
#include <fstream>
#include <functional>
#include <iostream>

#include <boost/program_options.hpp>

#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/BondIterators.h>


namespace RingDecomposerLib {
#include "RingDecomposerLib.h"
#include "RDLdimacs.h"
}

namespace options = boost::program_options;

options::variables_map parseCommandLine(int argc, char **argv) {
  options::options_description desc("Allowed options");
  desc.add_options()
      ("input,i", options::value<std::string>(),
          "Input file with molecules.")
      ("repeats,r", options::value<unsigned>(),
          "Number of benchmark repeats")
      ("rc", "calculate relevant cycles")
      ("debug", "Debug output")
      ("help", "show this help");

  options::variables_map parsedMap;
  options::store(options::command_line_parser(argc, argv).
                 options(desc).run(), parsedMap);
  options::notify(parsedMap);

  if (parsedMap.count("help")) {
    std::cout << desc << "\n";
    exit(0);
  }

  if (!parsedMap.count("input")) {
    std::cout << desc << "\n";
    exit(1);
  }

  if (!parsedMap.count("repeats")) {
    std::cout << desc << "\n";
    exit(1);
  }

  return parsedMap;
}


template <typename T>
struct graph_trait {};

// trait for converting RDKit molecule
template<>
struct graph_trait<const RDKit::ROMol> {
    static std::string get_name(const RDKit::ROMol& mol)
    {
      std::string name = "NO_NAME";

      if (mol.hasProp("_Name")) {
        mol.getProp("_Name", name);
      }

      return name;
    }

    static size_t get_nof_atoms(const RDKit::ROMol& mol)
    {
      return mol.getNumAtoms();
    }

    static size_t get_nof_bonds(const RDKit::ROMol& mol)
    {
      return mol.getNumBonds();
    }

    static RDKit::ConstBondIterator_ get_bonds_begin(const RDKit::ROMol& mol)
    {
      return mol.beginBonds();
    }

    static RDKit::ConstBondIterator_ get_bonds_end(const RDKit::ROMol& mol)
    {
      return mol.endBonds();
    }

    static unsigned get_first_id(const RDKit::Bond* b) {
      return b->getBeginAtomIdx();
    }

    static unsigned get_second_id(const RDKit::Bond* b) {
      return b->getEndAtomIdx();
    }
};

// trait for converting dimacs graph
template<>
struct graph_trait<const RingDecomposerLib::RDL_DimacsGraph> {
    static std::string get_name(const RingDecomposerLib::RDL_DimacsGraph& mol)
    {
      std::string name = "NO_NAME";

      if (mol.name) {
        name = mol.name;
      }

      return name;
    }

    static size_t get_nof_atoms(const RingDecomposerLib::RDL_DimacsGraph& mol)
    {
      return mol.nof_nodes;
    }

    static size_t get_nof_bonds(const RingDecomposerLib::RDL_DimacsGraph& mol)
    {
      return mol.nof_edges;
    }

    static RingDecomposerLib::RDL_edge* get_bonds_begin(const RingDecomposerLib::RDL_DimacsGraph& mol)
    {
      return mol.edges;
    }

    static RingDecomposerLib::RDL_edge* get_bonds_end(const RingDecomposerLib::RDL_DimacsGraph& mol)
    {
      return mol.edges + mol.nof_edges;
    }

    static unsigned get_first_id(const RingDecomposerLib::RDL_edge& b) {
      return b[0] - 1;
    }

    static unsigned get_second_id(const RingDecomposerLib::RDL_edge& b) {
      return b[1] - 1;
    }
};

// header for datapoint table
const std::string datapoint_header =
    "datapoint_header;"
    "name;"
    "calculation_time;"
    "nodes_edges_time;"
    "nof_rc_time;"
    "rcp_time;"
    "basis_time;"
    "rc_time";

// header for molecule characterization table
const std::string characterization_header =
    "characterization_header;"
    "name;"
    "nofAtoms;"
    "nofBonds;"
    "nofAtomsLargestRS;"
    "nofBondsLargestRS;"
    "nofUrf;"
    "nofRc";

// separator for results
const std::string separator = "####";

template <typename T>
bool analyze_molecule(
    T& mol,
    unsigned repeats,
    bool rc)
{
  std::string name = graph_trait<T>::get_name(mol);

  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;

  unsigned nof_urf = 0;
  double nof_rc = 0.0;

  std::cout << datapoint_header << std::endl;

  start = std::chrono::high_resolution_clock::now();
  // start timing construction & calculation
  for (unsigned r = 0; r < repeats; ++r) {
    RingDecomposerLib::RDL_graph* graph =
        RingDecomposerLib::RDL_initNewGraph(graph_trait<T>::get_nof_atoms(mol));

    for (auto it = graph_trait<T>::get_bonds_begin(mol);
          it != graph_trait<T>::get_bonds_end(mol); ++it) {
      RingDecomposerLib::RDL_addUEdge(graph, graph_trait<T>::get_first_id(*it),
          graph_trait<T>::get_second_id(*it));
    }
    RingDecomposerLib::RDL_data* data = RingDecomposerLib::RDL_calculate(graph);

    RingDecomposerLib::RDL_deleteData(data);
  }
  // end timing construction & calcuation
  end = std::chrono::high_resolution_clock::now();

  unsigned long long calculation_time =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
  double calculation_time_mean = static_cast<double>(calculation_time) / repeats;

  // construct new graph & calculate (for later calculations)
  RingDecomposerLib::RDL_graph* graph =
      RingDecomposerLib::RDL_initNewGraph(graph_trait<T>::get_nof_atoms(mol));

  for (auto it = graph_trait<T>::get_bonds_begin(mol);
        it != graph_trait<T>::get_bonds_end(mol); ++it) {
    unsigned add_result =
        RingDecomposerLib::RDL_addUEdge(graph, graph_trait<T>::get_first_id(*it),
        graph_trait<T>::get_second_id(*it));
    if (add_result == RingDecomposerLib::RDL_DUPLICATE_EDGE) {
      std::cerr << graph_trait<T>::get_first_id(*it) << " "
          << graph_trait<T>::get_second_id(*it) << std::endl;
      std::cerr << "WARNING: duplicate edge for " << name << std::endl;
    }
  }
  RingDecomposerLib::RDL_data* data = RingDecomposerLib::RDL_calculate(graph);
  if (!data) {
    std::cerr << "FATAL: calculation failed for " << name << std::endl;
    return false;
  }

  if (RingDecomposerLib::RDL_getNofURF(data) == RingDecomposerLib::RDL_INVALID_RESULT) {
    std::cerr << "FATAL: invalid result for " << name << std::endl;
    return false;
  }

  // start timing node and edge query
  start = std::chrono::high_resolution_clock::now();
  for (unsigned r = 0; r < repeats; ++r) {
    RingDecomposerLib::RDL_edge* edges;
    RingDecomposerLib::RDL_node* nodes;
    nof_urf = RingDecomposerLib::RDL_getNofURF(data);

    for (unsigned urfidx = 0; urfidx < nof_urf; ++urfidx) {
      RingDecomposerLib::RDL_getEdgesForURF(data, urfidx, &edges);
      RingDecomposerLib::RDL_getNodesForURF(data, urfidx, &nodes);
      free(edges);
      free(nodes);
    }
  }
  // end timing node and edge query
  end = std::chrono::high_resolution_clock::now();

  unsigned long long nodes_edges_time =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
  double nodes_edges_time_mean = static_cast<double>(nodes_edges_time) / repeats;

  if (RingDecomposerLib::RDL_getNofRC(data) == RingDecomposerLib::RDL_INVALID_RESULT) {
    std::cerr << "FATAL: invalid result for " << name << std::endl;
    return false;
  }

  // start timing calculate nofRC
  start = std::chrono::high_resolution_clock::now();
  for (unsigned r = 0; r < repeats; ++r) {
    nof_rc = RingDecomposerLib::RDL_getNofRC(data);
  }
  // end timing calculate nofRC
  end = std::chrono::high_resolution_clock::now();

  unsigned long long nof_rc_time =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
  double nof_rc_time_mean = static_cast<double>(nof_rc_time) / repeats;


  RingDecomposerLib::RDL_cycle** rcp_cycles;
  unsigned nof_rcp = RingDecomposerLib::RDL_getRCPrototypes(data, &rcp_cycles);
  if (nof_rcp == RingDecomposerLib::RDL_INVALID_RESULT) {
    std::cerr << "FATAL: invalid result for " << name << std::endl;
    return false;
  }
  RingDecomposerLib::RDL_deleteCycles(rcp_cycles, nof_rcp);

  // start timing get RCPs
  start = std::chrono::high_resolution_clock::now();
  for (unsigned r = 0; r < repeats; ++r) {
    RingDecomposerLib::RDL_cycle** cycles;
    unsigned nof_rcp = RingDecomposerLib::RDL_getRCPrototypes(data, &cycles);
    RingDecomposerLib::RDL_deleteCycles(cycles, nof_rcp);
  }
  // end timing get RCPs
  end = std::chrono::high_resolution_clock::now();

  unsigned long long rcp_time =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
  double rcp_time_mean = static_cast<double>(rcp_time) / repeats;

  RingDecomposerLib::RDL_cycle** basis_cycles;
  unsigned nof_basis = RingDecomposerLib::RDL_getSSSR(data, &basis_cycles);
  if (nof_basis == RingDecomposerLib::RDL_INVALID_RESULT) {
    std::cerr << "FATAL: invalid result for " << name << std::endl;
    return false;
  }
  RingDecomposerLib::RDL_deleteCycles(basis_cycles, nof_basis);

  // start timing get sssr
  start = std::chrono::high_resolution_clock::now();
  for (unsigned r = 0; r < repeats; ++r) {
    RingDecomposerLib::RDL_cycle** cycles;
    unsigned nof_basis = RingDecomposerLib::RDL_getSSSR(data, &cycles);
    RingDecomposerLib::RDL_deleteCycles(cycles, nof_basis);
  }
  // end timing get sssr
  end = std::chrono::high_resolution_clock::now();

  unsigned long long basis_time =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
  double basis_time_mean = static_cast<double>(basis_time) / repeats;


  unsigned long long rc_time = 0ul;
  double rc_time_mean;

  if (rc) {
    unsigned long long rc_count = 0ul;
    RingDecomposerLib::RDL_cycleIterator *it = RingDecomposerLib::RDL_getRCyclesIterator(data);

    if (it == NULL) {
      std::cerr << "FATAL: invalid result for " << name << std::endl;
      return false;
    }

    while (!RingDecomposerLib::RDL_cycleIteratorAtEnd(it)) {
      RingDecomposerLib::RDL_cycle* cycle = RingDecomposerLib::RDL_cycleIteratorGetCycle(it);
      RingDecomposerLib::RDL_deleteCycle(cycle);
      RingDecomposerLib::RDL_cycleIteratorNext(it);
      ++rc_count;
    }
    RingDecomposerLib::RDL_deleteCycleIterator(it);

    if (rc_count != nof_rc) {
      std::cerr << "FATAL: invalid result for " << name << std::endl;
      return false;
    }

    // start timing get RCs
    start = std::chrono::high_resolution_clock::now();
    for (unsigned r = 0; r < repeats; ++r) {
      RingDecomposerLib::RDL_cycleIterator *it = RingDecomposerLib::RDL_getRCyclesIterator(data);
      while (!RingDecomposerLib::RDL_cycleIteratorAtEnd(it)) {
        RingDecomposerLib::RDL_cycle* cycle = RingDecomposerLib::RDL_cycleIteratorGetCycle(it);
        RingDecomposerLib::RDL_deleteCycle(cycle);
        RingDecomposerLib::RDL_cycleIteratorNext(it);
      }
      RingDecomposerLib::RDL_deleteCycleIterator(it);
    }
    // end timing get RCs
    end = std::chrono::high_resolution_clock::now();

    rc_time =
        std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
  }
  rc_time_mean = static_cast<double>(rc_time) / repeats;

  std::cout << "datapoint;"
      << name << ";"
      << calculation_time_mean << ";"
      << nodes_edges_time_mean << ";"
      << nof_rc_time_mean << ";"
      << rcp_time_mean << ";"
      << basis_time_mean << ";"
      << rc_time_mean << std::endl;

  // find the BCC with largest z
  unsigned nof_atoms_largest_rs_by_z = 0;
  unsigned nof_bonds_largest_rs_by_z = 0;
  unsigned largest_z = 0;

  if (RingDecomposerLib::RDL_getNofRingsystems(data)
      == RingDecomposerLib::RDL_INVALID_RESULT) {
    std::cerr << "FATAL: invalid result for " << name << std::endl;
    return false;
  }

  for (unsigned i = 0; i < RingDecomposerLib::RDL_getNofRingsystems(data); ++i) {
    unsigned nof_atoms_rs =
        RingDecomposerLib::RDL_getNofNodesForRingsystem(data, i);
    unsigned nof_bonds_rs =
        RingDecomposerLib::RDL_getNofEdgesForRingsystem(data, i);
    unsigned z = nof_bonds_rs - nof_atoms_rs + 1;
    if (z > largest_z) {
      largest_z = z;
      nof_atoms_largest_rs_by_z = nof_atoms_rs;
      nof_bonds_largest_rs_by_z = nof_bonds_rs;
    }
    else if (z == largest_z && nof_bonds_rs > nof_bonds_largest_rs_by_z) {
      nof_atoms_largest_rs_by_z = nof_atoms_rs;
      nof_bonds_largest_rs_by_z = nof_bonds_rs;
    }
  }

  RingDecomposerLib::RDL_deleteData(data);
  std::cout << separator << std::endl;

  std::cout << characterization_header << std::endl;

  std::cout << "characterization;"
      << name << ";"
      << graph_trait<T>::get_nof_atoms(mol) << ";"
      << graph_trait<T>::get_nof_bonds(mol) << ";"
      << nof_atoms_largest_rs_by_z << ";"
      << nof_bonds_largest_rs_by_z << ";"
      << nof_urf << ";"
      << std::setprecision(0) << std::fixed << nof_rc << std::endl;
  std::cout << separator << std::endl;

  return true;
}

unsigned read_and_analyze_dimacsfile(
    const std::string& input,
    unsigned repeats,
    bool rc)
{
  // try to read dimacs graph
  RingDecomposerLib::RDL_DimacsGraph* graph =
      RingDecomposerLib::RDL_dimacsGraph_read(input.c_str());
  if (!graph || !analyze_molecule(
      *const_cast<const RingDecomposerLib::RDL_DimacsGraph*>(graph), repeats, rc)) {
    std::cout << datapoint_header << std::endl;
    std::cout << "datapoint;$ERROR" << std::endl;
    std::cout << separator << std::endl;

    std::cout << characterization_header << std::endl;
    std::cout << "characterization;$ERROR" << std::endl;
    std::cout << separator << std::endl;

    return EXIT_FAILURE;
  }

  RingDecomposerLib::RDL_dimacsGraph_delete(graph);

  return EXIT_SUCCESS;
}

unsigned read_and_analyze_molfile(
    const std::string& input,
    const std::string& extension,
    unsigned repeats,
    bool rc)
{
  // try reading RDKit molecule depending of file suffix
  RDKit::MolSupplier* supplier;
  if (extension == ".smi") {
    supplier = new RDKit::SmilesMolSupplier (input, " \t", 0, 1, false, true);
  }
  else {
    supplier = new RDKit::SDMolSupplier (input, true, true, true);
  }

  while (!supplier->atEnd()) {
    const RDKit::ROMol* mol = supplier->next();

    if (!mol || !analyze_molecule(*mol, repeats, rc)) {
      std::cout << datapoint_header << std::endl;
      std::cout << "datapoint;$ERROR" << std::endl;
      std::cout << separator << std::endl;

      std::cout << characterization_header << std::endl;
      std::cout << "characterization;$ERROR" << std::endl;
      std::cout << separator << std::endl;

      continue;
    }

    delete mol;
  }

  delete supplier;

  return EXIT_SUCCESS;
}

unsigned analyze_molecules(
    const std::string& input,
    unsigned repeats,
    bool rc)
{
  size_t point_pos = input.find_last_of(".");

  if (point_pos == std::string::npos) {
    std::cout << "could not determine file extension of: " << input << std::endl;
    return EXIT_FAILURE;
  }

  std::string extension = input.substr(point_pos);
  std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);

  if (extension == ".smi" || extension == ".sdf") {
    return read_and_analyze_molfile(input, extension, repeats, rc);
  }
  else if (extension == ".dimacs") {
    return read_and_analyze_dimacsfile(input, repeats, rc);
  }
  else {
    std::cout << "invalid file extension detected: " << extension << std::endl;
    return EXIT_FAILURE;
  }
}

int main(int argc, char *argv[])
{

  RingDecomposerLib::RDL_setOutputFunction(RingDecomposerLib::RDL_writeToStderr);

  options::variables_map parsedMap = parseCommandLine(argc, argv);
  std::string input(parsedMap["input"].as<std::string>());
  unsigned repeats = parsedMap["repeats"].as<unsigned>();
  bool calc_rc = parsedMap.count("rc");
  return analyze_molecules(input, repeats, calc_rc);
}
