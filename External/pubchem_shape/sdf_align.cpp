#include <iostream>
#include <sstream>
#include <stdexcept>

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Fingerprints/FingerprintUtil.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <RDGeneral/BoostEndInclude.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>

#include "PubChemShape.hpp"

using namespace std;

#define ERRORTHROW(msg_stream)     \
  do {                             \
    ostringstream os;              \
    os << msg_stream;              \
    throw runtime_error(os.str()); \
  } while (0)

// #define DEBUG_MSG(msg_stream) cout << msg_stream << '\n'
#define DEBUG_MSG(msg_stream)

using namespace RDKit;

int main(int argc, char **argv) {
  int status = -1;

  try {
    if (argc != 2)
      ERRORTHROW(
          "Usage: sdf_align <input.sdf>   The first molecule in the SDF "
          "is the reference");  //"opt_param
                                // max_preiters
                                // max_postiters");
    bool useColors = true;

    bool sanitize = true;
    bool removeHs = false;
    SDMolSupplier suppl(argv[1], sanitize, removeHs);
    SDWriter writer("sdf_align.out.sdf");
    unique_ptr<ROMol> ref{suppl[0]};
    if (!ref.get()) {
      ERRORTHROW("Failed to read ref conformer");
    }

    for (auto i = 1u; i < suppl.length(); ++i) {
      std::unique_ptr<ROMol> fit{suppl[i]};
      if (!fit) {
        continue;
      }
      vector<float> matrix(12, 0.0);
      int refConfId = -1;
      int fitConfId = -1;
      auto [nbr_st, nbr_ct] =
          AlignMolecule(*ref, *fit, matrix, refConfId, fitConfId, useColors);
      if (i == 1) {
        writer.write(*ref, refConfId);
      }
      writer.write(*fit, fitConfId);
    }

    status = 0;

  } catch (std::exception &e) {
    cerr << "Caught std::exception: " << e.what() << '\n';
  } catch (...) {
    cerr << "Caught unknown exception\n";
  }

  return status;
}
