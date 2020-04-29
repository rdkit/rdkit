#include <iostream>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

// fuzz_target.cc
extern "C" int LLVMFuzzerTestOneInput(const uint8_t *Data, size_t Size) {
  try {
    std::shared_ptr<RDKit::ROMol> mol3( RDKit::SmilesToMol(std::string (Data, Data + Size)) );
  } catch( RDKit::MolSanitizeException &e ) {
    std::cout << e.what() << std::endl;
  } catch( Invar::Invariant &e ) {
    std::cout << e.what() << std::endl;
  }
  return 0;  // Non-zero return values are reserved for future use.
}