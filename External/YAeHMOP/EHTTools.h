//
// Copyright (C) 2018 by Greg Landrum
//
#ifndef EHTTOOLS_H_20181226
#define EHTTOOLS_H_20181226

#include <string>

namespace RDKit {
class ROMol;
namespace EHTTools {
extern const std::string _EHTCharge;  // used to hold partial charges
extern const std::string _EHTMullikenOP; // used to hold overlap populations

bool runMol(const ROMol &mol, int confId=-1);

}
}
#endif
