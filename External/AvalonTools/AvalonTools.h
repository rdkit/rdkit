// $Id: AvalonTools.h 1897 2011-12-23 06:17:39Z glandrum $
//
// Created by Greg Landrum, July 2008
//
#ifndef __AVALONTOOLS_H__
#define __AVALONTOOLS_H__
#include <string>

#include <RDGeneral/export.h>
#include <GraphMol/RDKitBase.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/SparseIntVect.h>
#include <vector>
#include <boost/cstdint.hpp>

namespace AvalonTools {
const int avalonSSSBits = 0x007FFF;
const int avalonSimilarityBits = 0xF07FFF;

RDKIT_AVALONLIB_EXPORT std::string getCanonSmiles(RDKit::ROMol &mol,
                                                  int flags = -1);
RDKIT_AVALONLIB_EXPORT void getAvalonCountFP(
    const RDKit::ROMol &mol, RDKit::SparseIntVect<boost::uint32_t> &res,
    unsigned int nBits = 512, bool isQuery = false, bool resetVect = true,
    unsigned int bitFlags = avalonSSSBits);
RDKIT_AVALONLIB_EXPORT void getAvalonFP(const RDKit::ROMol &mol,
                                        ExplicitBitVect &res,
                                        unsigned int nBits = 512,
                                        bool isQuery = false,
                                        bool resetVect = true,
                                        unsigned int bitFlags = avalonSSSBits);
RDKIT_AVALONLIB_EXPORT void getAvalonFP(const RDKit::ROMol &mol,
                                        std::vector<boost::uint32_t> &res,
                                        unsigned int nBits = 512,
                                        bool isQuery = false,
                                        bool resetVect = true,
                                        unsigned int bitFlags = avalonSSSBits);
RDKIT_AVALONLIB_EXPORT unsigned int set2DCoords(RDKit::ROMol &mol,
                                                bool clearConfs = true);

RDKIT_AVALONLIB_EXPORT std::string getCanonSmiles(const std::string &data,
                                                  bool isSmiles,
                                                  int flags = -1);
RDKIT_AVALONLIB_EXPORT void getAvalonCountFP(
    const std::string &data, bool isSmiles,
    RDKit::SparseIntVect<boost::uint32_t> &res, unsigned int nBits = 512,
    bool isQuery = false, unsigned int bitFlags = avalonSSSBits);
RDKIT_AVALONLIB_EXPORT void getAvalonFP(const std::string &data, bool isSmiles,
                                        ExplicitBitVect &res,
                                        unsigned int nBits = 512,
                                        bool isQuery = false,
                                        bool resetVect = true,
                                        unsigned int bitFlags = avalonSSSBits);
RDKIT_AVALONLIB_EXPORT void getAvalonFP(const std::string &data, bool isSmiles,
                                        std::vector<boost::uint32_t> &res,
                                        unsigned int nBits = 512,
                                        bool isQuery = false,
                                        bool resetVect = true,
                                        unsigned int bitFlags = avalonSSSBits);

RDKIT_AVALONLIB_EXPORT std::string set2DCoords(const std::string &data,
                                               bool isSmiles);

RDKIT_AVALONLIB_EXPORT int initCheckMol(const std::string &optString);
RDKIT_AVALONLIB_EXPORT RDKit::ROMOL_SPTR checkMol(int &errors,
                                                  RDKit::ROMol &inMol);
RDKIT_AVALONLIB_EXPORT RDKit::ROMOL_SPTR checkMol(int &errors,
                                                  const std::string &data,
                                                  bool isSmiles);
RDKIT_AVALONLIB_EXPORT std::pair<std::string, int> checkMolString(
    const std::string &data, bool isSmiles);
RDKIT_AVALONLIB_EXPORT std::string getCheckMolLog();

RDKIT_AVALONLIB_EXPORT void closeCheckMolFiles();
}  // namespace AvalonTools
#endif
