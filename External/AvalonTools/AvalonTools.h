// $Id: AvalonTools.h 1897 2011-12-23 06:17:39Z glandrum $
//
// Created by Greg Landrum, July 2008
//
#ifndef __AVALONTOOLS_H__
#define __AVALONTOOLS_H__
#include <string>

#include <GraphMol/RDKitBase.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/SparseIntVect.h>
#include <vector>
#include <boost/cstdint.hpp>

namespace AvalonTools {
static int avalonSSSBits = 0x007FFF;
static int avalonSimilarityBits = 0xF07FFF;
std::string getCanonSmiles(RDKit::ROMol &mol, int flags = -1);
void getAvalonCountFP(const RDKit::ROMol &mol,
                      RDKit::SparseIntVect<boost::uint32_t> &res,
                      unsigned int nBits = 512, bool isQuery = false,
                      bool resetVect = true,
                      unsigned int bitFlags = avalonSSSBits);
void getAvalonFP(const RDKit::ROMol &mol, ExplicitBitVect &res,
                 unsigned int nBits = 512, bool isQuery = false,
                 bool resetVect = true, unsigned int bitFlags = avalonSSSBits);
void getAvalonFP(const RDKit::ROMol &mol, std::vector<boost::uint32_t> &res,
                 unsigned int nBits = 512, bool isQuery = false,
                 bool resetVect = true, unsigned int bitFlags = avalonSSSBits);
unsigned int set2DCoords(RDKit::ROMol &mol, bool clearConfs = true);

std::string getCanonSmiles(const std::string &data, bool isSmiles,
                           int flags = -1);
void getAvalonCountFP(const std::string &data, bool isSmiles,
                      RDKit::SparseIntVect<boost::uint32_t> &res,
                      unsigned int nBits = 512, bool isQuery = false,
                      unsigned int bitFlags = avalonSSSBits);
void getAvalonFP(const std::string &data, bool isSmiles, ExplicitBitVect &res,
                 unsigned int nBits = 512, bool isQuery = false,
                 bool resetVect = true, unsigned int bitFlags = avalonSSSBits);
void getAvalonFP(const std::string &data, bool isSmiles,
                 std::vector<boost::uint32_t> &res, unsigned int nBits = 512,
                 bool isQuery = false, bool resetVect = true,
                 unsigned int bitFlags = avalonSSSBits);

std::string set2DCoords(const std::string &data, bool isSmiles);

int initCheckMol(const std::string &optString);
RDKit::ROMOL_SPTR checkMol(int &errors, RDKit::ROMol &inMol);
RDKit::ROMOL_SPTR checkMol(int &errors, const std::string &data, bool isSmiles);
std::pair<std::string, int> checkMolString(const std::string &data,
                                           bool isSmiles);
std::string getCheckMolLog();

void closeCheckMolFiles();
}
#endif
