//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <cctype>
#include <memory.h>
#include <cerrno>
#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/property_tree/json_parser.hpp>
#include "../../RDGeneral/StreamOps.h"
#include "../FileParsers/FileParsers.h"
#include "../FileParsers/MolSupplier.h"  //SDF
#include "../SmilesParse/SmilesParse.h"
#include "StructCheckerOptions.h"

#include "AugmentedAtomData.cpp"

namespace RDKit {
namespace StructureCheck {

bool parseOptionsJSON(const std::string &json, StructCheckerOptions &op) {
  if (json.empty()) return false;
  try {
    std::istringstream ss;
    ss.str(json);
    boost::property_tree::ptree pt;
    boost::property_tree::read_json(ss, pt);
    op = StructCheckerOptions();  // reset to default values
    op.AcidityLimit = pt.get<double>("AcidityLimit", op.AcidityLimit);
    op.RemoveMinorFragments =
        pt.get<bool>("RemoveMinorFragments", op.RemoveMinorFragments);
    op.DesiredCharge = pt.get<int>("DesiredCharge", op.DesiredCharge);
    op.CheckCollisions = pt.get<bool>("CheckCollisions", op.CheckCollisions);
    op.CollisionLimitPercent =
        pt.get<int>("CollisionLimitPercent", op.CollisionLimitPercent);
    op.MaxMolSize = pt.get<int>("MaxMolSize", op.MaxMolSize);
    op.ConvertSText = pt.get<bool>("ConvertSText", op.ConvertSText);
    op.SqueezeIdentifiers =
        pt.get<bool>("SqueezeIdentifiers", op.SqueezeIdentifiers);
    op.StripZeros = pt.get<bool>("StripZeros", op.StripZeros);
    op.CheckStereo = pt.get<bool>("CheckStereo", op.CheckStereo);
    op.ConvertAtomTexts = pt.get<bool>("ConvertAtomTexts", op.ConvertAtomTexts);
    op.GroupsToSGroups = pt.get<bool>("GroupsToSGroups", op.GroupsToSGroups);
    op.Verbose = pt.get<bool>("Verbose", op.Verbose);
  } catch (boost::property_tree::json_parser_error &ex) {
    BOOST_LOG(rdErrorLog) << "JSON StructureCheck Options:" << ex.message()
                          << "\n";
    return false;
  } catch (std::exception &ex) {
    BOOST_LOG(rdErrorLog) << "JSON StructureCheck Options:" << ex.what()
                          << "\n";
    return false;
  } catch (...) {
    BOOST_LOG(rdErrorLog)
        << "JSON StructureCheck Options: Unknown exception.\n";
    return false;
  }
  return true;
}

bool loadOptionsFromFiles(
    StructCheckerOptions &op, const std::string &augmentedAtomTranslationsFile,
    const std::string &patternFile,        // file with clean patterns
    const std::string &rotatePatternFile,  // file with rotate patterns
    const std::string &stereoPatternFile,  // file with stereo patterns
    const std::string &tautomerFile) {
  bool res = true;
  if (!augmentedAtomTranslationsFile.empty())
    res &= op.loadAugmentedAtomTranslations(augmentedAtomTranslationsFile);
  if (!patternFile.empty()) res &= op.loadPatterns(patternFile);
  if (!rotatePatternFile.empty())
    res &= op.loadRotatePatterns(rotatePatternFile);
  if (!stereoPatternFile.empty())
    res &= op.loadStereoPatterns(stereoPatternFile);
  if (!tautomerFile.empty()) res &= op.loadTautomerData(tautomerFile);
  return res;
}

//=====================================================================
// File parsers helper functions:

static const char
    *bond_to_string[] =  // ordered in according with BondType values
    {"?", "-", "=", "#", "~", "-=", "-~", "=~", "*"};

// used in unit test
bool StringToAugmentedAtom(const char *str, AugmentedAtom &aa) {
  /*
   * The syntax of a augmented atom string is as follows:
   *
   *   <AAString> ::= ['@'|'!@'] <Atom Symbol> [<Charge>]
   *                 {'(' <Bond Symbol> <Atom Symbol> [<Charge>] ')'}.
   * Bond symbols and charge descriptors are defined in the symbol tables
   * 'bond_to_string' and 'charge_to_string'.
   * '@' means central atom is in ring, '!@' means central atom is not in ring,
   *omitting means any topography could match.
   */

  size_t i;
  aa.ShortName = str;

  if (str[0] == '@') {
    aa.Topology = RING;
    str++;
  } else if (str[0] == '!' && str[1] == '@') {
    aa.Topology = CHAIN;
    str++;
    str++;
  }
  aa.Ligands.clear();  // initially no ligands

  // fetch atom symbol
  if (!isalpha(*str)) {  // atom symbol is required
    BOOST_LOG(rdErrorLog) << "atom symbol is required '" << str << "'\n";
    return false;
  }
  i = 0;
  while (isalnum(str[i]) || str[i] == ',' || str[i] == '#') i++;
  aa.AtomSymbol = std::string(str, i);
  str += i;

  // charge definition
  if (str[0] == '+' && str[1] == '-') {  // "+-" - ANY_CHARGE
    aa.Charge = ANY_CHARGE;
    str += 2;
  } else if (*str == '-' || *str == '+') {
    int sign = (*str == '-') ? (-1) : (1);
    str++;
    if (*str >= '1' && *str <= '7') {
      aa.Charge = sign * (*str - '0');
      str++;
    } else {  // no match
      BOOST_LOG(rdErrorLog) << "syntax error '" << str - 1 << "'\n";
      return false;
    }
  } else
    aa.Charge = 0;  // NONE; default

  // radical definition
  if (str == strpbrk(str, "|.:") &&
      !isdigit(str[1])) {  // don't confuse with degree specification
    switch (*str) {
      case '|':
        aa.Radical = SINGLET;  // RadicalType::SINGLET;
        str++;
        break;
      case '.':
        aa.Radical = DOUBLET;
        str++;
        break;
      case ':':
        aa.Radical = TRIPLET;
        str++;
        break;
      default:
        return false;  // never
    };
  } else
    aa.Radical = RT_NONE;

  // read ligand descriptions
  while (*str == '(') {
    aa.Ligands.push_back(Ligand());
    str++;  // skip '('
    // read bond type
    bool found = false;
    for (i = 0; i < sizeof(bond_to_string) / sizeof(*bond_to_string); i++) {
      if (0 == memcmp(str, bond_to_string[i], isalpha(str[1]) ? 1 : 2)) {
        found = true;
        aa.Ligands.back().BondType = (AABondType)(i);
        str += strlen(bond_to_string[i]);
        break;
      }
    }
    if (!found) {
      BOOST_LOG(rdErrorLog) << "syntax error '" << str
                            << "'\nthere must be a bond type symbol.\n";
      return false;
    }
    if (!isalpha(*str)) { /* there must be an atom symbol */
      BOOST_LOG(rdErrorLog)
          << "syntax error '" << str << "'\nthere must be an atom symbol.\n";
      return false;
    }

    i = 0; /* fetch atom symbol */
    while (isalnum(str[i]) || str[i] == ',' || str[i] == '#') i++;
    aa.Ligands.back().AtomSymbol = std::string(str, i);
    str += i;

    // charge definition
    if (str[0] == '+' && str[1] == '-') {  // "+-" - ANY_CHARGE
      aa.Ligands.back().Charge = ANY_CHARGE;
      str += 2;
    } else if (*str == '-' || *str == '+') {
      int sign = (*str == '-') ? (-1) : (1);
      str++;
      if (*str >= '1' && *str <= '7') {
        aa.Ligands.back().Charge = sign * (*str - '0');
        str++;
      } else {  // no match
        BOOST_LOG(rdErrorLog) << "syntax error '" << str - 1 << "'\n";
        return false;
      }
    } else
      aa.Ligands.back().Charge = 0;  // NONE; default

    // radical definition
    if (str == strpbrk(str, "|.:") &&
        !isdigit(str[1])) {  // don't confuse with degree specification
      switch (*str) {
        case '|':
          aa.Ligands.back().Radical = SINGLET;
          str++;
          break;
        case '.':
          aa.Ligands.back().Radical = DOUBLET;
          str++;
          break;
        case ':':
          aa.Ligands.back().Radical = TRIPLET;
          str++;
          break;
        default:
          return false;  // never
      };
    } else
      aa.Ligands.back().Radical = RT_NONE;

    // substitution count descriptor
    if (str[0] == ':' && isdigit(str[1])) {
      aa.Ligands.back().SubstitutionCount = str[1] - '0';
      str += 2;
    } else
      aa.Ligands.back().SubstitutionCount = 0;

    if (*str == ')')  // check for close ')'
      str++;
    else {
      BOOST_LOG(rdErrorLog) << "unclosed ( '" << str << "'\n";
      return false;
    }
  }

  if (*str != '\0') {  // there are some extra symbols after scan (inside quoted
                       // AA string)
    BOOST_LOG(rdErrorLog) << "unexpected symbols after scan '" << str << "'\n";
    return false;
  }
  return true;
}

static bool ReadAugmentedAtoms(const std::string &path,
                               std::vector<AugmentedAtom> &atoms) {
  /*
   * '*.chk' and '*.aci' files loader
   * Only the portion of the strings between the two '"' characters is used,
   * other characters are regarded as comments.
   */
  FILE *fp = fopen(path.c_str(), "r");
  if (!fp) {
    BOOST_LOG(rdErrorLog) << "could not open file '" << path << "'\n"
                          << strerror(errno) << "\n";
    return false;
  }
  unsigned n = 0;

  char str[1024];
  while (fgets(str, sizeof(str), fp)) {  // comment block in '*.aci' file
    if (*str != '#') break;
  }
  sscanf(str, "%d", &n);

  atoms.clear();
  atoms.resize(n);
  /* read: "! _N10_" - unused version string      // '*.chk' file only
              char buffer[80];
              char *cp;
              int vers_len;
              fgets(buffer, 80, fp);
              cp = strchr(buffer, '_');
              if (cp && strchr(cp + 1, '_'))
              {
                  vers_len = (int)(strchr(cp + 1, '_') - cp) - 1;
                  strncpy(aa_check_version, cp + 1, vers_len);
                  aa_check_version[vers_len] = '\0';
              }

              if (log_file)
                  fprintf(log_file, "augmented atom check version = %s\n",
     aa_check_version);
  */
  for (unsigned i = 0; i < n && !feof(fp); i++) {
    str[0] = '\0';
    if (fgets(str, sizeof(str), fp)) {
      char *s = str;
      while (*s >= ' ' && *s != '"') s++;
      if (*s == '"') {
        s++;
        size_t len = 0;
        while (s[len] >= ' ' && s[len] != '"') len++;
        s[len] = '\0';  // == atom_string
        if (!StringToAugmentedAtom(s, atoms[i])) {
          BOOST_LOG(rdErrorLog)
              << "unsuccessful AA parsing of " << i + 1 << " " << str << "\n";
          fclose(fp);
          return false;
        }
      } else {  // incorrect text line
        BOOST_LOG(rdErrorLog)
            << "unsuccessful AA parsing of " << i + 1 << " " << str << "\n";
        fclose(fp);
        return false;
      }
    } else if (ferror(fp)) {
      BOOST_LOG(rdErrorLog)
          << "unsuccessful AA parsing of " << i + 1 << " " << str << "\n";
      fclose(fp);
      return false;
    }
  }
  fclose(fp);
  return true;
}

static bool ReadAAPairs(
    const std::string &path,
    std::vector<std::pair<AugmentedAtom, AugmentedAtom>> &trans_pairs) {
  /*
   * '*.trn' file loader
   * Reads file and constructs an array of pairs of augmented atom
   * descriptions. The first one being the search pattern and the second
   * one the target pattern.
   * The function expects the number of augmented atom description pairs
   * on the first line and the strings corresponding to the data structures
   * on the n following lines.
   * Only the portion of the strings between the two '"' characters is used,
   * other characters are regarded as comments.
   */
  FILE *fp = fopen(path.c_str(), "r");
  if (!fp) {
    BOOST_LOG(rdErrorLog) << "could not open file '" << path << "'\n"
                          << strerror(errno) << "\n";
    return false;
  }
  unsigned n = 0;
  int num_scan = fscanf(fp, "%d", &n);
  TEST_ASSERT(num_scan == 1);

  char buffer[80];

  // Reads the version from the first line in the file
  if (fgets(buffer, sizeof(buffer), fp)) {
    /*
                int  vers_len;
                char *cp;
                cp = strchr(buffer, '_');
                if (cp && strchr(cp + 1, '_'))
                {
                    vers_len = (int)(strchr(cp + 1, '_') - cp) - 1;
                    strncpy(aa_trans_version, cp + 1, vers_len);
                    aa_trans_version[vers_len] = '\0';
                }
                if (log_file)
                    fprintf(log_file, "augmented atom transformation version =
       %s\n", aa_trans_version);
    */
  }
  trans_pairs.clear();
  trans_pairs.resize(n);

  for (unsigned i = 0; i < n; i++) {
    char str[1024];
    str[0] = '\0';
    if (fgets(str, sizeof(str), fp)) {
      char *s = str;
      while (*s >= ' ' && *s != '"') s++;
      if (*s == '"') {
        s++;
        size_t len = 0;
        while (s[len] >= ' ' && s[len] != '"') len++;
        s[len] = '\0';  // == atom_string

        if (!StringToAugmentedAtom(s, trans_pairs[i].first)) {
          BOOST_LOG(rdErrorLog)
              << "unsuccessful translation of " << str << "\n";
          fclose(fp);
          return false;
        }

        s += len;
        s++;
        // second atom quoted string
        while (*s >= ' ' && *s != '"') s++;
        if (*s == '"') {
          s++;
          len = 0;
          while (s[len] >= ' ' && s[len] != '"') len++;
          s[len] = '\0';  // == atom_string

          if (!StringToAugmentedAtom(s, trans_pairs[i].second)) {
            BOOST_LOG(rdErrorLog)
                << "unsuccessful translation of " << str << "\n";
            fclose(fp);
            return false;
          }
        }
      }
    }
  }
  fclose(fp);
  if (trans_pairs.size() != n) {
    BOOST_LOG(rdErrorLog) << "syntax error in translation file " << path
                          << " pair " << trans_pairs.size() + 1 << "\n";
    return false;
  }
  return true;
}

//=====================================================================
static void loadDefaultAugmentedAtoms(StructCheckerOptions &struchkOpts) {
  std::vector<AugmentedAtom> good;
  good.reserve(sizeof(DefaultGoodAtoms) / sizeof(*DefaultGoodAtoms));

  for (const auto &DefaultGoodAtom : DefaultGoodAtoms) {
    good.push_back(AugmentedAtom());
    if (!StringToAugmentedAtom(DefaultGoodAtom.str, good.back())) {
      throw "INTERNAL ERROR in default data";
    }
  }

  std::vector<AugmentedAtom> acidic;
  acidic.reserve(sizeof(DefaultAcidicAtoms) / sizeof(*DefaultAcidicAtoms));

  for (const auto &DefaultAcidicAtom : DefaultAcidicAtoms) {
    acidic.push_back(AugmentedAtom());
    if (!StringToAugmentedAtom(DefaultAcidicAtom.str, acidic.back())) {
      throw "INTERNAL ERROR in default data";
    }
  }

  std::vector<std::pair<AugmentedAtom, AugmentedAtom>> trans_pairs;
  trans_pairs.resize(sizeof(DefaultAugmentedAtomTransforms) /
                     sizeof(*DefaultAugmentedAtomTransforms));
  for (size_t i = 0; i < sizeof(DefaultAugmentedAtomTransforms) /
                             sizeof(*DefaultAugmentedAtomTransforms);
       ++i) {
    if (!StringToAugmentedAtom(DefaultAugmentedAtomTransforms[i].from,
                               trans_pairs[i].first)) {
      throw "INTERNAL Error in default augmented atom transforms";
    }
    if (!StringToAugmentedAtom(DefaultAugmentedAtomTransforms[i].to,
                               trans_pairs[i].second)) {
      throw "INTERNAL Error in default augmented atom transforms";
    }
  }

  struchkOpts.setGoodAugmentedAtoms(good);
  struchkOpts.setAcidicAugmentedAtoms(acidic);
  struchkOpts.setAugmentedAtomTranslations(trans_pairs);
}
//=====================================================================

bool StructCheckerOptions::loadAugmentedAtomTranslations(
    const std::string &path) {
  AugmentedAtomPairs.clear();
  if (path.empty()) return false;
  return ReadAAPairs(path, AugmentedAtomPairs);
}

void StructCheckerOptions::setAugmentedAtomTranslations(
    const std::vector<std::pair<AugmentedAtom, AugmentedAtom>> &aaPairs) {
  AugmentedAtomPairs = aaPairs;
}

bool StructCheckerOptions::loadAcidicAugmentedAtoms(const std::string &path) {
  AcidicAtoms.clear();
  if (path.empty()) return false;
  return ReadAugmentedAtoms(path, AcidicAtoms);
}

void StructCheckerOptions::setAcidicAugmentedAtoms(
    const std::vector<AugmentedAtom> &atoms) {
  AcidicAtoms = atoms;
}

bool StructCheckerOptions::loadGoodAugmentedAtoms(const std::string &path) {
  GoodAtoms.clear();
  if (path.empty()) return false;
  return ReadAugmentedAtoms(path, GoodAtoms);
}
void StructCheckerOptions::setGoodAugmentedAtoms(
    const std::vector<AugmentedAtom> &atoms) {
  GoodAtoms = atoms;
}

static bool loadSDF(const std::string &path, std::vector<ROMOL_SPTR> &mols) {
  mols.clear();
  try {
    SDMolSupplier suppler(path);
    while (!suppler.atEnd()) {
      ROMol *m = suppler.next();
      if (m) {
        // ?? MakeHydrogensImplicit()
        mols.push_back(ROMOL_SPTR(m));
      }
    }
  } catch (std::exception &ex) {
    BOOST_LOG(rdErrorLog) << ex.what();  // log error
    return false;
  }
  return true;
}

bool StructCheckerOptions::loadPatterns(const std::string &path) {
  Patterns.clear();
  if (path.empty()) return false;
  return loadSDF(path, Patterns);
}

void StructCheckerOptions::setPatterns(const std::vector<ROMOL_SPTR> &p) {
  Patterns = p;
}

bool StructCheckerOptions::loadRotatePatterns(const std::string &path) {
  RotatePatterns.clear();
  if (path.empty()) return false;
  return loadSDF(path, RotatePatterns);
}

void StructCheckerOptions::setRotatePatterns(const std::vector<ROMOL_SPTR> &p) {
  RotatePatterns = p;
}

bool StructCheckerOptions::loadStereoPatterns(const std::string &path) {
  StereoPatterns.clear();
  if (path.empty()) return false;
  return loadSDF(path, StereoPatterns);
}

void StructCheckerOptions::setStereoPatterns(const std::vector<ROMOL_SPTR> &p) {
  StereoPatterns = p;
}

bool StructCheckerOptions::loadTautomerData(const std::string &path) {
  ToTautomer.clear();
  FromTautomer.clear();
  if (path.empty()) return false;
  try {
    std::string ext = path.substr(path.find_last_of(".") + 1);
    if (ext == "sdf" || ext == "SDF") {
      SDMolSupplier suppler(path);
      while (!suppler.atEnd()) {
        ROMol *m1 = suppler.next();
        if (suppler.atEnd()) break;
        ROMol *m2 = suppler.next();
        if (m1 && m2) {
          // ?? MakeHydrogensImplicit()
          FromTautomer.push_back(ROMOL_SPTR(m1));
          ToTautomer.push_back(ROMOL_SPTR(m2));
        }
      }
    } else if (ext == "rdf" ||
               ext == "RDF") {  // load RD-File, see CTfile format specification
      std::ifstream in(path.c_str());
      if (!in) {
        BOOST_LOG(rdErrorLog) << "could not open file '" << path << "'\n"
                              << strerror(errno) << "\n";
        return false;
      }
      std::string str;
      unsigned line = 0;
      do {
        str = getLine(in);
        line++;
        if (str.length() >= 4 && 0 == strcmp(str.c_str(), "$RXN")) {
          while (!in.eof()) {
            str = getLine(in);
            line++;
            if (str.length() >= 4 && 0 == strcmp(str.c_str(), "$MOL")) {
              RWMol *molFrom = MolDataStreamToMol(in, line);
              if (!molFrom) {
                BOOST_LOG(rdErrorLog) << "RD-file corrupted. " << path
                                      << " line " << line << "\n";  // log error
                return false;
              }
              while (molFrom && !in.eof()) {
                str = getLine(in);
                line++;
                if (str.length() >= 4 && 0 == strcmp(str.c_str(), "$MOL")) {
                  RWMol *molTo = MolDataStreamToMol(in, line);
                  if (!molTo) {
                    BOOST_LOG(rdErrorLog)
                        << "RD-file corrupted. " << path << " line " << line
                        << "\n";  // log error
                    return false;
                  }
                  // ?? MakeHydrogensImplicit()
                  FromTautomer.push_back(ROMOL_SPTR(molFrom));
                  ToTautomer.push_back(ROMOL_SPTR(molTo));
                }
              }
            }
          }
        }
      } while (!in.eof());
    } else {
      BOOST_LOG(rdErrorLog)
          << "unsupported file type. " << path << "\n";  // log error
      return false;
    }
  } catch (std::exception &ex) {
    BOOST_LOG(rdErrorLog) << ex.what() << "\n";  // log error
    return false;
  }
  return true;
}

void StructCheckerOptions::setTautomerData(const std::vector<ROMOL_SPTR> &from,
                                           const std::vector<ROMOL_SPTR> &to) {
  FromTautomer = from;
  ToTautomer = to;
}

static void parseSMARTS(std::vector<ROMOL_SPTR> &mols,
                        const std::vector<std::string> &smarts) {
  mols.clear();
  for (const auto &patt : smarts) {
    ROMol *m = SmartsToMol(patt);
    if (m) mols.push_back(ROMOL_SPTR(m));
  }
}

void StructCheckerOptions::parsePatterns(
    const std::vector<std::string> &smarts) {
  parseSMARTS(Patterns, smarts);
}

void StructCheckerOptions::parseRotatePatterns(
    const std::vector<std::string> &smarts) {
  parseSMARTS(RotatePatterns, smarts);
}

void StructCheckerOptions::parseStereoPatterns(
    const std::vector<std::string> &smarts) {
  parseSMARTS(StereoPatterns, smarts);
}

void StructCheckerOptions::parseTautomerData(
    const std::vector<std::string> &smartsFrom,
    const std::vector<std::string> &smartsTo) {
  parseSMARTS(FromTautomer, smartsFrom);
  parseSMARTS(ToTautomer, smartsTo);
}
//--------------------------------------------------------------------------

StructCheckerOptions::StructCheckerOptions()
    : AcidityLimit(0.0),
      RemoveMinorFragments(false),
      DesiredCharge(0),
      CheckCollisions(false),
      CollisionLimitPercent(15),
      MaxMolSize(255),
      ConvertSText(false),
      SqueezeIdentifiers(false),
      StripZeros(false),
      CheckStereo(false),
      ConvertAtomTexts(false),
      GroupsToSGroups(false),
      Verbose(false),
      Elneg0(0.0)  // elneg_table[0].value;
      ,
      Alpha(0.0),
      Beta(0.0) {
  loadDefaultAugmentedAtoms(*this);
}

}  // namespace StructureCheck
}  // namespace RDKit
