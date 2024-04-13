//
//  Copyright (C) 2022-2023 Tad Hurst, Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/RDLog.h>
#include "MarvinDefs.h"

namespace RDKit {

std::string getDoubleAsText(double val, unsigned int precision = 6) {
  // this is left here for backwards compatibility
  if (precision == 6 && fabs(val) < 0.00001) {
    return "0.00000";
  }

  std::ostringstream valstr;
  valstr << std::setprecision(precision) << val;
  return valstr.str();
}

std::string MarvinArrow::toString() const {
  std::ostringstream out;
  out << "<arrow type=\"" << type << "\" x1=\"" << x1 << "\" y1=\"" << y1
      << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\"/>";

  return out.str();
}

ptree MarvinMolBase::toPtree() const {
  ptree out;
  out.put("<xmlattr>.molID", molID);

  if (parent != nullptr) {
    out.put("<xmlattr>.id", id);
    out.put("<xmlattr>.role", role());
  }
  if (this->hasAtomBondBlocks()) {
    for (auto atom : atoms) {
      out.add_child("atomArray.atom", atom->toPtree(this->coordinatePrecision));
    }

    for (auto bond : bonds) {
      out.add_child("bondArray.bond", bond->toPtree());
    }
  }

  return out;
}

void MarvinMolBase::addSgroupsToPtree(ptree &out) const {
  for (auto &sgroup : sgroups) {
    out.add_child("molecule", sgroup->toPtree());
  }
}

template <typename T>
bool getCleanNumber(std::string strToParse, T &outVal) {
  if (boost::algorithm::trim_copy(strToParse) !=
      strToParse) {  // should be no white space
    return false;
  }
  try {
    outVal = boost::lexical_cast<T>(strToParse);
  } catch (const std::exception &e) {
    return false;
  }

  return true;
}

void MarvinMolBase::parseAtomsAndBonds(ptree &molTree) {
    boost::property_tree::ptree atomArray = molTree.get_child("atomArray");

    // there are two types of atom arrays:
    // <atomArray atomID="a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11" elementType="C
    // C C C C C Cl C N O O" formalCharge="0 0 0 0 0 0 0 0 1 0 -1"
    // lonePair="0 0 0 0 0 0 3 0 0 2 3" x2="-4.3334 -5.6670 -5.6670 -4.3334
    // -2.9997 -2.9997 -4.3335 -1.6660 -7.0007 -1.6660 -0.3323" y2="1.7693
    // 0.9993 -0.5409 -1.3109 -0.5409 0.9993 3.3093 -1.3109 -1.3109 -2.8509
    // -0.5410"></atomArray>
    //  AND
    // <atomArray>
    //       <atom id="a1" elementType="C" x2="-9.4583" y2="1.9358"
    //       mrvStereoGroup="and1"/> <atom id="a2" elementType="C"
    //       x2="-10.7921" y2="1.1658"/> <atom id="a3" elementType="C"
    //       x2="-10.7921" y2="-0.3744"/> <atom id="a8" elementType="O"
    //       x2="-12.1257" y2="-1.1444" lonePair="2"/>
    //   </atomArray>

    // See which one we have

    std::string atomID =
        atomArray.get<std::string>("<xmlattr>.atomID", "UseLongForm");
    if (atomID == "UseLongForm") {
      // long form - each atom on a line

      for (auto &v : molTree.get_child("atomArray")) {
        if (v.first != "atom") {
          continue;
        }

        auto *mrvAtom = new MarvinAtom();
        this->pushOwnedAtom(mrvAtom);
        this->atoms.push_back(mrvAtom);

        mrvAtom->id = v.second.get<std::string>("<xmlattr>.id", "");
        mrvAtom->elementType =
            v.second.get<std::string>("<xmlattr>.elementType", "");

        if (mrvAtom->id == "" || mrvAtom->elementType == "") {
          throw FileParseException(
              "Expected id, elementType for an atom definition in MRV file");
        }

        std::string x2 = v.second.get<std::string>("<xmlattr>.x2", "");
        std::string y2 = v.second.get<std::string>("<xmlattr>.y2", "");

        // x2 and y2 are doubles

        if (x2 != "" && y2 != "" &&
            (!getCleanNumber(x2, mrvAtom->x2) ||
             !getCleanNumber(y2, mrvAtom->y2))) {
          throw FileParseException(
              "The values x2 and y2 must be large floats in MRV file");
        }

        std::string x3 = v.second.get<std::string>("<xmlattr>.x3", "");
        std::string y3 = v.second.get<std::string>("<xmlattr>.y3", "");
        std::string z3 = v.second.get<std::string>("<xmlattr>.z3", "");

        // x3, y3, and z3 are doubles

        if (x3 != "" && y3 != "" && z3 != "" &&
            (!getCleanNumber(x3, mrvAtom->x3) ||
             !getCleanNumber(y3, mrvAtom->y3) ||
             !getCleanNumber(z3, mrvAtom->z3))) {
          throw FileParseException(
              "The values x3, y3,  and z2 must be large floats in MRV file");
        }

        std::string formalCharge =
            v.second.get<std::string>("<xmlattr>.formalCharge", "");
        if (formalCharge != "") {
          if (!getCleanNumber(formalCharge, mrvAtom->formalCharge)) {
            throw FileParseException(
                "The value for formalCharge must be an integer in MRV file");
          }
        } else {
          mrvAtom->formalCharge = 0;
        }

        mrvAtom->radical = v.second.get<std::string>("<xmlattr>.radical", "");
        if (mrvAtom->radical != "") {
          if (!boost::algorithm::contains(
                  marvinRadicalVals,
                  std::vector<std::string>{mrvAtom->radical})) {
            std::ostringstream err;
            err << "The value for radical must be one of "
                << boost::algorithm::join(marvinRadicalVals, ", ")
                << " in MRV file";
            throw FileParseException(err.str());
          }
        } else {
          mrvAtom->radical = "";
        }

        std::string isotopeStr =
            v.second.get<std::string>("<xmlattr>.isotope", "");
        if (isotopeStr != "") {
          if (!getCleanNumber(isotopeStr, mrvAtom->isotope) ||
              mrvAtom->isotope <= 0) {
            throw FileParseException(
                "The value for isotope must be a positive number in MRV file");
          }
        } else {
          mrvAtom->isotope = 0;
        }

        std::string valenceStr =
            v.second.get<std::string>("<xmlattr>.mrvValence", "");
        if (valenceStr != "") {
          if (!getCleanNumber(valenceStr, mrvAtom->mrvValence) ||
              mrvAtom->mrvValence < 0) {
            throw FileParseException(
                "The value for mrvValence must be a positive number in MRV file");
          }
        } else {
          mrvAtom->mrvValence = -1;
        }

        std::string hCountStr =
            v.second.get<std::string>("<xmlattr>.hydrogenCount", "");
        if (hCountStr != "") {
          if (!getCleanNumber(hCountStr, mrvAtom->hydrogenCount) ||
              mrvAtom->hydrogenCount < 0) {
            throw FileParseException(
                "The value for hydrogenCount must be a non-negative number in MRV file");
          }
        } else {
          mrvAtom->hydrogenCount = -1;
        }

        mrvAtom->mrvAlias = v.second.get<std::string>("<xmlattr>.mrvAlias", "");

        mrvAtom->rgroupRef = v.second.get<int>("<xmlattr>.rgroupRef", -1);

        mrvAtom->mrvStereoGroup =
            v.second.get<std::string>("<xmlattr>.mrvStereoGroup", "");
        if (mrvAtom->mrvStereoGroup == "0") {
          mrvAtom->mrvStereoGroup = "";
        }

        std::string mrvMap = v.second.get<std::string>("<xmlattr>.mrvMap", "");
        if (mrvMap != "") {
          if (!getCleanNumber(mrvMap, mrvAtom->mrvMap) ||
              mrvAtom->mrvMap <= 0) {
            throw FileParseException(
                "The value for mrvMap must be an non-=negative integer in MRV file");
          }
        } else {
          mrvAtom->mrvMap = 0;
        }

        mrvAtom->sgroupRef =
            v.second.get<std::string>("<xmlattr>.sgroupRef", "");

        mrvAtom->sgroupAttachmentPoint =
            v.second.get<std::string>("<xmlattr>.sgroupAttachmentPoint", "");
      }
    } else  // single line form of atoms
    {
      // <atomArray atomID="a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11"
      // elementType="C C C C C C Cl C N O O" formalCharge="0 0 0 0 0 0 0 0
      // 1 0 -1" lonePair="0 0 0 0 0 0 3 0 0 2 3" x2="-4.3334 -5.6670
      // -5.6670 -4.3334 -2.9997 -2.9997 -4.3335 -1.6660 -7.0007 -1.6660
      // -0.3323" y2="1.7693 0.9993 -0.5409 -1.3109 -0.5409 0.9993 3.3093
      // -1.3109 -1.3109 -2.8509 -0.5410"></atomArray>

      std::vector<std::string> atomIds;
      size_t atomCount;
      if (atomID == "") {
        atomCount = 0;
      } else {
        boost::algorithm::split(atomIds, atomID, boost::algorithm::is_space());
        atomCount = atomIds.size();
      }

      std::vector<std::string> elementTypes;
      std::string elementType =
          atomArray.get<std::string>("<xmlattr>.elementType", "");
      boost::algorithm::split(elementTypes, elementType,
                              boost::algorithm::is_space());

      std::vector<std::string> x2s;
      std::string x2 = atomArray.get<std::string>("<xmlattr>.x2", "");
      boost::algorithm::split(x2s, x2, boost::algorithm::is_space());

      std::vector<std::string> y2s;
      std::string y2 = atomArray.get<std::string>("<xmlattr>.y2", "");
      boost::algorithm::split(y2s, y2, boost::algorithm::is_space());

      std::vector<std::string> formalCharges;
      std::string formalCharge =
          atomArray.get<std::string>("<xmlattr>.formalCharge", "");
      boost::algorithm::split(formalCharges, formalCharge,
                              boost::algorithm::is_space());

      std::vector<std::string> isotopes;
      std::string isotope = atomArray.get<std::string>("<xmlattr>.isotope", "");
      boost::algorithm::split(isotopes, isotope, boost::algorithm::is_space());

      std::vector<std::string> radicals;
      std::string radical = atomArray.get<std::string>("<xmlattr>.radical", "");
      boost::algorithm::split(radicals, radical, boost::algorithm::is_space());

      std::vector<std::string> hydrogenCounts;
      std::string hydrogenCount =
          atomArray.get<std::string>("<xmlattr>.hydrogenCount", "");
      boost::algorithm::split(hydrogenCounts, hydrogenCount,
                              boost::algorithm::is_space());

      std::vector<std::string> mrvValences;
      std::string mrvValence =
          atomArray.get<std::string>("<xmlattr>.mrvValence", "");
      boost::algorithm::split(mrvValences, mrvValence,
                              boost::algorithm::is_space());

      std::vector<std::string> mrvAliases;
      std::string mrvAlias =
          atomArray.get<std::string>("<xmlattr>.mrvAlias", "");
      boost::algorithm::split(mrvAliases, mrvAlias,
                              boost::algorithm::is_space());

      std::vector<std::string> rgroupRefs;
      std::string rgroupRef =
          atomArray.get<std::string>("<xmlattr>.rgroupRef", "");
      boost::algorithm::split(rgroupRefs, rgroupRef,
                              boost::algorithm::is_space());

      std::vector<std::string> mrvStereoGroups;
      std::string mrvStereoGroup =
          atomArray.get<std::string>("<xmlattr>.mrvStereoGroup", "");
      boost::algorithm::split(mrvStereoGroups, mrvStereoGroup,
                              boost::algorithm::is_space());

      std::vector<std::string> mrvMaps;
      std::string mrvMap = atomArray.get<std::string>("<xmlattr>.mrvMap", "");
      boost::algorithm::split(mrvMaps, mrvMap, boost::algorithm::is_space());

      std::vector<std::string> sgroupRefs;
      std::string sgroupRef =
          atomArray.get<std::string>("<xmlattr>.sgroupRef", "");
      boost::algorithm::split(sgroupRefs, sgroupRef,
                              boost::algorithm::is_space());

      std::vector<std::string> sgroupAttachmentPoints;
      std::string sgroupAttachmentPoint =
          atomArray.get<std::string>("<xmlattr>.sgroupAttachmentPoint", "");
      boost::algorithm::split(sgroupAttachmentPoints, sgroupAttachmentPoint,
                              boost::algorithm::is_space());

      if (atomID != "") {
        if (elementType == "") {
          throw FileParseException(
              "Expected an elementType array for an atomArray definition in MRV file");
        }
        if (elementTypes.size() < atomCount) {
          throw FileParseException(
              "There must be an element type for each atom id");
        }
      }

      for (size_t i = 0; i < atomCount; ++i) {
        auto *mrvAtom = new MarvinAtom();
        this->pushOwnedAtom(mrvAtom);

        this->atoms.push_back(mrvAtom);

        mrvAtom->id = atomIds[i];

        mrvAtom->elementType = elementTypes[i];

        if (x2 != "" && y2 != "" && x2s.size() > i && y2s.size() > i) {
          if (!getCleanNumber(x2s[i], mrvAtom->x2) ||
              !getCleanNumber(y2s[i], mrvAtom->y2)) {
            throw FileParseException(
                "The values x2 and y2 must be large floats in MRV file");
          }
        }

        if (formalCharge != "" && formalCharges.size() > i) {
          if (!getCleanNumber(formalCharges[i], mrvAtom->formalCharge)) {
            throw FileParseException(
                "The value for formalCharge must be an integer in MRV file");
          }
        } else {
          mrvAtom->formalCharge = 0;
        }

        if (isotope != "" && isotopes.size() > i) {
          if (!getCleanNumber(isotopes[i], mrvAtom->isotope)) {
            throw FileParseException(
                "The value for formalCharge must be an integer in MRV file");
          }
        } else {
          mrvAtom->isotope = 0;
        }

        if (mrvValence != "" && mrvValences.size() > i) {
          if (mrvValences[i] == "-") {
            mrvAtom->mrvValence = -1;
          } else if (!getCleanNumber(mrvValences[i], mrvAtom->mrvValence)) {
            throw FileParseException(
                "The value for mrvValences must be an integer in MRV file");
          }
        } else {
          mrvAtom->mrvValence = -1;
        }

        if (hydrogenCount != "" && hydrogenCounts.size() > i) {
          if (hydrogenCounts[i] != "-" &&
              !getCleanNumber(hydrogenCounts[i], mrvAtom->hydrogenCount)) {
            throw FileParseException(
                "The value for hydrogenCount must be an integer in MRV file");
          }
        } else {
          mrvAtom->hydrogenCount = -1;
        }

        if (radical != "" && radicals.size() > i && radicals[i] != "0") {
          mrvAtom->radical = radicals[i];
          if (!boost::algorithm::contains(
                  marvinRadicalVals,
                  std::vector<std::string>{mrvAtom->radical})) {
            std::ostringstream err;
            err << "The value for radical must be one of "
                << boost::algorithm::join(marvinRadicalVals, ", ")
                << " in MRV file";
            throw FileParseException(err.str());
          }
        } else {
          mrvAtom->radical = "";
        }

        if (mrvAlias != "" && mrvAliases.size() > i) {
          mrvAtom->mrvAlias = mrvAliases[i];
        } else {
          mrvAtom->mrvAlias = "";
        }

        if (rgroupRef != "" && rgroupRefs.size() > i) {
          if (!getCleanNumber(rgroupRefs[i], mrvAtom->rgroupRef)) {
            throw FileParseException(
                "rgroupRef value must be an integer in MRV file");
          }
        } else {
          mrvAtom->rgroupRef = -1;
        }

        if (mrvStereoGroup != "" && mrvStereoGroups.size() > i &&
            mrvStereoGroups[i] != "0")  // "0" is NOT a stereo group
        {
          mrvAtom->mrvStereoGroup = mrvStereoGroups[i];
        } else {
          mrvAtom->mrvStereoGroup = "";
        }

        if (mrvMap != "" && mrvMaps.size() > i) {
          if (!getCleanNumber(mrvMaps[i], mrvAtom->mrvMap) ||
              mrvAtom->mrvMap < 0) {
            throw FileParseException(
                "The value for mrvMap must be an non-negative integer in MRV file");
          }
        } else {
          mrvAtom->mrvMap = 0;
        }

        if (sgroupRef != "" && sgroupRefs.size() > i && sgroupRefs[i] != "0") {
          mrvAtom->sgroupRef = sgroupRefs[i];
        } else {
          mrvAtom->sgroupRef = "";
        }

        if (sgroupAttachmentPoint != "" && sgroupAttachmentPoints.size() > i &&
            sgroupAttachmentPoints[i] != "0") {
          mrvAtom->sgroupAttachmentPoint = sgroupAttachmentPoints[i];
        } else {
          mrvAtom->sgroupAttachmentPoint = "";
        }
      }
    }

    auto bondArray = molTree.get_child_optional("bondArray");
    if (bondArray) {
      for (auto &v : molTree.get_child("bondArray")) {
        if (v.first != "bond") {
          continue;
        }

        auto *mrvBond = new MarvinBond();
        this->pushOwnedBond(mrvBond);
        this->bonds.push_back(mrvBond);

        mrvBond->id = v.second.get<std::string>("<xmlattr>.id", "");
        if (mrvBond->id == "") {
          throw FileParseException(
              "Expected id for an bond definition in MRV file");
        }

        std::string atomRefs2 =
            v.second.get<std::string>("<xmlattr>.atomRefs2", "");

        std::vector<std::string> atomRefs2s;
        boost::algorithm::split(atomRefs2s, atomRefs2,
                                boost::algorithm::is_space());
        mrvBond->atomRefs2[0] = atomRefs2s[0];
        mrvBond->atomRefs2[1] = atomRefs2s[1];
        if (atomRefs2s.size() != 2 ||
            !boost::algorithm::contains(
                this->atoms, std::vector<std::string>{mrvBond->atomRefs2[0]},
              atomRefInAtoms) ||
            !boost::algorithm::contains(
                this->atoms, std::vector<std::string>{mrvBond->atomRefs2[1]},
              atomRefInAtoms)) {
          throw FileParseException(
              "atomRefs2 must contain two atom refs that must appear in the atoms array in MRV file");
        }

        mrvBond->order = v.second.get<std::string>("<xmlattr>.order", "");
        if (mrvBond->order != "") {
          if (!boost::algorithm::contains(
                  marvinBondOrders, std::vector<std::string>{mrvBond->order})) {
            std::ostringstream err;
            err << "Expected one of  "
                << boost::algorithm::join(marvinBondOrders, ", ")
                << " for order for an bond definition in MRV file";
            throw FileParseException(err.str());
          }
        }

        mrvBond->queryType =
            v.second.get<std::string>("<xmlattr>.queryType", "");
        if (mrvBond->queryType != "") {
          if (!boost::algorithm::contains(
                  marvinQueryBondsTypes,
                  std::vector<std::string>{mrvBond->queryType})) {
            std::ostringstream err;
            err << "Expected one of  "
                << boost::algorithm::join(marvinQueryBondsTypes, ", ")
                << " for queryType for an bond definition in MRV file";
            throw FileParseException(err.str());
          }
        }

        mrvBond->convention =
            v.second.get<std::string>("<xmlattr>.convention", "");
        if (mrvBond->convention != "") {
          if (!boost::algorithm::contains(
                  marvinConventionTypes,
                  std::vector<std::string>{mrvBond->convention})) {
            std::ostringstream err;
            err << "Expected one of  "
                << boost::algorithm::join(marvinConventionTypes, ", ")
                << " for convention for an bond definition in MRV file";
            throw FileParseException(err.str());
          }
        }

        int bondStereoDeclCount = 0;
        mrvBond->bondStereo.value = v.second.get<std::string>("bondStereo", "");
        if (mrvBond->bondStereo.value != "") {
          bondStereoDeclCount++;
        if (boost::algorithm::to_lower_copy(mrvBond->bondStereo.value) == "w" ||
            boost::algorithm::to_lower_copy(mrvBond->bondStereo.value) == "h") {
          // do nothing  - this is OK
        } else if (boost::algorithm::to_lower_copy(mrvBond->bondStereo.value) ==
                       "c" ||
              boost::algorithm::to_lower_copy(mrvBond->bondStereo.value) ==
                       "t") {
            mrvBond->bondStereo.value = "";  // cis and trans are ignored
          } else {
            throw FileParseException(
                "The value for bondStereo must be \"H\", \"W\", \"C\" or \"T\" in MRV file (\"C\" and \"T\" are ignored)");
          }
        }

        // see if bondstereo has a dictRef or convention

        auto bondStereoItem = v.second.get_child_optional("bondStereo");
        if (bondStereoItem) {
          for (auto &ww : v.second) {
            mrvBond->bondStereo.convention =
                ww.second.get<std::string>("<xmlattr>.convention", "");
            if (mrvBond->bondStereo.convention != "") {
              bondStereoDeclCount++;
              if (mrvBond->bondStereo.convention != "MDL") {
                throw FileParseException(
                    "Expected MDL as value for the bond convention attribute");
              }
              mrvBond->bondStereo.conventionValue =
                  ww.second.get<std::string>("<xmlattr>.conventionValue", "");
              if (!boost::algorithm::contains(
                      marvinStereoConventionTypes,
                      std::vector<std::string>{
                          mrvBond->bondStereo.conventionValue})) {
                std::ostringstream err;
                err << "Expected one of  "
                    << boost::algorithm::join(marvinStereoConventionTypes, ", ")
                    << " for a bond convention for an bond stereo def";
                throw FileParseException(err.str());
              }
            }

            mrvBond->bondStereo.dictRef =
                ww.second.get<std::string>("<xmlattr>.dictRef", "");
            if (mrvBond->bondStereo.dictRef != "") {
              bondStereoDeclCount++;
              if (!boost::algorithm::contains(
                      marvinStereoDictRefTypes,
                      std::vector<std::string>{mrvBond->bondStereo.dictRef})) {
                std::ostringstream err;
                err << "Expected one of  "
                    << boost::algorithm::join(marvinStereoDictRefTypes, ", ")
                    << " for a discRef value for an bond stereo def";
                throw FileParseException(err.str());
              }
            }
          }
        }

        // check that there were not too many different declarations

        if (bondStereoDeclCount > 1) {
          throw FileParseException(
              "bondStereo can either have only one of: a value, dictRef, or convention value");
        }
      }
    }
  }

ptree MarvinArrow::toPtree() const {
  ptree out;
  out.put("<xmlattr>.type", type);
  out.put("<xmlattr>.x1", getDoubleAsText(x1));
  out.put("<xmlattr>.y1", getDoubleAsText(y1));
  out.put("<xmlattr>.x2", getDoubleAsText(x2));
  out.put("<xmlattr>.y2", getDoubleAsText(y2));

  return out;
}

std::string MarvinPlus::toString() const {
  std::ostringstream out;

  out << "<MReactionSign id=\"" << id
      << "\" toptions=\"NOROT\" fontScale=\"14.0\" halign=\"CENTER\" valign=\"CENTER\" autoSize=\"true\">"
         "<Field name=\"text\"><![CDATA[ {D font=SansSerif,size=18,bold}+]]></Field>"
         "<MPoint x=\""
      << x1 << "\" y=\"" << y1
      << "\"/>"
         "<MPoint x=\""
      << x2 << "\" y=\"" << y1
      << "\"/>"
         "<MPoint x=\""
      << x2 << "\" y=\"" << y2
      << "\"/>"
         "<MPoint x=\""
      << x1 << "\" y=\"" << y2
      << "\"/>"
         "</MReactionSign>";

  return out.str();
}

ptree MarvinPlus::toPtree() const {
  ptree out;

  out.put("<xmlattr>.id", id);
  out.put("<xmlattr>.toptions", "NOROT");
  out.put("<xmlattr>.fontScale", "14.0");
  out.put("<xmlattr>.halign", "CENTER");
  out.put("<xmlattr>.valign", "CENTER");
  out.put("<xmlattr>.autoSize", "true");
  ptree field;
  field.put("<xmlattr>.name", "text");

  field.put_value("{D font=SansSerif,size=18,bold}+");

  out.put_child("Field", field);

  ptree p1, p2, p3, p4;
  p1.put("<xmlattr>.x", getDoubleAsText(x1));
  p1.put("<xmlattr>.y", getDoubleAsText(y1));
  out.add_child("MPoint", p1);

  p2.put("<xmlattr>.x", getDoubleAsText(x2));
  p2.put("<xmlattr>.y", getDoubleAsText(y1));
  out.add_child("MPoint", p2);

  p3.put("<xmlattr>.x", getDoubleAsText(x2));
  p3.put("<xmlattr>.y", getDoubleAsText(y2));
  out.add_child("MPoint", p3);

  p4.put("<xmlattr>.x", getDoubleAsText(x1));
  p4.put("<xmlattr>.y", getDoubleAsText(y2));
  out.add_child("MPoint", p4);

  return out;
}

std::string MarvinCondition::toString() const {
  std::ostringstream out;

  out << "<MTextBox id=\"" << id << "\" toption=\"NOROT\" halign=\"" << halign
      << "\" valign=\"" << valign << "\" autoSize=\"true\"";
  if (fontScale > 0) {
    out << " fontScale=\"" << fontScale << "\"";
  }

  out << ">"
         "<Field name=\"text\">"
      << text
      << "</Field>"
         "<MPoint x=\""
      << x << "\" y=\"" << y
      << "\"/>"
         "<MPoint x=\""
      << x << "\" y=\"" << y
      << "\"/>"
         "<MPoint x=\""
      << x << "\" y=\"" << y
      << "\"/>"
         "<MPoint x=\""
      << x << "\" y=\"" << y
      << "\"/>"
         "</MTextBox>";

  return out.str();
}

ptree MarvinCondition::toPtree() const {
  ptree out;

  out.put("<xmlattr>.id", id);
  if (fontScale > 0) {
    out.put("<xmlattr>.fontScale", fontScale);
  }
  out.put("<xmlattr>.toption", "NOROT");
  out.put("<xmlattr>.halign", halign);
  out.put("<xmlattr>.valign", valign);
  out.put("<xmlattr>.autoSize", "true");
  ptree field;
  field.put("<xmlattr>.name", "text");

  field.put_value("![CDATA[ {D font=SansSerif,size=18,bold}+]]");

  out.put_child("Field", field);

  // add four points

  std::ostringstream xstr, ystr;
  xstr << x;
  ystr << y;

  for (int i = 0; i < 4; ++i) {
    ptree p1;
    p1.put("<xmlattr>.x", xstr.str());
    p1.put("<xmlattr>.y", ystr.str());
    out.add_child("MPoint", p1);
  }

  return out;
}

std::string MarvinAttachmentPoint::toString() const {
  std::ostringstream out;

  out << "<attachmentPoint atom=\"" << atom << "\" order=\"" << order
      << "\" bond=\"" << bond << "\"/>";

  return out.str();
}

ptree MarvinAttachmentPoint::toPtree() const {
  ptree out;

  out.put("<xmlattr>.atom", atom);
  out.put("<xmlattr>.order", order);
  out.put("<xmlattr>.bond", bond);

  return out;
}

MarvinAtom::MarvinAtom()
    : x2(DBL_MAX),
      y2(DBL_MAX),
      x3(DBL_MAX),
      y3(DBL_MAX),
      z3(DBL_MAX),
      formalCharge(0),
      mrvValence(-1),
      hydrogenCount(-1),
      mrvMap(0),
      sGroupRefIsSuperatom(false),
      rgroupRef(-1)  // indicates that it was not specified

{}

MarvinAtom::MarvinAtom(const MarvinAtom &atomToCopy, std::string newId)
    : id(newId),
      elementType(atomToCopy.elementType),
      x2(atomToCopy.x2),
      y2(atomToCopy.y2),
      x3(atomToCopy.x3),
      y3(atomToCopy.y3),
      z3(atomToCopy.z3),
      formalCharge(atomToCopy.formalCharge),
      radical(atomToCopy.radical),
      isotope(atomToCopy.isotope),
      mrvValence(atomToCopy.mrvValence),
      hydrogenCount(atomToCopy.hydrogenCount),
      mrvAlias(atomToCopy.mrvAlias),
      mrvStereoGroup(atomToCopy.mrvStereoGroup),
      mrvMap(atomToCopy.mrvMap),
      sgroupRef(atomToCopy.sgroupRef),
      sGroupRefIsSuperatom(atomToCopy.sGroupRefIsSuperatom),
      sgroupAttachmentPoint(atomToCopy.sgroupAttachmentPoint),
      rgroupRef(atomToCopy.rgroupRef)  // indicates that it was not specified
{}

bool MarvinAtom::operator==(const MarvinAtom &rhs) const {
  return this->id == rhs.id;
}

bool MarvinAtom::operator==(const MarvinAtom *rhs) const {
  return this->id == rhs->id;
}

bool MarvinAtom::isElement() const {
  return this->elementType != "R" && this->elementType != "*";
}

std::string MarvinAtom::toString() const {
  // <atom id="a7" elementType="C" x2="15.225" y2="-8.3972"
  // sgroupAttachmentPoint="1"/>

  std::ostringstream out;
  out << "<atom id=\"" << id << "\" elementType=\"" << elementType << "\"";

  if (x2 != DBL_MAX && y2 != DBL_MAX) {
    out << " x2=\"" << x2 << "\" y2=\"" << y2 << "\"";
  }

  if (x3 != DBL_MAX && y3 != DBL_MAX && z3 != DBL_MAX) {
    out << " x3=\"" << x3 << "\" y3=\"" << y3 << "\" z3=\"" << z3 << "\"";
  }

  if (formalCharge != 0) {
    out << " formalCharge=\"" << formalCharge << "\"";
  }

  if (radical != "") {
    out << " radical=\"" << radical << "\"";
  }

  if (isElement() && isotope != 0) {
    out << " isotope=\"" << isotope << "\"";
  }

  if (mrvAlias != "") {
    out << " mrvAlias=\"" << mrvAlias << "\"";
  }

  if (mrvValence >= 0) {
    out << " mrvValence=\"" << mrvValence << "\"";
  }

  if (hydrogenCount > 0) {
    out << " hydrogenCount=\"" << hydrogenCount << "\"";
  }

  if (mrvStereoGroup != "") {
    out << " mrvStereoGroup=\"" << mrvStereoGroup << "\"";
  }

  if (mrvMap != 0) {
    out << " mrvMap=\"" << mrvMap << "\"";
  }

  if (sgroupRef != "") {
    out << " sgroupRef=\"" << sgroupRef << "\"";
  }

  if (rgroupRef >= 0) {
    out << " rgroupRef=\"" << rgroupRef << "\"";
  }

  if (sgroupAttachmentPoint != "") {
    out << " sgroupAttachmentPoint=\"" << sgroupAttachmentPoint << "\"";
  }

  out << "/>";

  return out.str();
}

ptree MarvinAtom::toPtree(unsigned int coordinatePrecision) const {
  ptree out;

  out.put("<xmlattr>.id", id);
  out.put("<xmlattr>.elementType", elementType);

  if (x2 != DBL_MAX && y2 != DBL_MAX) {
    out.put("<xmlattr>.x2", getDoubleAsText(x2, coordinatePrecision));
    out.put("<xmlattr>.y2", getDoubleAsText(y2, coordinatePrecision));
  }

  if (x3 != DBL_MAX && y3 != DBL_MAX && z3 != DBL_MAX) {
    out.put("<xmlattr>.x3", getDoubleAsText(x3, coordinatePrecision));
    out.put("<xmlattr>.y3", getDoubleAsText(y3, coordinatePrecision));
    out.put("<xmlattr>.z3", getDoubleAsText(z3, coordinatePrecision));
  }

  if (formalCharge != 0) {
    out.put("<xmlattr>.formalCharge", formalCharge);
  }

  if (radical != "") {
    out.put("<xmlattr>.radical", radical);
  }

  if (isElement() && isotope != 0) {
    out.put("<xmlattr>.isotope", isotope);
  }

  if (mrvAlias != "") {
    out.put("<xmlattr>.mrvAlias", mrvAlias);
  }

  if (mrvValence >= 0) {
    out.put("<xmlattr>.mrvValence", mrvValence);
  }

  if (hydrogenCount > 0) {
    out.put("<xmlattr>.hydrogenCount", hydrogenCount);
  }

  if (mrvStereoGroup != "") {
    out.put("<xmlattr>.mrvStereoGroup", mrvStereoGroup);
  }

  if (mrvMap != 0) {
    out.put("<xmlattr>.mrvMap", mrvMap);
  }

  if (sgroupRef != "") {
    out.put("<xmlattr>.sgroupRef", sgroupRef);
  }

  if (rgroupRef >= 0) {
    out.put("<xmlattr>.rgroupRef", rgroupRef);
  }

  if (sgroupAttachmentPoint != "") {
    out.put("<xmlattr>.sgroupAttachmentPoint", sgroupAttachmentPoint);
  }

  return out;
}

MarvinBond::MarvinBond(const MarvinBond &bondToCopy, std::string newId,
                       std::string atomRef1, std::string atomRef2)
    : id(newId),
      order(bondToCopy.order),
      bondStereo(bondToCopy.bondStereo),
      queryType(bondToCopy.queryType),
      convention(bondToCopy.convention) {
  atomRefs2[0] = atomRef1;
  atomRefs2[1] = atomRef2;
}

const std::string MarvinBond::getBondType() const {
  std::string tempQueryType = boost::algorithm::to_upper_copy(queryType);
  std::string tempOrder = boost::algorithm::to_upper_copy(order);
  std::string tempConvention = boost::algorithm::to_upper_copy(convention);

  if (tempQueryType != "") {
    if (tempQueryType == "SD" || tempQueryType == "SA" ||
        tempQueryType == "DA" || tempQueryType == "ANY") {
      return tempQueryType;
    } else {
      std::ostringstream err;
      err << "unrecognized query bond type " << queryType << " in MRV File ";
      throw FileParseException(err.str());
    }
  } else if (tempConvention != "")  // if no query type, check for convention
  {
    if (tempConvention == "CXN:COORD") {
      return "DATIVE";
    } else {
      std::ostringstream err;
      err << "unrecognized convention " << convention << " in MRV File ";
      throw FileParseException(err.str());
    }
  } else if (tempOrder !=
             "")  // if no query type not convention,  so check for order
  {
    if (tempOrder == "1" || tempOrder == "2" || tempOrder == "3" ||
        tempOrder == "A") {
      return tempOrder;
    } else {
      std::ostringstream err;
      err << "unrecognized bond type " << order << " in MRV File ";
      throw FileParseException(err.str());
    }
  } else {
    throw FileParseException(
        "bond must have one of:  order, queryType, or convention in MRV File ");
  }
}

bool MarvinBond::isEqual(const MarvinAtom &other) const {
  return this->id == other.id;
}

bool MarvinBond::operator==(const MarvinAtom &rhs) const {
  return this->isEqual(rhs);
}

std::string MarvinBond::toString() const {
  // <bond id="b8" atomRefs2="a1 a7" order="1">

  std::ostringstream out;

  out << "<bond id=\"" << id << "\" atomRefs2=\"" << atomRefs2[0] << " "
      << atomRefs2[1] << "\"";

  if (order != "") {
    out << " order=\"" << order << "\"";
  }

  if (queryType != "") {
    out << " queryType=\"" << queryType << "\"";
  }

  if (convention != "") {
    out << " convention=\"" << convention << "\"";
  }

  if (bondStereo.value != "" || bondStereo.dictRef != "" ||
      (bondStereo.convention != "" && bondStereo.conventionValue != "")) {
    out << "><bondStereo";
    if (bondStereo.value != "") {
      out << ">" << bondStereo.value << "</bondStereo>";
    } else if (bondStereo.dictRef != "") {
      out << " dictRef=\"" << bondStereo.dictRef << "\"/>";
    } else {
      out << " convention=\"" << bondStereo.convention
          << "\" conventionValue=\"" << bondStereo.conventionValue << "\"/>";
    }
    out << "</bond>";
  } else {
    out << "/>";  // just end the bond
  }
  return out.str();
}

ptree MarvinBond::toPtree() const {
  ptree out;

  out.put("<xmlattr>.id", id);
  out.put("<xmlattr>.atomRefs2", atomRefs2[0] + " " + atomRefs2[1]);

  if (order != "") {
    out.put("<xmlattr>.order", order);
  }

  if (queryType != "") {
    out.put("<xmlattr>.queryType", queryType);
  }

  if (convention != "") {
    out.put("<xmlattr>.convention", convention);
  }

  if (bondStereo.value != "" || bondStereo.dictRef != "" ||
      (bondStereo.convention != "" && bondStereo.conventionValue != "")) {
    ptree bondStereoPt;

    if (bondStereo.value != "") {
      bondStereoPt.put_value(bondStereo.value);
    } else if (bondStereo.dictRef != "") {
      bondStereoPt.put("<xmlattr>.dictRef", bondStereo.dictRef);
    } else {
      bondStereoPt.put("<xmlattr>.convention", bondStereo.convention);
      bondStereoPt.put("<xmlattr>.conventionValue", bondStereo.conventionValue);
    }

    out.put_child("bondStereo", bondStereoPt);
  }
  return out;
}

MarvinMolBase::~MarvinMolBase() {}

int MarvinMolBase::getAtomIndex(std::string id) const {
  auto atomIter =
      find_if(atoms.begin(), atoms.end(),
              [id](const MarvinAtom *arg) { return arg->id == id; });
  if (atomIter != atoms.end()) {
    return atomIter - atoms.begin();
  } else {
    return -1;
  }
}

void MarvinMolBase::pushOwnedAtom(MarvinAtom *atom) {
  PRECONDITION(atom, "bad atom");
  PRECONDITION(this->parent, "only sgroups should call the base class version");
  this->parent->pushOwnedAtom(atom);
}

void MarvinMolBase::pushOwnedBond(MarvinBond *bond) {
  PRECONDITION(bond, "bad bond");
  PRECONDITION(this->parent, "only sgroups should call the base class version");
  this->parent->pushOwnedBond(std::move(bond));
}

void MarvinMolBase::removeOwnedAtom(MarvinAtom *atom) {
  PRECONDITION(this->parent, "only sgroups should call the base class version");
  this->parent->removeOwnedAtom(atom);
}

void MarvinMolBase::removeOwnedBond(MarvinBond *bond) {
  PRECONDITION(this->parent, "only sgroups should call the base class version");
  this->parent->removeOwnedBond(bond);
}

int MarvinMolBase::getBondIndex(std::string id) const {
  auto bondIter =
      find_if(bonds.begin(), bonds.end(),
              [id](const MarvinBond *arg) { return arg->id == id; });
  if (bondIter != bonds.end()) {
    return bondIter - bonds.begin();
  } else {
    return -1;
  }
}

MarvinAtom *MarvinMolBase::findAtomByRef(std::string atomId) {
  auto atomIter =
      find_if(this->atoms.begin(), this->atoms.end(),
              [atomId](const MarvinAtom *arg) { return arg->id == atomId; });
  if (atomIter != this->atoms.end()) {
    return *atomIter;
  }

  // try the parent if there is one

  if (this->parent != nullptr) {
    return this->parent->findAtomByRef(atomId);
  }

  return nullptr;
}

MarvinBond *MarvinMolBase::findBondByRef(std::string bondId) {
  auto bondIter =
      find_if(this->bonds.begin(), this->bonds.end(),
              [bondId](const MarvinBond *arg) { return arg->id == bondId; });
  if (bondIter != this->bonds.end()) {
    return *bondIter;
  }

  // try the parent if there is one

  if (this->parent != nullptr) {
    return this->parent->findBondByRef(bondId);
  }

  return nullptr;
}

const std::vector<std::string> MarvinMolBase::getBondList() const {
  std::vector<std::string> bondList;
  for (auto bond : bonds) {
    bondList.push_back(bond->id);
  }

  return bondList;
}

const std::vector<std::string> MarvinMolBase::getAtomList() const {
  std::vector<std::string> atomList;
  for (auto atom : atoms) {
    atomList.push_back(atom->id);
  }

  return atomList;
}

bool MarvinMolBase::has2dCoords() const {
  for (auto atom : atoms) {
    if (atom->x2 == DBL_MAX || atom->y2 == DBL_MAX) {
      return false;
    }
  }

  return true;
}

bool MarvinMolBase::has3dCoords() const {
  for (auto atom : atoms) {
    if (atom->x3 == DBL_MAX || atom->y3 == DBL_MAX || atom->z3 == DBL_MAX) {
      return false;
    }
  }

  return true;
}

bool MarvinMolBase::hasAny2dCoords() const {
  for (auto atom : atoms) {
    if (atom->x2 != DBL_MAX && atom->y2 != DBL_MAX) {
      return true;
    }
  }

  return false;
}

bool MarvinMolBase::hasAny3dCoords() const {
  for (auto atom : atoms) {
    if (atom->x3 != DBL_MAX && atom->y3 != DBL_MAX && atom->z3 != DBL_MAX) {
      return true;
    }
  }

  return false;
}

bool MarvinMolBase::hasCoords() const { return has2dCoords() || has3dCoords(); }

void MarvinMolBase::removeCoords() {
  for (auto atom : atoms) {
    atom->x2 = DBL_MAX;
    atom->y2 = DBL_MAX;
  }
}

void MarvinMolBase::setPrecision(unsigned int precision) {
  this->coordinatePrecision = precision;
}

int MarvinMolBase::getExplicitValence(const MarvinAtom &marvinAtom) const {
  unsigned int resTimes10 = 0;  // calculated as 10 * the actual value so we
                                // can use int match, and have 1.5 order bonds

  for (auto bondPtr : bonds) {
    if (bondPtr->atomRefs2[0] != marvinAtom.id &&
        bondPtr->atomRefs2[1] != marvinAtom.id) {
      continue;  // this bond is NOT to the atom
    }

    std::string marvinBondType = bondPtr->getBondType();

    if (marvinBondType == "SD" || marvinBondType == "SA" ||
        marvinBondType == "DA") {
      resTimes10 += 15;  // really 1.5 order bond
    } else if (marvinBondType == "ANY") {
      resTimes10 +=
          10;  // no good answer for Any bonds - treat as a single bond
    } else if (marvinBondType == "DATIVE") {
      //
      //  the following code has been removed because we really need both ends
      //  of the dative bond to be counted as zero
      //
      // if (bondPtr->atomRefs2[1] == marvinAtom.id) //second atom of dative
      // bond count as 1 (first atom is zero)
      //   resTimes10 += 10;  // really 1 order bond
    } else if (marvinBondType == "1") {
      resTimes10 += 10;  // really 1 order bond
    } else if (marvinBondType == "2") {
      resTimes10 += 20;  // really 2 order bond
    } else if (marvinBondType == "3") {
      resTimes10 += 30;  // really 3 order bond
    } else if (marvinBondType == "A") {
      resTimes10 += 15;  // really 1.5 order bond
    }
  }

  return resTimes10 / 10;
}

MarvinSruCoModSgroup::MarvinSruCoModSgroup(std::string roleNameInit,
                                           MarvinMolBase *parentInit) {
  PRECONDITION(parentInit != nullptr, "parentInit cannot be null");

  parent = parentInit;

  if (roleNameInit != "SruSgroup" && roleNameInit != "CopolymerSgroup" &&
      roleNameInit != "ModificationSgroup") {
    throw FileParseException(
        "A MarvinSruCoModSgroup type must be one of \"SruSgroup\". \"CopolymerSgroup\", and \"ModificationSgroup\"");
  }

  this->roleName = roleNameInit;
}

MarvinSruCoModSgroup::MarvinSruCoModSgroup(MarvinMolBase *parentInit,
                                           std::string roleNameInit,
                                           ptree &molTree) {
  PRECONDITION(parentInit != nullptr, "parentInit cannot be null");

  //      <molecule molID="m2" id="sg1" role="SruSgroup" atomRefs="a1 a3
  //      a4" title="n" connect="ht" correspondence="" bondList=""/>
  //      <molecule molID="m3" id="sg2" role="SruSgroup" atomRefs="a5 a6
  //      a7 a8" title="n" connect="hh" correspondence="" bondList=""/>
  //      <molecule molID="m4" id="sg3" role="SruSgroup" atomRefs="a10
  //      a11" title="n" connect="eu" correspondence=""
  //      bondList=""/></molecule>

  parent = parentInit;
  this->id = molTree.get<std::string>("<xmlattr>.id", "");
  if (this->id == "") {
    throw FileParseException("Expected a iD in MRV file");
  }
  this->molID = molTree.get<std::string>("<xmlattr>.molID", "");
  if (this->molID == "") {
    throw FileParseException("Expected a molID in MRV file");
  }

  this->roleName = roleNameInit;
  std::string atomRefsStr = molTree.get<std::string>("<xmlattr>.atomRefs", "");
  if (atomRefsStr == "") {
    throw FileParseException(
        "Expected  atomRefs for a SruSgroup definition in MRV file");
  }

  std::vector<std::string> atomList;
  boost::algorithm::split(atomList, atomRefsStr, boost::algorithm::is_space());
  for (auto &it : atomList) {
    auto atomPtr = this->parent->findAtomByRef(it);

    if (atomPtr == nullptr) {
      throw FileParseException(
          "AtomRef specification for an SRU, Copolymer, or Modification group definition was not found in any parent");
    }
    this->atoms.push_back(atomPtr);
  }

  this->title = molTree.get<std::string>("<xmlattr>.title", "");
  if (this->title == "") {
    throw FileParseException(
        "Expected title for a SruSgroup definition in MRV file");
  }

  // the title for an SRU group must be a lower case letter, or a range
  // of two positive ints (4-6)

  this->connect = molTree.get<std::string>("<xmlattr>.connect", "");
  if (this->connect == "") {
    this->connect = "ht";
  }
  if (!boost::algorithm::contains(sruSgroupConnectChoices,
                                  std::vector<std::string>{this->connect})) {
    std::ostringstream err;
    err << "Expected a connect  string of \""
        << boost::algorithm::join(sruSgroupConnectChoices, ", ")
        << "\" for a SruSgroup definition in MRV file";
    throw FileParseException(err.str());
  }

  this->correspondence =
      molTree.get<std::string>("<xmlattr>.correspondence", "");

  std::string bondListStr = molTree.get<std::string>("<xmlattr>.bondList", "");
  if (bondListStr != "") {
    std::vector<std::string> bondList;
    boost::algorithm::split(bondList, bondListStr,
                            boost::algorithm::is_space());

    for (auto &it : bondList) {
      auto bondPtr = this->parent->findBondByRef(it);

      if (bondPtr == nullptr) {
        throw FileParseException(
            "BondList specification for an SRU, Copolymer, or Modification group definition was not found in any parent");
      }
      this->bonds.push_back(bondPtr);
    }
  }
}

std::string MarvinSruCoModSgroup::toString() const {
  std::ostringstream out;

  out << "<molecule molID=\"" << molID << "\" id=\"" << id << "\" role=\""
      << roleName << "\" atomRefs=\""
      << boost::algorithm::join(getAtomList(), " ") << "\" title=\"" << title
      << "\" connect=\"" << connect << "\" correspondence=\"" << correspondence
      << "\" bondList=\"" << boost::algorithm::join(getBondList(), " ")
      << "\"/>";

  return out.str();
}

ptree MarvinSruCoModSgroup::toPtree() const {
  ptree out = MarvinMolBase::toPtree();

  out.put("<xmlattr>.atomRefs", boost::algorithm::join(getAtomList(), " "));
  out.put("<xmlattr>.title", title);
  out.put("<xmlattr>.connect", connect);
  out.put("<xmlattr>.correspondence", correspondence);
  out.put("<xmlattr>.bondList", boost::algorithm::join(getBondList(), " "));

  return out;
}

MarvinMolBase *getActualParent(const MarvinMolBase *child) {
  PRECONDITION(child, "child cannot be null");
  for (auto actualParent = child->parent;;
       actualParent = actualParent->parent) {
    TEST_ASSERT(actualParent);

    if (actualParent->role() != "MultipleSgroup") {
      return actualParent;
      ;
    }
  }

  TEST_ASSERT(false);  // should not be reachable
}

MarvinMolBase *MarvinSruCoModSgroup::copyMol(
    const std::string &idAppendage) const {
  auto outSgroup = new MarvinSruCoModSgroup(this->roleName, this->parent);
  this->parent->sgroups.emplace_back(outSgroup);

  outSgroup->molID = this->molID + idAppendage;
  outSgroup->title = this->title;
  outSgroup->connect = this->connect;
  outSgroup->correspondence = this->correspondence;

  auto actualParent = getActualParent(this);

  // the only time this is to be called is when a multiple group above it is
  // being expanded and the new atoms and bonds have already been made
  //  here we just need to find them and add refs to the new multiple sgroup

  for (auto atomToCopy : atoms) {
    auto atomIndex = actualParent->getAtomIndex(atomToCopy->id + idAppendage);
    if (atomIndex < 0) {
      throw FileParseException(
          "unexpected error - new atom not found in parent");
    }

    outSgroup->atoms.push_back(actualParent->atoms[atomIndex]);
  }

  for (auto bondToCopy : bonds) {
    auto bondIndex = actualParent->getBondIndex(bondToCopy->id + idAppendage);
    if (bondIndex < 0) {
      throw FileParseException(
          "unexpected error - new bond not found in parent");
    }

    outSgroup->bonds.push_back(actualParent->bonds[bondIndex]);
  }

  return outSgroup;
}

std::string MarvinSruCoModSgroup::role() const { return roleName; }

bool MarvinSruCoModSgroup::hasAtomBondBlocks() const { return false; }

MarvinDataSgroup::MarvinDataSgroup(MarvinMolBase *parentInit) {
  parent = parentInit;
}

MarvinDataSgroup::MarvinDataSgroup(MarvinMolBase *parentInit, ptree &molTree) {
  PRECONDITION(parentInit != nullptr, "parentInit cannot be null");
  parent = parentInit;

  this->id = molTree.get<std::string>("<xmlattr>.id", "");
  if (this->id == "") {
    throw FileParseException("Expected a iD in MRV file");
  }
  this->molID = molTree.get<std::string>("<xmlattr>.molID", "");
  if (this->molID == "") {
    throw FileParseException("Expected a molID in MRV file");
  }

  std::string atomRefsStr = molTree.get<std::string>("<xmlattr>.atomRefs", "");
  if (atomRefsStr == "") {
    throw FileParseException(
        "Expected  atomRefs for a DataSgroup definition in MRV file");
  }

  std::vector<std::string> atomList;
  boost::algorithm::split(atomList, atomRefsStr, boost::algorithm::is_space());
  for (auto &it : atomList) {
    auto atomPtr = this->parent->findAtomByRef(it);

    if (atomPtr == nullptr) {
      throw FileParseException(
          "AtomRef specification for a DataSgroup definition was not found in any parent");
    }

    this->atoms.push_back(atomPtr);
  }

  this->context = molTree.get<std::string>("<xmlattr>.context", "");

  this->fieldName = molTree.get<std::string>("<xmlattr>.fieldName", "");
  this->placement = molTree.get<std::string>("<xmlattr>.placement", "");

  this->unitsDisplayed = molTree.get<std::string>("<xmlattr>.unitsDisplayed",
                                                  "Unit not displayed");
  if (this->unitsDisplayed == "") {
    this->unitsDisplayed = "Unit not displayed";
  }
  std::string unitsDisplayed =
      boost::algorithm::to_lower_copy(this->unitsDisplayed);

  if (unitsDisplayed != "unit displayed" &&
      unitsDisplayed != "unit not displayed") {
    throw FileParseException(
        "Expected unitsDisplayed to be either \"Unit displayed\" or \"Unit not displayed\" for a DataSgroup definition in MRV file");
  }

  this->queryType = molTree.get<std::string>("<xmlattr>.queryType", "");
  this->queryOp = molTree.get<std::string>("<xmlattr>.queryOp", "");
  this->units = molTree.get<std::string>("<xmlattr>.units", "");

  this->fieldData = molTree.get<std::string>("<xmlattr>.fieldData", "");

  std::string x = molTree.get<std::string>("<xmlattr>.x", "0.0");
  std::string y = molTree.get<std::string>("<xmlattr>.y", "0.0");

  if (x == "") {
    throw FileParseException(
        "Expected x for a DataSgroup definition in MRV file");
  }
  if (!getCleanNumber(x, this->x)) {
    throw FileParseException(
        "The value for x must be a floating point value in MRV file");
  }
  if (y == "") {
    throw FileParseException(
        "Expected y for a DataSgroup definition in MRV file");
  }
  if (!getCleanNumber(y, this->y)) {
    throw FileParseException(
        "The value for y must be a floating point value in MRV file");
  }
}

MarvinMolBase *MarvinDataSgroup::copyMol(const std::string &idAppendage) const {
  auto outSgroup = new MarvinDataSgroup(this->parent);
  this->parent->sgroups.emplace_back(outSgroup);

  outSgroup->molID = this->molID + idAppendage;
  outSgroup->id = this->id + idAppendage;
  outSgroup->context = this->context;
  outSgroup->fieldName = this->fieldName;
  outSgroup->placement = this->placement;
  outSgroup->unitsDisplayed = this->unitsDisplayed;
  outSgroup->units = this->units;
  outSgroup->fieldData = this->fieldData;
  outSgroup->queryType = this->queryType;
  outSgroup->queryOp = this->queryOp;
  outSgroup->x = this->x;
  outSgroup->y = this->y;

  // the only time this is to be called is when a multiple group above it is
  // being expanded and the new atoms and bonds have already been made
  //  here we just need to find them and add refs to the new multiple sgroup

  for (auto atomToCopy : atoms) {
    auto atomIndex = this->parent->getAtomIndex(atomToCopy->id + idAppendage);
    if (atomIndex < 0) {
      throw FileParseException(
          "unexpected error - new atom not found in parent");
    }

    outSgroup->atoms.push_back(this->parent->atoms[atomIndex]);
  }

  return outSgroup;
}

std::string MarvinDataSgroup::toString() const {
  std::ostringstream out;

  out << "<molecule molID=\"" << molID << "\" id=\"" << id
      << "\" role=\"DataSgroup\" atomRefs=\""
      << boost::algorithm::join(getAtomList(), " ") << "\" context=\""
      << context << "\" fieldName=\"" << fieldName << "\" placement=\""
      << placement << "\" unitsDisplayed=\"" << unitsDisplayed
      << "\" fieldData=\"" << fieldData;

  if (units != "") {
    out << "\" units=\"" << units;
  }

  if (queryType != "" && queryOp != "") {
    out << "\" queryType=\"" << queryType << "\" queryOp=\""
        << boost::property_tree::xml_parser::encode_char_entities(queryOp);
  }

  out << "\" x=\"" << x << "\" y=\"" << y << "\"/>";

  return out.str();
}

ptree MarvinDataSgroup::toPtree() const {
  ptree out = MarvinMolBase::toPtree();

  out.put("<xmlattr>.atomRefs", boost::algorithm::join(getAtomList(), " "));
  out.put("<xmlattr>.context", context);
  out.put("<xmlattr>.fieldName", fieldName);
  out.put("<xmlattr>.placement", placement);
  out.put("<xmlattr>.unitsDisplayed", unitsDisplayed);
  out.put("<xmlattr>.fieldData", fieldData);

  if (units != "") {
    out.put("<xmlattr>.units", units);
  }

  if (queryType != "" && queryOp != "") {
    out.put("<xmlattr>.queryType", queryType);
    out.put("<xmlattr>.queryOp", queryOp);
  }
  std::ostringstream xstr, ystr;
  xstr << x;
  ystr << y;
  out.put("<xmlattr>.x", xstr.str());
  out.put("<xmlattr>.y", ystr.str());

  return out;
}

std::string MarvinDataSgroup::role() const { return "DataSgroup"; }

bool MarvinDataSgroup::hasAtomBondBlocks() const { return false; }

MarvinMultipleSgroup::MarvinMultipleSgroup(MarvinMolBase *parentInit) {
  parent = parentInit;
}

MarvinMultipleSgroup::MarvinMultipleSgroup(MarvinMolBase *parentInit,
                                           ptree &molTree) {
  PRECONDITION(parentInit != nullptr, "parentInit cannot be null");
  parent = parentInit;
  this->id = molTree.get<std::string>("<xmlattr>.id", "");
  if (this->id == "") {
    throw FileParseException("Expected a iD in MRV file");
  }
  this->molID = molTree.get<std::string>("<xmlattr>.molID", "");
  if (this->molID == "") {
    throw FileParseException("Expected a molID in MRV file");
  }

  // first make sure this can be made into a proper MultipleSgroup.
  // To be valid, the title must be a positive integer and there must be
  // exactly two bonds the connect to the repeating group

  std::string atomRefsStr = molTree.get<std::string>("<xmlattr>.atomRefs", "");
  if (atomRefsStr == "") {
    throw FileParseException(
        "Expected  atomRefs for a MultipleSgroup definition in MRV file");
  }

  std::vector<std::string> atomList;
  boost::algorithm::split(atomList, atomRefsStr, boost::algorithm::is_space());
  for (auto &it : atomList) {
    auto atomPtr = this->parent->findAtomByRef(it);

    if (atomPtr == nullptr) {
      throw FileParseException(
          "AtomRef specification for a MultipleSgroup definition was not found in any parent");
    }

    this->atoms.push_back(atomPtr);
    this->parentAtoms.push_back(atomPtr);
  }

  this->title = molTree.get<std::string>("<xmlattr>.title", "");
  int testInt;
  if (this->title == "" || !getCleanNumber(this->title, testInt) ||
      testInt <= 0) {
    throw FileParseException(
        "Expected a positive integer title for a MultipleSgroup definition in MRV file");
  }
}

std::string MarvinMultipleSgroup::role() const { return "MultipleSgroup"; }

bool MarvinMultipleSgroup::hasAtomBondBlocks() const { return false; }

MarvinMolBase *MarvinMultipleSgroup::copyMol(
    const std::string &idAppendage) const {
  auto outSgroup = new MarvinMultipleSgroup(this->parent);
  this->parent->sgroups.emplace_back(outSgroup);

  outSgroup->molID = this->molID + idAppendage;
  outSgroup->id = this->id + idAppendage;
  outSgroup->title = this->title;

  // the only time this is to be called is when a multiple group above it is
  // being expanded and the new atoms and bonds have already been made
  //  here we just need to find them and add refs to the new multiple sgroup

  for (auto atomToCopy : atoms) {
    auto atomIndex = this->parent->getAtomIndex(atomToCopy->id + idAppendage);
    if (atomIndex < 0) {
      throw FileParseException(
          "unexpected error - new atom not found in parent");
    }

    outSgroup->atoms.push_back(this->parent->atoms[atomIndex]);
  }

  for (auto atomToCopy : parentAtoms) {
    auto atomIndex = this->parent->getAtomIndex(atomToCopy->id + idAppendage);
    if (atomIndex < 0) {
      throw FileParseException(
          "unexpected error - new atom not found in parent");
    }

    outSgroup->parentAtoms.push_back(this->parent->atoms[atomIndex]);
  }

  for (auto bondToCopy : bonds) {
    auto bondIndex = this->parent->getBondIndex(bondToCopy->id + idAppendage);
    if (bondIndex < 0) {
      throw FileParseException(
          "unexpected error - new bond not found in parent");
    }

    outSgroup->bonds.push_back(this->parent->bonds[bondIndex]);
  }

  outSgroup->isExpanded = isExpanded;

  return outSgroup;
}

std::string MarvinMultipleSgroup::toString() const {
  std::ostringstream out;

  out << "<molecule molID=\"" << molID << "\" id=\"" << id
      << "\" role=\"MultipleSgroup\" atomRefs=\""
      << boost::algorithm::join(getAtomList(), " ") << "\" title=\"" << title
      << "\">";

  for (auto &sgroup : sgroups) {
    out << sgroup->toString();
  }

  out << "</molecule>";

  return out.str();
}

ptree MarvinMultipleSgroup::toPtree() const {
  ptree out = MarvinMolBase::toPtree();

  out.put("<xmlattr>.atomRefs", boost::algorithm::join(getAtomList(), " "));
  out.put("<xmlattr>.title", title);

  MarvinMolBase::addSgroupsToPtree(out);

  return out;
}

MarvinMulticenterSgroup::MarvinMulticenterSgroup(MarvinMolBase *parentInit) {
  PRECONDITION(parentInit != nullptr, "parentInit cannot be null");
  parent = parentInit;
}

MarvinMulticenterSgroup::MarvinMulticenterSgroup(MarvinMolBase *parentInit,
                                                 ptree &molTree) {
  PRECONDITION(parentInit != nullptr, "parentInit cannot be null");
  parent = parentInit;
  this->id = molTree.get<std::string>("<xmlattr>.id", "");
  if (this->id == "") {
    throw FileParseException("Expected a iD in MRV file");
  }
  this->molID = molTree.get<std::string>("<xmlattr>.molID", "");
  if (this->molID == "") {
    throw FileParseException("Expected a molID in MRV file");
  }

  std::string atomRefsStr = molTree.get<std::string>("<xmlattr>.atomRefs", "");
  if (atomRefsStr == "") {
    throw FileParseException(
        "Expected  atomRefs for a MulticenterSgroup definition in MRV file");
  }

  std::vector<std::string> atomList;
  boost::algorithm::split(atomList, atomRefsStr, boost::algorithm::is_space());
  for (auto &it : atomList) {
    auto atomPtr = this->parent->findAtomByRef(it);

    if (atomPtr == nullptr) {
      throw FileParseException(
          "AtomRef specification for a MulticenterSgroup definition was not found in any parent");
    }

    this->atoms.push_back(atomPtr);
  }

  std::string centerId = molTree.get<std::string>("<xmlattr>.center", "");

  auto atomPtr = this->parent->findAtomByRef(centerId);

  if (atomPtr == nullptr) {
    throw FileParseException(
        "Center specification for a MulticenterSgroup definition was not found in any parent");
  }

  this->center = atomPtr;
}

std::string MarvinMulticenterSgroup::role() const {
  return "MulticenterSgroup";
}

bool MarvinMulticenterSgroup::hasAtomBondBlocks() const { return false; }

MarvinMolBase *MarvinMulticenterSgroup::copyMol(
    const std::string &idAppendage) const {
  auto outSgroup = new MarvinMulticenterSgroup(this->parent);
  this->parent->sgroups.emplace_back(outSgroup);

  outSgroup->molID = this->molID + idAppendage;
  outSgroup->id = this->id;

  // the only time this is to be called is when a multiple group above it is
  // being expanded and the new atoms and bonds have already been made
  //  here we just need to find them and add refs to the new multiple sgroup

  for (auto atomToCopy : atoms) {
    auto atomIndex = this->parent->getAtomIndex(atomToCopy->id + idAppendage);
    if (atomIndex < 0) {
      throw FileParseException(
          "unexpected error - new atom not found in parent");
    }

    outSgroup->atoms.push_back(this->parent->atoms[atomIndex]);
  }

  auto centerIndex = this->parent->getAtomIndex(this->center->id + idAppendage);
  if (centerIndex < 0) {
    throw FileParseException("unexpected error - new atom not found in parent");
  }
  outSgroup->center = this->parent->atoms[centerIndex];

  return outSgroup;
}

std::string MarvinMulticenterSgroup::toString() const {
  // <molecule molID="m2" id="sg1" role="MulticenterSgroup" atomRefs="a2 a6 a5
  // a4 a3" center="a18"/>
  std::ostringstream out;

  out << "<molecule molID=\"" << molID << "\" id=\"" << id
      << "\" role=\"MulticenterSgroup\" atomRefs=\""
      << boost::algorithm::join(getAtomList(), " ") << "\" center=\""
      << center->id << "\"/>";

  return out.str();
}

ptree MarvinMulticenterSgroup::toPtree() const {
  ptree out = MarvinMolBase::toPtree();

  out.put("<xmlattr>.atomRefs", boost::algorithm::join(getAtomList(), " "));
  out.put("<xmlattr>.center", center->id);

  return out;
}

MarvinGenericSgroup::MarvinGenericSgroup(MarvinMolBase *parentInit) {
  PRECONDITION(parentInit != nullptr, "parentInit cannot be null");
  parent = parentInit;
}

MarvinGenericSgroup::MarvinGenericSgroup(MarvinMolBase *parentInit,
                                         ptree &molTree) {
  PRECONDITION(parentInit != nullptr, "parentInit cannot be null");
  parent = parentInit;
  this->id = molTree.get<std::string>("<xmlattr>.id", "");
  if (this->id == "") {
    throw FileParseException("Expected a iD in MRV file");
  }
  this->molID = molTree.get<std::string>("<xmlattr>.molID", "");
  if (this->molID == "") {
    throw FileParseException("Expected a molID in MRV file");
  }

  std::string atomRefsStr = molTree.get<std::string>("<xmlattr>.atomRefs", "");
  if (atomRefsStr == "") {
    throw FileParseException(
        "Expected  atomRefs for a GenericSgroup definition in MRV file");
  }

  std::vector<std::string> atomList;
  boost::algorithm::split(atomList, atomRefsStr, boost::algorithm::is_space());
  for (auto &it : atomList) {
    auto atomPtr = this->parent->findAtomByRef(it);

    if (atomPtr == nullptr) {
      throw FileParseException(
          "AtomRef specification for a GenericSgroup definition was not found in any parent");
    }

    this->atoms.push_back(atomPtr);
  }

  this->charge = molTree.get<std::string>("<xmlattr>.charge", "");
  if (this->charge != "onAtoms" && this->charge != "onBracket") {
    throw FileParseException(
        "Expected  omAtoms or onBracket for a charge attr of a GenericSgroup definition in MRV file");
  }
}

std::string MarvinGenericSgroup::role() const { return "GenericSgroup"; }

bool MarvinGenericSgroup::hasAtomBondBlocks() const { return false; }

MarvinMolBase *MarvinGenericSgroup::copyMol(
    const std::string &idAppendage) const {
  auto outSgroup = new MarvinGenericSgroup(this->parent);
  this->parent->sgroups.emplace_back(outSgroup);

  outSgroup->molID = this->molID + idAppendage;
  outSgroup->id = this->id + idAppendage;

  // the only time this is to be called is when a multiple group above it is
  // being expanded and the new atoms and bonds have already been made
  //  here we just need to find them and add refs to the new multiple sgroup

  for (auto atomToCopy : atoms) {
    auto atomIndex = this->parent->getAtomIndex(atomToCopy->id + idAppendage);
    if (atomIndex < 0) {
      throw FileParseException(
          "unexpected error - new atom not found in parent");
    }

    outSgroup->atoms.push_back(this->parent->atoms[atomIndex]);
  }

  return outSgroup;
}

std::string MarvinGenericSgroup::toString() const {
  // <molecule molID="m2" id="sg1" role="MulticenterSgroup" atomRefs="a2 a6 a5
  // a4 a3" center="a18"/>
  std::ostringstream out;

  out << "<molecule molID=\"" << molID << "\" id=\"" << id
      << "\" role=\"GenericSgroup\" atomRefs=\""
      << boost::algorithm::join(getAtomList(), " ") << "\"/>";

  return out.str();
}

ptree MarvinGenericSgroup::toPtree() const {
  ptree out = MarvinMolBase::toPtree();

  out.put("<xmlattr>.atomRefs", boost::algorithm::join(getAtomList(), " "));

  return out;
}

MarvinMonomerSgroup::MarvinMonomerSgroup(MarvinMolBase *parentInit) {
  PRECONDITION(parentInit != nullptr, "parentInit cannot be null");
  parent = parentInit;
}

MarvinMonomerSgroup::MarvinMonomerSgroup(MarvinMolBase *parentInit,
                                         ptree &molTree) {
  PRECONDITION(parentInit != nullptr, "parentInit cannot be null");
  parent = parentInit;
  this->id = molTree.get<std::string>("<xmlattr>.id", "");
  if (this->id == "") {
    throw FileParseException("Expected a iD in MRV file");
  }
  this->molID = molTree.get<std::string>("<xmlattr>.molID", "");
  if (this->molID == "") {
    throw FileParseException("Expected a molID in MRV file");
  }

  std::string atomRefsStr = molTree.get<std::string>("<xmlattr>.atomRefs", "");
  if (atomRefsStr == "") {
    throw FileParseException(
        "Expected  atomRefs for a MonomerSgroup definition in MRV file");
  }

  std::vector<std::string> atomList;
  boost::algorithm::split(atomList, atomRefsStr, boost::algorithm::is_space());
  for (auto &it : atomList) {
    auto atomPtr = this->parent->findAtomByRef(it);

    if (atomPtr == nullptr) {
      throw FileParseException(
          "AtomRef specification for a MonomerSgroup definition was not found in any parent");
    }

    this->atoms.push_back(atomPtr);
  }

  this->title = molTree.get<std::string>("<xmlattr>.title", "");
  if (this->title == "") {
    throw FileParseException(
        "Expected  title for a MonomerSgroup definition in MRV file");
  }
  this->charge = molTree.get<std::string>("<xmlattr>.charge", "");
  if (this->charge == "") {
    throw FileParseException(
        "Expected  omAtoms or onBracket for a charge attr of a MonomerSgroup definition in MRV file");
  }
}

std::string MarvinMonomerSgroup::role() const { return "MonomerSgroup"; }

bool MarvinMonomerSgroup::hasAtomBondBlocks() const { return false; }

MarvinMolBase *MarvinMonomerSgroup::copyMol(
    const std::string &idAppendage) const {
  auto outSgroup = new MarvinMonomerSgroup(this->parent);
  this->parent->sgroups.emplace_back(outSgroup);

  outSgroup->molID = this->molID + idAppendage;
  outSgroup->id = this->id + idAppendage;
  outSgroup->title = this->title;
  outSgroup->charge = this->charge;

  // the only time this is to be called is when a multiple group above it is
  // being expanded and the new atoms and bonds have already been made
  //  here we just need to find them and add refs to the new multiple sgroup

  for (auto atomToCopy : atoms) {
    auto atomIndex = this->parent->getAtomIndex(atomToCopy->id + idAppendage);
    if (atomIndex < 0) {
      throw FileParseException(
          "unexpected error - new atom not found in parent");
    }

    outSgroup->atoms.push_back(this->parent->atoms[atomIndex]);
  }

  return outSgroup;
}

std::string MarvinMonomerSgroup::toString() const {
  // <molecule id="sg1" role="MonomerSgroup" title="mon" charge="onAtoms"
  // molID="m2" atomRefs="a2 a1 a3 a4">
  // </molecule>

  std::ostringstream out;

  out << "<molecule molID=\"" << molID << "\" id=\"" << id
      << "\" role=\"MonomerSgroup\" atomRefs=\""
      << boost::algorithm::join(getAtomList(), " ") << "\" title=\"" << title
      << "\" charge=\"" << charge << "\">"
      << "</molecule>";

  return out.str();
}

ptree MarvinMonomerSgroup::toPtree() const {
  ptree out = MarvinMolBase::toPtree();

  out.put("<xmlattr>.atomRefs", boost::algorithm::join(getAtomList(), " "));
  out.put("<xmlattr>.title", title);
  out.put("<xmlattr>.charge", charge);

  return out;
}

MarvinSuperatomSgroupExpanded::MarvinSuperatomSgroupExpanded(
    MarvinMolBase *parentInit) {
  PRECONDITION(parentInit != nullptr, "parentInit cannot be null");
  parent = parentInit;
}

MarvinSuperatomSgroupExpanded::MarvinSuperatomSgroupExpanded(
    MarvinMolBase *parentInit, ptree &molTree) {
  PRECONDITION(parentInit != nullptr, "parentInit cannot be null");
  parent = parentInit;
  this->id = molTree.get<std::string>("<xmlattr>.id", "");
  if (this->id == "") {
    throw FileParseException("Expected a iD in MRV file");
  }
  this->molID = molTree.get<std::string>("<xmlattr>.molID", "");
  if (this->molID == "") {
    throw FileParseException("Expected a molID in MRV file");
  }

  std::string atomRefs = molTree.get<std::string>("<xmlattr>.atomRefs", "");
  if (atomRefs != "") {
    std::vector<std::string> atomList;
    boost::algorithm::split(atomList, atomRefs, boost::algorithm::is_space());
    for (auto &it : atomList) {
      auto atomPtr = this->parent->findAtomByRef(it);

      if (atomPtr == nullptr) {
        throw FileParseException(
            "AtomRef specification for an SuperatomSgroupExpanded group definition was not found in any parent");
      }
      this->atoms.push_back(atomPtr);
    }
  } else  // must have no atomRefs nor AtomArray - get the atoms from
          // the parent block - the ones that reference this
          // superSgroup
  {
    for (MarvinAtom *atom : parent->atoms) {
      if (atom->sgroupRef == this->id) {
        this->atoms.push_back(atom);
      }
    }
  }

  this->title = molTree.get<std::string>("<xmlattr>.title", "");
  if (this->title == "") {
    throw FileParseException(
        "Expected  title for a SuperatomSgroupExpanded definition in MRV file");
  }
}

MarvinSuperatomSgroupExpanded::~MarvinSuperatomSgroupExpanded() {}

MarvinMolBase *MarvinSuperatomSgroupExpanded::copyMol(
    const std::string &idAppendage) const {
  auto outSgroup = new MarvinSuperatomSgroupExpanded(this->parent);
  this->parent->sgroups.emplace_back(outSgroup);

  outSgroup->molID = this->molID + idAppendage;
  outSgroup->id = this->id + idAppendage;
  outSgroup->title = this->title;

  auto actualParent = getActualParent(this);

  // the only time this is to be called is when a multiple group above it is
  // being expanded and the new atoms and bonds have already been made
  //  here we just need to find them and add refs to the new multiple sgroup

  for (auto atomToCopy : atoms) {
    auto atomIndex = actualParent->getAtomIndex(atomToCopy->id + idAppendage);
    if (atomIndex < 0) {
      throw FileParseException(
          "unexpected error - new atom not found in parent");
    }

    outSgroup->atoms.push_back(actualParent->atoms[atomIndex]);
  }

  for (auto bondToCopy : bonds) {
    auto bondIndex = actualParent->getBondIndex(bondToCopy->id + idAppendage);
    if (bondIndex < 0) {
      throw FileParseException(
          "unexpected error - new bond not found in parent");
    }

    outSgroup->bonds.push_back(actualParent->bonds[bondIndex]);
  }

  return outSgroup;
}

std::string MarvinSuperatomSgroupExpanded::toString() const {
  std::ostringstream out;

  out << "<molecule molID=\"" << molID << "\" id=\"" << id
      << "\" role=\"SuperatomSgroup\" atomRefs=\""
      << boost::algorithm::join(getAtomList(), " ") << "\" title=\"" << title
      << "\"/>";

  return out.str();
}

ptree MarvinSuperatomSgroupExpanded::toPtree() const {
  ptree out = MarvinMolBase::toPtree();

  out.put("<xmlattr>.atomRefs", boost::algorithm::join(getAtomList(), " "));
  out.put("<xmlattr>.title", title);
  return out;
}

std::string MarvinSuperatomSgroupExpanded::role() const {
  return "SuperatomSgroupExpanded";
}

bool MarvinSuperatomSgroupExpanded::hasAtomBondBlocks() const { return false; }

MarvinSuperatomSgroup::MarvinSuperatomSgroup(MarvinMolBase *parentInit) {
  PRECONDITION(parentInit != nullptr, "parentInit cannot be null");
  parent = parentInit;
}

MarvinSuperatomSgroup::MarvinSuperatomSgroup(MarvinMolBase *parentInit,
                                             ptree &molTree) {
  PRECONDITION(parentInit != nullptr, "parentInit cannot be null");
  parent = parentInit;
  this->id = molTree.get<std::string>("<xmlattr>.id", "");
  if (this->id == "") {
    throw FileParseException("Expected a iD in MRV file");
  }
  this->molID = molTree.get<std::string>("<xmlattr>.molID", "");
  if (this->molID == "") {
    throw FileParseException("Expected a molID in MRV file");
  }

  this->title = molTree.get<std::string>("<xmlattr>.title", "");
  if (this->title == "") {
    throw FileParseException(
        "Expected  title for a SuperatomSgroup definition in MRV file");
  }

  this->parseAtomsAndBonds(molTree);

  bool found;

  try {
    boost::property_tree::ptree AttachmentPointArrayTree =
        molTree.get_child("AttachmentPointArray");
    found = true;
  } catch (const std::exception &e) {
    found = false;
  }

  if (found) {
    for (auto &v : molTree.get_child("AttachmentPointArray")) {
      std::string bondId = v.second.get<std::string>("<xmlattr>.bond", "");
      if (bondId == "") {  // this can happen if the attachment point is
                           // not actually used - as in Amino acids that
                           // have non-used crosslink atoms
        continue;
      }

      auto &marvinAttachmentPoint =
          this->attachmentPoints.emplace_back(new MarvinAttachmentPoint);

      marvinAttachmentPoint->atom =
          v.second.get<std::string>("<xmlattr>.atom", "");
      marvinAttachmentPoint->order =
          v.second.get<std::string>("<xmlattr>.order", "");
      marvinAttachmentPoint->bond = bondId;

      if (marvinAttachmentPoint->atom == "" ||
          marvinAttachmentPoint->order == "") {
        throw FileParseException(
            "Expected atom, order and bond, for an AttachmentPoint definition in MRV file");
      }

      if (!boost::algorithm::contains(
              this->atoms,
              std::vector<std::string>{marvinAttachmentPoint->atom},
              MarvinMol::atomRefInAtoms)) {
        throw FileParseException(
            "Atom specification for an AttachmentPoint definition must be in the parent's atom array in MRV file");
      }

      // bond must be found in the  bonds vector

      if (this->parent->findBondByRef(marvinAttachmentPoint->bond) == nullptr) {
        throw FileParseException(
            "Bond specification for an AttachmentPoint definition must be in the bond array in MRV file");
      }

      // the order must be an integer

      int orderInt;
      if (!getCleanNumber(marvinAttachmentPoint->order, orderInt)) {
        throw FileParseException(
            "Order for an AttachmentPoint definition must be an integer in MRV file");
      }
    }
  }
}

MarvinSuperatomSgroup::~MarvinSuperatomSgroup() {
  this->attachmentPoints.clear();
}

std::string MarvinSuperatomSgroup::role() const { return "SuperatomSgroup"; }

bool MarvinSuperatomSgroup::hasAtomBondBlocks() const { return true; }

MarvinMolBase *MarvinSuperatomSgroup::copyMol(const std::string &) const {
  PRECONDITION(0, "Internal error:  copying a SuperatomSgroup");
}

std::string MarvinSuperatomSgroup::toString() const {
  std::ostringstream out;

  out << "<molecule molID=\"" << molID << "\" id=\"" << id
      << "\" role=\"SuperatomSgroup\" title=\"" << title << "\">";

  out << "<atomArray>";
  for (auto atom : atoms) {
    out << atom->toString();
  }
  out << "</atomArray>";

  out << "<bondArray>";
  for (auto bond : bonds) {
    out << bond->toString();
  }
  out << "</bondArray>";

  if (attachmentPoints.size() > 0) {
    out << "<AttachmentPointArray>";
    for (auto &attachmentPoint : attachmentPoints) {
      out << attachmentPoint->toString();
    }
    out << "</AttachmentPointArray>";
  }

  for (auto &sgroup : sgroups) {
    out << sgroup->toString();
  }

  out << "</molecule>";

  return out.str();
}

ptree MarvinSuperatomSgroup::toPtree() const {
  ptree out = MarvinMolBase::toPtree();

  out.put("<xmlattr>.title", title);

  if (attachmentPoints.size() > 0) {
    ptree attachmentpointArray;
    for (auto &attachmentPoint : attachmentPoints) {
      attachmentpointArray.add_child("attachmentPoint",
                                     attachmentPoint->toPtree());
    }
    out.put_child("AttachmentPointArray", attachmentpointArray);
  }
  MarvinMolBase::addSgroupsToPtree(out);

  return out;
}

MarvinMol::MarvinMol() { parent = nullptr; }
MarvinMol::MarvinMol(ptree &molTree) {
  parent = nullptr;

  // get atoms if this mol is supposed to have them

  this->parseAtomsAndBonds(molTree);
}

MarvinMol::~MarvinMol() {}

std::string MarvinMol::role() const { return "MarvinMol"; }

bool MarvinMol::hasAtomBondBlocks() const { return true; }

void MarvinMol::pushOwnedAtom(MarvinAtom *atom) {
  PRECONDITION(atom, "bad atom");
  ownedAtoms.emplace_back(atom);
}

void MarvinMol::pushOwnedBond(MarvinBond *bond) {
  PRECONDITION(bond, "bad bond");
  ownedBonds.emplace_back(bond);
}

void MarvinMol::removeOwnedAtom(MarvinAtom *atom) {
  PRECONDITION(atom != nullptr, "atom cannot be null");
  eraseUniquePtr<MarvinAtom>(ownedAtoms, atom);
}

void MarvinMol::removeOwnedBond(MarvinBond *bond) {
  PRECONDITION(bond != nullptr, "bond cannot be null");
  eraseUniquePtr<MarvinBond>(ownedBonds, bond);
}

bool MarvinMolBase::atomRefInAtoms(MarvinAtom *a, std::string b) {
  PRECONDITION(a != nullptr, "a cannot be null");
  return a->id == b;
}

bool MarvinMolBase::molIDInSgroups(std::string a, std::string b) {
  return a == b;
}

bool MarvinMolBase::bondRefInBonds(MarvinBond *a, std::string b) {
  PRECONDITION(a != nullptr, "a cannot be null");
  return a->id == b;
}

void MarvinMolBase::cleanUpSgNumbering(
    int &sgCount, std::map<std::string, std::string> &sgMap) {
  for (auto &sgroup : this->sgroups) {
    std::string newId = "sg" + std::to_string(++sgCount);
    sgMap[sgroup->id] = newId;
    sgroup->id = newId;

    sgroup->cleanUpSgNumbering(sgCount, sgMap);
  }
}

void MarvinSuperatomSgroup::cleanUpNumberingMolsAtomsBonds(
    int &molCount,   // this is the starting mol count, and receives the ending
                     // mol count - THis is used when
                     // MarvinMol->convertToSuperAtoms is called multiple
                     // times from a RXN
    int &atomCount,  // starting and ending atom count
    int &bondCount,  // starting and ending bond count
    std::map<std::string, std::string> &sgMap,
    std::map<std::string, std::string> &atomMap,
    std::map<std::string, std::string> &bondMap) {
  // the common part for all classes
  MarvinMolBase::cleanUpNumberingMolsAtomsBonds(molCount, atomCount, bondCount,
                                                sgMap, atomMap, bondMap);

  // the specific part for this class
  auto marvinSuperatomSgroup = (MarvinSuperatomSgroup *)this;
  for (auto &attachmentPoint : marvinSuperatomSgroup->attachmentPoints) {
    attachmentPoint->atom = atomMap[attachmentPoint->atom];
    attachmentPoint->bond =
        bondMap[attachmentPoint->bond];  // bond is actually in the parent
  }
}

void MarvinMolBase::cleanUpNumberingMolsAtomsBonds(
    int &molCount,   // this is the starting mol count, and receives the ending
                     // mol count - THis is used when
                     // MarvinMol->convertToSuperAtoms is called multiple
                     // times from a RXN
    int &atomCount,  // starting and ending atom count
    int &bondCount,  // starting and ending bond count)
    std::map<std::string, std::string> &sgMap,
    std::map<std::string, std::string> &atomMap,
    std::map<std::string, std::string> &bondMap) {
  // clean up the mol ids, the atomIds and bondIds.  make  map of the old to
  // new atom ids and bond ids

  this->molID = "m" + std::to_string(++molCount);

  if (this->hasAtomBondBlocks()) {
    for (auto atomPtr : this->atoms) {
      std::string newId = "a" + std::to_string(++atomCount);
      atomMap[atomPtr->id] = newId;
      atomPtr->id = newId;

      // fix the sgroupRef
      if (atomPtr->sgroupRef != "") {
        atomPtr->sgroupRef = sgMap[atomPtr->sgroupRef];
      }
    }

    for (auto bondPtr : this->bonds) {
      std::string newId = "b" + std::to_string(++bondCount);
      bondMap[bondPtr->id] = newId;
      bondPtr->id = newId;

      // fix the bond's references to atoms

      bondPtr->atomRefs2[0] = atomMap[bondPtr->atomRefs2[0]];
      bondPtr->atomRefs2[1] = atomMap[bondPtr->atomRefs2[1]];
    }
  }

  // Now the sgroups

  for (auto &sgroup : sgroups) {
    sgroup->cleanUpNumberingMolsAtomsBonds(molCount, atomCount, bondCount,
                                           sgMap, atomMap, bondMap);
  }
}

void MarvinMolBase::cleanUpNumbering(
    int &molCount  // this is the starting mol count, and receives the ending
                   // mol count - THis is used when
                   // MarvinMol->convertToSuperAtoms is called multiple times
                   // from a RXN
    ,
    int &atomCount  // starting and ending atom count
    ,
    int &bondCount  // starting and ending bond count
    ,
    int &sgCount  // starting and ending sg  count)
    ,
    std::map<std::string, std::string>
        &sgMap  // map from old sg number to new sg number
    ,
    std::map<std::string, std::string>
        &atomMap  // map from old atom number to new atom number
    ,
    std::map<std::string, std::string>
        &bondMap  // map from old bond number to new bond number
)

{
  this->cleanUpSgNumbering(sgCount, sgMap);
  this->cleanUpNumberingMolsAtomsBonds(molCount, atomCount, bondCount, sgMap,
                                       atomMap, bondMap);
}

bool MarvinMolBase::AnyOverLappingAtoms(const MarvinMolBase *otherMol) const {
  PRECONDITION(otherMol != nullptr, "otherMol cannot be null");
  std::vector<std::string> firstList = getAtomList();
  std::sort(firstList.begin(), firstList.end());
  std::vector<std::string> secondList = otherMol->getAtomList();
  std::sort(secondList.begin(), secondList.end());
  std::vector<std::string> intersection;
  std::set_intersection(firstList.begin(), firstList.end(), secondList.begin(),
                        secondList.end(), back_inserter(intersection));

  return intersection.size() > 0;
}

// the following routine determines if a sgroup role is passive when expanding
// in preparation for creating an RDKit mol.  sgroups are NOT passive if they
// create atoms or bonds in the parent when expanding

bool MarvinMolBase::isPassiveRoleForExpansion() const {
  // some derived classes override this and return false
  return true;
}

bool MarvinSuperatomSgroup::isPassiveRoleForExpansion() const { return false; }

bool MarvinMultipleSgroup::isPassiveRoleForExpansion() const { return false; }

// this routine determines if the sgroup is passive for contraction, whern
// coming from an RDKit mol to a Marvin Mol
//  It is NOT passive if atoms in the parent are removed by contracting it.

bool MarvinMolBase::isPassiveRoleForContraction() const {
  // some derived classes override this and return false

  return true;
}

bool MarvinMol::isPassiveRoleForContraction() const { return false; }

bool MarvinSuperatomSgroupExpanded::isPassiveRoleForContraction() const {
  return false;
}

bool MarvinMultipleSgroup::isPassiveRoleForContraction() const { return false; }

void MarvinMolBase::parseMoleculeSpecific(RDKit::RWMol *mol,
                                          std::unique_ptr<SubstanceGroup> &,
                                          int) {
  PRECONDITION(mol != nullptr, "mol cannot be null");
  // default is nothing - most derived classes override this
};

void MarvinSuperatomSgroupExpanded::parseMoleculeSpecific(
    RDKit::RWMol *mol, std::unique_ptr<SubstanceGroup> &sgroup,
    int sequenceId) {
  PRECONDITION(mol != nullptr, "mol cannot be null");

  std::string typ = "SUP";
  sgroup.reset(new SubstanceGroup(mol, typ));
  sgroup->setProp<unsigned int>("index", sequenceId);

  for (auto atomPtr : this->atoms) {
    sgroup->addAtomWithIdx(this->parent->getAtomIndex(atomPtr->id));
  }

  for (auto bondPtr : this->bonds) {
    sgroup->addBondWithIdx(this->parent->getBondIndex(bondPtr->id));
  }

  sgroup->setProp("LABEL", this->title);
}
void MarvinMultipleSgroup::parseMoleculeSpecific(
    RDKit::RWMol *mol, std::unique_ptr<SubstanceGroup> &sgroup,
    int sequenceId) {
  PRECONDITION(mol != nullptr, "mol cannot be null");

  // now the MultipleSgroups
  // note: sequence continues counting from the loop above

  std::string typ = "MUL";
  sgroup.reset(new SubstanceGroup(mol, typ));
  sgroup->setProp<unsigned int>("index", sequenceId);

  for (auto atomPtr : this->atoms) {
    sgroup->addAtomWithIdx(this->parent->getAtomIndex(atomPtr->id));
  }

  for (auto atomPtr : this->parentAtoms) {
    sgroup->addParentAtomWithIdx(this->parent->getAtomIndex(atomPtr->id));
  }

  // the connection bonds

  for (auto bondptr : this->bondsToAtomsNotInExpandedGroup) {
    sgroup->addBondWithIdx(this->parent->getBondIndex(bondptr->id));
  }

  sgroup->setProp("MULT", this->title);
}
void MarvinSruCoModSgroup::parseMoleculeSpecific(
    RDKit::RWMol *mol, std::unique_ptr<SubstanceGroup> &sgroup,
    int sequenceId) {
  PRECONDITION(mol != nullptr, "mol cannot be null");

  std::string typ;
  if (this->role() == "SruSgroup") {
    typ = "SRU";
  } else if (this->role() == "CopolymerSgroup") {
    typ = "COP";
  } else if (this->role() == "ModificationSgroup") {
    typ = "MOD";
  } else {
    throw FileParseException(
        "Internal error: unrecognized role in a MarvinSruCoPolSgroup");
  }

  sgroup.reset(new SubstanceGroup(mol, typ));
  sgroup->setProp<unsigned int>("index", sequenceId);

  sgroup->setProp("CONNECT", this->connect);

  for (auto atomPtr : this->atoms) {
    sgroup->addAtomWithIdx(this->parent->getAtomIndex(atomPtr->id));
  }

  for (auto bondPtr : this->bonds) {
    sgroup->addBondWithIdx(this->parent->getBondIndex(bondPtr->id));
  }

  sgroup->setProp("LABEL", this->title);
}

void MarvinDataSgroup::parseMoleculeSpecific(
    RDKit::RWMol *mol, std::unique_ptr<SubstanceGroup> &sgroup,
    int sequenceId) {
  PRECONDITION(mol != nullptr, "mol cannot be null");
  // Now the data groups

  std::string typ = "DAT";
  sgroup.reset(new SubstanceGroup(mol, typ));
  sgroup->setProp<unsigned int>("index", sequenceId);

  for (auto atomPtr : this->atoms) {
    sgroup->addAtomWithIdx(this->parent->getAtomIndex(atomPtr->id));
  }

  sgroup->setProp("FIELDNAME", this->fieldName);

  if (this->queryType != "") {
    sgroup->setProp("QUERYTYPE", this->queryType);
  }
  if (this->queryOp != "") {
    sgroup->setProp("QUERYOP", this->queryOp);
  }

  std::ostringstream out;
  out << std::fixed << std::setw(10) << std::setprecision(4) << this->x
      << std::fixed << std::setw(10) << std::setprecision(4) << this->y
      << "    DRU   ALL  0       0";

  sgroup->setProp("FIELDDISP", out.str());  // really not used by RDKIT

  std::vector<std::string> fieldDatas;
  fieldDatas.push_back(this->fieldData);
  sgroup->setProp("DATAFIELDS", fieldDatas);

  // The following props are not part of the RDKit structure for MOL
  // files, but we save them so that we can round-trip the MRV

  sgroup->setProp("UNITS", this->units);
  sgroup->setProp("UNITSDISPLAYED", this->unitsDisplayed);
  sgroup->setProp("CONTEXT", this->context);
  sgroup->setProp("PLACEMENT", this->placement);
  sgroup->setProp("X", this->x);
  sgroup->setProp("Y", this->y);
}

void MarvinMulticenterSgroup::parseMoleculeSpecific(
    RDKit::RWMol *, std::unique_ptr<SubstanceGroup> &, int) {
  // the MultiCenter Sgroups

  // There is really no place to put these in RDKit.  We should have
  // removed these already
  PRECONDITION(
      0, "Internal error:   a MarvinMulticenterSgroup has not been removed");
}

// the Generic groups

void MarvinGenericSgroup::parseMoleculeSpecific(
    RDKit::RWMol *mol, std::unique_ptr<SubstanceGroup> &sgroup,
    int sequenceId) {
  PRECONDITION(mol != nullptr, "mol cannot be null");

  std::string typ = "GEN";
  sgroup.reset(new SubstanceGroup(mol, typ));
  sgroup->setProp<unsigned int>("index", sequenceId);

  for (auto atomPtr : this->atoms) {
    sgroup->addAtomWithIdx(this->parent->getAtomIndex(atomPtr->id));
  }

  // note: there is no place to put the change="onAtoms" flag
}

// now the MonomerSgroups
// note: sequence continues counting from the loop above

void MarvinMonomerSgroup::parseMoleculeSpecific(
    RDKit::RWMol *mol, std::unique_ptr<SubstanceGroup> &sgroup,
    int sequenceId) {
  PRECONDITION(mol != nullptr, "mol cannot be null");

  std::string typ = "MON";
  sgroup.reset(new SubstanceGroup(mol, typ));
  sgroup->setProp<unsigned int>("index", sequenceId);

  for (auto atomPtr : this->parent->atoms) {
    int atomIndex = this->getAtomIndex(atomPtr->id);
    sgroup->addAtomWithIdx(atomIndex);
  }
  sgroup->setProp("LABEL", this->title);

  // Note: RDKit does not have a place for the Bracket information nor
  // the charge="onAtoms" attr
}

void moveSgroup(MarvinMolBase *sgroupToMove, MarvinMolBase *newParent,
                bool removeFromOldParent = true) {
  PRECONDITION(sgroupToMove, "bad sgroupToMove pointer");

  MarvinMolBase *parent = sgroupToMove->parent;
  PRECONDITION(parent, "bad sgroupToMove parent pointer");
  auto moveUniqIter =
      find_if(sgroupToMove->parent->sgroups.begin(),
              sgroupToMove->parent->sgroups.end(),
              [sgroupToMove](std::unique_ptr<MarvinMolBase> &sgptr) {
                return sgptr->id == sgroupToMove->id;
              });

  if (moveUniqIter == sgroupToMove->parent->sgroups.end()) {
    throw FileParseException(
        "Unexpected error - sgroup not found in its parent");
  }

  newParent->sgroups.push_back(std::move(*moveUniqIter));
  if (removeFromOldParent) {
    parent->sgroups.erase(moveUniqIter);
  }
  sgroupToMove->parent = newParent;
}

void promoteChild(MarvinMolBase *molToMove) {
  PRECONDITION(molToMove, "bad molToMove pointer");
  moveSgroup(molToMove, molToMove->parent->parent);
}

void MarvinMultipleSgroup::expandOneMultipleSgroup() {
  // Mulitplesgroups are handled differently in Marvin and RDKit/mol file
  // format
  //
  // In Marvin, the atoms of the group are indicated, and the "title" contains
  // the number of replicates In mol files, the atoms of the group are
  // actually replicated in the atomBlock, and the MUL constructs indicates the
  // list of ALL such atoms and the list of the original parent atoms from
  // which the others were replicated.  It also indicate the two bonds from
  // the group atoms to atoms NOT in the group (although these could be
  // derived
  //
  //  This routine takes a MarvinMol and  prepares it for generation of an
  //  RDKit mol.  Each MultipleSgroup is expanded by:
  //    1) expanding the atom block - the new replicate sets of atoms are
  //    places just after the last atom in the original (parent) set 2)
  //    records the two bonds to atoms NOT in the expanded group - one is from
  //    a "parent" atom and one from a replicated atom. 3) bonds the parent
  //    and replicated sets to each other in a head to tail fashion.

  // make sure it is not already expanded - this should not happen

  if (this->isExpanded) {
    return;
  }

  // find the two bonds (or zero bonds) to atoms NOT in the group

  std::vector<std::string>
      originalConnectorAtoms;  // the bonds in the parent that were to outside
                               // atoms BEFORE replication
  std::string connectorAtoms[2];  // working atoms to connect the replicates
  std::string outsideConnectorAtom;
  std::vector<MarvinBond *> bondsInGroup;
  this->bondsToAtomsNotInExpandedGroup.clear();

  // find the actual parent - multiple groups can be nested, so the atoms and
  // bonds of the lower one are really in the grandparent or higher)

  auto actualParent = getActualParent(this);

  for (MarvinBond *bondPtr : actualParent->bonds) {
    bool atom1InSet = boost::algorithm::contains(
        this->atoms, std::vector<std::string>{bondPtr->atomRefs2[0]},
        atomRefInAtoms);
    bool atom2InSet = boost::algorithm::contains(
        this->atoms, std::vector<std::string>{bondPtr->atomRefs2[1]},
        atomRefInAtoms);
    if (atom1InSet && atom2InSet) {
      bondsInGroup.push_back(bondPtr);
    } else if (atom1InSet) {
      originalConnectorAtoms.push_back(bondPtr->atomRefs2[0]);
      this->bondsToAtomsNotInExpandedGroup.push_back(bondPtr);
      outsideConnectorAtom =
          bondPtr
              ->atomRefs2[1];  // last one set wins - when done should be the
                               // outside atom for the 2nd bond to the outside
    } else if (atom2InSet) {
      originalConnectorAtoms.push_back(bondPtr->atomRefs2[1]);
      this->bondsToAtomsNotInExpandedGroup.push_back(bondPtr);
      outsideConnectorAtom =
          bondPtr
              ->atomRefs2[0];  // last one set wins - when done should be the
                               // outside atom for the 2nd bond to the outside
    }
  }

  // should be two bonds to the outside (or zero)
  if (this->bondsToAtomsNotInExpandedGroup.size() != 2 &&
      this->bondsToAtomsNotInExpandedGroup.size() != 0) {
    throw FileParseException(
        "Error - there must be zero or two bonds from the group to atoms not in the group");
  }

  // find the hightest id of an atom in the Group - we will insert replicated
  // atoms after that one

  int highestIndex = 0;
  for (MarvinAtom *atomPtr : this->atoms) {
    int thisIndex = actualParent->getAtomIndex(atomPtr->id);
    if (thisIndex > highestIndex) {
      highestIndex = thisIndex;
    }
  }
  MarvinAtom *lastAtomInGroupPtr = actualParent->atoms[highestIndex];

  // ready to make copies of the group

  int copyCount = std::stoi(this->title);

  if (this->bondsToAtomsNotInExpandedGroup.size() == 2) {
    connectorAtoms[1] =
        originalConnectorAtoms[1];  // this one does NOT have an _# on it
  }

  for (int copyIndex = 2; copyIndex <= copyCount;
       ++copyIndex)  // start at 2, we already have the first copy
  {
    std::string idAppendage = "_" + this->id + ":" + std::to_string(copyIndex);
    for (auto parentAtomPtr : this->parentAtoms) {
      auto copyAtom =
          new MarvinAtom(*parentAtomPtr, parentAtomPtr->id + idAppendage);
      actualParent->pushOwnedAtom(copyAtom);
      if (parentAtomPtr->sgroupRef != "" && !copyAtom->sGroupRefIsSuperatom) {
        copyAtom->sgroupRef = parentAtomPtr->sgroupRef + idAppendage;
      }

      // copy atoms to all parents up to the actual parent (multipleSgroups in
      // multipleSgroups)

      for (auto thisParent = this->parent;; thisParent = thisParent->parent) {
        TEST_ASSERT(thisParent != nullptr);
        auto insertItr = find(thisParent->atoms.begin(),
                              thisParent->atoms.end(), lastAtomInGroupPtr);
        if (insertItr != thisParent->atoms.end()) {
          ++insertItr;  // insert after the last one
        }
        thisParent->atoms.insert(insertItr, copyAtom);

        // for multiple groups, also insert in the parentAtoms array - that
        // multiple sgroup has NOT been done yet

        if (thisParent->role() == "MultipleSgroup") {
          auto thisMultipleParent = (MarvinMultipleSgroup *)thisParent;
          insertItr =
              find(thisMultipleParent->parentAtoms.begin(),
                   thisMultipleParent->parentAtoms.end(), lastAtomInGroupPtr);
          if (insertItr != thisMultipleParent->parentAtoms.end()) {
            ++insertItr;  // insert after the last one
          }
          thisMultipleParent->parentAtoms.insert(insertItr, copyAtom);
        }

        if (thisParent == actualParent) {
          break;
        }
      }
      lastAtomInGroupPtr = copyAtom;
      this->atoms.push_back(copyAtom);
    }

    // add the bonds for this copy

    for (auto parentBondPtr : bondsInGroup) {
      auto copyBond =
          new MarvinBond(*parentBondPtr, parentBondPtr->id + idAppendage,
                         parentBondPtr->atomRefs2[0] + idAppendage,
                         parentBondPtr->atomRefs2[1] + idAppendage);
      actualParent->pushOwnedBond(copyBond);
      actualParent->bonds.push_back(copyBond);
    }

    for (int i = this->sgroups.size() - 1; i >= 0; --i) {
      auto copiedMol = this->sgroups[i]->copyMol(idAppendage);
      moveSgroup(copiedMol, this->parent);
    }

    // add the bond from the last group to the new one

    if (this->bondsToAtomsNotInExpandedGroup.size() == 2) {
      connectorAtoms[0] = originalConnectorAtoms[0] + idAppendage;
      auto connectorBond = new MarvinBond(
          *this->bondsToAtomsNotInExpandedGroup[0],
          this->bondsToAtomsNotInExpandedGroup[0]->id + idAppendage,
          connectorAtoms[0], connectorAtoms[1]);
      actualParent->pushOwnedBond(connectorBond);
      actualParent->bonds.push_back(connectorBond);

      // update connectorAtom[1] for the next pass (or the final bond to the
      // outside atom)

      connectorAtoms[1] = originalConnectorAtoms[1] + idAppendage;
    }
  }

  // fix the bond to the second outside atom to reflect that it comes from the
  // last replicate

  if (this->bondsToAtomsNotInExpandedGroup.size() == 2) {
    this->bondsToAtomsNotInExpandedGroup[1]->atomRefs2[0] =
        outsideConnectorAtom;
    this->bondsToAtomsNotInExpandedGroup[1]->atomRefs2[1] = connectorAtoms[1];
  }

  // move any sgroups below this one up to the parent - we have already made
  // copies of it

  for (int i = this->sgroups.size() - 1; i >= 0; --i) {
    promoteChild(sgroups[i].get());
  }

  this->isExpanded = true;
}

void MarvinSuperatomSgroup::convertFromOneSuperAtom() {
  // the mol-style super atoms are significanTly different than the Marvin
  // super atoms.  The mol style has all atoms and bonds in the main mol, and
  // parameter lines that indicate the name of the super atom and the atom
  // indices affected.
  //
  //  The SuperatomSgroup form can appear in contracted or expanded form in
  //  MRV blocks:
  //
  // In the contracted form, the MRV block can have one or more dummy atoms
  // with the super atom name in the main molecule, and a separate sub-mol
  // with the atoms and bonds of the superatom. It also can have one or more
  // separate records called attachmentPoints that specify which atom(s) in
  // the super atom sub-mol that replaces the dummy atom(s), and also a bond
  // pointer (in the parent mol).
  //
  // In the expanded form, all of the atoms are in the parent molecule, and
  // the sub-mol only refers to them.  The attachment points refer to the
  // atoms in the parent mol. The MultipleSgroup and SruGroup seem very much
  // alike.   In both, all atoms are in the parent mol, and the sub-mol refers
  // to a group of atoms in that parent. The SruGroup specifies the name and
  // also the connection information (head-to-tail, head-head, and unspecified).
  // The Multiple S-group specifies the report count for the group
  //
  // This routine deals with only the contracted form of the SuperatomSgroups,
  // and copies the atoms and bonds from the sub=mol to the parent mol, and
  // deletes the dummy atom form the parent mol.
  //  It also saves the info creates a replacement
  //  MarvinSuperatomSgroupExpanded item.
  try {
    // save the name of the superatom

    auto superatomSgroupExpanded =
        new MarvinSuperatomSgroupExpanded(this->parent);
    this->parent->sgroups.push_back(
        std::unique_ptr<MarvinMolBase>(superatomSgroupExpanded));

    superatomSgroupExpanded->title = this->title;
    superatomSgroupExpanded->id =
        this->id;  // the expanded sgroup will replace the contracted one
    superatomSgroupExpanded->molID =
        this->molID;  // the expanded sgroup will replace the contracted one

    // find the actual parent - multiple groups can be nested, so the atoms
    // and bonds of the lower one are really in the grandparent or higher)

    auto actualParent = getActualParent(this);

    //  remove and delete the dummy atom from the parent.

    bool coordExist = this->has2dCoords();
    auto dummyAtomIter = find_if(
        actualParent->atoms.begin(), actualParent->atoms.end(),
        [this](const MarvinAtom *arg) { return arg->sgroupRef == this->id; });
    if (dummyAtomIter == actualParent->atoms.end()) {
      throw FileParseException(
          "No contracted atom found for a superatomSgroup");
    }
    auto dummyAtomPtr = *dummyAtomIter;

    // before we move or delete any atoms from the superatomSgroup,  Make sure
    // that the  sibling sgroups are handled if a sibling contains the dummy R
    // atom, then the atoms to be promoted from this sgroup should be member
    // of the sibling too

    for (auto &sibling : parent->sgroups) {
      if (sibling.get() == this) {
        continue;  // not a sibling
      }

      auto dummyInSibling =
          find_if(sibling->atoms.begin(), sibling->atoms.end(),
                  [dummyAtomPtr](const MarvinAtom *arg) {
                    return arg == dummyAtomPtr;
                  });

      if (dummyInSibling != sibling->atoms.end()) {
        sibling->atoms.erase(dummyInSibling);

        for (auto atom : this->atoms) {
          sibling->atoms.push_back(atom);
        }
      }
    }

    // get a list of the bonds that will need to be fixed after the copy of
    // atoms and bonds

    std::vector<MarvinAttachmentPoint *> orphanedBonds;
    RDGeom::Point3D centerOfAttachmentPoints;
    RDGeom::Point3D centerOfGroup;
    for (auto orphanedBond : actualParent->bonds) {
      std::string attachedAtomId = "";
      if (orphanedBond->atomRefs2[0] == dummyAtomPtr->id) {
        attachedAtomId = orphanedBond->atomRefs2[1];
      } else if (orphanedBond->atomRefs2[1] == dummyAtomPtr->id) {
        attachedAtomId = orphanedBond->atomRefs2[0];
      }

      if (attachedAtomId != "") {
        auto attachmentPoint = find_if(
            this->attachmentPoints.begin(), this->attachmentPoints.end(),
            [orphanedBond](const std::unique_ptr<MarvinAttachmentPoint> &arg) {
              return arg->bond == orphanedBond->id;
            });
        if (attachmentPoint == this->attachmentPoints.end()) {
          throw FileParseException(
              "No attachment point found for bond to the condensed atom in a superatomSgroup");
        }
        orphanedBonds.push_back(attachmentPoint->get());

        auto subMolAtomIndex = this->getAtomIndex((*attachmentPoint)->atom);
        if (subMolAtomIndex < 0) {
          throw FileParseException(
              "Attachment atom not found in the superatomsGroup");
        }
        MarvinAtom *attachedAtom = this->atoms[subMolAtomIndex];
        if (coordExist) {
          centerOfAttachmentPoints.x += attachedAtom->x2;
          centerOfAttachmentPoints.y += attachedAtom->y2;
        }
      }
    }

    if (coordExist) {
      RDGeom::Point3D offset;

      if (orphanedBonds.size() > 0) {
        centerOfAttachmentPoints.x /= orphanedBonds.size();
        centerOfAttachmentPoints.y /= orphanedBonds.size();
        offset.x = dummyAtomPtr->x2 - centerOfAttachmentPoints.x;
        offset.y = dummyAtomPtr->y2 - centerOfAttachmentPoints.y;
      } else  // use the center of all atoms in the group
      {
        RDGeom::Point3D centerOfGroup;

        for (auto atom : this->atoms) {
          centerOfGroup.x += atom->x2;
          centerOfGroup.y += atom->y2;
        }
        centerOfGroup.x /= this->atoms.size();
        centerOfGroup.y /= this->atoms.size();
        offset.x = dummyAtomPtr->x2 - centerOfGroup.x;
        offset.y = dummyAtomPtr->y2 - centerOfGroup.y;
      }

      for (auto atom : this->atoms) {
        atom->x2 += offset.x;
        atom->y2 += offset.y;
      }
    }

    for (auto thisParent = parent;; thisParent = thisParent->parent) {
      TEST_ASSERT(thisParent != nullptr);
      auto dummyInParent =
          find_if(thisParent->atoms.begin(), thisParent->atoms.end(),
                  [dummyAtomPtr](const MarvinAtom *arg) {
                    return arg == dummyAtomPtr;
                  });
      TEST_ASSERT(dummyInParent != thisParent->atoms.end());
      thisParent->atoms.erase(dummyInParent);  // get rid of the atoms pointer
                                               // to the old dummy atom

      // for multiple groups, also delete it from the parentAtoms array - that
      // multiple sgroup has NOT been done yet

      if (thisParent->role() == "MultipleSgroup") {
        auto thisMultipleParent = (MarvinMultipleSgroup *)thisParent;

        auto deleteIter =
            find(thisMultipleParent->parentAtoms.begin(),
                 thisMultipleParent->parentAtoms.end(), dummyAtomPtr);
        TEST_ASSERT(deleteIter != thisMultipleParent->parentAtoms.end());
        thisMultipleParent->parentAtoms.erase(deleteIter);
      }

      if (thisParent == actualParent) {
        break;
      }
    }

    removeOwnedAtom(dummyAtomPtr);

    // add the atoms and bonds from the super group to the parent

    for (auto subAtomPtr : this->atoms) {
      for (auto thisParent = parent;; thisParent = thisParent->parent) {
        TEST_ASSERT(thisParent != nullptr);
        thisParent->atoms.push_back(subAtomPtr);

        if (thisParent->role() == "MultipleSgroup") {
          auto thisMultipleParent = (MarvinMultipleSgroup *)thisParent;
          thisMultipleParent->parentAtoms.push_back(subAtomPtr);
        }

        if (thisParent == actualParent) {
          break;
        }
      }

      superatomSgroupExpanded->atoms.push_back(subAtomPtr);

      // remove the sgroupRef from the atom if it has one
      (subAtomPtr)->sgroupAttachmentPoint = "";
    }

    for (auto &bond : this->bonds) {
      actualParent->bonds.push_back(bond);
    }

    // process the attachment points - fix the bond that was made wrong by
    // deleting the dummy atom(s)

    for (auto &attachmentPoint : this->attachmentPoints) {
      // find the bond in the parent

      auto bondIter =
          find_if(actualParent->bonds.begin(), actualParent->bonds.end(),
                  [&attachmentPoint](const MarvinBond *arg) {
                    return arg->id == attachmentPoint->bond;
                  });
      if (bondIter == actualParent->bonds.end()) {
        throw FileParseException(
            "Bond specification for an AttachmentPoint definition was not found in the bond array in MRV file");
      }

      superatomSgroupExpanded->bonds.push_back(*bondIter);

      // one of the two atoms in the bond is NOT in the mol - we deleted the
      // dummy atom.

      int atomIndex;
      for (atomIndex = 0; atomIndex < 2; ++atomIndex) {
        if (!boost::algorithm::contains(
                actualParent->atoms,
                std::vector<std::string>{(*bondIter)->atomRefs2[atomIndex]},
                atomRefInAtoms)) {
          (*bondIter)->atomRefs2[atomIndex] =
              attachmentPoint->atom;  // the attach atom
          break;
        }
      }
      if (atomIndex == 2)  // not found?
      {
        std::ostringstream err;
        err << "Bond " << attachmentPoint->bond.c_str()
            << " from attachment point did not have a missing atom in MRV file";
        throw FileParseException(err.str());
      }
    }

    this->attachmentPoints.clear();
    this->atoms.clear();
    this->bonds.clear();

    // promote any children sgroups

    for (auto &childSgroup : this->sgroups) {
      moveSgroup(childSgroup.get(), actualParent, false);
    }
    this->sgroups.clear();
    auto sgroupItr =
        std::find_if(parent->sgroups.begin(), parent->sgroups.end(),
                     [this](std::unique_ptr<MarvinMolBase> &sgptr) {
                       return sgptr->id == this->id;
                     });
    TEST_ASSERT(sgroupItr != parent->sgroups.end());

    parent->sgroups.erase(sgroupItr);

  } catch (const std::exception &e) {
    throw;
  }
}

void MarvinMulticenterSgroup::processOneMulticenterSgroup() {
  auto actualParent = getActualParent(this);

  // delete the bonds to the dummy atom
  std::vector<MarvinBond *> orphanedBonds;  // list of bonds to delete
  for (auto orphanedBond : actualParent->bonds) {
    std::string attachedAtomId = "";
    if (orphanedBond->atomRefs2[0] == this->center->id ||
        orphanedBond->atomRefs2[1] == this->center->id) {
      orphanedBonds.push_back(orphanedBond);
    }
  }
  for (auto orphanedBond : orphanedBonds) {
    auto orphanedBondIter = std::find(actualParent->bonds.begin(),
                                      actualParent->bonds.end(), orphanedBond);
    TEST_ASSERT(orphanedBondIter != actualParent->bonds.end());
    actualParent->bonds.erase(orphanedBondIter);
    actualParent->removeOwnedBond(orphanedBond);
  }

  auto centerIter = std::find(actualParent->atoms.begin(),
                              actualParent->atoms.end(), this->center);
  if (centerIter !=
      actualParent->atoms.end())  // it might already have been deleted by
                                  // another multicenter group
  {
    auto centerPtr = *centerIter;
    for (auto thisParent = parent;; thisParent = thisParent->parent) {
      TEST_ASSERT(thisParent != nullptr);
      auto parentAtomIter = find(
          thisParent->atoms.begin(), thisParent->atoms.end(),
          centerPtr);  // get rid of the atoms pointer to the old dummy atom
      TEST_ASSERT(parentAtomIter != thisParent->atoms.end());
      thisParent->atoms.erase(parentAtomIter);  // get rid of the atoms pointer
                                                // to the old dummy atom

      if (thisParent->role() == "MultipleSgroup") {
        auto thisMultipleParent = (MarvinMultipleSgroup *)thisParent;
        auto deleteIter =
            find(thisMultipleParent->parentAtoms.begin(),
                 thisMultipleParent->parentAtoms.end(), centerPtr);
        TEST_ASSERT(deleteIter != thisMultipleParent->parentAtoms.end());
        thisMultipleParent->parentAtoms.erase(deleteIter);
      }

      if (thisParent == actualParent) {
        break;
      }
    }

    removeOwnedAtom(centerPtr);
  }

  // erase the multicenter group from its parent

  for (auto thisParent = parent;; thisParent = thisParent->parent) {
    TEST_ASSERT(thisParent != nullptr);
    auto sgrupIter =
        std::find_if(thisParent->sgroups.begin(), thisParent->sgroups.end(),
                     [this](std::unique_ptr<MarvinMolBase> &sgptr) {
                       return sgptr->id == this->id;
                     });
    TEST_ASSERT(sgrupIter != thisParent->sgroups.end());
    thisParent->sgroups.erase(sgrupIter);
    if (thisParent == actualParent) {
      break;
    }
  }
}

void MarvinMolBase::prepSgroupsForRDKit() {
  // this routine recursively fixes the hierarchy of sgroups - some may NOT
  // actually belong underneath their parent and can be promoted to the
  // grandparent

  // first fix all the children

  std::vector<std::string> sgroupsMolIdsDone;
  for (bool allDone = false;
       !allDone;)  // until all at this level are done - but fixing one child
                   // can add sgroups to this level
  {
    allDone = true;  // until it is not true in the loop below
    for (auto &childSgroup : this->sgroups) {
      if (!boost::algorithm::contains(
              sgroupsMolIdsDone, std::vector<std::string>{childSgroup->molID},
              molIDInSgroups)) {
        sgroupsMolIdsDone.push_back(childSgroup->molID);
        childSgroup->prepSgroupsForRDKit();
        allDone = false;
        break;  // have to start the loop over - we might have changed the
                // vector
      }
    }
  }
  // if this is top level, quit

  if (this->parent == nullptr) {
    return;
  }

  // check to see if we can move this sgroup up a level. We can if the parent
  // is not the root of the tree, and its parent is passive (it does not
  // change the atoms of the parent when it is processed or expanded, or if
  // all of its atoms do not belong to the the non-passive parent

  // if this is already a child of the root, it cannot be moved up
  if (this->parent->parent != nullptr) {
    // check this one

    if (this->parent->isPassiveRoleForExpansion()) {
      // move the child up

      promoteChild(this);
      return;
    }

    // check to see if all of its atomRefs of the child are in the parent.
    //  there are 3 possibilities:
    //  1) all atoms are in the parent - do NOT move the child - it belongs
    //  under the parent 2) NO atoms are in the parent - move the child up 3)
    //  SOME atoms are in the parent - throw an error

    auto result = this->isSgroupInSetOfAtoms(this->parent->atoms);
    switch (result) {
      case SgroupNotInAtomSet:
        promoteChild(this);
        return;

      case SgroupInAtomSet:
        // do noting here
        break;

      case SgroupBothInAndNotInAtomSet:
        throw FileParseException(
            "Child sGroup has atoms both in the parent and NOT in the parent");
      default:
        throw FileParseException(
            "Unexpected error: unrecognized result from isSgroupInSetOfAtoms");
    }
  }

  // made to here, so the child belongs in the parent.
  //  do any sgroup specific processing

  processSpecialSgroups();
  return;
}

void MarvinMolBase::processSpecialSgroups() {
  // by default, nothing is done here.  Some subclasses will override this.
}

void MarvinSuperatomSgroup::processSpecialSgroups() {
  convertFromOneSuperAtom();
}

void MarvinMulticenterSgroup::processSpecialSgroups() {
  processOneMulticenterSgroup();
}
void MarvinMultipleSgroup::processSpecialSgroups() {
  expandOneMultipleSgroup();
}

MarvinMolBase *MarvinSuperatomSgroupExpanded::convertToOneSuperAtom() {
  // the mol-style super atoms are significantly different than the Marvin
  // super atoms.  The mol style has all atoms and bonds in the main mol, and
  // parameter lines that indicate the name of the super atom and the atom
  // indices affected.
  //
  // The Marvin has a dummy atom with the super atom name in the main
  // molecule, and a separate sub-mol with the atoms and bonds of the
  // superatom.  It also has a separate record called attachmentPoint that
  // atom in the super atom sub=mol that is replaces the dummy atom, and also
  // a bond pointer and bond order in case the bond order changes.
  //
  //  Takes information from a MarvinSuperatomSgroupExpanded and converts the
  //  MarvinSuperatomSgroup structure

  PRECONDITION(this->parent, "invalid parent");
  auto actualParent = getActualParent(this);

  // make a new sub mol

  auto marvinSuperatomSgroup = new MarvinSuperatomSgroup(this->parent);
  this->parent->sgroups.push_back(
      std::unique_ptr<MarvinMolBase>(marvinSuperatomSgroup));
  std::string newAtomName = "NA_" + this->molID;

  marvinSuperatomSgroup->molID = this->molID;
  marvinSuperatomSgroup->title = this->title;
  marvinSuperatomSgroup->id = this->id;

  bool coordsExist = has2dCoords();  // we have to check before we add the dummy
                                     // atom  - it will not have coords (yet)

  // add the dummy atom into the parent

  auto dummyParentAtom = new MarvinAtom();
  actualParent->pushOwnedAtom(dummyParentAtom);
  dummyParentAtom->elementType = "R";
  dummyParentAtom->id = newAtomName;
  dummyParentAtom->sgroupRef = marvinSuperatomSgroup->id;
  dummyParentAtom->sGroupRefIsSuperatom = true;

  for (auto thisParent = this->parent;;
       thisParent = thisParent->parent)  // do all parents and grandparents
                                         // ... to to the actual parent
  {
    thisParent->atoms.push_back(dummyParentAtom);

    if (thisParent->role() == "MultipleSgroup") {
      ((MarvinMultipleSgroup *)thisParent)
          ->parentAtoms.push_back(dummyParentAtom);
    }

    if (thisParent == actualParent) {
      break;
    }
  }

  RDGeom::Point3D centerOfGroup;

  for (auto atom : this->atoms) {
    atom->sgroupRef = "";
    marvinSuperatomSgroup->atoms.push_back(atom);

    for (auto thisParent = this->parent;;
         thisParent = thisParent->parent)  // do all parents and grandparents
                                           // ... to to the actual parent
    {
      auto deleteIter =
          find(thisParent->atoms.begin(), thisParent->atoms.end(), atom);
      TEST_ASSERT(deleteIter != thisParent->atoms.end());
      thisParent->atoms.erase(deleteIter);

      if (thisParent->role() == "MultipleSgroup") {
        auto marvinMultipleSgroup = (MarvinMultipleSgroup *)thisParent;
        auto sgroupIter = find(marvinMultipleSgroup->parentAtoms.begin(),
                               marvinMultipleSgroup->parentAtoms.end(), atom);
        TEST_ASSERT(sgroupIter != marvinMultipleSgroup->parentAtoms.end());
        marvinMultipleSgroup->parentAtoms.erase(sgroupIter);
      }

      if (thisParent == actualParent) {
        break;
      }
    }

    if (coordsExist)  // get the center of all atoms in the group - we might
                      // use this if there are no attachment points
    {
      centerOfGroup.x += atom->x2;
      centerOfGroup.y += atom->y2;
    }
  }

  // move the bonds of the group

  int attachmentPointsAdded = 0;
  MarvinAtom *atomPtr = nullptr;
  RDGeom::Point3D centerOfAttachmentPoints;

  for (auto bond : actualParent->bonds) {
    bool atom1InGroup = boost::algorithm::contains(
        marvinSuperatomSgroup->atoms,
        std::vector<std::string>{(bond)->atomRefs2[0]}, atomRefInAtoms);
    bool atom2InGroup = boost::algorithm::contains(
        marvinSuperatomSgroup->atoms,
        std::vector<std::string>{(bond)->atomRefs2[1]}, atomRefInAtoms);
    if (atom1InGroup && atom2InGroup) {  // both are in, so copy the bond
      marvinSuperatomSgroup->bonds.push_back(bond);

    } else if (atom1InGroup ||
               atom2InGroup)  // one is in so this is a connection point
    {
      // fix the bonds to the dummy atom and add an attachment point

      if (atom1InGroup) {
        atomPtr =
            marvinSuperatomSgroup->atoms[marvinSuperatomSgroup->getAtomIndex(
                (bond)->atomRefs2[0])];
        (bond)->atomRefs2[0] = newAtomName;

      } else {
        atomPtr =
            marvinSuperatomSgroup->atoms[marvinSuperatomSgroup->getAtomIndex(
                (bond)->atomRefs2[1])];
        (bond)->atomRefs2[1] = newAtomName;
      }

      atomPtr->sgroupAttachmentPoint = std::to_string(++attachmentPointsAdded);

      // add an attachmentPoint structure

      auto &marvinAttachmentPoint =
          marvinSuperatomSgroup->attachmentPoints.emplace_back(
              new MarvinAttachmentPoint);
      marvinAttachmentPoint->atom = atomPtr->id;
      marvinAttachmentPoint->bond = bond->id;

      marvinAttachmentPoint->order = std::to_string(attachmentPointsAdded);

      if (coordsExist) {
        centerOfAttachmentPoints.x += atomPtr->x2;
        centerOfAttachmentPoints.y += atomPtr->y2;
      }
    }
  }

  // now remove the bonds that were moved to the superGroup from the parent

  for (auto bond : marvinSuperatomSgroup->bonds) {
    int index = actualParent->getBondIndex(bond->id);
    TEST_ASSERT(index >= 0 &&
                static_cast<unsigned int>(index) < actualParent->bonds.size());
    actualParent->bonds.erase(actualParent->bonds.begin() + index);
  }

  if (coordsExist) {
    if (marvinSuperatomSgroup->attachmentPoints.size() >
        0)  // Any attachment points?  if so we use the center of the attached
            // atoms in the supergroup
    {
      // put the new dummy atom at the center of the removed group

      dummyParentAtom->x2 = centerOfAttachmentPoints.x /
                            marvinSuperatomSgroup->attachmentPoints.size();
      dummyParentAtom->y2 = centerOfAttachmentPoints.y /
                            marvinSuperatomSgroup->attachmentPoints.size();
    } else {
      // No attachments to the supergroup - (probably all atoms in the mol are
      // in the supergroup) -use the center of all atoms in the group
      dummyParentAtom->x2 =
          centerOfGroup.x / marvinSuperatomSgroup->atoms.size();
      dummyParentAtom->y2 =
          centerOfGroup.y / marvinSuperatomSgroup->atoms.size();
    }
  }

  // move any children from the old expanded group to the new superatomSgroup

  for (auto &sgroupPtr : this->sgroups) {
    moveSgroup(sgroupPtr.get(), marvinSuperatomSgroup, false);
  }
  this->sgroups.clear();

  // remove the old expanded superatom sgroup
  auto sgroupIter = std::find_if(parent->sgroups.begin(), parent->sgroups.end(),
                                 [this](std::unique_ptr<MarvinMolBase> &sgptr) {
                                   return sgptr->id == this->id;
                                 });
  TEST_ASSERT(sgroupIter != parent->sgroups.end());
  parent->sgroups.erase(sgroupIter);

  // fix up any siblings that contain all the atoms of this group.

  for (auto &sibling : marvinSuperatomSgroup->parent->sgroups) {
    if (sibling.get() == marvinSuperatomSgroup) {
      continue;  // not a sibling
    }

    // see if the new superatom group has all of its atoms in the sibling.  If
    // to add the dummy atom to the sibling and remove the atoms (there must
    // be at least one atom in the sibling that is not in the new superatom,
    // or it would be a child of the superatom already

    auto isSgroupInSetOfAtomsResult =
        marvinSuperatomSgroup->isSgroupInSetOfAtoms(sibling->atoms);

    switch (isSgroupInSetOfAtomsResult) {
      case SgroupInAtomSet:
        // remove the group atoms from the sibling and add the dummy atom

        sibling->atoms.push_back(dummyParentAtom);
        for (auto atomToRemove : marvinSuperatomSgroup->atoms) {
          int index = sibling->getAtomIndex(atomToRemove->id);
          TEST_ASSERT(index >= 0 &&
                      static_cast<unsigned int>(index) < sibling->atoms.size());
          sibling->atoms.erase(sibling->atoms.begin() + index);
        }

        break;

      case SgroupNotInAtomSet:
        // do nothing
        break;

      case SgroupBothInAndNotInAtomSet:
        throw FileParseException(
            "a subgroup contains some of the atoms of a superatomSgroup - this should not happen");
    }
  }

  return marvinSuperatomSgroup;
}

int MarvinMultipleSgroup::getMatchedOrphanBondIndex(
    std::string atomIdToCheck, std::vector<MarvinBond *> &bondsToTry,
    std::vector<MarvinBond *> &orphanedBonds) const {
  for (auto testBond = bondsToTry.begin(); testBond != bondsToTry.end();
       ++testBond) {
    if (*testBond == nullptr) {
      continue;
    }
    std::string otherAtomId;
    if ((*testBond)->atomRefs2[0] == atomIdToCheck) {
      otherAtomId = (*testBond)->atomRefs2[1];
    } else if ((*testBond)->atomRefs2[1] == atomIdToCheck) {
      otherAtomId = (*testBond)->atomRefs2[0];
    } else {
      continue;  // bond does not have the atomId
    }

    auto testBondIter =
        find(orphanedBonds.begin(), orphanedBonds.end(), *testBond);
    if (testBondIter != orphanedBonds.end()) {
      return testBondIter - orphanedBonds.begin();
    }

    *testBond = nullptr;

    // try the children

    int childId =
        getMatchedOrphanBondIndex(otherAtomId, bondsToTry, orphanedBonds);
    if (childId >= 0) {
      return childId;
    }
  }

  return -1;
}

void MarvinMultipleSgroup::contractOneMultipleSgroup() {
  // this routine takes the expanded MultipleSgroup (which comes from the
  // RDKit version, or is ready to produce the RDKit version), and contracts
  // it
  //  to the Marvin format version.
  //  the replicates are deleted (atoms and bonds), and the two orphaned bonds
  //  are  repaired

  // make sure it is not already contracted - this should not happen

  if (!this->isExpanded) {
    return;
  }

  // find the actual parent - multiple groups can be nested, so the atoms and
  // bonds of the lower one are really in the grandparent or higher)

  auto actualParent = getActualParent(this);

  // find the atoms to be deleted

  std::vector<MarvinAtom *> atomsToDelete;
  std::vector<MarvinBond *> bondsToDelete;
  std::vector<MarvinBond *> orphanedBonds;
  std::vector<MarvinMolBase *> sgroupsToDelete;

  for (MarvinAtom *atomPtr : this->atoms) {
    // Atoms in the group but NOT in the parent part are to be deleted
    if (!boost::algorithm::contains(this->parentAtoms,
                                    std::vector<std::string>{atomPtr->id},
                                    atomRefInAtoms)) {
      atomsToDelete.push_back(atomPtr);
    }
  }

  for (MarvinBond *bondPtr : actualParent->bonds) {
    bool atom1ToBeDeleted = boost::algorithm::contains(
        atomsToDelete, std::vector<std::string>{bondPtr->atomRefs2[0]},
        atomRefInAtoms);
    bool atom2ToBeDeleted = boost::algorithm::contains(
        atomsToDelete, std::vector<std::string>{bondPtr->atomRefs2[1]},
        atomRefInAtoms);
    if (atom1ToBeDeleted && atom2ToBeDeleted) {
      bondsToDelete.push_back(bondPtr);
    } else if (atom1ToBeDeleted) {
      orphanedBonds.push_back(bondPtr);
      swap(bondPtr->atomRefs2[0],
           bondPtr->atomRefs2[1]);  // make sure the first ref atom is the
                                    // orphaned one (still in the mol), and
                                    // the second is the atom to be deleted
      bondsToDelete.push_back(
          bondPtr);  // one of each pair of orphaned bonds will be UNdeleted
    } else if (atom2ToBeDeleted) {
      orphanedBonds.push_back(bondPtr);
      bondsToDelete.push_back(
          bondPtr);  // one of each pair of orphaned bonds will be UNdeleted
    }
  }

  // there must be zero two or four orphaned bonds
  // zero happens when the whole group is replicated and not chained
  // two is most common - the replicates being deleted are connected to each
  // other 4 happens when the parent replicate (not deleted) is in between two
  // replicates that are deleted.
  //  there are two orphaned bonds on the parent replicate, and two on the
  //  rest of the molecule

  if (orphanedBonds.size() != 0 && orphanedBonds.size() != 2 &&
      orphanedBonds.size() != 4) {
    throw FileParseException(
        "Error: there should be zero, two or four orphaned bonds while contracting a MultipleSgroup");
  }

  // delete any of this group's children that were using the atoms to be
  // deleted

  for (auto &childSgroup : sgroups) {
    if (childSgroup->isSgroupInSetOfAtoms(atomsToDelete) !=
        SgroupNotInAtomSet) {
      sgroupsToDelete.push_back(childSgroup.get());
    }
  }

  // now fix the orphaned bonds - The first gets the atoms from the matched
  // second orphan bond  and was NOT removed. the matched second orphaned bond
  // is deleted

  while (orphanedBonds.size() > 0) {
    int matchedOrphanBondIndex = -1;
    auto orphanedBondToFix = orphanedBonds[0];
    orphanedBonds.erase(orphanedBonds.begin());

    if (orphanedBonds.size() == 3)  // was 4 but the first one was erased
    {
      std::vector<MarvinBond *> bondsToTry =
          bondsToDelete;  // copy of bonds to delete
      auto deleteIter =
          find(bondsToTry.begin(), bondsToTry.end(), orphanedBondToFix);
      TEST_ASSERT(deleteIter != bondsToTry.end());
      bondsToTry.erase(deleteIter);

      matchedOrphanBondIndex = this->getMatchedOrphanBondIndex(
          orphanedBondToFix->atomRefs2[1], bondsToTry, orphanedBonds);
      if (matchedOrphanBondIndex < 0) {
        throw FileParseException(
            "For a MultipleSgroup expansion, the matching orphaned bond was not found");
      }
    } else if (orphanedBonds.size() ==
               1) {                // was 2 but the first one was erased
      matchedOrphanBondIndex = 0;  // only one left - it must be the match
    } else {
      throw FileParseException(
          "For a MultipleSgroup contraction, the orphaned bonds must be 0,2 or 4");
    }

    orphanedBondToFix->atomRefs2[1] =
        orphanedBonds[matchedOrphanBondIndex]
            ->atomRefs2[0];  // [0] is the orphaned atom (still in the mol)

    // undelete the bond which has been fixed (really remove it from the list
    // of bonds to be delete)

    auto undeleteIter =
        find(bondsToDelete.begin(), bondsToDelete.end(), orphanedBondToFix);
    TEST_ASSERT(undeleteIter != bondsToDelete.end());
    bondsToDelete.erase(undeleteIter);

    // any siblings or children that reference the matched Orphaned bond must
    // be fixed to reference the retained (fixed) orphan bond

    for (auto &siblingSgroup :
         this->parent->sgroups)  // note: includes this sgroups, which is
                                 // technically NOT a sibling
    {
      for (auto &bond : siblingSgroup->bonds) {
        if (bond->id == orphanedBonds[matchedOrphanBondIndex]->id) {
          bond = orphanedBondToFix;
        }
      }
    }

    // it COULD happen that children referenced the matched orphan bond

    for (auto &childSgroup : this->sgroups) {
      for (auto childBond : childSgroup->bonds) {
        if (childBond->id == orphanedBonds[matchedOrphanBondIndex]->id) {
          childBond = orphanedBondToFix;
        }
      }
    }
    TEST_ASSERT(matchedOrphanBondIndex >= 0 &&
                static_cast<unsigned int>(matchedOrphanBondIndex) <
                    orphanedBonds.size());
    orphanedBonds.erase(orphanedBonds.begin() + matchedOrphanBondIndex);
  };

  // if we get here, all the changes to be made have been determined. so
  // make them

  for (auto childSgroup : sgroupsToDelete) {
    auto sgroupIter =
        std::find_if(sgroups.begin(), sgroups.end(),
                     [childSgroup](std::unique_ptr<MarvinMolBase> &sgptr) {
                       return sgptr->id == childSgroup->id;
                     });
    TEST_ASSERT(sgroupIter != sgroups.end());
    sgroups.erase(sgroupIter);
  }

  // remove the atoms

  for (auto atomPtr : atomsToDelete) {
    for (auto thisParent = this->parent;;
         thisParent = thisParent->parent) {  // do all parents and grandparents
                                             // ... to to the actual parent

      auto deleteIter = std::find(thisParent->atoms.begin(),
                                  thisParent->atoms.end(), atomPtr);
      TEST_ASSERT(deleteIter != thisParent->atoms.end());
      thisParent->atoms.erase(deleteIter);

      if (thisParent == actualParent) {
        break;
      }
    }

    // also delete the atoms from this multiple group AND any sibling sgroups

    for (auto &siblingSgroup : this->parent->sgroups) {
      auto siblingAtomPtr = std::find(siblingSgroup->atoms.begin(),
                                      siblingSgroup->atoms.end(), atomPtr);
      if (siblingAtomPtr != siblingSgroup->atoms.end()) {
        siblingSgroup->atoms.erase(siblingAtomPtr);
      }
    }

    removeOwnedAtom(atomPtr);
  }

  atomsToDelete.clear();

  // remove the bonds

  for (MarvinBond *bondPtr : bondsToDelete) {
    auto deleteIter = std::find(actualParent->bonds.begin(),
                                actualParent->bonds.end(), bondPtr);
    TEST_ASSERT(deleteIter != actualParent->bonds.end());
    actualParent->bonds.erase(deleteIter);
    actualParent->removeOwnedBond(bondPtr);
  }
  bondsToDelete.clear();

  this->isExpanded = false;
}

IsSgroupInAtomSetResult MarvinMolBase::isSgroupInSetOfAtoms(
    const std::vector<MarvinAtom *> &setOfAtoms) const {
  // superatom sgroups are different - it has an override for this call
  // if not overridden,  types are in the group based on their atoms

  bool isInGroup = false;
  bool isNotInGroup = false;
  for (auto oneAtom : this->atoms) {
    if (boost::algorithm::contains(setOfAtoms,
                                   std::vector<std::string>{oneAtom->id},
                                   atomRefInAtoms)) {
      isInGroup = true;
    } else {
      isNotInGroup = true;
    }
  }

  if (isInGroup && isNotInGroup) {
    return SgroupBothInAndNotInAtomSet;
  }

  if (isInGroup) {
    return SgroupInAtomSet;
  }

  return SgroupNotInAtomSet;
}

IsSgroupInAtomSetResult MarvinSuperatomSgroup::isSgroupInSetOfAtoms(
    const std::vector<MarvinAtom *> &setOfAtoms) const {
  // superatom sgroups are different - they are in the set if one of the
  // condensed atom in the parent is in group

  auto dummyAtomIter = find_if(
      parent->atoms.begin(), parent->atoms.end(),
      [this](const MarvinAtom *arg) { return arg->sgroupRef == this->id; });
  if (dummyAtomIter == parent->atoms.end()) {
    throw FileParseException("No contracted atom found for a superatomSgroup");
  }

  if (boost::algorithm::contains(setOfAtoms,
                                 std::vector<std::string>{(*dummyAtomIter)->id},
                                 atomRefInAtoms)) {
    return SgroupInAtomSet;
  }
  return SgroupNotInAtomSet;
}

void MarvinMolBase::processSgroupsFromRDKit() {
  // this routine recursively fixes sgroups to be in the heirarchy expected by
  // Marvin format

  if (this->isPassiveRoleForContraction()) {
    return;  // passive mols are only leaves of the tree - nothing to do
  }

  // see if any siblings should be children of this one

  MarvinMolBase *molToProcess =
      this;  // the processing MIGHT change the mol from one kind to another
             // (e.g. MarvinSuperatomSgroupExpanded -> MarvinSuperatomSgroup)

  if (this->parent != nullptr)  // NOT the top level mol
  {
    std::vector<std::string> sgroupsMolIdsDone;
    for (bool allDone = false; !allDone;) {
      // until all at this level are done - but fixing one
      // child can add sgroups to this level
      allDone = true;  // until it is not true in the loop below
      for (auto &sibling : this->parent->sgroups) {
        if (boost::algorithm::contains(sgroupsMolIdsDone,
                                       std::vector<std::string>{sibling->molID},
                                       molIDInSgroups)) {
          continue;  // look at the next one to see if it is done
        }

        allDone = false;
        sgroupsMolIdsDone.push_back(sibling->molID);

        if (sibling.get() == this) {
          continue;  // cannot be a child of itself
        }

        //  there are 2 possibilities:
        //  1) all atoms are in this one - make the sibling into a child of
        //  this one 2) SOME atoms are NOT this one - it remains a sibling

        if (sibling->isSgroupInSetOfAtoms(this->atoms) != SgroupInAtomSet) {
          continue;
        }

        moveSgroup(sibling.get(), this);
        break;  // have to re=start the loop over sgroups - it has just
                // changed
      }
    }

    std::string thisSgroupRole = this->role();
    if (thisSgroupRole == "SuperatomSgroupExpanded") {
      molToProcess =
          ((MarvinSuperatomSgroupExpanded *)this)->convertToOneSuperAtom();
    } else if (thisSgroupRole == "MultipleSgroup" &&
               ((MarvinMultipleSgroup *)this)->isExpanded == true) {
      try {
        ((MarvinMultipleSgroup *)this)->contractOneMultipleSgroup();
      } catch (FileParseException &e) {
        BOOST_LOG(rdErrorLog)
            << e.what() << std::endl
            << "The multiple sgroup will be ignored" << std::endl;

        auto sgroupIter = std::find_if(
            this->parent->sgroups.begin(), this->parent->sgroups.end(),
            [this](std::unique_ptr<MarvinMolBase> &sgptr) {
              return sgptr->id == this->id;
            });
        TEST_ASSERT(sgroupIter != this->parent->sgroups.end());
        this->parent->sgroups.erase(sgroupIter);
        return;
      }
    }
  }

  // now fix all this groups children

  std::vector<std::string> childrenSgroupsMolIdsDone;
  for (bool allDone = false;
       !allDone;)  // until all at this level are done - but fixing one child
                   // can add sgroups to this level
  {
    allDone = true;  // until it is not true in the loop below
    for (auto &childSgroup : molToProcess->sgroups) {
      if (boost::algorithm::contains(
              childrenSgroupsMolIdsDone,
              std::vector<std::string>{childSgroup->molID}, molIDInSgroups)) {
        continue;  // this one is done, look at the next one
      }

      childrenSgroupsMolIdsDone.push_back(childSgroup->molID);
      childSgroup->processSgroupsFromRDKit();
      allDone = false;
      break;  // have to start the loop over - we might have changed the
              // vector
    }
  }

  return;
}

MarvinMolBase *MarvinMol::copyMol(const std::string &) const {
  PRECONDITION(0, "Internal error:  copying a MarvinMol");
}

std::string MarvinMol::toString() const {
  std::ostringstream out;

  out << "<molecule molID=\"" << molID << "\">";

  out << "<atomArray>";
  for (auto atom : atoms) {
    out << atom->toString();
  }
  out << "</atomArray>";

  out << "<bondArray>";
  for (auto bond : bonds) {
    out << bond->toString();
  }
  out << "</bondArray>";

  for (auto &sgroup : sgroups) {
    out << sgroup->toString();
  }

  out << "</molecule>";

  return out.str();
}

ptree MarvinMol::toPtree() const {
  ptree out = MarvinMolBase::toPtree();

  MarvinMolBase::addSgroupsToPtree(out);

  return out;
}

std::string MarvinMol::generateMolString() {
  std::ostringstream out;

  out << "<cml xmlns=\"http://www.chemaxon.com\"  xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
         "xsi:schemaLocation=\"http://www.chemaxon.com http://www.chemaxon.com/marvin/schema/mrvSchema_20_20_0.xsd\">"
         "<MDocument><MChemicalStruct>";

  out << toString();

  out << "</MChemicalStruct></MDocument></cml>";

  return out.str();
}

ptree MarvinMol::toMolPtree() const {
  ptree out;

  out.put("cml.<xmlattr>.xmlns", "http://www.chemaxon.com");
  out.put("cml.<xmlattr>.xmlns:xsi",
          "http://www.w3.org/2001/XMLSchema-instance");
  out.put(
      "cml.<xmlattr>.xsi:schemaLocation",
      "http://www.chemaxon.com http://www.chemaxon.com/marvin/schema/mrvSchema_20_20_0.xsd");
  out.put_child("cml.MDocument.MChemicalStruct.molecule", toPtree());
  return out;
}

MarvinReaction::~MarvinReaction() {}

void MarvinReaction::prepSgroupsForRDKit() {
  // This routine converts all the mols in the rxn to be ready for conversion
  // to RDKIT mols
  int molCount = 0, atomCount = 0, bondCount = 0, sgCount = 0;

  std::map<std::string, std::string> sgMap;
  std::map<std::string, std::string> atomMap;
  std::map<std::string, std::string> bondMap;

  for (auto &reactant : reactants) {
    reactant->prepSgroupsForRDKit();

    reactant->cleanUpNumbering(molCount, atomCount, bondCount, sgCount, sgMap,
                               atomMap, bondMap);
  }
  for (auto &agent : agents) {
    agent->prepSgroupsForRDKit();
    agent->cleanUpNumbering(molCount, atomCount, bondCount, sgCount, sgMap,
                            atomMap, bondMap);
  }
  for (auto &product : products) {
    product->prepSgroupsForRDKit();
    product->cleanUpNumbering(molCount, atomCount, bondCount, sgCount, sgMap,
                              atomMap, bondMap);
  }
}

std::string MarvinReaction::toString() {
  std::ostringstream out;

  out << "<cml xmlns=\"http://www.chemaxon.com\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\""
         " xsi:schemaLocation=\"http://www.chemaxon.com http://www.chemaxon.com/marvin/schema/mrvSchema_20_20_0.xsd\">"
         "<MDocument><MChemicalStruct><reaction>";

  out << "<reactantList>";
  for (auto &reactant : reactants) {
    out << reactant->toString();
  }
  out << "</reactantList>";
  out << "<agentList>";
  for (auto &agent : agents) {
    out << agent->toString();
  }
  out << "</agentList>";
  out << "<productList>";
  for (auto &product : products) {
    out << product->toString();
  }
  out << "</productList>";

  out << arrow.toString();

  out << "</reaction></MChemicalStruct>";

  for (auto &pluse : pluses) {
    out << pluse->toString();
  }
  for (auto &condition : conditions) {
    out << condition->toString();
  }

  out << "</MDocument></cml>";
  return out.str();
}

ptree MarvinReaction::toPtree() const {
  ptree out;

  out.put("cml.<xmlattr>.xmlns", "http://www.chemaxon.com");
  out.put("cml.<xmlattr>.xmlns:xsi",
          "http://www.w3.org/2001/XMLSchema-instance");
  out.put(
      "cml.<xmlattr>.xsi:schemaLocation",
      "http://www.chemaxon.com http://www.chemaxon.com/marvin/schema/mrvSchema_20_20_0.xsd");

  for (auto &reactant : reactants) {
    out.add_child(
        "cml.MDocument.MChemicalStruct.reaction.reactantList.molecule",
        reactant->toPtree());
  }

  for (auto &agent : agents) {
    out.add_child("cml.MDocument.MChemicalStruct.reaction.agentList.molecule",
                  agent->toPtree());
  }

  for (auto &product : products) {
    out.add_child("cml.MDocument.MChemicalStruct.reaction.productList.molecule",
                  product->toPtree());
  }

  out.put_child("cml.MDocument.MChemicalStruct.reaction.arrow",
                arrow.toPtree());

  for (auto &plus : pluses) {
    out.add_child("cml.MDocument.MReactionSign", plus->toPtree());
  }

  for (auto &condition : conditions) {
    out.add_child("cml.MDocument.MTextBox", condition->toPtree());
  }

  return out;
}

MarvinRectangle::MarvinRectangle(double left, double right, double top,
                                 double bottom) {
  upperLeft.x = left;
  upperLeft.y = top;
  lowerRight.x = right;
  lowerRight.y = bottom;
  centerIsStale = true;
}

MarvinRectangle::MarvinRectangle(const RDGeom::Point3D &upperLeftInit,
                                 const RDGeom::Point3D &lowerRightInit) {
  upperLeft = upperLeftInit;
  lowerRight = lowerRightInit;
  centerIsStale = true;
}

MarvinRectangle::MarvinRectangle(const std::vector<MarvinAtom *> &atoms) {
  centerIsStale = true;

  if (atoms.size() == 0) {
    return;
  }
  upperLeft.x = DBL_MAX;
  upperLeft.y = -DBL_MAX;
  lowerRight.x = -DBL_MAX;
  lowerRight.y = DBL_MAX;

  for (auto atom : atoms) {
    upperLeft.x = std::min(upperLeft.x, atom->x2);
    lowerRight.x = std::max(lowerRight.x, atom->x2);
    upperLeft.y = std::max(upperLeft.y, atom->y2);
    lowerRight.y = std::min(lowerRight.y, atom->y2);
  }
}

MarvinRectangle::MarvinRectangle(const std::vector<MarvinRectangle> &rects) {
  centerIsStale = true;

  if (rects.size() == 0) {
    return;
  }
  upperLeft.x = DBL_MAX;
  upperLeft.y = -DBL_MAX;
  lowerRight.x = -DBL_MAX;
  lowerRight.y = DBL_MAX;

  for (auto rect : rects) {
    this->extend(rect);
  }
}

void MarvinRectangle::extend(const MarvinRectangle &otherRectangle) {
  if (otherRectangle.upperLeft.x < upperLeft.x) {
    upperLeft.x = otherRectangle.upperLeft.x;
  }
  if (otherRectangle.lowerRight.x > lowerRight.x) {
    lowerRight.x = otherRectangle.lowerRight.x;
  }

  if (otherRectangle.upperLeft.y > upperLeft.y) {
    upperLeft.y = otherRectangle.upperLeft.y;
  }
  if (otherRectangle.lowerRight.y < lowerRight.y) {
    lowerRight.y = otherRectangle.lowerRight.y;
  }

  centerIsStale = true;
}

RDGeom::Point3D &MarvinRectangle::getCenter() {
  if (centerIsStale) {
    center.x = (lowerRight.x + upperLeft.x) / 2.0;
    center.y = (lowerRight.y + upperLeft.y) / 2.0;
    centerIsStale = false;
  }
  return center;
}

bool MarvinRectangle::overlapsVertically(
    const MarvinRectangle &otherRectangle) const {
  if (otherRectangle.upperLeft.y < lowerRight.y ||
      otherRectangle.lowerRight.y > upperLeft.y) {
    return false;
  }
  return true;
}

bool MarvinRectangle::overlapsVHorizontally(
    const MarvinRectangle &otherRectangle) const {
  if (otherRectangle.upperLeft.x > lowerRight.x ||
      otherRectangle.lowerRight.x < upperLeft.x) {
    return false;
  }
  return true;
}

bool MarvinRectangle::compareRectanglesByX(MarvinRectangle &r1,
                                           MarvinRectangle &r2) {
  return (r1.getCenter().x < r2.getCenter().x);
}

bool MarvinRectangle::compareRectanglesByYReverse(MarvinRectangle &r1,
                                                  MarvinRectangle &r2) {
  return (r1.getCenter().y > r2.getCenter().y);
}

MarvinStereoGroup::MarvinStereoGroup(StereoGroupType grouptypeInit,
                                     int groupNumberInit) {
  groupType = grouptypeInit;
  groupNumber = groupNumberInit;
}
}  // namespace RDKit
