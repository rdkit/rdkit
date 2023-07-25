//
//  Copyright (C) 2002-2022 Greg Landrum and other RDKit contributors
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

thread_local std::map<std::string, std::string> sgMap;
thread_local std::map<std::string, std::string> atomMap;
thread_local std::map<std::string, std::string> bondMap;

std::string getDoubleAsText(double val) {
  if (fabs(val) < 0.00001) {
    return "0.00000";
  }

  std::ostringstream valstr;
  valstr << val;
  return valstr.str();
}

void MarvinMolBase::clearMaps() {
  sgMap.clear();
  atomMap.clear();
  bondMap.clear();
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
      out.add_child("atomArray.atom", atom->toPtree());
    }

    for (auto bond : bonds) {
      out.add_child("bondArray.bond", bond->toPtree());
    }
  }

  return out;
}

void MarvinMolBase::addSgroupsToPtree(ptree &out) const {
  for (auto sgroup : sgroups) {
    out.add_child("molecule", sgroup->toPtree());
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

ptree MarvinAtom::toPtree() const {
  ptree out;

  out.put("<xmlattr>.id", id);
  out.put("<xmlattr>.elementType", elementType);

  if (x2 != DBL_MAX && y2 != DBL_MAX) {
    out.put("<xmlattr>.x2", getDoubleAsText(x2));
    out.put("<xmlattr>.y2", getDoubleAsText(y2));
  }

  if (x3 != DBL_MAX && y3 != DBL_MAX && z3 != DBL_MAX) {
    out.put("<xmlattr>.x3", getDoubleAsText(x3));
    out.put("<xmlattr>.y3", getDoubleAsText(y3));
    out.put("<xmlattr>.z3", getDoubleAsText(z3));
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
             "")  // if no query type not conventtion,  so check for order
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
    std::ostringstream err;
    err << "bond must have one of:  order, queryType, or convention in MRV File ";
    throw FileParseException(err.str());
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

MarvinMolBase::~MarvinMolBase() {
  for (auto sgroup : sgroups) {
    delete sgroup;
  }
}

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

bool MarvinMolBase::hasCoords() const { return has3dCoords() || has3dCoords(); }

void MarvinMolBase::removeCoords() {
  for (auto atom : atoms) {
    atom->x2 = DBL_MAX;
    atom->y2 = DBL_MAX;
  }
}

int MarvinMolBase::getExplicitValence(const MarvinAtom &marvinAtom) const {
  unsigned int resTimes10 = 0;  // calculated as 10 * the actual value so we
                                // can use int match, and have 1.5 order bonds

  for (auto bondPtr : bonds) {
    if (bondPtr->atomRefs2[0] != marvinAtom.id &&
        bondPtr->atomRefs2[1] != marvinAtom.id) {
      continue;  // this bond is NOT to the atom
    }

    std::string tempConvention =
        boost::algorithm::to_upper_copy(bondPtr->convention);
    std::string marvinBondType = bondPtr->getBondType();

    if (marvinBondType == "SD" || marvinBondType == "SA" ||
        marvinBondType == "DA") {
      resTimes10 += 15;  // really 1.5 order bond
    } else if (marvinBondType == "ANY") {
      resTimes10 +=
          10;  // no good answer for Any bonds - treat as a single bond
    } else if (marvinBondType == "DATIVE") {
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
  parent = parentInit;

  if (roleNameInit != "SruSgroup" && roleNameInit != "CopolymerSgroup" &&
      roleNameInit != "ModificationSgroup") {
    throw FileParseException(
        "A MarvinSruCoModSgroup type must be one of \"SruSgroup\". \"CopolymerSgroup\", and \"ModificationSgroup\"");
  }

  this->roleName = roleNameInit;
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

MarvinMolBase *MarvinSruCoModSgroup::copyMol(std::string idAppendage) const {
  auto outSgroup = new MarvinSruCoModSgroup(this->roleName, this->parent);
  this->parent->sgroups.push_back(outSgroup);

  outSgroup->molID = this->molID + idAppendage;
  outSgroup->title = this->title;
  outSgroup->connect = this->connect;
  outSgroup->correspondence = this->correspondence;

  auto actualParent = this->parent;
  while (actualParent->role() == "MultipleSgroup") {
    actualParent = actualParent->parent;
  }

  // the only time this is to be called is when a mutliple group above it is
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

std::string MarvinSruCoModSgroup::role() const { return std::string(roleName); }

bool MarvinSruCoModSgroup::hasAtomBondBlocks() const { return false; }

MarvinDataSgroup::MarvinDataSgroup(MarvinMolBase *parentInit) {
  parent = parentInit;
}

MarvinMolBase *MarvinDataSgroup::copyMol(std::string idAppendage) const {
  auto outSgroup = new MarvinDataSgroup(this->parent);
  this->parent->sgroups.push_back(outSgroup);

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

  // the only time this is to be called is when a mutliple group above it is
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

std::string MarvinDataSgroup::role() const { return std::string("DataSgroup"); }

bool MarvinDataSgroup::hasAtomBondBlocks() const { return false; }

MarvinMultipleSgroup::MarvinMultipleSgroup(MarvinMolBase *parentInit) {
  parent = parentInit;
}

std::string MarvinMultipleSgroup::role() const {
  return std::string("MultipleSgroup");
}

bool MarvinMultipleSgroup::hasAtomBondBlocks() const { return false; }

MarvinMolBase *MarvinMultipleSgroup::copyMol(std::string idAppendage) const {
  auto outSgroup = new MarvinMultipleSgroup(this->parent);
  this->parent->sgroups.push_back(outSgroup);

  outSgroup->molID = this->molID + idAppendage;
  outSgroup->id = this->id + idAppendage;
  outSgroup->title = this->title;

  // the only time this is to be called is when a mutliple group above it is
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

  for (auto sgroup : sgroups) {
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
  parent = parentInit;
}

std::string MarvinMulticenterSgroup::role() const {
  return std::string("MulticenterSgroup");
}

bool MarvinMulticenterSgroup::hasAtomBondBlocks() const { return false; }

MarvinMolBase *MarvinMulticenterSgroup::copyMol(std::string idAppendage) const {
  auto outSgroup = new MarvinMulticenterSgroup(this->parent);
  this->parent->sgroups.push_back(outSgroup);

  outSgroup->molID = this->molID + idAppendage;
  outSgroup->id = this->id;

  // the only time this is to be called is when a mutliple group above it is
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
  parent = parentInit;
}

std::string MarvinGenericSgroup::role() const {
  return std::string("GenericSgroup");
}

bool MarvinGenericSgroup::hasAtomBondBlocks() const { return false; }

MarvinMolBase *MarvinGenericSgroup::copyMol(std::string idAppendage) const {
  auto outSgroup = new MarvinGenericSgroup(this->parent);
  this->parent->sgroups.push_back(outSgroup);

  outSgroup->molID = this->molID + idAppendage;
  outSgroup->id = this->id + idAppendage;

  // the only time this is to be called is when a mutliple group above it is
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
  parent = parentInit;
}

std::string MarvinMonomerSgroup::role() const {
  return std::string("MonomerSgroup");
}

bool MarvinMonomerSgroup::hasAtomBondBlocks() const { return false; }

MarvinMolBase *MarvinMonomerSgroup::copyMol(std::string idAppendage) const {
  auto outSgroup = new MarvinMonomerSgroup(this->parent);
  this->parent->sgroups.push_back(outSgroup);

  outSgroup->molID = this->molID + idAppendage;
  outSgroup->id = this->id + idAppendage;
  outSgroup->title = this->title;
  outSgroup->charge = this->charge;

  // the only time this is to be called is when a mutliple group above it is
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
  parent = parentInit;
}

MarvinSuperatomSgroupExpanded::~MarvinSuperatomSgroupExpanded() {}

MarvinMolBase *MarvinSuperatomSgroupExpanded::copyMol(
    std::string idAppendage) const {
  auto outSgroup = new MarvinSuperatomSgroupExpanded(this->parent);
  this->parent->sgroups.push_back(outSgroup);

  outSgroup->molID = this->molID + idAppendage;
  outSgroup->id = this->id + idAppendage;
  outSgroup->title = this->title;

  auto actualParent = this->parent;
  while (actualParent->role() == "MultipleSgroup") {
    actualParent = actualParent->parent;
  }

  // the only time this is to be called is when a mutliple group above it is
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
  return std::string("SuperatomSgroupExpanded");
}

bool MarvinSuperatomSgroupExpanded::hasAtomBondBlocks() const { return false; }

MarvinSuperatomSgroup::MarvinSuperatomSgroup(MarvinMolBase *parentInit) {
  parent = parentInit;
}

MarvinSuperatomSgroup::~MarvinSuperatomSgroup() {
  for (auto &attachmentPoint : attachmentPoints) {
    delete attachmentPoint;
  }

  for (auto &atom : atoms) {
    delete atom;
    atom = nullptr;
  }
  for (auto &bond : bonds) {
    delete bond;
    bond = nullptr;
  }
}

std::string MarvinSuperatomSgroup::role() const {
  return std::string("SuperatomSgroup");
}

bool MarvinSuperatomSgroup::hasAtomBondBlocks() const { return true; }

MarvinMolBase *MarvinSuperatomSgroup::copyMol(std::string idAppendage) const {
  throw FileParseException("Internal error:  copying a SuperatomSgroup");
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
    for (auto attachmentPoint : attachmentPoints) {
      out << attachmentPoint->toString();
    }
    out << "</AttachmentPointArray>";
  }

  for (auto sgroup : sgroups) {
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
    for (auto attachmentPoint : attachmentPoints) {
      attachmentpointArray.add_child("attachmentPoint",
                                     attachmentPoint->toPtree());
    }
    out.put_child("AttachmentPointArray", attachmentpointArray);
  }
  MarvinMolBase::addSgroupsToPtree(out);

  return out;
}

MarvinMol::MarvinMol() { parent = nullptr; }

MarvinMol::~MarvinMol() {
  for (auto &atom : atoms) {
    delete atom;
    atom = nullptr;
  }
  for (auto &bond : bonds) {
    delete bond;
    bond = nullptr;
  }
}

std::string MarvinMol::role() const { return std::string("MarvinMol"); }

bool MarvinMol::hasAtomBondBlocks() const { return true; }

bool MarvinMolBase::atomRefInAtoms(MarvinAtom *a, std::string b) {
  return a->id == b;
}

bool MarvinMolBase::molIDInSgroups(std::string a, std::string b) {
  return a == b;
}

bool MarvinMolBase::bondRefInBonds(MarvinBond *a, std::string b) {
  return a->id == b;
}

void MarvinMolBase::cleanUpSgNumbering(int &sgCount) {
  for (auto *sgroup : this->sgroups) {
    std::string newId = "sg" + std::to_string(++sgCount);
    sgMap[sgroup->id] = newId;
    sgroup->id = newId;

    sgroup->cleanUpSgNumbering(sgCount);
  }
}

void MarvinMolBase::cleanUpNumberingMolsAtomsBonds(
    int &molCount    // this is the starting mol count, and receives the ending
                     // mol count - THis is used when
                     // MarvinMol->convertToSuperAtaoms is called multiple times
                     // from a RXN
    ,
    int &atomCount   // starting and ending atom count
    ,
    int &bondCount)  // starting and ending bond count)
{
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

  if (this->role() == "SuperatomSgroup") {
    auto marvinSuperatomSgroup = (MarvinSuperatomSgroup *)this;
    for (auto attachmentPoint : marvinSuperatomSgroup->attachmentPoints) {
      attachmentPoint->atom = atomMap[attachmentPoint->atom];
      attachmentPoint->bond =
          bondMap[attachmentPoint->bond];  // bond is actually in the parent
    }
  }

  // //Now the sgroups

  for (auto sgroup : sgroups) {
    sgroup->cleanUpNumberingMolsAtomsBonds(molCount, atomCount, bondCount);
  }
}

void MarvinMolBase::cleanUpNumbering(
    int &molCount   // this is the starting mol count, and receives the ending
                    // mol count - THis is used when
                    // MarvinMol->convertToSuperAtaoms is called multiple times
                    // from a RXN
    ,
    int &atomCount  // starting and ending atom count
    ,
    int &bondCount  // starting and ending bond count
    ,
    int &sgCount)   // starting and ending sg  count)
{
  this->cleanUpSgNumbering(sgCount);
  this->cleanUpNumberingMolsAtomsBonds(molCount, atomCount, bondCount);
}

bool MarvinMolBase::AnyOverLappingAtoms(const MarvinMolBase *otherMol) const {
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
  if (this->role() == "SuperatomSgroup" ||
      (this->role() == "MultipleSgroup" &&
       ((MarvinMultipleSgroup *)this)->isExpanded == false)) {
    return false;
  }
  return true;
}

// this routine determines if the sgroup is passrive for contraction, whern
// coming from an RDKit mol to a Marvin Mol
//  It is NOT passive if atoms in the parent are removed by contracting it.

bool MarvinMolBase::isPassiveRoleForContraction() const {
  if (this->role() == "MarvinMol" ||
      this->role() == "SuperatomSgroupExpanded" ||
      (this->role() == "MultipleSgroup" &&
       ((MarvinMultipleSgroup *)this)->isExpanded == true)) {
    return false;
  }
  return true;
}

void moveSgroup(MarvinMolBase *sgroupToMove, MarvinMolBase *newParent,
                bool removeFromOldParent = true) {
  newParent->sgroups.push_back(sgroupToMove);
  if (removeFromOldParent) {
    auto deleteIter = find(sgroupToMove->parent->sgroups.begin(),
                           sgroupToMove->parent->sgroups.end(), sgroupToMove);
    if (deleteIter == sgroupToMove->parent->sgroups.end()) {
      throw FileParseException(
          "Unexpected error - sgroup not found in its parent");
    }
    sgroupToMove->parent->sgroups.erase(
        find(sgroupToMove->parent->sgroups.begin(),
             sgroupToMove->parent->sgroups.end(), sgroupToMove));
  }
  sgroupToMove->parent = newParent;
}

void promoteChild(MarvinMolBase *molToMove) {
  moveSgroup(molToMove, molToMove->parent->parent);
}

void MarvinMultipleSgroup::expandOneMultipleSgroup() {
  // Mulitplesgroups are handled differently in Marvin and RDKit/mol file
  // format
  //
  // In Marvin, the atoms of the group are indicated, and the "title" contains
  // the number of replicates In mol files, the atoms of the group are
  // actually replicated in the atomBlock, and the MUL contructs indicates the
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
      originalConnectorAtoms;     // the bonds inthe parent that were to outside
                                  // atoms BEFORE replication
  std::string connectorAtoms[2];  // working atoms to connect the replicates
  std::string outsideConnectorAtom;
  std::vector<MarvinBond *> bondsInGroup;
  this->bondsToAtomsNotInExpandedGroup.clear();

  // find the actual parent - multiple groups can be nested, so the atoms and
  // bonds of the lower one are really in the grandparent or higher)

  MarvinMolBase *actualParent = this->parent;
  while (actualParent->role() == "MultipleSgroup") {
    actualParent = actualParent->parent;
  }

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

  std::vector<MarvinMolBase *> sgroupsToAdd;
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
      if (parentAtomPtr->sgroupRef != "" && !copyAtom->sGroupRefIsSuperatom) {
        copyAtom->sgroupRef = parentAtomPtr->sgroupRef + idAppendage;
      }

      // copy atoms to all parents up to the actual parent (multipleSgrups in
      // multipleSgroups)

      for (auto thisParent = this->parent;; thisParent = thisParent->parent) {
        auto insertItr = find(thisParent->atoms.begin(),
                              thisParent->atoms.end(), lastAtomInGroupPtr);
        if (insertItr != thisParent->atoms.end()) {
          insertItr++;  // insert after the last one
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
            insertItr++;  // insert after the last one
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
      actualParent->bonds.push_back(copyBond);
    }

    for (int i = this->sgroups.size() - 1; i >= 0; --i) {
      auto copiedMol = sgroups[i]->copyMol(idAppendage);
      moveSgroup(copiedMol, this->parent);
    }

    // add the bond from the last group to the new one

    if (this->bondsToAtomsNotInExpandedGroup.size() == 2) {
      connectorAtoms[0] = originalConnectorAtoms[0] + idAppendage;
      auto connectorBond = new MarvinBond(
          *this->bondsToAtomsNotInExpandedGroup[0],
          this->bondsToAtomsNotInExpandedGroup[0]->id + idAppendage,
          connectorAtoms[0], connectorAtoms[1]);
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
    promoteChild(sgroups[i]);
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
  // separate records called attachmentPoints that specifiy which atom(s) in
  // the super atom sub-mol that replaces the dummy atom(s), and also a bond
  // pointer (in the parent mol).
  //
  // In the expanded form, all of the atoms are in the parent molecule, and
  // the sub-mol only refers to them.  The attachement points refer to the
  // atoms in the parent mol. The MultipleSgroup and SruGroup seem very much
  // alike.   In both, all atoms are in the parent mol, and the sub-mol refers
  // to a group of atoms in that parent. The SruGroup specifies the name and
  // also the connection informat (head-to-tail, head-head, and unspecified).
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
    this->parent->sgroups.push_back(superatomSgroupExpanded);

    superatomSgroupExpanded->title = this->title;
    superatomSgroupExpanded->id =
        this->id;     // the expanded sgroup will replace the contracted one
    superatomSgroupExpanded->molID =
        this->molID;  // the expanded sgroup will replace the contracted one

    // find the actual parent - multiple groups can be nested, so the atoms
    // and bonds of the lower one are really in the grandparent or higher)

    MarvinMolBase *actualParent = this->parent;
    while (actualParent->role() == "MultipleSgroup") {
      actualParent = actualParent->parent;
    }

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
    // of the sibing too

    for (auto sibling : parent->sgroups) {
      if (sibling == this) {
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
            [orphanedBond](const MarvinAttachmentPoint *arg) {
              return arg->bond == orphanedBond->id;
            });
        if (attachmentPoint == this->attachmentPoints.end()) {
          throw FileParseException(
              "No attachment point found for bond to the condensed atom in a superatomSgroup");
        }
        orphanedBonds.push_back(*attachmentPoint);

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
      auto dummyInParent =
          find_if(thisParent->atoms.begin(), thisParent->atoms.end(),
                  [dummyAtomPtr](const MarvinAtom *arg) {
                    return arg == dummyAtomPtr;
                  });
      thisParent->atoms.erase(dummyInParent);  // get rid of the atoms pointer
                                               // to the old dummy atom

      // for multiple groups, also delete it from the parentAtoms array - that
      // multiple sgroup has NOT been done yet

      if (thisParent->role() == "MultipleSgroup") {
        auto thisMultipleParent = (MarvinMultipleSgroup *)thisParent;
        auto deleteIter =
            find(thisMultipleParent->parentAtoms.begin(),
                 thisMultipleParent->parentAtoms.end(), dummyAtomPtr);
        thisMultipleParent->parentAtoms.erase(deleteIter);
      }

      if (thisParent == actualParent) {
        break;
      }
    }

    delete dummyAtomPtr;  // get rid of the MolAtom

    // add the atoms and bonds from the super group to the parent

    for (auto subAtomPtr : this->atoms) {
      for (auto thisParent = parent;; thisParent = thisParent->parent) {
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

      subAtomPtr = nullptr;  // so that it is not in two places  - prevents
                             // double deleting
    }
    for (auto &bond : this->bonds) {
      actualParent->bonds.push_back(bond);
      bond = nullptr;  // so that it is not in two places  - prevents double
                       // deleting
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

      delete attachmentPoint;
    }

    this->attachmentPoints.clear();
    this->atoms.clear();
    this->bonds.clear();

    // promate any children sgroups

    for (auto childSgroup : this->sgroups) {
      moveSgroup(childSgroup, actualParent, false);
    }
    this->sgroups.clear();

    parent->sgroups.erase(
        std::find(parent->sgroups.begin(), parent->sgroups.end(), this));
    delete this;
  } catch (const std::exception &e) {
    throw;
  }
}

void MarvinMulticenterSgroup::processOneMulticenterSgroup() {
  MarvinMolBase *actualParent = this->parent;
  while (actualParent->role() == "MultipleSgroup") {
    actualParent = actualParent->parent;
  }

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
    delete orphanedBond;
    auto orphanedBondIter = std::find(actualParent->bonds.begin(),
                                      actualParent->bonds.end(), orphanedBond);
    actualParent->bonds.erase(orphanedBondIter);
  }

  auto centerIter = std::find(actualParent->atoms.begin(),
                              actualParent->atoms.end(), this->center);
  if (centerIter !=
      actualParent->atoms.end())  // it might already have been deleted by
                                  // another multicenter group
  {
    auto centerPtr = *centerIter;
    for (auto thisParent = parent;; thisParent = thisParent->parent) {
      thisParent->atoms.erase(find(
          thisParent->atoms.begin(), thisParent->atoms.end(),
          centerPtr));  // get rid of the atoms pointer to the old dummy atom
      if (thisParent->role() == "MultipleSgroup") {
        auto thisMultipleParent = (MarvinMultipleSgroup *)thisParent;
        auto deleteIter =
            find(thisMultipleParent->parentAtoms.begin(),
                 thisMultipleParent->parentAtoms.end(), centerPtr);
        thisMultipleParent->parentAtoms.erase(deleteIter);
      }

      if (thisParent == actualParent) {
        break;
      }
    }
    delete centerPtr;  // get rid of the MolAtom
  }

  // erase the multicenter group from its parent

  for (auto thisParent = parent;; thisParent = thisParent->parent) {
    thisParent->sgroups.erase(std::find(thisParent->sgroups.begin(),
                                        thisParent->sgroups.end(), this));

    if (thisParent == actualParent) {
      break;
    }
  }

  delete this;
}

void MarvinMolBase::prepSgroupsForRDKit() {
  // this routine recursively fixes the hierarchy of sgroups - some may NOT
  // actaully belong underneath their parent and can be promnoted to the
  // grandparent

  // first fix all the children

  std::vector<std::string> sgroupsMolIdsDone;
  for (bool allDone = false;
       !allDone;)    // until all at this level are done - but fixing one child
                     // can add sgroups to this level
  {
    allDone = true;  // unitl it is not true in the loop below
    for (auto childSgroup : this->sgroups) {
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

  std::string parentRole = this->parent->role();
  std::string childRole = this->role();

  if (parentRole != "MarvinMol") {
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
    }
  }

  // made to here, so the child belongs in the parent.

  if (childRole == "SuperatomSgroup") {
    ((MarvinSuperatomSgroup *)this)->convertFromOneSuperAtom();
  } else if (childRole == "MulticenterSgroup") {
    ((MarvinMulticenterSgroup *)this)->processOneMulticenterSgroup();
  } else if (childRole == "MultipleSgroup") {
    ((MarvinMultipleSgroup *)this)->expandOneMultipleSgroup();
  }

  return;
}

MarvinMolBase *MarvinSuperatomSgroupExpanded::convertToOneSuperAtom() {
  // the mol-style super atoms are significatnly different than the Marvin
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

  MarvinMolBase *actualParent = this->parent;
  while (actualParent->role() == "MultipleSgroup") {
    actualParent = actualParent->parent;
  }

  // make a new sub mol

  auto marvinSuperatomSgroup = new MarvinSuperatomSgroup(this->parent);
  this->parent->sgroups.push_back(marvinSuperatomSgroup);
  std::string newAtomName = "NA_" + this->molID;

  marvinSuperatomSgroup->molID = this->molID;
  marvinSuperatomSgroup->title = this->title;
  marvinSuperatomSgroup->id = this->id;

  bool coordsExist = has2dCoords();  // we have to check before we add the dummy
                                     // atom  - it will not have coords (yet)

  // add the dummy atom into the parent

  auto dummyParentAtom = new MarvinAtom();
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
      thisParent->atoms.erase(
          find(thisParent->atoms.begin(), thisParent->atoms.end(), atom));

      if (thisParent->role() == "MultipleSgroup") {
        auto marvinMultipleSgroup = (MarvinMultipleSgroup *)thisParent;
        marvinMultipleSgroup->parentAtoms.erase(
            find(marvinMultipleSgroup->parentAtoms.begin(),
                 marvinMultipleSgroup->parentAtoms.end(), atom));
      }

      if (thisParent == actualParent) {
        break;
      }
    }

    if (coordsExist)  // get the center of all atoms in the group - we might
                      // use this if there are no attachement points
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

      // add an attachentPoint structure

      auto marvinAttachmentPoint = new MarvinAttachmentPoint();
      marvinSuperatomSgroup->attachmentPoints.push_back(marvinAttachmentPoint);
      marvinAttachmentPoint->atom = atomPtr->id;
      marvinAttachmentPoint->bond = bond->id;
      ;
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

  for (auto sgroupPtr : this->sgroups) {
    moveSgroup(sgroupPtr, marvinSuperatomSgroup, false);
  }
  this->sgroups.clear();

  // remove the old expanded superatom sgroup

  parent->sgroups.erase(
      std::find(parent->sgroups.begin(), parent->sgroups.end(), this));
  delete this;

  // fix up any siblings that contain all the atoms of this group.

  for (auto sibling : marvinSuperatomSgroup->parent->sgroups) {
    if (sibling == marvinSuperatomSgroup) {
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
        // remove the group atoms from the sibling and add the dumm atom

        sibling->atoms.push_back(dummyParentAtom);
        for (auto atomToRemove : marvinSuperatomSgroup->atoms) {
          int index = sibling->getAtomIndex(atomToRemove->id);
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
    if (*testBond == NULL) {
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

    *testBond = NULL;

    // try the children

    int childId =
        getMatchedOrphanBondIndex(otherAtomId, bondsToTry, orphanedBonds);
    if (childId >= 0) {
      return childId;
    }
  }

  return (-1);
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

  MarvinMolBase *actualParent = this->parent;
  while (actualParent->role() == "MultipleSgroup") {
    actualParent = actualParent->parent;
  }

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

  for (auto childSgroup : sgroups) {
    if (childSgroup->isSgroupInSetOfAtoms(atomsToDelete) !=
        SgroupNotInAtomSet) {
      sgroupsToDelete.push_back(childSgroup);
    }
  }

  // now fix the orphaned bonds - The first gets the atoms from the matched
  // second orphan bond  and was NOT removed. the matched second orphaned bond
  // is deleted

  while (orphanedBonds.size() > 0) {
    int matchedOrphanBondIndex = (-1);
    auto orphanedBondToFix = orphanedBonds[0];
    orphanedBonds.erase(orphanedBonds.begin());

    if (orphanedBonds.size() == 3)  // was 4 but the first one was erased
    {
      std::vector<MarvinBond *> bondsToTry =
          bondsToDelete;  // copy of bonds to delete
      bondsToTry.erase(
          find(bondsToTry.begin(), bondsToTry.end(), orphanedBondToFix));
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
            ->atomRefs2[0];  // [0] is the orpahened atom (still in the mol)

    // undelete the bond which has been fixed (really remove it from the list
    // of bonds to be delete)

    bondsToDelete.erase(
        find(bondsToDelete.begin(), bondsToDelete.end(), orphanedBondToFix));

    // any siblings or children that reference the matched Orphaned bond must
    // be fixed to reference the retained (fixed) orphan bond

    for (auto siblingSgroup :
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

    for (auto childSgroup : this->sgroups) {
      for (auto childBond : childSgroup->bonds) {
        if (childBond->id == orphanedBonds[matchedOrphanBondIndex]->id) {
          childBond = orphanedBondToFix;
        }
      }
    }

    orphanedBonds.erase(orphanedBonds.begin() + matchedOrphanBondIndex);
  };

  // if we get here, all the changes to be made have been determined. so
  // make them

  for (auto childSgroup : sgroupsToDelete) {
    sgroups.erase(std::find(sgroups.begin(), sgroups.end(), childSgroup));
    delete childSgroup;
  }

  // remove the atoms

  for (auto atomPtr : atomsToDelete) {
    for (auto thisParent = this->parent;;
         thisParent = thisParent->parent) {  // do all parents and grandparents
                                             // ... to to the actual parent

      thisParent->atoms.erase(std::find(thisParent->atoms.begin(),
                                        thisParent->atoms.end(), atomPtr));

      if (thisParent == actualParent) {
        break;
      }
    }

    // also delete the atoms from this multiple group AND any sibling sgroups

    for (auto siblingSgroup : this->parent->sgroups) {
      auto siblingAtomPtr = std::find(siblingSgroup->atoms.begin(),
                                      siblingSgroup->atoms.end(), atomPtr);
      if (siblingAtomPtr != siblingSgroup->atoms.end()) {
        siblingSgroup->atoms.erase(siblingAtomPtr);
      }
    }

    delete atomPtr;
  }

  atomsToDelete.clear();

  // remove the bonds

  for (MarvinBond *bondPtr : bondsToDelete) {
    actualParent->bonds.erase(std::find(actualParent->bonds.begin(),
                                        actualParent->bonds.end(), bondPtr));
    delete bondPtr;
  }
  bondsToDelete.clear();

  this->isExpanded = false;
}

IsSgroupInAtomSetResult MarvinMolBase::isSgroupInSetOfAtoms(
    std::vector<MarvinAtom *> &setOfAtoms) const {
  // superatom sgrups are different - they are in the set if one of the
  // condensed atom in the parent is in group

  if (this->role() == "SuperatomSgroup") {
    auto dummyAtomIter = find_if(
        parent->atoms.begin(), parent->atoms.end(),
        [this](const MarvinAtom *arg) { return arg->sgroupRef == this->id; });
    if (dummyAtomIter == parent->atoms.end()) {
      throw FileParseException(
          "No contracted atom found for a superatomSgroup");
    }

    if (boost::algorithm::contains(
            setOfAtoms, std::vector<std::string>{(*dummyAtomIter)->id},
            atomRefInAtoms)) {
      return SgroupInAtomSet;
    } else {
      return SgroupNotInAtomSet;
    }
  }

  // all other types are in the group based on their atoms

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

void MarvinMolBase::processSgroupsFromRDKit() {
  // this routine recursively fixes sgroups to be in the heirarchy expected by
  // Marvin format

  if (this->isPassiveRoleForContraction()) {
    return;  // passive mols are only leaves of the tree - nothing to do
  }

  // see if any siblings should be chilren of this one

  MarvinMolBase *molToProcess =
      this;  // the processing MIGHT change the mol from one kind to another
             // (e.g. MarvinSuperatomSgroupExpanded -> MarvinSuperatomSgroup)

  if (this->parent != nullptr)  // NOT the top level mol
  {
    std::vector<std::string> sgroupsMolIdsDone;
    for (bool allDone = false;
         !allDone;)    // until all at this level are done - but fixing one
                       // child can add sgroups to this level
    {
      allDone = true;  // unitl it is not true in the loop below
      for (auto sibling : this->parent->sgroups) {
        if (boost::algorithm::contains(sgroupsMolIdsDone,
                                       std::vector<std::string>{sibling->molID},
                                       molIDInSgroups)) {
          continue;  // look at the next one to see if it is done
        }

        allDone = false;
        sgroupsMolIdsDone.push_back(sibling->molID);

        if (sibling == this) {
          continue;  // cannot be a child of itself
        }

        //  there are 2 possibilities:
        //  1) all atoms are in this one - make the sibling into a child of
        //  this one 2) SOME atoms are NOT this one - it remains a sibling

        if (sibling->isSgroupInSetOfAtoms(this->atoms) != SgroupInAtomSet) {
          continue;
        }

        moveSgroup(sibling, this);
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
            << "The Mutliple sgroup will be ignored" << std::endl;

        this->parent->sgroups.erase(std::find(
            this->parent->sgroups.begin(), this->parent->sgroups.end(), this));
        delete this;
        return;
      }
    }
  }

  // now fix all this groups children

  std::vector<std::string> chilrenSgroupsMolIdsDone;
  for (bool allDone = false;
       !allDone;)    // until all at this level are done - but fixing one child
                     // can add sgroups to this level
  {
    allDone = true;  // unitl it is not true in the loop below
    for (auto childSgroup : molToProcess->sgroups) {
      if (boost::algorithm::contains(
              chilrenSgroupsMolIdsDone,
              std::vector<std::string>{childSgroup->molID}, molIDInSgroups)) {
        continue;  // this one is done, look at the next one
      }

      chilrenSgroupsMolIdsDone.push_back(childSgroup->molID);
      childSgroup->processSgroupsFromRDKit();
      allDone = false;
      break;  // have to start the loop over - we might have changed the
              // vector
    }
  }

  return;
}

MarvinMolBase *MarvinMol::copyMol(std::string idAppendage) const {
  throw FileParseException("Internal error:  copying a MarvinMol");
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

  for (auto sgroup : sgroups) {
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

MarvinReaction::~MarvinReaction() {
  for (auto &reactant : reactants) {
    delete reactant;
  }
  for (auto &agent : agents) {
    delete agent;
  }
  for (auto &product : products) {
    delete product;
  }
  for (auto &pluse : pluses) {
    delete pluse;
  }
  for (auto &condition : conditions) {
    delete condition;
  }
}

void MarvinReaction::prepSgroupsForRDKit() {
  // This routine converts all the mols in the rxn to be ready for conversion
  // to RDKIT mols
  int molCount = 0, atomCount = 0, bondCount = 0, sgCount = 0;
  sgMap.clear();
  atomMap.clear();
  bondMap.clear();

  for (auto &reactant : reactants) {
    reactant->prepSgroupsForRDKit();
    // clean and renumber the sgroups
    reactant->clearMaps();
    reactant->cleanUpNumbering(molCount, atomCount, bondCount, sgCount);
  }
  for (auto &agent : agents) {
    agent->prepSgroupsForRDKit();
    agent->cleanUpNumbering(molCount, atomCount, bondCount, sgCount);
  }
  for (auto &product : products) {
    product->prepSgroupsForRDKit();
    product->cleanUpNumbering(molCount, atomCount, bondCount, sgCount);
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

  for (auto &pluse : pluses) {
    out.add_child("cml.MDocument.MReactionSign", pluse->toPtree());
  }

  for (auto &condition : conditions) {
    out.add_child("cml.MDocument.MTextBox", condition->toPtree());
  }

  return out;
}

MarvinRectangle::MarvinRectangle() {
  upperLeft.x = 0.0;
  upperLeft.y = 0.0;
  lowerRight.x = 0.0;
  lowerRight.y = 0.0;
  centerIsStale = true;
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

MarvinRectangle::MarvinRectangle(const std::vector<MarvinAtom *> atoms) {
  centerIsStale = true;

  if (atoms.size() == 0) {
    return;
  }
  upperLeft.x = DBL_MAX;
  upperLeft.y = -DBL_MAX;
  lowerRight.x = -DBL_MAX;
  lowerRight.y = DBL_MAX;

  for (auto atom : atoms) {
    if (atom->x2 < upperLeft.x) {
      upperLeft.x = atom->x2;
    }
    if (atom->x2 > lowerRight.x) {
      lowerRight.x = atom->x2;
    }

    if (atom->y2 > upperLeft.y) {
      upperLeft.y = atom->y2;
    }
    if (atom->y2 < lowerRight.y) {
      lowerRight.y = atom->y2;
    }
  }
}

MarvinRectangle::MarvinRectangle(const std::vector<MarvinRectangle> rects) {
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
