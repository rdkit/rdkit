#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include "FileParsers/FileParsers.h"
#include "FileParsers/FileParserUtils.h"

#include "MACROMol.h"
#include "Atom.h"

namespace RDKit {

void RDKit::MACROMolTemplate::findMainSgroupForTemplate(
    std::string className, std::string templateName) const {
  p_mainSgroupIdx = UINT_MAX;

  auto sgroups = RDKit::getSubstanceGroups(*this);

  for (unsigned int sgroupIdx = 0; sgroupIdx < sgroups.size(); ++sgroupIdx) {
    auto &sgroupToTest = sgroups[sgroupIdx];
    std::string sup;
    std::string sgroupAtomClass;
    if (sgroupToTest.getPropIfPresent("TYPE", sup) && sup == "SUP" &&
        sgroupToTest.getPropIfPresent("CLASS", sgroupAtomClass) &&
        sgroupAtomClass == className) {
      p_mainSgroupIdx = sgroupIdx;
      break;
    }
  }

  if (p_mainSgroupIdx == UINT_MAX) {
    std::ostringstream errout;
    errout << "No main sgroup found for template: " << className << "/"
           << templateName;
    throw RDKit::FileParseException(errout.str());
  }
}

void RDKit::MACROMolTemplate::init(std::string className,
                                   std::vector<std::string> templateNames,
                                   std::vector<std::string> templateAttrs) {
  PRECONDITION(className.empty() == false, "no className for template");
  PRECONDITION(templateNames.size() > 0, "no template names for template");

  this->setProp(RDKit::common_properties::molAtomClass, className);
  this->setProp(RDKit::common_properties::templateNames, templateNames);

  for (auto templateAttr : templateAttrs) {
    std::vector<std::string> subTokens;
    boost::algorithm::split(subTokens, templateAttr,
                            boost::algorithm::is_any_of("="));
    if (subTokens.size() != 2) {
      std::ostringstream errout;
      errout << "Attribute string is not of the form \"AttrName=value\": "
             << templateAttr;
      throw RDKit::FileParseException(errout.str());
    }
    this->setProp(subTokens[0], subTokens[1]);
  }
  p_mainSgroupIdx = UINT_MAX;
}

MACROMolTemplate::MACROMolTemplate(std::unique_ptr<RWMol> &mol,
                                   std::string className,
                                   std::vector<std::string> templateNames,
                                   std::vector<std::string> templateAttrs)
    : RWMol(std::move(*mol)) {
  init(className, templateNames, templateAttrs);
}

MACROMolTemplate::MACROMolTemplate(const MACROMolTemplate &other)
    : RWMol(other) {
  // p_mol.reset(new RWMol(*other.p_mol));

  p_mainSgroupIdx = UINT_MAX;
}

MACROMolTemplate::MACROMolTemplate(std::unique_ptr<RWMol> &mol,
                                   std::string className,
                                   std::string templateName,
                                   std::vector<std::string> templateAttrs)
    : RWMol(std::move(*mol)) {
  PRECONDITION(!templateName.empty(), "no name for template");

  std::vector<std::string> templateNames(1, templateName);
  MACROMolTemplate::init(className, templateNames, templateAttrs);
}

unsigned int RDKit::MACROMol::addMacroAtom(std::string className,
                                           std::string templateName) {
  auto atom = new Atom(0);
  atom->setAtomicNum(0);

  atom->setProp(common_properties::dummyLabel, templateName);
  atom->setProp(common_properties::molAtomClass, className);
  p_atomIdxToTemplateIdxIsStale = true;
  return this->addAtom(atom, false, true);
}

void RDKit::MACROMol::addMacroBond(unsigned int fromAtomIdx,
                                   unsigned int toAtomIdx,
                                   Bond::BondType bondType,
                                   std::string fromConnectionPoint,
                                   std::string toConnectionPoint) {
  auto bondIdx = this->addBond(fromAtomIdx, toAtomIdx, bondType);
  auto bond = this->getBondWithIdx(bondIdx);
  this->setBondBookmark(bond, bondIdx);

  if (!toConnectionPoint.empty()) {
    bond->setProp(common_properties::_MolFileBondAttachPt2, toConnectionPoint);
  }

  if (!fromConnectionPoint.empty()) {
    bond->setProp(common_properties::_MolFileBondAttachPt1,
                  fromConnectionPoint);
  }
}

unsigned int RDKit::MACROMol::atomIdxToTemplateIdx(unsigned int atomIdx) {
  if (p_atomIdxToTemplateIdxIsStale) {
    // rebuild the map
    p_atomIdxToTemplateIdx.clear();

    for (auto atom : this->atoms()) {
      std::string atomClass;
      std::string dummyLabel = "";

      if (!atom->getPropIfPresent(common_properties::dummyLabel, dummyLabel) ||
          dummyLabel == "" ||
          !atom->getPropIfPresent(common_properties::molAtomClass, atomClass) ||
          atomClass == "") {
        p_atomIdxToTemplateIdx[atom->getIdx()] = UINT_MAX;
        continue;
      }

      for (unsigned int tIdx = 0; tIdx < getTemplateCount(); ++tIdx) {
        const MACROMolTemplate *templateMol = getTemplate(tIdx);
        std::string templateAtomClass;
        std::vector<std::string> templateNames;
        templateMol->getPropIfPresent<std::string>(
            common_properties::molAtomClass, templateAtomClass);
        templateMol->getPropIfPresent<std::vector<std::string>>(
            common_properties::templateNames, templateNames);
        if (templateAtomClass == atomClass &&
            std::find(templateNames.begin(), templateNames.end(), dummyLabel) !=
                templateNames.end()) {
          p_atomIdxToTemplateIdx[atom->getIdx()] = tIdx;
          break;
        }
      }

      if (p_atomIdxToTemplateIdx.contains(atom->getIdx()) == false) {
        throw FileParseException("Template for macro atom not found");
      }
    }

    p_atomIdxToTemplateIdxIsStale = false;
  }
  return p_atomIdxToTemplateIdx[atomIdx];
}

MACROMolTemplate *RDKit::MACROMol::atomIdxToMACROMolTemplate(
    unsigned int atomIdx) {
  return this->getTemplate(this->atomIdxToTemplateIdx(atomIdx));
}

}  // namespace RDKit
