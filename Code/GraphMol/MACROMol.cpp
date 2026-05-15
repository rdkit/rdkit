#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include "FileParsers/FileParsers.h"
#include "FileParsers/FileParserUtils.h"


#include "MACROMol.h"
#include "Atom.h"

namespace RDKit {

void RDKit::MACROMolTemplate::findMainSgroupForTemplate(
    std::string className, std::string templateName) const {
  d_mainSgroupIdx = UINT_MAX;

  auto sgroups = RDKit::getSubstanceGroups(*this);

  for (unsigned int sgroupIdx = 0; sgroupIdx < sgroups.size(); ++sgroupIdx) {
    auto &sgroupToTest = sgroups[sgroupIdx];
    std::string sup;
    std::string sgroupAtomClass;
    if (sgroupToTest.getPropIfPresent("TYPE", sup) && sup == "SUP" &&
        sgroupToTest.getPropIfPresent("CLASS", sgroupAtomClass) &&
        sgroupAtomClass == className) {
      d_mainSgroupIdx = sgroupIdx;
      break;
    }
  }

  if (d_mainSgroupIdx == UINT_MAX) {
    std::ostringstream errout;
    errout << "No main sgroup found for template: " << className << "/"
           << templateName;
    throw RDKit::FileParseException(errout.str());
  }
}

void RDKit::MACROMolTemplate::init(std::string className,
                                   std::vector<std::string> templateNames,
                                   std::vector<std::pair<std::string,std::string>> templateAttrs) {
  PRECONDITION(className.empty() == false, "no className for template");
  PRECONDITION(templateNames.size() > 0, "no template names for template");

  this->setProp(RDKit::common_properties::molAtomClass, className);
  this->setProp(RDKit::common_properties::templateNames, templateNames);

  for (auto templateAttr : templateAttrs) {
    this->setProp(templateAttr.first, templateAttr.second);

  }
  d_mainSgroupIdx = UINT_MAX;
}

MACROMolTemplate::MACROMolTemplate(std::unique_ptr<RWMol> &mol,
                                   std::string className,
                                   std::vector<std::string> templateNames,
                                   std::vector<std::pair<std::string,std::string>> templateAttrs)
    : RWMol(std::move(*mol)) {
  init(className, templateNames, templateAttrs);
}

MACROMolTemplate::MACROMolTemplate(const MACROMolTemplate &other)
    : RWMol(other) {
  d_mainSgroupIdx = UINT_MAX;
}

MACROMolTemplate::MACROMolTemplate(std::unique_ptr<RWMol> &mol,
                                   std::string className,
                                   std::string templateName,
                                   std::vector<std::pair<std::string,std::string>> templateAttrs)
    : RWMol(std::move(*mol)) {
  PRECONDITION(!templateName.empty(), "no name for template");

  std::vector<std::string> templateNames(1, templateName);
  MACROMolTemplate::init(className, templateNames, templateAttrs);
}

 RDKit::SubstanceGroup *MACROMolTemplate::getMainSgroup() {
    
    std::call_once(d_mainSgroupIdxOnceFlag, [this]() {
      std::string className = "";
      std::vector<std::string> templateNames;
      if (!getPropIfPresent(RDKit::common_properties::molAtomClass,
                                  className) ||
          !getPropIfPresent(RDKit::common_properties::templateNames,
                                  templateNames)) {
        std::ostringstream errout;
        errout << "Template molecule is missing required properties: "
               << RDKit::common_properties::molAtomClass << " and/or "
               << RDKit::common_properties::templateNames;
        throw RDKit::FileParseException(errout.str());
      }
   

      findMainSgroupForTemplate(className, templateNames[0]);
    });
    return &RDKit::getSubstanceGroups(*this)[d_mainSgroupIdx];
  }



unsigned int RDKit::MACROMol::addMacroAtom(std::string className,
                                           std::string templateName) {
  auto atom = new Atom(0);
  atom->setAtomicNum(0);

  atom->setProp(common_properties::dummyLabel, templateName);
  atom->setProp(common_properties::molAtomClass, className);
  d_atomIdxToTemplatePtrIsStale = true;
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

 MACROMolTemplate *RDKit::MACROMol::atomIdxToTemplatePtr(unsigned int atomIdx) const {
  if (d_atomIdxToTemplatePtrIsStale) {
    // rebuild the map
    d_atomIdxToTemplatePtr.clear();

    for (auto atom : this->atoms()) {
      std::string atomClass;
      std::string dummyLabel = "";

      d_atomIdxToTemplatePtr[atom->getIdx()] = nullptr;

      if (!atom->getPropIfPresent(common_properties::dummyLabel, dummyLabel) ||
          dummyLabel == "" ||
          !atom->getPropIfPresent(common_properties::molAtomClass, atomClass) ||
          atomClass == "") {
        continue;
      }

      // first look in the local temmplates, if any
      bool foundTemplate = false;
      if (d_templateLibrary.libContains(atomClass, dummyLabel)) {
        d_atomIdxToTemplatePtr[atom->getIdx()]  = this->d_templateLibrary.getTemplate(atomClass, dummyLabel);
        foundTemplate = true;
      } else {
        // try all the external libs

        for (const auto extLib :  d_externalTemplateLibs) {
          if( extLib->libContains(atomClass, dummyLabel)) {
            d_atomIdxToTemplatePtr[atom->getIdx()]  = extLib->getTemplate(atomClass, dummyLabel);
            foundTemplate = true;
            break;
          }
        }
      }
      if (!foundTemplate) {
          std::ostringstream errout;
          errout << "Template not found for atom " << dummyLabel;
          throw RDKit::FileParseException(errout.str());
      }
    }


    d_atomIdxToTemplatePtrIsStale = false;
  }
  return d_atomIdxToTemplatePtr[atomIdx];
}

void MACROMolTemplateLib::addTemplate(std::unique_ptr<MACROMolTemplate> &templateMolToAdd) {
    PRECONDITION(templateMolToAdd, "bad template molecule");

    this->push_back(std::move(templateMolToAdd));

    auto templateMol = this->back().get();
    std::string templateClass;
    std::vector<std::string> templateNames;
    templateMol->getPropIfPresent<std::string>(
        common_properties::molAtomClass, templateClass);
    templateMol->getPropIfPresent<std::vector<std::string>>(
        common_properties::templateNames, templateNames);

    for (auto templateName : templateNames){
      MACROMolTemplateKey key(std::pair(templateClass, templateName));

      d_keyToIndex[key] = this->size() - 1;
    }
  }

  void MACROMolTemplateLib::addTemplateLib(MACROMolTemplateLib &libToAdd) {
    for (auto &libTemplate : libToAdd) {
      addTemplate(libTemplate);
    }
  }

  void MACROMolTemplateLib::copyTemplateLib(const MACROMolTemplateLib &libToCopy){
    this->clearTemplateLib();

  
    
    for (auto &templateToCopy : libToCopy) {
      // make a copy of the template

      auto templateCopy = std::unique_ptr<MACROMolTemplate>(new MACROMolTemplate(*(templateToCopy.get())));
      this->addTemplate(templateCopy);

    }

  }
}  // namespace RDKit
