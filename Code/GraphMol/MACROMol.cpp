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
                                   std::vector<std::pair<std::string,std::string>> templateAttrs) {
  PRECONDITION(className.empty() == false, "no className for template");
  PRECONDITION(templateNames.size() > 0, "no template names for template");

  this->setProp(RDKit::common_properties::molAtomClass, className);
  this->setProp(RDKit::common_properties::templateNames, templateNames);

  for (auto templateAttr : templateAttrs) {
    this->setProp(templateAttr.first, templateAttr.second);

  }
  p_mainSgroupIdx = UINT_MAX;
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
  p_mainSgroupIdx = UINT_MAX;
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
    
    std::call_once(p_mainSgroupIdxOnceFlag, [this]() {
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
    return &RDKit::getSubstanceGroups(*this)[p_mainSgroupIdx];
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

      p_atomIdxToTemplateIdx[atom->getIdx()]  = this->d_templateLibrary.getMACROMolTemplateIndex(atomClass, dummyLabel);
  
    }

    p_atomIdxToTemplateIdxIsStale = false;
  }
  return p_atomIdxToTemplateIdx[atomIdx];
}

MACROMolTemplate *RDKit::MACROMol::atomIdxToMACROMolTemplate(
    unsigned int atomIdx) {
  return this->getTemplate(this->atomIdxToTemplateIdx(atomIdx));
}

 void MACROMolTemplateLib::addTemplate(std::unique_ptr<MACROMolTemplate> &templateMol, bool takeOwnership) {
    PRECONDITION(templateMol, "bad template molecule");

    addTemplate(templateMol.get());

    if (takeOwnership) {
      this->ownedTemplates.emplace_back(std::move(templateMol));
    }
 }


  void MACROMolTemplateLib::addTemplate(MACROMolTemplate *templateMol) {
    PRECONDITION(templateMol, "bad template molecule");

    this->push_back(templateMol);


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

  void MACROMolTemplateLib::addTemplateLib(MACROMolTemplateLib &libToAdd, bool takeOwnership) {
    for (auto libTemplate : libToAdd) {
      addTemplate(libTemplate);
    }

    if (takeOwnership) {
      for (auto &libTemplate : libToAdd.ownedTemplates) {
        ownedTemplates.emplace_back(std::move(libTemplate));
      }
      libToAdd.ownedTemplates.clear();
    }
  }

  void MACROMolTemplateLib::copyTemplateLib(const MACROMolTemplateLib &libToCopy){
    this->clearTemplateLib();

    // copy and add the ownedTemplates
    for (const auto &templateToCopy :libToCopy.ownedTemplates) {
      // make a copy of the template

      auto templateCopy = std::unique_ptr<MACROMolTemplate>(new MACROMolTemplate(*(templateToCopy.get())));
      this->addTemplate(templateCopy, true /* take ownership*/);

    }

    //Copy the non-owned ones (from a global or external library)

    for (auto templateToCopy : libToCopy) {
      if (libToCopy.isTemplateOwnedByLib(templateToCopy)) {
        auto templateCopy = std::unique_ptr<MACROMolTemplate>(templateToCopy);
        this->addTemplate(templateCopy, true /* take ownership*/);
      } else {
        this->addTemplate(templateToCopy);
      }
    }

  }

}  // namespace RDKit
