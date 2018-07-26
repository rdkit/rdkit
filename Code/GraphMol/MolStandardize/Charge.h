#ifndef __RD_CHARGE_H__
#define __RD_CHARGE_H__

#include <Catalogs/Catalog.h>
#include <GraphMol/MolStandardize/AcidBaseCatalog/AcidBaseCatalogEntry.h>
#include <GraphMol/MolStandardize/AcidBaseCatalog/AcidBaseCatalogParams.h>

namespace RDKit {
class RWMol;
class ROMol;

namespace MolStandardize {

typedef RDCatalog::HierarchCatalog<AcidBaseCatalogEntry, AcidBaseCatalogParams,
                                   int>
    AcidBaseCatalog;

struct ChargeCorrection {
  std::string Name;
  std::string Smarts;
  int Charge;

  ChargeCorrection(std::string name, std::string smarts, int charge)
      : Name(name), Smarts(smarts), Charge(charge) {}
};

// The default list of ChargeCorrections.
extern std::vector<ChargeCorrection> CHARGE_CORRECTIONS;

class Reionizer {
  // A class to fix charges and reionize a molecule such that the strongest
  // acids
  //  ionize first.

 public:
  ROMol *reionize(const ROMol &mol, AcidBaseCatalog *abcat,
                  std::vector<ChargeCorrection> ccs = CHARGE_CORRECTIONS);

 private:
  std::pair<unsigned int, std::vector<unsigned int>> *strongestProtonated(
      const ROMol &mol,
      const std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> &abpairs);
  std::pair<unsigned int, std::vector<unsigned int>> *weakestIonized(
      const ROMol &mol,
      const std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> &abpairs);

};  // Reionizer class

class Uncharger {
  // Class for neutralizing ionized acids and bases.

 public:
  Uncharger();
  Uncharger(const Uncharger &other);
  ~Uncharger();

  ROMol *uncharge(const ROMol &mol);

 private:
  std::shared_ptr<ROMol> pos_h;
  std::shared_ptr<ROMol> pos_quat;
  std::shared_ptr<ROMol> neg;
  std::shared_ptr<ROMol> neg_acid;
};  // Uncharger class

}  // namespace MolStandardize
}  // namespace RDKit
#endif
