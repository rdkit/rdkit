#include "Charge.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <boost/range/adaptor/reversed.hpp>

namespace RDKit {
namespace MolStandardize {

// The default list of ChargeCorrections.
std::vector<ChargeCorrection> CHARGE_CORRECTIONS = 
	{ ChargeCorrection("[Li,Na,K]", "[Li,Na,K;X0+0]", 1),
	  ChargeCorrection("[Mg,Ca]", "[Mg,Ca;X0+0]", 2),
	  ChargeCorrection("[Cl]", "[Cl;X0+0]", -1)
	};

ROMol* Reionizer::reionize(const ROMol &mol, AcidBaseCatalog *abcat, 
								std::vector<ChargeCorrection> ccs) {
	PRECONDITION(abcat, "");
	const AcidBaseCatalogParams *abparams = abcat->getCatalogParams();

	PRECONDITION(abparams, "");
	const std::vector<std::pair<ROMol*, ROMol*>> &abpairs = abparams->getPairs();
	const std::vector<std::pair<ROMol*, ROMol*>>* p_abpairs = &abpairs;

	
	ROMol* omol = new ROMol(mol);
	int start_charge = MolOps::getFormalCharge(*omol);

	for (const auto &cc : ccs) {
		RDKit::MatchVectType res;
		std::cout << cc.Smarts << std::endl;
		unsigned int matches = SubstructMatch(*omol, *(SmartsToMol(cc.Smarts)), res);
		if ( matches != 0 ) {
			for (auto &pair : res) {
				// std::cout << pair.first << ", " << pair.second << std::endl; 
				Atom* atom = omol->getAtomWithIdx(pair.second);
				// std::cout << atom->getSymbol() << std::endl;
				std::cout << "Applying charge correction " << cc.Name << " " << 
								atom->getSymbol() << " " << cc.Charge << std::endl;
				atom->setFormalCharge(cc.Charge);
			}
		}
	}
	std::cout << MolToSmiles(*omol) << std::endl;
	int current_charge = MolOps::getFormalCharge(*omol);
	int charge_diff = current_charge - start_charge;
	std::cout << "Current charge: " << current_charge << std::endl;
	std::cout << "Charge diff: " << charge_diff << std::endl;

	// If molecule is now neutral, assume everything is now fixed
	// But otherwise, if charge has become more positive, 
	// look for additional protonated acid groups to ionize
	
	if (current_charge != 0) {
		while (charge_diff > 0) {
			// returns the acid strength ranking (ppos) 
			// and the substruct match (poccur) in a pair
			std::pair<unsigned int, std::vector<unsigned int>>* res = 
							this->strongestProtonated(mol, p_abpairs);
			if (res == nullptr) {break;}
			else {
				unsigned int ppos = res->first;
				std::vector<unsigned int> poccur = res->second;
				std::string abname;
				std::pair<ROMol*, ROMol*> abpair = (*p_abpairs)[ppos];
				(abpair.first)->getProp(common_properties::_Name, abname);
				std::cout << "Ionizing " << abname << 
								" to balance previous charge corrections" << std::endl;
				Atom* patom = omol->getAtomWithIdx(poccur.back());
//				std::cout << patom->getSymbol() << std::endl;
				patom->setFormalCharge( patom->getFormalCharge() - 1 );
				
				if (patom->getNumExplicitHs() > 0) {
					patom->setNumExplicitHs( patom->getNumExplicitHs() - 1 );
				}

				patom->updatePropertyCache();
				--charge_diff;
			}
		}
	}

	std::cout << MolToSmiles(*omol) << std::endl;
	std::cout << "Charge diff: " << charge_diff << std::endl;

	std::set< std::vector<int> > already_moved;
	while (true) {
		
	std::pair<unsigned int, std::vector<unsigned int>>* sp_res =
              this->strongestProtonated(*omol, p_abpairs);
	std::pair<unsigned int, std::vector<unsigned int>>* wi_res =
              this->weakestIonized(*omol, p_abpairs);
	if ( sp_res != nullptr && wi_res != nullptr ) {
		unsigned int ppos = sp_res->first;
		unsigned int ipos = wi_res->first;
		std::vector<unsigned int> poccur = sp_res->second;
		std::vector<unsigned int> ioccur = wi_res->second;
		std::cout << "ppos: " << ppos << std::endl;
		std::cout << "ipos: " << ipos << std::endl;
		if ( ppos < ipos ) {
			if ( poccur.back() == ioccur.back() ) {
				// Bad! H wouldn't be moved, resulting in infinite loop.
				std::cout << "Aborted reionization due to unexpected situation" << std::endl;
				break;
			}

			std::vector<int> key = {poccur.back(), ioccur.back()};
			std::sort(key.begin(), key.end());
			
			already_moved.insert(key);

			std::string prot_name, ionized_name;
			std::pair<ROMol*, ROMol*> prot_pair = (*p_abpairs)[ppos];
			std::pair<ROMol*, ROMol*> ionized_pair = (*p_abpairs)[ipos];
			(prot_pair.first)->getProp(common_properties::_Name, prot_name);
			(ionized_pair.first)->getProp(common_properties::_Name, ionized_name);

			std::cout << "Moved proton from " << prot_name << " to " 
							<< ionized_name << std::endl;
			// Remove hydrogen from strongest protonated
			Atom* patom = omol->getAtomWithIdx(poccur.back());
			patom->setFormalCharge( patom->getFormalCharge() - 1 );
			// If no implicit Hs to autoremove, and at least 1 explicit H to remove,
			//  reduce explicit count by 1
			if ( patom->getNumImplicitHs() == 0 && patom->getNumExplicitHs() > 0 ) {
				patom->setNumExplicitHs( patom->getNumExplicitHs() - 1 );
				// TODO: Remove any chiral label on patom?
			}
			patom->updatePropertyCache();

			// Add hydrogen to weakest ionized
			Atom* iatom = omol->getAtomWithIdx(ioccur.back());
			iatom->setFormalCharge( iatom->getFormalCharge() + 1 );
			// Increase explicit H count if no implicit, or aromatic N or P, 
			// or non default valence state
			const PeriodicTable *table = PeriodicTable::getTable();
			INT_VECT valence_list = table->getValenceList(iatom->getAtomicNum());
			bool found = (std::find(valence_list.begin(), valence_list.end(), 
															iatom->getTotalValence()) != valence_list.end());
			if ( iatom->getNoImplicit() || 
											(( patom->getAtomicNum() == 7 || patom->getAtomicNum() == 15 )
											 && patom->getIsAromatic() ) || !found ) {
				iatom->setNumExplicitHs( iatom->getNumExplicitHs() + 1 );
			}
			iatom->updatePropertyCache();
		} else {break;}
	} else {break;}
	} // while loop

	std::cout << MolToSmiles(*omol) << std::endl;
	MolOps::sanitizeMol(* (static_cast<RWMol*>(omol)) );

	return omol;
}

std::pair<unsigned int, std::vector<unsigned int>>* Reionizer::strongestProtonated(
								const ROMol &mol, 
								const std::vector<std::pair<ROMol*, ROMol*>> *abpairs) {

	// position is the position in the acid list. 
	unsigned int position = 0;
	for (auto &abpair: *abpairs) {
		RDKit::MatchVectType res;
	//	std::cout << MolToSmiles(*(abpair.first)) << std::endl;
		unsigned int matches = SubstructMatch(mol, *(abpair.first), res);
		if (matches > 0) {
			std::vector<unsigned int> occurence;
				for (const auto &pair : res) {
					occurence.push_back(pair.second);
	//				std::cout << pair.second << std::endl;
				}
			return new std::pair<unsigned int, std::vector<unsigned int>>(position, occurence);
		}
		++position;
	}
	return nullptr;
}

std::pair<unsigned int, std::vector<unsigned int>>* Reionizer::weakestIonized(
								const ROMol &mol, 
								const std::vector<std::pair<ROMol*, ROMol*>> *abpairs) {

	// position is the position in the acid list. 
	unsigned int position = 0;
	for ( auto &abpair: boost::adaptors::reverse(*abpairs) ) {
		RDKit::MatchVectType res;
	//	std::cout << MolToSmiles(*(abpair.first)) << std::endl;
		unsigned int matches = SubstructMatch(mol, *(abpair.second), res);
		if (matches > 0) {
			std::vector<unsigned int> occurence;
				for (const auto &pair : res) {
					occurence.push_back(pair.second);
	//				std::cout << pair.second << std::endl;
				}
			return new std::pair<unsigned int, std::vector<unsigned int>>(
											(abpairs->size() - position - 1), occurence);
		}
		++position;
	}
	return nullptr;
}

Uncharger::Uncharger()
	: pos_h(SmartsToMol("[+!H0!$(*~[-])]")),
		pos_quat(SmartsToMol("[+H0!$(*~[-])]")),
		neg(SmartsToMol("[-!$(*~[+H0])]")),
		neg_acid(SmartsToMol("[$([O-][C,P,S]=O),$([n-]1nnnc1),$(n1[n-]nnc1)]")) {};

Uncharger::Uncharger(const Uncharger &other) {
	pos_h = other.pos_h;
	pos_quat = other.pos_quat;
	neg = other.neg;
	neg_acid = other.neg_acid;
};

Uncharger::~Uncharger() {};

ROMol* Uncharger::uncharge(const ROMol &mol) {

	ROMol* omol = new ROMol(mol);

	std::vector<MatchVectType> p_matches;
	std::vector<MatchVectType> q_matches;
	std::vector<MatchVectType> n_matches;
	std::vector<MatchVectType> a_matches;

	// Get atom ids for matches
	unsigned int p_matched = SubstructMatch( *omol, *(this->pos_h), p_matches );
	unsigned int q_matched = SubstructMatch( *omol, *(this->pos_quat), q_matches );
	unsigned int n_matched = SubstructMatch( *omol, *(this->neg), n_matches );
	unsigned int a_matched = SubstructMatch( *omol, *(this->neg_acid), a_matches );

	// trying to understand how to use n_matches as a vector...
//	std::cout << "Size " << n_matches.size() << std::endl;
//	std::cout << n_matches[0][0].second << std::endl;
//	std::cout << n_matches[1][0].second << std::endl;
//	for (auto &i : n_matches) {
//		for (auto &j : i) {
//			std::cout << j.second << std::endl;
//		}
//	}
	//
		
	// Neutralize negative charges
	if (q_matched > 0) {
		// Surplus negative charges more than non-neutralizable positive charges
		int neg_surplus = n_matched - q_matched;
		if (a_matched > 0 && neg_surplus > 0) {
			// zwitterion with more negative charges than quaternary positive centres
			while (neg_surplus > 0 && a_matched > 0) {
				// Add hydrogen to first negative acid atom, increase formal charge
				// Until quaternary positive == negative total or no more negative acid
				Atom* atom = omol->getAtomWithIdx(a_matches[0][0].second);
				a_matches.erase(a_matches.begin());
				atom->setNumExplicitHs(atom->getNumExplicitHs() + 1);
				atom->setFormalCharge(atom->getFormalCharge() + 1);
				--neg_surplus;
				std::cout << "Removed negative charge." << std::endl;
//				std::cout << MolToSmiles(*omol) << std::endl;
			}
		}
	} else {
			std::vector<unsigned int> n_idx_matches;
			std::cout << "Negative matches " << std::endl;
			for (const auto &match : n_matches) {
				for (const auto &pair : match) {
					n_idx_matches.push_back(pair.second);
				//std::cout << pair.second << std::endl;
				}
			}
			for (const auto &idx : n_idx_matches) {
				Atom *atom = omol->getAtomWithIdx(idx);
				while (atom->getFormalCharge() < 0) {
					atom->setNumExplicitHs(atom->getNumExplicitHs() + 1);
					atom->setFormalCharge(atom->getFormalCharge() + 1);
					std::cout << "Removed negative charge." << std::endl;
//					std::cout << MolToSmiles(*omol) << std::endl;
				}
			}	
	}
	// Neutralize positive charges
	std::vector<unsigned int> p_idx_matches;
	for (const auto &match : p_matches) {
		for (const auto &pair : match) {
			p_idx_matches.push_back(pair.second);
		}
	}
	for (const auto &idx : p_idx_matches) {
		Atom *atom = omol->getAtomWithIdx(idx);
		while (atom->getFormalCharge() > 0 && atom->getNumExplicitHs() > 0) {
			atom->setNumExplicitHs(atom->getNumExplicitHs() - 1);
			atom->setFormalCharge(atom->getFormalCharge() - 1);
			std::cout << "Removed positive charge." << std::endl;
		}
	}
	return omol;
}

} // namesepace MolStandardize
} // namespace RDKit
