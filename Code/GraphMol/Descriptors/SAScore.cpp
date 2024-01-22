#include "SAScore.h"

namespace RDKit::Descriptors::SAScore {

double FeatureLibrary::GetFeatureScore(std::uint32_t feature) const {
	FeatureScores::const_iterator it = feature_scores.find(feature);
	if (it != feature_scores.end()) {
		return it->second;
	};
	return -4.0;
};

std::uint64_t FeatureLibrary::FeatureCountPercentile(double percentile) const {
	std::uint64_t total_count = 0;
	std::vector<std::uint32_t> features;
	features.reserve(feature_counts.size());
	for (const auto [feature, count] : feature_counts) {
		features.push_back(feature);
		total_count += count;
	};
	std::sort(features.begin(), features.end(),
		[&](std::uint32_t feature1, std::uint32_t feature2) {
			return feature_counts.at(feature1) > feature_counts.at(feature2);
		});
	double target_count = total_count * percentile / 100.0;
	std::uint64_t cumulative_count = 0, pct = 0;
	for (std::uint32_t feature : features) {
		std::uint64_t feature_count = feature_counts.at(feature);
		cumulative_count += feature_count;
		if (cumulative_count >= target_count) {
			pct = feature_count;
			break;
		};
	};
	return pct;
};

void FeatureLibrary::AddMolecule(const RDKit::ROMol& molecule) {
	const RDKit::SparseIntVect<std::uint32_t>* fingerprint =
		RDKit::MorganFingerprints::getFingerprint(molecule, feature_radius);
	for (const auto [feature, count] : fingerprint->getNonzeroElements()) {
		auto [it, emplaced] = feature_counts.try_emplace(feature, count);
		if (!emplaced) {
			it->second += count;
		};
	};
	delete fingerprint;
};

void FeatureLibrary::CalcFeatureScores() {
	feature_scores.clear();
	double percentile = FeatureCountPercentile();
	for (const auto [feature, count] : feature_counts) {
		double score = log10(count / percentile);
		feature_scores.emplace(feature, score);
	};
};

double FeatureLibrary::FeatureScore(
	const RDKit::SparseIntVect<std::uint32_t>& fingerprint,
	unsigned n_atoms) const {
	double score = 0.0;
	double n_features = 0.0;
	double n_unique_features = 0.0;
	for (const auto [feature, count] : fingerprint.getNonzeroElements()) {
		score += GetFeatureScore(feature) * count;
		n_features += count;
		++n_unique_features;
	};
	if (n_features == 0.0) {
		n_features = 1.0;
		n_unique_features = 1.0;
	};
	score /= n_features;
	// Calculate a correction term for the fingerprint density that deems
	// symmetrical molecules easier to synthesize and viceversa.
	// NOTE: This term isn't part of the original publication, but was added by
	// the paper's author to the Python RDKit implementation.
	double symmetry_correction = 0.0;
	if (n_atoms > n_unique_features) {
		symmetry_correction = log(n_atoms / n_unique_features) * 0.5;
	};
	return score + symmetry_correction;
};

double FeatureLibrary::FeatureScore(const RDKit::ROMol& molecule) const {
	const RDKit::SparseIntVect<std::uint32_t>* fingerprint =
		RDKit::MorganFingerprints::getFingerprint(molecule, feature_radius);
	double feature_score = FeatureScore(*fingerprint, molecule.getNumAtoms());
	delete fingerprint;
	return feature_score;
};

Components::Components(
	const FeatureLibrary& feature_library,
	const RDKit::ROMol& molecule) {
	unsigned n_atoms = molecule.getNumAtoms();
	size_penalty = std::pow(n_atoms, 1.005) - n_atoms;
	stereo_penalty = log10(CalcNumChiralCenters(molecule) + 1);
	spiro_penalty = log10(RDKit::Descriptors::calcNumSpiroAtoms(molecule) + 1);
	bridgehead_penalty = log10(
		RDKit::Descriptors::calcNumBridgeheadAtoms(molecule) + 1);
	macrocycle_penalty = 0.0;
	// This differs from the paper, which defines the macrocycle penalty to be:
	// log10(n_macrocycles + 1)
	// This form generates better results when 2 or more macrocycles are present.
	if (HasMacrocycle(molecule)) {
		macrocycle_penalty = 0.30103; // log10(2)
	};
	complexity_score = 0.0 - size_penalty - stereo_penalty - spiro_penalty -
		bridgehead_penalty - macrocycle_penalty;
	feature_score = feature_library.FeatureScore(molecule);
	raw_sascore = complexity_score + feature_score;
	sascore = ScaleSAScore(raw_sascore);
};

double Components::ScaleSAScore(double raw_sascore) const {
	double scaled_sascore = 11.0 - (raw_sascore + 5.0) / 6.5 * 9.0;
	if (scaled_sascore > 8.0) {
		scaled_sascore = 8.0 + log(scaled_sascore - 8.0);
	};
	if (scaled_sascore > 10.0) {
		scaled_sascore = 10.0;
	} else if (scaled_sascore < 1.0) {
		scaled_sascore = 1.0;
	};
	return scaled_sascore;
};

unsigned CalcNumChiralCenters(const RDKit::ROMol& molecule) {
	unsigned n_chiral_centers = 0;
	for (const RDKit::Atom* atom : molecule.atoms()) {
		if (atom->getChiralTag() != RDKit::Atom::CHI_UNSPECIFIED) {
			++n_chiral_centers;
		};
	};
	return n_chiral_centers;
};

bool HasMacrocycle(const RDKit::ROMol& molecule) {
	const RDKit::RingInfo* ring_info = molecule.getRingInfo();
	const std::vector<std::vector<int>>& rings = ring_info->atomRings();
	for (const auto& ring : rings) {
		if (ring.size() > 8) {
			return true;
		};
	};
	return false;
};

}; // ! namespace RDKit::Descriptors::SAScore

namespace RDKit::Descriptors {

double calcSAScore(
	const SAScore::FeatureLibrary& feature_library,
	const RDKit::ROMol& molecule) {
	return SAScore::Components(feature_library, molecule).sascore;
};

}; // ! namespace RDKit::Descriptors
