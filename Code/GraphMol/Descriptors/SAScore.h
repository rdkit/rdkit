#pragma once
#ifndef RDKIT_SASCORE_H
#define RDKIT_SASCORE_H

#include <RDGeneral/export.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <boost/serialization/unordered_map.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <math.h>

namespace RDKit::Descriptors::SAScore {

static const unsigned feature_library_version = 20230224;

//! Class to calculate and store SAScore feature contributions to SAScore
/*!
	For information on SAScore see the corresponding paper:
	P. Ertl and A. Schuffenhauer. Journal of Cheminformatics 1, 8 (2009)
*/
class FeatureLibrary {
	typedef std::unordered_map<std::uint32_t, std::uint64_t> FeatureCounts;
	typedef std::unordered_map<std::uint32_t, double> FeatureScores;

	unsigned feature_radius = 2;
	FeatureCounts feature_counts;
	FeatureScores feature_scores;

private:
	std::uint64_t FeatureCountPercentile(
		double percentile = 80.0) const;

	double GetFeatureScore(std::uint32_t feature) const;

	template <class Archive>
	void serialize(
		Archive& ar,
		const unsigned version = feature_library_version) {
		ar & feature_radius & feature_counts & feature_scores;
	};

	friend class boost::serialization::access;

public:
	FeatureLibrary(unsigned feature_radius = 2) :
		feature_radius(feature_radius) {};

	void AddMolecule(const RDKit::ROMol& molecule);

	void CalcFeatureScores();

	double FeatureScore(
		const RDKit::SparseIntVect<std::uint32_t>& fingerprint,
		unsigned n_atoms) const;
	double FeatureScore(const RDKit::ROMol& molecule) const;
};

//! Class to calculate and store the components used to calculate SAScore
struct Components {
	double size_penalty = 0.0;
	double stereo_penalty = 0.0;
	double spiro_penalty = 0.0;
	double bridgehead_penalty = 0.0;
	double macrocycle_penalty = 0.0;
	double complexity_score = 0.0;
	double feature_score = 0.0;
	double raw_sascore = 0.0;
	double sascore = 0.0;

	//! Calculate the individual components of a molecule's SAScore
	Components(
		const FeatureLibrary& feature_library,
		const RDKit::ROMol& molecule);

	double ScaleSAScore(double raw_sascore) const;
};

unsigned CalcNumChiralCenters(const RDKit::ROMol& molecule);
bool HasMacrocycle(const RDKit::ROMol& molecule);

}; // ! namespace RDKit::Descriptors::SAScore

namespace RDKit::Descriptors {

//! Calculate the SAScore of a molecule
/*!
	For information on SAScore see the corresponding paper:
	P. Ertl and A. Schuffenhauer. Journal of Cheminformatics 1, 8 (2009)
  \param feature_library  Library of feature contributions to SAScore
  \param molecule         The molecule of interest
*/
RDKIT_DESCRIPTORS_EXPORT
double calcSAScore(
	const SAScore::FeatureLibrary& feature_library,
	const RDKit::ROMol& molecule);

}; // ! namespace RDKit::Descriptors

#endif // ! RDKIT_SASCORE_H
