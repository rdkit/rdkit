#define BOOST_TEST_MODULE Test_Coordgen

#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>
#include <unordered_set>

#include "../CoordgenFragmenter.h"
#include "../sketcherMinimizer.h"
#include "../sketcherMinimizerMaths.h"
#include "../sketcherMinimizerStretchInteraction.h"
#include "../sketcherMinimizerBendInteraction.h"
#include "../sketcherMaeReading.h"
#include "coordgenBasicSMILES.h"

#include "maeparser/MaeConstants.hpp"
#include "maeparser/Reader.hpp"

using std::unordered_set;
using namespace schrodinger;

const boost::filesystem::path test_samples_path(TEST_SAMPLES_PATH);

namespace
{
std::map<sketcherMinimizerAtom*, int>
getReportingIndices(sketcherMinimizerMolecule& mol)
{
    std::map<sketcherMinimizerAtom*, int> fakeIndices;
    int index = 0;
    for (auto& atom : mol.getAtoms()) {
        fakeIndices.emplace(atom, ++index);
    }
    return fakeIndices;
}

bool areBondsNearIdeal(sketcherMinimizerMolecule& mol,
                       std::map<sketcherMinimizerAtom*, int>& indices,
                       std::set<std::pair <int, int> > skip = std::set<std::pair <int, int> > ())
{
    const float targetBondLength = BONDLENGTH * BONDLENGTH;
    const auto tolerance = static_cast<float>(targetBondLength * 0.1);

    bool passed = true;
    for (auto& bond : mol.getBonds()) {
        auto bondPair = std::pair<int, int> (indices[bond->getStartAtom()], indices[bond->getEndAtom()]);
        if (skip.find(bondPair) != skip.end()) {
            continue;
        }
        auto& startCoordinates = bond->getStartAtom()->getCoordinates();
        auto& endCoordinates = bond->getEndAtom()->getCoordinates();

        const auto sq_distance = sketcherMinimizerMaths::squaredDistance(
                                                                         startCoordinates, endCoordinates);
        const auto deviation = sq_distance - targetBondLength;
        if (deviation < -tolerance || deviation > tolerance) {
            std::cerr << "Bond" << indices[bond->getStartAtom()] << '-'
            << indices[bond->getEndAtom()] << " has length "
            << sq_distance << " (" << targetBondLength << ")\n";
            passed = false;
        }
    }
    return passed;
}

bool noCrossingBonds(sketcherMinimizerMolecule& mol,
                         std::map<sketcherMinimizerAtom*, int>& indices)
{
    bool passed = true;
    for (auto& bond : mol.getBonds()) {
        for (auto& bond2 : mol.getBonds()) {
            if (bond == bond2) continue;
            if (bond->getStartAtom() == bond2->getStartAtom()) continue;
            if (bond->getStartAtom() == bond2->getEndAtom()) continue;
            if (bond->getEndAtom() == bond2->getStartAtom()) continue;
            if (bond->getEndAtom() == bond2->getEndAtom()) continue;

            auto& startCoordinates1 = bond->getStartAtom()->getCoordinates();
            auto& endCoordinates1 = bond->getEndAtom()->getCoordinates();
            auto& startCoordinates2 = bond2->getStartAtom()->getCoordinates();
            auto& endCoordinates2 = bond2->getEndAtom()->getCoordinates();

            if (sketcherMinimizerMaths::intersectionOfSegments(startCoordinates1,
                                                               endCoordinates1,
                                                               startCoordinates2,
                                                               endCoordinates2)) {
                std::cerr << "Bond" << indices[bond->getStartAtom()] << '-'
                << indices[bond->getEndAtom()] << " intersects bond "
                << indices[bond2->getStartAtom()]<< '-'
                << indices[bond2->getEndAtom()]<<")\n";
                passed = false;
            }
        }
    }
    return passed;
}

} // namespace


static sketcherMinimizerMolecule* operator"" _smiles(const char * smiles, size_t len)
{
    return approxSmilesParse({smiles, len});
}

BOOST_AUTO_TEST_CASE(SampleTest)
{
    // A small sample test showing how to import a molecule from a .mae file and
    // initialize the minimizer with it.

    const std::string testfile = (test_samples_path / "test.mae").string();

    mae::Reader r(testfile);
    auto block = r.next(mae::CT_BLOCK);
    BOOST_REQUIRE(block != nullptr);

    auto* mol = mol_from_mae_block(*block);
    BOOST_REQUIRE(mol != nullptr);
    BOOST_REQUIRE_EQUAL(mol->getAtoms().size(), 26);
    BOOST_REQUIRE_EQUAL(mol->getBonds().size(), 26);

    sketcherMinimizer minimizer;
    minimizer.initialize(mol); // minimizer takes ownership of mol
    minimizer.runGenerateCoordinates();

    for (auto& atom : mol->getAtoms()) {
        auto c = atom->getCoordinates();

        // This is best we can do with checking: there are precision and
        // rounding issues depending on platform and environment.
        BOOST_CHECK(c.x() != 0 || c.y() != 0);
    }

    auto indices = getReportingIndices(*mol);
    BOOST_CHECK(areBondsNearIdeal(*mol, indices));
}


BOOST_AUTO_TEST_CASE(TemplateTest)
{
    ///
    // Do the structures in the templates file get the same coordinates that
    // were supplied in the templates file?

    const boost::filesystem::path source_dir(SOURCE_DIR);
    const std::string templates_file = (source_dir / "templates.mae").string();

    // Known issues. See issue #52
    const unordered_set<size_t> no_match = {1, 8, 19, 20, 22, 32, 43, 53, 65, 66, 67};
    // 32 is odd. minimization removes atoms?? But it matches??
    const unordered_set<size_t> match_incorrectly = {18, 27};

    mae::Reader r(templates_file);
    std::shared_ptr<mae::Block> b;
    size_t template_index = 0;
    while ((b = r.next(mae::CT_BLOCK)) != nullptr) {
        if (no_match.count(template_index) > 0) {
            ++template_index;
            continue;
        }

        auto* mol = mol_from_mae_block(*b);
        BOOST_REQUIRE(mol != nullptr);
        const auto original_atom_count = mol->getAtoms().size();

        sketcherMinimizer minimizer;
        minimizer.setTemplateFileDir(source_dir.string());

        minimizer.initialize(mol); // minimizer takes ownership of mol
        minimizer.runGenerateCoordinates();

        BOOST_CHECK_EQUAL(original_atom_count, mol->getAtoms().size());

        bool any_rigid = false;
        bool all_rigid = true;
        for (auto a: mol->getAtoms()) {
            if (a->rigid) {
                any_rigid = true;
            } else {
                all_rigid = false;
            }
        }
        const bool matches_incorrectly = match_incorrectly.count(template_index) > 0;
        if (!any_rigid) {
            BOOST_CHECK_MESSAGE(any_rigid, "No template found for " << template_index);
        } else if (!matches_incorrectly) {
            BOOST_CHECK_MESSAGE(all_rigid, "Not all atoms templated for " << template_index);
        }

        ++template_index;
    }
}

BOOST_AUTO_TEST_CASE(ClearWedgesTest)
{
    // test that when writing stereochemistry we first clear the existing one

    const std::string testfile = (test_samples_path / "testChirality.mae").string();

    mae::Reader r(testfile);
    auto block = r.next(mae::CT_BLOCK);
    BOOST_REQUIRE(block != nullptr);

    auto* mol = mol_from_mae_block(*block);
    BOOST_REQUIRE(mol != nullptr);

     /*set wedges on all bonds*/
    mol->getBonds().at(0)->hasStereochemistryDisplay = true;
    mol->getBonds().at(1)->hasStereochemistryDisplay = true;
    mol->getBonds().at(2)->hasStereochemistryDisplay = true;
    mol->getBonds().at(3)->hasStereochemistryDisplay = true;

    /*set chirality on the atom*/
    auto carbon = mol->getAtoms().at(0);
    BOOST_REQUIRE_EQUAL(carbon->atomicNumber, 6);
    carbon->hasStereochemistrySet = true;
    carbon->isR = true;

    /*run minimization*/
    sketcherMinimizer minimizer;
    minimizer.initialize(mol); // minimizer takes ownership of mol
    minimizer.runGenerateCoordinates();

    /*make sure that previous wedges are reset and only 2 are now used*/
    BOOST_REQUIRE_EQUAL(mol->getBonds().at(0)->hasStereochemistryDisplay, false);
    BOOST_REQUIRE_EQUAL(mol->getBonds().at(1)->hasStereochemistryDisplay, true);
    BOOST_REQUIRE_EQUAL(mol->getBonds().at(2)->hasStereochemistryDisplay, false);
    BOOST_REQUIRE_EQUAL(mol->getBonds().at(3)->hasStereochemistryDisplay, true);
}


BOOST_AUTO_TEST_CASE(DisableMetalZOBs)
{
    // Load a molecule with a single bond to a metal. Make sure that disabling the automatic conversion
    // to zob behaves as expected

    const std::string testfile = (test_samples_path / "nonterminalMetalZobs.mae").string();

    mae::Reader r(testfile);
    auto block = r.next(mae::CT_BLOCK);
    BOOST_REQUIRE(block != nullptr);

    auto* mol = mol_from_mae_block(*block);
    BOOST_REQUIRE(mol != nullptr);

    sketcherMinimizer minimizer;
    minimizer.setTreatNonterminalBondsToMetalAsZOBs(false);
    auto Al = mol->getAtoms().at(0);
    auto N = mol->getAtoms().at(1);
    //make sure we got the right atoms
    BOOST_REQUIRE_EQUAL(Al->atomicNumber, 13);
    BOOST_REQUIRE_EQUAL(N->atomicNumber, 7);

    minimizer.initialize(mol); // minimizer takes ownership of mol
    minimizer.runGenerateCoordinates();
    auto bondLength = (Al->coordinates - N->coordinates).length();
    auto expectedLength = 50.f;
    auto tolerance = 2.f;
    BOOST_REQUIRE ((bondLength > expectedLength-tolerance) && (bondLength < expectedLength+tolerance));
    auto indices = getReportingIndices(*mol);
    BOOST_CHECK(areBondsNearIdeal(*mol, indices));
}


BOOST_AUTO_TEST_CASE(terminalMetalZOBs)
{
    // Load a molecule with a single bond to a terminal metal. Make sure the bondlength is consistant with a terminal bond

    const std::string testfile = (test_samples_path / "metalZobs.mae").string();

    mae::Reader r(testfile);
    auto block = r.next(mae::CT_BLOCK);
    BOOST_REQUIRE(block != nullptr);

    auto* mol = mol_from_mae_block(*block);
    BOOST_REQUIRE(mol != nullptr);

    sketcherMinimizer minimizer;
    auto Al = mol->getAtoms().at(0);
    auto N = mol->getAtoms().at(1);
    //make sure we got the right atoms
    BOOST_REQUIRE_EQUAL(Al->atomicNumber, 13);
    BOOST_REQUIRE_EQUAL(N->atomicNumber, 7);

    minimizer.initialize(mol); // minimizer takes ownership of mol
    minimizer.runGenerateCoordinates();
    auto bondLength = (Al->coordinates - N->coordinates).length();
    auto expectedLength = 50.f;
    auto tolerance = 2.f;
    BOOST_REQUIRE ((bondLength > expectedLength-tolerance) && (bondLength < expectedLength+tolerance));
    auto indices = getReportingIndices(*mol);
    BOOST_CHECK(areBondsNearIdeal(*mol, indices));
}

BOOST_AUTO_TEST_CASE(testMinimizedRingsShape)
{
    // minimize a complex macrocycle. Make sure rings have the shape of regular polygons

    const std::string testfile = (test_samples_path / "macrocycle.mae").string();

    mae::Reader r(testfile);
    auto block = r.next(mae::CT_BLOCK);
    BOOST_REQUIRE(block != nullptr);

    auto* mol = mol_from_mae_block(*block);
    BOOST_REQUIRE(mol != nullptr);

    sketcherMinimizer minimizer;


    minimizer.initialize(mol); // minimizer takes ownership of mol
    minimizer.runGenerateCoordinates();
    //check the length of every non macrocycle-ring bond
    int bondsN = 0;
    for (auto interaction : minimizer.getStretchInteractions()) {
        auto ring = sketcherMinimizer::sameRing(interaction->atom1, interaction->atom2);
        if (ring && !ring->isMacrocycle()) {
            auto expectedLength = 50.f;
            auto tolerance = 2.f;
            auto bondLength = (interaction->atom1->coordinates - interaction->atom2->coordinates).length();
            BOOST_REQUIRE ((bondLength > expectedLength-tolerance) && (bondLength < expectedLength+tolerance));
            bondsN++;
        }
    }
//    check the angles
    int anglesN = 0;
    for (auto interaction : minimizer.getBendInteractions()) {
        if (interaction->isRing) {
            auto ring = sketcherMinimizer::sameRing(interaction->atom1, interaction->atom2, interaction->atom3);
            BOOST_REQUIRE (ring != nullptr);
            BOOST_REQUIRE (!ring->isMacrocycle());
            auto expectedValue = interaction->restV;
            auto tolerance = 2.f;
            auto angle = interaction->angle();
            BOOST_REQUIRE ((angle > expectedValue-tolerance) && (angle < expectedValue+tolerance));
            anglesN++;
        }
    }
    BOOST_REQUIRE (anglesN  == 32);
}


BOOST_AUTO_TEST_CASE(testPolyominoCoordinatesOfSubstituent)
{
    Polyomino p;
    p.addHex(hexCoords(0, 0));
    auto substCoords =
        p.coordinatesOfSubstituent(vertexCoords(1, 0, 0));
    BOOST_REQUIRE(substCoords == vertexCoords(1, -1, -1));

    p.addHex(hexCoords(1, 0));
    substCoords = p.coordinatesOfSubstituent(vertexCoords(1, 0, 0));
    BOOST_REQUIRE(substCoords == vertexCoords(0, 0, -1));
}

BOOST_AUTO_TEST_CASE(testPolyominoSameAs)
{
    // identity
    Polyomino p1;
    p1.addHex(hexCoords(0, 0));
    p1.addHex(hexCoords(1, 0));
    p1.addHex(hexCoords(2, 0));
    p1.addHex(hexCoords(0, 1));
    BOOST_REQUIRE(p1.isTheSameAs(p1));

    // order-independence
    Polyomino p2;
    p2.addHex(hexCoords(2, 0));
    p2.addHex(hexCoords(0, 0));
    p2.addHex(hexCoords(0, 1));
    p2.addHex(hexCoords(1, 0));
    BOOST_REQUIRE(p1.isTheSameAs(p2));
    BOOST_REQUIRE(p2.isTheSameAs(p1));

    // translation
    Polyomino p3;
    p3.addHex(hexCoords(4, 2));
    p3.addHex(hexCoords(5, 2));
    p3.addHex(hexCoords(6, 2));
    p3.addHex(hexCoords(4, 3));
    BOOST_REQUIRE(p1.isTheSameAs(p3));
    BOOST_REQUIRE(p3.isTheSameAs(p1));

    // rotation
    Polyomino p4;
    p4.addHex(hexCoords(0, 0));
    p4.addHex(hexCoords(-1, 0));
    p4.addHex(hexCoords(-2, 0));
    p4.addHex(hexCoords(0, -1));
    BOOST_REQUIRE(p1.isTheSameAs(p4));
    BOOST_REQUIRE(p4.isTheSameAs(p1));

    // symmetry
    Polyomino p5;
    p5.addHex(hexCoords(0, 0));
    p5.addHex(hexCoords(0, 1));
    p5.addHex(hexCoords(0, 2));
    p5.addHex(hexCoords(1, 0));
    BOOST_REQUIRE(!p1.isTheSameAs(p5));
    BOOST_REQUIRE(!p5.isTheSameAs(p1));

    // different number of points
    Polyomino p6;
    p6.addHex(hexCoords(1, 0));
    p6.addHex(hexCoords(2, 0));
    p6.addHex(hexCoords(0, 1));
    BOOST_REQUIRE(!p1.isTheSameAs(p6));
    BOOST_REQUIRE(!p6.isTheSameAs(p1));
}

BOOST_AUTO_TEST_CASE(testGetDoubleBondConstraints)
{
    /*
     test that getDoubleBondConstraints behaves as expected:
     notably don't set up interactions for double bonds that are also part of
     a small ring.
     CRDGEN-160
     */
    sketcherMinimizer min;
    CoordgenFragmentBuilder fragmentBuilder;
    CoordgenMacrocycleBuilder macrocycleBuilder;

    const std::string testfile = (test_samples_path / "test_mol.mae").string();
    mae::Reader r(testfile);
    auto block = r.next(mae::CT_BLOCK);
    BOOST_REQUIRE(block != nullptr);

    auto* mol = mol_from_mae_block(*block);
    BOOST_REQUIRE(mol != nullptr);

    min.initialize(mol); // minimizer takes ownership of mol
    for (auto molecule : min.getMolecules()) {
        for (auto ring : molecule->getRings()) {
            std::vector<sketcherMinimizerAtom*> atoms =
                fragmentBuilder.orderRingAtoms(ring);
            std::vector<doubleBondConstraint> constraints =
                macrocycleBuilder.getDoubleBondConstraints(atoms);
            BOOST_REQUIRE(constraints.empty());
        }
    }
}

BOOST_AUTO_TEST_CASE(testClockwiseOrderedSubstituents)
{
    auto mol = "CN(C)C"_smiles;

    sketcherMinimizer minimizer;
    minimizer.initialize(mol); // minimizer takes ownership of mol
    minimizer.runGenerateCoordinates();

    const auto& atoms = minimizer.getAtoms();
    sketcherMinimizerAtom* center = atoms.at(0);
    sketcherMinimizerAtom* neigh1 = atoms.at(1);
    sketcherMinimizerAtom* neigh2 = atoms.at(2);
    sketcherMinimizerAtom* neigh3 = atoms.at(3);
    BOOST_REQUIRE_EQUAL(center->getAtomicNumber(), 7);

    sketcherMinimizerPointF origin(0, 0);
    sketcherMinimizerPointF above(0, 50);
    sketcherMinimizerPointF left(-50, 0);
    sketcherMinimizerPointF right(50, 0);

    center->coordinates = origin;
    neigh1->coordinates = above;
    neigh2->coordinates = left;
    neigh3->coordinates = right;

    auto orderedNeighbors =
        center->clockwiseOrderedNeighbors();

    BOOST_REQUIRE((orderedNeighbors[0]->coordinates - above).length() == 0);
    BOOST_REQUIRE((orderedNeighbors[1]->coordinates - left).length() == 0);
    BOOST_REQUIRE((orderedNeighbors[2]->coordinates - right).length() == 0);

    BOOST_REQUIRE_EQUAL(orderedNeighbors[0], neigh1);
    BOOST_REQUIRE_EQUAL(orderedNeighbors[1], neigh2);
    BOOST_REQUIRE_EQUAL(orderedNeighbors[2], neigh3);

    neigh1->coordinates = above;
    neigh3->coordinates = left;
    neigh2->coordinates = right;

    orderedNeighbors = center->clockwiseOrderedNeighbors();

    BOOST_REQUIRE((orderedNeighbors[0]->coordinates - above).length() == 0);
    BOOST_REQUIRE((orderedNeighbors[1]->coordinates - left).length() == 0);
    BOOST_REQUIRE((orderedNeighbors[2]->coordinates - right).length() == 0);

    BOOST_REQUIRE_EQUAL(orderedNeighbors[0], neigh1);
    BOOST_REQUIRE_EQUAL(orderedNeighbors[1], neigh3);
    BOOST_REQUIRE_EQUAL(orderedNeighbors[2], neigh2);
}

BOOST_AUTO_TEST_CASE(testClockwiseOrderedNaN)
{
    std::unique_ptr<sketcherMinimizerMolecule> mol("CN(C)C"_smiles);
    auto& atoms = mol->getAtoms();
    sketcherMinimizerAtom* center = atoms.at(0);
    sketcherMinimizerAtom* neigh1 = atoms.at(1);
    neigh1->coordinates = sketcherMinimizerPointF(std::nanf("name"), std::nanf("name"));
    const auto orderedNeighbors = center->clockwiseOrderedNeighbors();

    // We usually allow the sketcher minimizer to clean up bonds & atoms,
    // but we don't have a minimizer here
    for (auto& bond : mol->getBonds()) {
        delete bond;
    }
    for (auto& atom : atoms) {
        delete atom;
    }
}


BOOST_AUTO_TEST_CASE(testbicyclopentane)
{
    /*
     test that bicyclo(1,1,1)pentane is rendered without clashes between the bridge carbons
     CRDGEN-270
     */

    auto mol = "C1C2CC1C2"_smiles;
    sketcherMinimizer minimizer;
    minimizer.initialize(mol); // minimizer takes ownership of mol
    minimizer.runGenerateCoordinates();

    const auto& atoms = minimizer.getAtoms();
    auto bridgeAtom1 = atoms.at(1);
    auto bridgeAtom2 = atoms.at(2);
    auto bridgeAtom3 = atoms.at(3);
    BOOST_TEST(bridgeAtom1->neighbors.size() == 2);
    BOOST_TEST(bridgeAtom2->neighbors.size() == 2);
    BOOST_TEST(bridgeAtom3->neighbors.size() == 2);
    auto distance1 = (bridgeAtom1->getCoordinates() - bridgeAtom2->getCoordinates()).length();
    auto distance2 = (bridgeAtom1->getCoordinates() - bridgeAtom3->getCoordinates()).length();
    auto distance3 = (bridgeAtom2->getCoordinates() - bridgeAtom3->getCoordinates()).length();
    auto minimumDistance = 15.f;
    BOOST_TEST(distance1 > minimumDistance);
    BOOST_TEST(distance2 > minimumDistance);
    BOOST_TEST(distance3 > minimumDistance);
}

BOOST_AUTO_TEST_CASE(testFusedRings)
{
    /*
     CRDGEN271, CRDGEN-272
     */

    std::vector<std::string> smiles {"C1CCC23CCCCC2CC3C1",
        "C1=CC2C3CC4C(CC3NC3CCCC(C32)N1)NC1CCCC2C1C4CCN2"};
    for (auto smile : smiles) {
        auto mol = approxSmilesParse(smile);
        sketcherMinimizer minimizer;
        minimizer.initialize(mol); // minimizer takes ownership of mol
        minimizer.runGenerateCoordinates();
        auto indices = getReportingIndices(*mol);
        BOOST_CHECK(areBondsNearIdeal(*mol, indices));
        BOOST_CHECK(noCrossingBonds(*mol, indices));
    }
}


BOOST_AUTO_TEST_CASE(testTemplates)
{
    auto mol = "C12CC3CC(CC3C1)C2"_smiles;
    sketcherMinimizer minimizer;
    minimizer.initialize(mol); // minimizer takes ownership of mol
    minimizer.runGenerateCoordinates();
    auto indices = getReportingIndices(*mol);
    //template has two stretched bonds
    std::set<std::pair<int, int> > skip;
    skip.insert(std::pair<int, int> (2, 5));
    skip.insert(std::pair<int, int> (6, 2));
    BOOST_CHECK(areBondsNearIdeal(*mol, indices, skip));
}


BOOST_AUTO_TEST_CASE(testRingComplex)
{
    auto mol = "CC1CC2CCCC(C3CC4CCC(C3)C4C)C2O1"_smiles;
    sketcherMinimizer minimizer;
    minimizer.initialize(mol); // minimizer takes ownership of mol
    minimizer.runGenerateCoordinates();
    auto indices = getReportingIndices(*mol);
    BOOST_CHECK(noCrossingBonds(*mol, indices));
}

BOOST_AUTO_TEST_CASE(testCoordgenFragmenter)
{
    /*
     Test that a molecule is fragmented as expected and fragments are given the
     correct flags. Atoms 3-8 and 11 are constrained.
                            6
                          /   \
     1 -- 2 -- 3 -- 4 -- 5     7 -- 8 -- 9 -- 10
                          \   /
                            11
    */
    auto mol = "CCCCC1CC(CCC)C1"_smiles;
    auto atoms = mol->getAtoms();
    for (int i = 3; i <= 8; ++i) {
        atoms[i-1]->constrained = true;
    }
    atoms[10]->constrained = true;
    auto atom_map = getReportingIndices(*mol);

    sketcherMinimizer minimizer;
    minimizer.initialize(mol);
    CoordgenFragmenter::splitIntoFragments(mol);

    std::vector<std::set<int>> expected_fragments = {{1, 2}, {3}, {4}, {5, 6, 7, 11}, {8}, {9, 10}};
    std::vector<std::set<int>> actual_fragments;
    for (auto fragment : mol->_fragments) {
        std::set<int> fragment_idices;
        for (auto at : fragment->getAtoms()) {
            fragment_idices.insert(atom_map[at]);
        }
        actual_fragments.push_back(fragment_idices);
    }
    std::sort(actual_fragments.begin(), actual_fragments.end());

    // Check that the fragmenting is correct
    BOOST_REQUIRE_EQUAL(actual_fragments.size(), expected_fragments.size());
    for (size_t i = 0; i < actual_fragments.size(); ++i) {
        BOOST_CHECK_EQUAL_COLLECTIONS(expected_fragments[i].begin(), expected_fragments[i].end(), actual_fragments[i].begin(), actual_fragments[i].end());
    }

    // Fragment containing atoms (1, 2)
    BOOST_TEST(atoms[0]->fragment->constrained == false);
    BOOST_TEST(atoms[0]->fragment->constrainedFlip == false);

    // Fragment containing atom 3. Flip should not be constrained since
    // the fragment does not have any constrained child fragments
    BOOST_TEST(atoms[2]->fragment->constrained == true);
    BOOST_TEST(atoms[2]->fragment->constrainedFlip == false);

    // Fragment containing atom 4. Flip should be constrained since fragment
    // has a child fragment that is constrained
    BOOST_TEST(atoms[3]->fragment->constrained == true);
    BOOST_TEST(atoms[3]->fragment->constrainedFlip == true);

    // Fragment containing atoms (5, 6, 7, 11)
    BOOST_TEST(atoms[4]->fragment->constrained == true);
    BOOST_TEST(atoms[4]->fragment->constrainedFlip == true);

    // Fragment containing atom 8
    BOOST_TEST(atoms[7]->fragment->constrained == true);
    BOOST_TEST(atoms[7]->fragment->constrainedFlip == false);

    // Fragment containing atoms (9, 10)
    BOOST_TEST(atoms[8]->fragment->constrained == false);
    BOOST_TEST(atoms[8]->fragment->constrainedFlip == false);

    for (auto& fragment : mol->_fragments) {
        delete fragment;
    }
}
