#define BOOST_TEST_MODULE Test_smilesparser

#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include "../sketcherMinimizer.h"
#include "coordgenBasicSMILES.h"

using namespace schrodinger;


static std::unique_ptr<sketcherMinimizerMolecule> operator"" _smiles(const char * smiles, size_t len)
{
    return std::unique_ptr<sketcherMinimizerMolecule>(approxSmilesParse({smiles, len}));
}

BOOST_AUTO_TEST_CASE(Basics)
{
    auto mol = "CCCC"_smiles;
    BOOST_TEST(mol->getAtoms().size() == 4);

    mol = "CNO"_smiles;
    auto& atoms = mol->getAtoms();
    auto c = atoms[0];
    auto n = atoms[1];
    auto o = atoms[2];
    BOOST_TEST(c->getAtomicNumber() == 6);
    BOOST_TEST(n->getAtomicNumber() == 7);
    BOOST_TEST(o->getAtomicNumber() == 8);
    BOOST_TEST(c->isNeighborOf(n));
    BOOST_TEST(n->isNeighborOf(o));
    BOOST_TEST(!c->isNeighborOf(o));
}

BOOST_AUTO_TEST_CASE(Rings)
{
    // 4 membered ring with a tail
    auto mol = "C1CCC1C"_smiles;
    BOOST_TEST(mol->getAtoms().size() == 5);
    BOOST_TEST(mol->getBonds().size() == 5);
    BOOST_TEST(mol->getAtoms()[0]->isNeighborOf(mol->getAtoms()[3]));
}

BOOST_AUTO_TEST_CASE(Branching)
{
    auto mol = "CC(C)(C)C"_smiles;
    BOOST_TEST(mol->getAtoms()[1]->getBonds().size() == 4);

    mol = "CC(C)(CC)C"_smiles;
    BOOST_TEST(mol->getAtoms()[3]->getBonds().size() == 2);

}

BOOST_AUTO_TEST_CASE(BondOrder)
{
    auto mol = "C=C"_smiles;
    BOOST_TEST(mol->getBonds()[0]->getBondOrder() == 2);
}