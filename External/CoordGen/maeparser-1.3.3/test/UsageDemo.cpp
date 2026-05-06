///
// Demonstrate reading a schrodinger/Maestro formatted file and gleaning the
// bonding information, coordinates, and atomic number.
//
// Note that Maestro format uses 1-based indices when referring to indexed
// data in blocks. (for example, the atom numbers in bond blocks are 1 indexed.)
//
// Maestro "structures" may contain multiple non-bonded molecules in a
// coherent environment. For instance, both a ligand and a receptor may exist
// in a single f_m_ct block.
//
#include <array>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <vector>

#include "MaeConstants.hpp"
#include "Reader.hpp"

#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

using namespace schrodinger::mae;

// These classes are not intended for production use. The are merely intended
// to illustrate where data is stored in the "block" data structures.
class Bond
{
  public:
    Bond(int atom0, int atom1, int bond_order)
        : atom0(atom0), atom1(atom1), bond_order(bond_order)
    {
    }

    const int atom0;
    const int atom1;
    const int bond_order;
};

class Structure
{
  public:
    std::string title;
    std::vector<int> atomic_numbers;
    std::vector<std::array<double, 3>> coordinates;

    std::vector<Bond> bonds;

    // A "property" that some atoms have (others may not have this property)
    std::unordered_map<size_t, int> demo_property;
};

const boost::filesystem::path test_samples_path(TEST_SAMPLES_PATH);
const std::string compressed_sample =
    (test_samples_path / "test2.maegz").string();

BOOST_AUTO_TEST_SUITE(DemoSuite)

// Reads all atom and bond information from test.mae, which is a standard
// mae formatted file. Only accesses properties that are gauranteed to
// exist in every f_m_ct block.
BOOST_AUTO_TEST_CASE(maeBlock)
{
    schrodinger::mae::Reader r(compressed_sample);

    std::vector<std::shared_ptr<Structure>> structures;
    std::shared_ptr<Block> b;
    while ((b = r.next(CT_BLOCK)) != nullptr) {
        auto st = std::make_shared<Structure>();
        st->title = b->getStringProperty(CT_TITLE);

        // Atom data is in the m_atom indexed block
        {
            const auto atom_data = b->getIndexedBlock(ATOM_BLOCK);
            // All atoms are gauranteed to have these three field names:
            const auto atomic_numbers =
                atom_data->getIntProperty(ATOM_ATOMIC_NUM);
            const auto xs = atom_data->getRealProperty(ATOM_X_COORD);
            const auto ys = atom_data->getRealProperty(ATOM_Y_COORD);
            const auto zs = atom_data->getRealProperty(ATOM_Z_COORD);
            const auto size = atomic_numbers->size();
            BOOST_REQUIRE_EQUAL(size, xs->size());
            BOOST_REQUIRE_EQUAL(size, ys->size());
            BOOST_REQUIRE_EQUAL(size, zs->size());

            // atomic numbers, and x, y, and z coordinates
            for (size_t i = 0; i < size; ++i) {
                st->atomic_numbers.push_back(atomic_numbers->at(i));
                st->coordinates.push_back({{xs->at(i), ys->at(i), zs->at(i)}});
            }

            // Other properties could fail, because not all properties have
            // values for all atoms. The last atom of the first structure does
            // not have the "i_m_template_index" property.
            const auto template_indices =
                atom_data->getIntProperty("i_m_template_index");
            for (size_t i = 0; i < size; ++i) {
                if (template_indices->isDefined(i)) {
                    st->demo_property[i] = template_indices->at(i);
                }
            }
        }

        // Bond data is in the m_bond indexed block
        {
            const auto bond_data = b->getIndexedBlock(BOND_BLOCK);
            // All bonds are gauranteed to have these three field names:
            auto bond_atom_1s = bond_data->getIntProperty(BOND_ATOM_1);
            auto bond_atom_2s = bond_data->getIntProperty(BOND_ATOM_2);
            auto orders = bond_data->getIntProperty(BOND_ORDER);
            const auto size = bond_atom_1s->size();

            for (size_t i = 0; i < size; ++i) {
                // Atom indices in the bond data structure are 1 indexed!
                const auto bond_atom_1 = bond_atom_1s->at(i) - 1;
                const auto bond_atom_2 = bond_atom_2s->at(i) - 1;
                const auto order = orders->at(i);

                // Only one direction of the bond is recorded in the file
                st->bonds.emplace_back(bond_atom_1, bond_atom_2, order);
                st->bonds.emplace_back(bond_atom_2, bond_atom_1, order);
            }
        }

        structures.push_back(st);
    }

    // Check that all three f_m_ct blocks were read.
    BOOST_CHECK_EQUAL(structures.size(), 3u);
}

BOOST_AUTO_TEST_SUITE_END()
