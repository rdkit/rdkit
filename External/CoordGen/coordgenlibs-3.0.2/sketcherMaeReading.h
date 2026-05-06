#pragma once

///
// Shim for creating sketcherMolecules from .mae files.

#include "sketcherMinimizerMolecule.h"
#include "maeparser/MaeConstants.hpp"
#include "maeparser/Reader.hpp"

///
// A very simple utility function to parse a mae::Block into a 2D
// sketcherMinimizerMolecule. Anything beyond atomic number, x and y coordinates
// and bond orders will be ignored (i.e. no chiralities or stereo bonds will be
// parsed).
//
sketcherMinimizerMolecule* mol_from_mae_block(schrodinger::mae::Block& block)
{
    auto molecule = new sketcherMinimizerMolecule();
    // Atom data is in the m_atom indexed block
    {
        const auto atom_data =
            block.getIndexedBlock(schrodinger::mae::ATOM_BLOCK);
        // All atoms are gauranteed to have these three field names:
        const auto atomic_numbers =
            atom_data->getIntProperty(schrodinger::mae::ATOM_ATOMIC_NUM);
        const auto xs =
            atom_data->getRealProperty(schrodinger::mae::ATOM_X_COORD);
        const auto ys =
            atom_data->getRealProperty(schrodinger::mae::ATOM_Y_COORD);
        const auto size = atomic_numbers->size();

        // atomic numbers, and x, y, and z coordinates
        for (size_t i = 0; i < size; ++i) {
            auto atom = molecule->addNewAtom();
            atom->setAtomicNumber(atomic_numbers->at(i));
            atom->setCoordinates(sketcherMinimizerPointF(
                static_cast<float>(xs->at(i)), static_cast<float>(ys->at(i))));
        }
    }

    // Bond data is in the m_bond indexed block
    try {
        const auto bond_data =
            block.getIndexedBlock(schrodinger::mae::BOND_BLOCK);
        // All bonds are gauranteed to have these three field names:
        auto from_atoms =
            bond_data->getIntProperty(schrodinger::mae::BOND_ATOM_1);
        auto to_atoms =
            bond_data->getIntProperty(schrodinger::mae::BOND_ATOM_2);
        auto orders = bond_data->getIntProperty(schrodinger::mae::BOND_ORDER);
        const auto size = from_atoms->size();

        for (size_t i = 0; i < size; ++i) {
            // Maestro atoms are 1 indexed!
            auto* from_atom = molecule->getAtoms().at(from_atoms->at(i) - 1);
            auto* to_atom = molecule->getAtoms().at(to_atoms->at(i) - 1);
            auto bond = molecule->addNewBond(from_atom, to_atom);
            bond->setBondOrder(orders->at(i));
        }
    } catch (const std::out_of_range&) {
        // no bonds.
    }

    return molecule;
}
