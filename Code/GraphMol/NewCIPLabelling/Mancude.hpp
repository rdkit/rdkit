//
//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once

#include <queue>
#include <vector>

namespace RDKit
{

namespace NewCIPLabelling
{

template <typename A, typename B> class BaseMol;

namespace Mancude
{

class Fraction
{
  private:
    template <typename A, typename B>
    friend std::vector<Fraction> calcFracAtomNums(const BaseMol<A, B>* mol);
    int num = 0;
    int den = 1;

  public:
    static int compare(int anum, int aden, int bnum, int bden)
    {
        if (anum * bden == aden * bnum) {
            return 0;
        } else {
            return static_cast<double>(anum) / aden -
                   static_cast<double>(bnum) / bden;
        }
    }

    Fraction() = default;
    Fraction(int num, int den) : num{num}, den{den} {}

    int getNum() const { return num; }

    int getDen() const { return den; }

    int compareTo(const Fraction& o)
    {
        return compare(this->num, this->den, o.num, o.den);
    }
};

enum class Type {
    Cv4D4,      // =CH-
    Nv3D2,      // =N-
    Nv4D3Plus,  // =[N+]<
    Nv2D2Minus, // -[N-]-
    Cv3D3Minus, // -[CH-]-
    Ov3D2Plus,  // -[O+]=
    Other
};

namespace
{

template <typename A, typename B>
bool SeedTypes(std::vector<Type>& types, const BaseMol<A, B>* mol)
{
    bool result = false;
    for (A atom : mol->atoms()) {
        const int aidx = mol->getAtomIdx(atom);
        types[aidx] = Type::Other;

        // check ring
        int btypes = mol->getNumHydrogens(atom);
        bool ring = false;
        for (B bond : mol->getBonds(atom)) {
            switch (mol->getBondOrder(bond)) {
            case 1:
                btypes += 0x00000001;
                break;
            case 2:
                btypes += 0x00000100;
                break;
            case 3:
                btypes += 0x00010000;
                break;
            default:
                btypes += 0x01000000;
                break;
            }
            if (mol->isInRing(bond)) {
                ring = true;
            }
        }
        if (ring) {
            int q = mol->getCharge(atom);
            switch (mol->getAtomicNum(atom)) {
            case 6:  // C
            case 14: // Si
            case 32: // Ge
                if (q == 0 && btypes == 0x0102)
                    types[aidx] = Type::Cv4D4;
                else if (q == -1 && btypes == 0x0003) {
                    types[aidx] = Type::Cv3D3Minus;
                    result = true;
                }
                break;
            case 7:  // N
            case 15: // P
            case 33: // As
                if (q == 0 && btypes == 0x0101) {
                    types[aidx] = Type::Nv3D2;
                    result = true;
                } else if (q == -1 && btypes == 0x0002) {
                    types[aidx] = Type::Nv2D2Minus;
                    result = true;
                } else if (q == +1 && btypes == 0x0102) {
                    types[aidx] = Type::Nv4D3Plus;
                    result = true;
                }
                break;
            case 8: // O
                if (q == 1 && btypes == 0x0101) {
                    types[aidx] = Type::Ov3D2Plus;
                    result = true;
                }
                break;
            }
        }
    }
    return result;
}

template <typename A, typename B>
void RelaxTypes(std::vector<Type>& types, const BaseMol<A, B>* mol)
{
    auto counts = std::vector<int>(mol->getNumAtoms());
    auto queue = std::queue<A>();

    for (A atom : mol->atoms()) {
        int aidx = mol->getAtomIdx(atom);
        for (B bond : mol->getBonds(atom)) {
            A nbr = mol->getOther(bond, atom);
            if (types[mol->getAtomIdx(nbr)] != Type::Other) {
                ++counts[aidx];
            }
        }

        if (counts[aidx] == 1) {
            queue.push(atom);
        }
    }

    while (!queue.empty()) {
        A atom = queue.front();
        queue.pop();

        int aidx = mol->getAtomIdx(atom);
        if (types[aidx] != Type::Other) {
            types[aidx] = Type::Other;

            for (B bond : mol->getBonds(atom)) {
                A nbr = mol->getOther(bond, atom);
                int nbridx = nbr->getIdx();
                --counts[nbridx];
                if (counts[nbridx] == 1) {
                    queue.push(nbr);
                }
            }
        }
    }
}

template <typename A, typename B>
void VisitPart(std::vector<int>& parts, const std::vector<Type>& types,
               int part, A atom, const BaseMol<A, B>* mol)
{
    A next;
    do {
        next = nullptr;
        for (B bond : mol->getBonds(atom)) {
            if (!mol->isInRing(bond)) {
                continue;
            }

            A nbr = mol->getOther(bond, atom);
            int aidx = mol->getAtomIdx(nbr);

            if (parts[aidx] == 0 && types[aidx] != Type::Other) {
                parts[aidx] = part;
                if (next != nullptr) {
                    VisitPart(parts, types, part, nbr, mol);
                } else {
                    next = nbr;
                }
            }
        }
        atom = next;
    } while (atom != nullptr);
}

template <typename A, typename B>
int VisitParts(std::vector<int>& parts, const std::vector<Type>& types,
               const BaseMol<A, B>* mol)
{
    int numparts = 0;
    for (A atom : mol->atoms()) {
        int aidx = mol->getAtomIdx(atom);
        if (parts[aidx] == 0 && types[aidx] != Type::Other) {
            parts[aidx] = ++numparts;
            VisitPart(parts, types, parts[aidx], atom, mol);
        }
    }
    return numparts;
}

} // namespace

template <typename A, typename B>
std::vector<Fraction> calcFracAtomNums(const BaseMol<A, B>* mol)
{
    auto num_atoms = mol->getNumAtoms();
    auto fractions = std::vector<Fraction>(num_atoms);

    for (int aidx = 0; aidx < num_atoms; ++aidx) {
        fractions[aidx] = Fraction(mol->getAtomicNum(mol->getAtom(aidx)), 1);
    }

    auto types = std::vector<Type>(num_atoms, Type::Other);
    if (SeedTypes(types, mol)) {
        RelaxTypes(types, mol);

        auto parts = std::vector<int>(num_atoms);
        int numparts = VisitParts(parts, types, mol);

        auto resparts = std::vector<int>(numparts);
        int numres = 1;

        if (numparts > 0) {
            for (int i = 0; i < num_atoms; ++i) {
                if (parts[i] == 0) {
                    continue;
                }
                A atom = mol->getAtom(i);

                if (types[i] == Type::Cv3D3Minus ||
                    types[i] == Type::Nv2D2Minus) {
                    int j = 0;
                    for (; j < numres; ++j)
                        if (resparts[j] == parts[i]) {
                            break;
                        }
                    if (j >= numres) {
                        resparts[numres] = parts[i];
                        ++numres;
                    }
                }

                fractions[i].num = 0;
                fractions[i].den = 0;

                for (B bond : mol->getBonds(atom)) {
                    A nbr = mol->getOther(bond, atom);
                    if (parts[mol->getAtomIdx(nbr)] == parts[i]) {
                        fractions[i].num += mol->getAtomicNum(nbr);
                        ++fractions[i].den;
                    }
                }
            }
        }

        if (numres > 0) {
            for (int j = 0; j < numres; ++j) {
                auto frac = Fraction(0, 0);
                int part = resparts[j];
                for (int i = 0; i < num_atoms; ++i) {
                    if (parts[i] == part) {
                        fractions[i] = frac;
                        ++frac.den;
                        A atom = mol->getAtom(i);

                        for (B bond : mol->getBonds(atom)) {
                            A nbr = mol->getOther(bond, atom);
                            int bord = mol->getBondOrder(bond);
                            if (bord > 1 &&
                                parts[mol->getAtomIdx(nbr)] == part) {
                                frac.num += (bord - 1) * mol->getAtomicNum(nbr);
                            }
                        }
                    }
                }
            }
        }
    }
    return fractions;
}

} // namespace Mancude

} // namespace NewCIPLabelling
} // namespace RDKit
