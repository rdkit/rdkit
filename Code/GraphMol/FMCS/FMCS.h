//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once
#include <vector>
#include <string>
#include <stdexcept>
#include "../RDKitBase.h"

namespace RDKit {
    struct MCSAtomCompareParameters {
        bool    MatchValences;
    public:
        MCSAtomCompareParameters() : MatchValences(false) {}
    };

    struct MCSBondCompareParameters {
        bool    RingMatchesRingOnly;
        bool    CompleteRingsOnly;
    public:
        MCSBondCompareParameters() : RingMatchesRingOnly(false), CompleteRingsOnly(false) {}
    };

    typedef bool (*MCSAtomCompareFunction)(const MCSAtomCompareParameters& p, const ROMol& mol1, unsigned int atom1, const ROMol& mol2, unsigned int atom2, void* userData);
    typedef bool (*MCSBondCompareFunction)(const MCSBondCompareParameters& p, const ROMol& mol1, unsigned int bond1, const ROMol& mol2, unsigned int bond2, void* userData);

    // Some predefined functors:
    bool MCSAtomCompareAny      (const MCSAtomCompareParameters& p, const ROMol& mol1, unsigned int atom1, const ROMol& mol2, unsigned int atom2, void* userData);
    bool MCSAtomCompareElements (const MCSAtomCompareParameters& p, const ROMol& mol1, unsigned int atom1, const ROMol& mol2, unsigned int atom2, void* userData);
    bool MCSAtomCompareIsotopes (const MCSAtomCompareParameters& p, const ROMol& mol1, unsigned int atom1, const ROMol& mol2, unsigned int atom2, void* userData);

    bool MCSBondCompareAny      (const MCSBondCompareParameters& p, const ROMol& mol1, unsigned int bond1, const ROMol& mol2, unsigned int bond2, void* userData);
    bool MCSBondCompareOrder    (const MCSBondCompareParameters& p, const ROMol& mol1, unsigned int bond1, const ROMol& mol2, unsigned int bond2, void* userData); // ignore Aromatization
    bool MCSBondCompareOrderExact(const MCSBondCompareParameters& p, const ROMol& mol1, unsigned int bond1, const ROMol& mol2, unsigned int bond2, void* userData);

    struct MCSProgressData {
        unsigned NumAtoms;
        unsigned NumBonds;
        unsigned SeedProcessed;
    public:
        MCSProgressData() : NumAtoms(0), NumBonds(0), SeedProcessed(0) {}
    };
    struct MCSParameters;
    typedef bool (*MCSProgressCallback)(const MCSProgressData& stat, const MCSParameters &params, void* userData);
    bool MCSProgressCallbackTimeout(const MCSProgressData& stat, const MCSParameters &params, void* userData);

    struct MCSParameters {
        bool     MaximizeBonds;
        double   Threshold;
        unsigned Timeout;   // in seconds
        bool     Verbose;
        MCSAtomCompareParameters    AtomCompareParameters;
        MCSBondCompareParameters    BondCompareParameters;
        MCSAtomCompareFunction      AtomTyper;
        MCSBondCompareFunction      BondTyper;
        void*                       CompareFunctionsUserData;
        MCSProgressCallback         ProgressCallback;       // return false to interrupt execution
        void*                       ProgressCallbackUserData;
    public:
        MCSParameters(): MaximizeBonds(true)
            , Threshold(1.0)    // match to all
            , Timeout(-1)
            , Verbose(false)
            , AtomTyper(MCSAtomCompareElements)
            , BondTyper(MCSBondCompareOrder)
            , CompareFunctionsUserData(0)
            , ProgressCallback(MCSProgressCallbackTimeout)
            , ProgressCallbackUserData(0)
        {}
    };

    struct MCSResult {
        unsigned    NumAtoms;
        unsigned    NumBonds;
        std::string SmartsString;
        bool        Canceled;   // interrupted by timeout or user defined progress callback. Contains valid current MCS !
    public:
        MCSResult() : NumAtoms(0), NumBonds(0), Canceled(false) {}
        bool isCompleted()const {
            return !Canceled;
        }
    };

    MCSResult findMCS (const std::vector<ROMOL_SPTR>& mols, const MCSParameters* params=0);

} // namespace RDKit

