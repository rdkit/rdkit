//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once
#ifndef _RD_STRUCTCHECKER_H__
#define _RD_STRUCTCHECKER_H__

#include <string>
#include <vector>
#include "../RDKitBase.h"

/* Example of Usage
1)  StructChecker chk;
    int flags = StructureCheck::checkMolStructure( mol ); // use defaults
 or
2)
StructureCheck::StructCheckerOptions options;   // use defaults
    // To use external data
    StructureCheck::loadOptionsFromFiles(options, file1, file2, …);
    StructChecker chk(options);

    for( mol in mols ) {
        int flags = StructureCheck::checkMolStructure( mol, &options);
        if (0!=(flags & StructureCheck::StructureFlags::BAD_SET)) {
        // write to error file
        } else if (0!=(flags & StructureCheck::StructureFlags::TRANSFORMED_SET)) {
        // input molecule was transformed
        } else { // flag == NO_CHANGE 
        // no change
        }
    }
*/

namespace RDKit {
 namespace StructureCheck {

    // Flags for the return values of the StructureChecker
    // TypeDefs for translating augmented atom pairs

    struct Ligand {
        std::string AtomSymbol;
        int         Charge;
        int         Radical;
        int         SubstitutionCount; // substitution count 0 = don't care
        Bond::BondType BondType; // was: short 
    };

    struct AugmentedAtom {
        std::string AtomSymbol;
        std::string ShortName;
        int         Charge;
        int         Radical;
        int         Topography; // 1=ring 2=chain 0 don't care
        std::vector<Ligand> Ligands;
    };

    ////////////////////////////////////////////////////////////////////////////
    // Structure Check Options
    // Holds all the user options for the StructureChecking.
    // Can be initialized from factory functions, perhaps serialized

    struct StructCheckerOptions {
        double  AcidityLimit;
        bool    RemoveMinorFragments;
        int     DesiredCharge;
        bool    CheckCollisions;
        int     CollisionLimitPercent;
        int     MaxMolSize;
        bool    ConvertSText;
        bool    SqueezeIdentifiers;
        bool    StripZeros;
        bool    CheckStereo;
        bool    ConvertAtomTexts;
        bool    GroupsToSGroups;
        bool    Verbose;
        // Internal data for struchk
        std::vector<std::pair<AugmentedAtom, AugmentedAtom> > AugmentedAtomPairs;
        std::vector<AugmentedAtom> AcidicAtoms;
        std::vector<AugmentedAtom> GoodAtoms;
        std::vector<ROMOL_SPTR> Patterns;
        std::vector<ROMOL_SPTR> RotatePatterns;
        std::vector<ROMOL_SPTR> StereoPatterns;
        std::vector<ROMOL_SPTR> FromTautomer;
        std::vector<ROMOL_SPTR> ToTautomer;
    public:
        StructCheckerOptions() : AcidityLimit(0.0)
            , RemoveMinorFragments(false)
            , DesiredCharge(0)
            , CheckCollisions(false)
            , CollisionLimitPercent(0)
            , MaxMolSize(255)           //// Really NEED it ???
            , ConvertSText(false)
            , SqueezeIdentifiers(false)
            , StripZeros(false)
            , CheckStereo(false)
            , ConvertAtomTexts(false)
            , GroupsToSGroups(false)
            , Verbose(false)
        {}

        void clear() { *this = StructCheckerOptions(); }
    };

    bool parseOptionsJSON(const std::string &json, StructCheckerOptions &op);

    bool loadOptionsFromFiles(StructCheckerOptions &op,
        const std::string &augmentedAtomTranslationsFile = "",
        const std::string &patternFile       = "",       // file with clean patterns
        const std::string &rotatePatternFile = "", // file with rotate patterns
        const std::string &stereoPatternFile = "", // file with stereo patterns
        const std::string &tautomerFile      = "");

    class StructChecker
    {
    public:
        typedef enum StructureFlags {
            NO_CHANGE = 0,
            BAD_MOLECULE            = 0x0001,
            ALIAS_CONVERSION_FAILED = 0x0002,
            STEREO_ERROR            = 0x0004,
            STEREO_FORCED_BAD       = 0x0008,
            ATOM_CLASH              = 0x0010,
            ATOM_CHECK_FAILED       = 0x0020,
            SIZE_CHECK_FAILED       = 0x0040,
            // reserved error = 0x0080,
            TRANSFORMED             = 0x0100,
            FRAGMENTS_FOUND         = 0x0200,
            EITHER_WARNING          = 0x0400,
            DUBIOUS_STEREO_REMOVED  = 0x0800,
            RECHARGED               = 0x1000,
            STEREO_TRANSFORMED      = 0x2000,
            TEMPLATE_TRANSFORMED    = 0x4000,
            TAUTOMER_TRANSFORMED    = 0x8000,
            // mask:
            BAD_SET = (BAD_MOLECULE |
                ALIAS_CONVERSION_FAILED |
                STEREO_ERROR |
                STEREO_FORCED_BAD |
                ATOM_CLASH |
                ATOM_CHECK_FAILED |
                SIZE_CHECK_FAILED),

            TRANSFORMED_SET = (TRANSFORMED |
                FRAGMENTS_FOUND |
                EITHER_WARNING |
                DUBIOUS_STEREO_REMOVED |
                STEREO_TRANSFORMED |
                TEMPLATE_TRANSFORMED |
                TAUTOMER_TRANSFORMED |
                RECHARGED),
        } StructureFlags;
    // attributes:
    private:
        const StructCheckerOptions Options;
    public:
        inline StructChecker() {}
        inline StructChecker(const StructCheckerOptions &options) : Options(options) {}

        // Check structure and return a set of StructureFlags
        unsigned checkMolStructure(RWMol &mol) const;

        // an instance independed helper methods:
        // Converts structure property flags to a comma seperated string
        static std::string StructureFlagsToString(unsigned flags);
        // Converts a comma seperated string to a StructureFlag unsigned integer
        static unsigned StringToStructureFlags(const std::string& str);
    // internal implementation:
    private:
    };
 }
}
#endif
