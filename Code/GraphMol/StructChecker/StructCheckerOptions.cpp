//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <iostream>
#include <sstream>
#include <boost/property_tree/json_parser.hpp>
#include "StructChecker.h"

namespace RDKit {
 namespace StructureCheck {

     bool parseOptionsJSON(const std::string &json, StructCheckerOptions &op) {
        if (json.empty())
            return false;
        try {
            std::istringstream ss;
            ss.str(json);
            std::istream& iss = ss;
            boost::property_tree::ptree pt;
            boost::property_tree::read_json(ss, pt);
            op = StructCheckerOptions(); //reset to default values
            op.AcidityLimit = pt.get<double>("AcidityLimit", op.AcidityLimit);
            op.RemoveMinorFragments = pt.get<bool>("RemoveMinorFragments", op.RemoveMinorFragments);
            op.DesiredCharge = pt.get<int >("DesiredCharge", op.DesiredCharge);
            op.CheckCollisions = pt.get<bool>("CheckCollisions", op.CheckCollisions);
            op.CollisionLimitPercent = pt.get<int >("CollisionLimitPercent", op.CollisionLimitPercent);
            op.MaxMolSize = pt.get<int >("MaxMolSize", op.MaxMolSize);
            op.ConvertSText = pt.get<bool>("ConvertSText", op.ConvertSText);
            op.SqueezeIdentifiers = pt.get<bool>("SqueezeIdentifiers", op.SqueezeIdentifiers);
            op.StripZeros = pt.get<bool>("StripZeros", op.StripZeros);
            op.CheckStereo = pt.get<bool>("CheckStereo", op.CheckStereo);
            op.ConvertAtomTexts = pt.get<bool>("ConvertAtomTexts", op.ConvertAtomTexts);
            op.GroupsToSGroups = pt.get<bool>("GroupsToSGroups", op.GroupsToSGroups);
            op.Verbose = pt.get<bool>("Verbose", op.Verbose);
        }
        catch (boost::property_tree::json_parser_error &ex) {
            std::cerr<<"JSON:"<< ex.message() <<"\n";
            return false;
        }
        catch (std::exception &ex) {
            std::cerr<<"JSON:"<< ex.what() << "\n";
            return false;
        }
        catch (...) {
            std::cerr << "JSON: Unknown exception.\n";
            return false;
        }
        return true;
    }

    bool loadOptionsFromFiles(StructCheckerOptions &op,
        const std::string &augmentedAtomTranslationsFile,
        const std::string &patternFile,       // file with clean patterns
        const std::string &rotatePatternFile, // file with rotate patterns
        const std::string &stereoPatternFile, // file with stereo patterns
        const std::string &tautomerFile)    {

//TODO: ...

        return true;
    }

 }// namespace StructureCheck
} // namespace RDKit
