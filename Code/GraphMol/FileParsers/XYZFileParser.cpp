#include <boost/lexical_cast.hpp>

#include "FileParsers.h"
#include "FileParserUtils.h"
#include <RDGeneral/StreamOps.h>

#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <exception>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>
#include <cstdlib>
#include <cstdio>

namespace RDKit {

void ParseExtraLine(std::string_view extraLine) {
    std::string_view whitespace{" \v\t"};
    if (extraLine.find_first_not_of(whitespace) != std::string_view::npos) {
        std::ostringstream errout;
        errout << "More lines than expected" << std::endl;
        throw FileParseException(errout.str());
    }
}

Atom *ParseXYZFileAtomLine(std::string_view atomLine, RDGeom::Point3D &pos, unsigned int line) {
    // find text
    std::string_view whitespace{" \v\t"};
    size_t delim0 = atomLine.find_first_not_of(whitespace);
    size_t delim1 = atomLine.find_first_of(whitespace, delim0);
    if (delim1 == std::string_view::npos) {
        std::ostringstream errout;
        errout << "Missing coordinates on line " << line << std::endl;
        throw FileParseException(errout.str());
    }
    size_t delim2 = atomLine.find_first_not_of(whitespace, delim1);
    if (delim2 == std::string_view::npos) {
        std::ostringstream errout;
        errout << "Missing coordinates on line " << line << std::endl;
        throw FileParseException(errout.str());
    }
    size_t delim3 = atomLine.find_first_of(whitespace, delim2);
    if (delim3 == std::string_view::npos) {
        std::ostringstream errout;
        errout << "Missing coordinates on line " << line << std::endl;
        throw FileParseException(errout.str());
    }
    size_t delim4 = atomLine.find_first_not_of(whitespace, delim3);
    if (delim4 == std::string_view::npos) {
        std::ostringstream errout;
        errout << "Missing coordinates on line " << line << std::endl;
        throw FileParseException(errout.str());
    }
    size_t delim5 = atomLine.find_first_of(whitespace, delim4);
    if (delim5 == std::string_view::npos) {
        std::ostringstream errout;
        errout << "Missing coordinates on line " << line << std::endl;
        throw FileParseException(errout.str());
    }
    size_t delim6 = atomLine.find_first_not_of(whitespace, delim5);
    if (delim6 == std::string_view::npos) {
        std::ostringstream errout;
        errout << "Missing coordinates on line " << line << std::endl;
        throw FileParseException(errout.str());
    }
    size_t delim7 = atomLine.find_last_not_of(whitespace) + 1;
    
    // set conformer
    try {
        pos.x = FileParserUtils::toDouble(atomLine.substr(delim2, delim3 - delim2), false);
    } catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert '" << atomLine.substr(delim2, delim3 - delim2) << "' to double on line " << line << std::endl;
        throw FileParseException(errout.str());
    }

    try {
        pos.y = FileParserUtils::toDouble(atomLine.substr(delim4, delim5 - delim4), false);
    } catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert '" << atomLine.substr(delim4, delim5 - delim4) << "' to double on line " << line << std::endl;
        throw FileParseException(errout.str());
    }
    
    try {
        pos.z = FileParserUtils::toDouble(atomLine.substr(delim6, delim7 - delim6), false);
    } catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert '" << atomLine.substr(delim6, delim7 - delim6) << "' to double on line " << line << std::endl;
        throw FileParseException(errout.str());
    }
    
    // create atom
    Atom *atom = new Atom;
    
    std::string symb{atomLine.substr(delim0, delim1 - delim0)};
    if (symb.size() == 2 && symb[1] >= 'A' && symb[1] <= 'Z') {
        symb[1] = static_cast<char>(tolower(symb[1]));
    }
    
    try {
      atom->setAtomicNum(PeriodicTable::getTable()->getAtomicNumber(symb));
    } catch (const Invar::Invariant &e) {
      delete atom;
      throw FileParseException(e.what());
    }
    
    return atom;
}


RWMol *XYZDataStreamToMol(std::istream &inStream) {
    PRECONDITION(inStream, "no stream");
    unsigned int numAtoms = 0;
    
    std::string_view num{getLine(inStream)};
    try {
        numAtoms = FileParserUtils::toUnsigned(num);
    } catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Unable to recognize the number of atoms: cannot convert '" << num << "' to unsigned int on line 0" << std::endl;
        throw FileParseException(errout.str());
    }
    
    std::string comment{getLine(inStream)};
    
    RWMol *mol = nullptr;
    if (numAtoms != 0) {
        mol = new RWMol();
        Conformer *conf = new Conformer(numAtoms);
        if (comment != "") {
            mol->setProp("_FileComments", comment);
        }
        for (unsigned int i = 0; i < numAtoms; i++) {
            if (inStream.eof()) {
                throw FileParseException("EOF hit while reading atoms");
            }
            RDGeom::Point3D pos;
            std::string_view atomLine{getLine(inStream)};
            Atom *atom = ParseXYZFileAtomLine(atomLine, pos, i + 2);
            unsigned int idx = mol->addAtom(atom, false, true);
            conf->setAtomPos(idx, pos);
            mol->setAtomBookmark(atom, idx);
        }
        mol->addConformer(conf);
    }
    
    while (!inStream.eof()) {
        std::string_view extraLine{getLine(inStream)};
        ParseExtraLine(extraLine);
    }
    
    return mol;
}

RWMol *XYZFileToMol(const std::string &fName, int charge=0) {
    
    std::ifstream xyzFile(fName);
    if (!xyzFile || (xyzFile.bad())) {
      std::ostringstream errout;
      errout << "Bad input file " << fName;
      throw BadFileException(errout.str());
    }
    
    RWMol *mol = nullptr;
    if (!xyzFile.eof()) {
        mol = XYZDataStreamToMol(xyzFile);
    }

    return mol;
}

} // namespace RDKit

