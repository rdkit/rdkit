#include "TautomerCatalogUtils.h"
#include <RDGeneral/BadFileException.h>
#include <boost/tokenizer.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
#include <fstream>
#include <string>

namespace RDKit {
namespace {
MolStandardize::TautomerTransform* getTautomer(const std::string &tmpStr) {
MolStandardize::TautomerTransform* transform = nullptr;
	if (tmpStr.length() == 0) {
    // empty line
    return transform;
  }
  if (tmpStr.substr(0, 2) == "//") {
    // comment line
    return transform;
  }
  boost::char_separator<char> tabSep("\t");
  tokenizer tokens(tmpStr, tabSep);
	std::vector<std::string> result(tokens.begin(), tokens.end());

	// tautomer information to collect from each line
	std::string name = "";
	std::string smarts = "";
	std::string bonds = ""; 
	std::string charges = ""; 

	// line must have at least two tab separated values
	if ( result.size() < 2 ) {
		std::cout << "Invalid line." << std::endl;
		return transform;
	}
	// line only has name and smarts
	if ( result.size() == 2 ) {
		name = result[0];
		smarts = result[1];
	}
	// line has name, smarts, bonds
	if ( result.size() == 3 ) {
		name = result[0];
		smarts = result[1];
		bonds = result[2];
	}
	// line has name, smarts, bonds, charges
	if ( result.size() == 4 ) {
		name = result[0];
		smarts = result[1];
		bonds = result[2];
		charges = result[3];
	}

  boost::erase_all(smarts, " ");
  boost::erase_all(name, " ");
  boost::erase_all(bonds, " ");
  boost::erase_all(charges, " ");

	ROMol* tautomer = SmartsToMol(smarts);
  CHECK_INVARIANT(tautomer, smarts);
  tautomer->setProp(common_properties::_Name, name);
	transform = new MolStandardize::TautomerTransform( tautomer, bonds, charges );
	
	return transform;
}
} // end of local utility namespace
	


namespace MolStandardize {

std::vector<TautomerTransform> readTautomers(std::string fileName) {
	std::ifstream inStream(fileName.c_str());
	if ((!inStream) || (inStream.bad())) {
		std::ostringstream errout;
		errout << "Bad input file " << fileName;
		throw BadFileException(errout.str());
	}
	std::vector<TautomerTransform> tautomers = readTautomers(inStream);
	return tautomers;
}

std::vector<TautomerTransform> readTautomers(std::istream &inStream, int nToRead) {
	std::vector<TautomerTransform> tautomers;
	tautomers.clear();
	if (inStream.bad()) {
		throw BadFileException("Bad stream contents.");
	}
	const int MAX_LINE_LEN = 512;
	char inLine[MAX_LINE_LEN];
	std::string tmpstr;
	int nRead = 0;
	while (!inStream.eof() && (nToRead < 0 || nRead < nToRead)) {
		inStream.getline(inLine, MAX_LINE_LEN, '\n');
		tmpstr = inLine;
		// parse the molpair on this line (if there is one)
		TautomerTransform* transform = getTautomer(tmpstr) ;
		if (transform != nullptr) {
//			std::cout << MolToSmiles(*(transform->Mol) ) << std::endl;
			tautomers.push_back(*transform);
			nRead++;
		}
	}
	return tautomers;
}

} // namespace MolStandardize
} // namespace RDKit
