//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "Reporter.h"
#include "Util.h"

namespace GarethUtil {

using namespace std;

ostream * Reporter::fileStream = &cout;
Reporter::ReportingLevel Reporter::minReportingLevel = NORMAL;

std::ostringstream& Reporter::get(ReportingLevel level) {
	os << currentTime();
	os << " " << levelToString(level) << ": ";
	//os << std::string(level > DETAIL ? 0 : DETAIL-level, '\t');
	reportingLevel = level;
	return os;
}

Reporter::~Reporter() {
	if (reportingLevel >= minReportingLevel) {
		os << std::endl;
		*fileStream << os.str();
		fileStream->flush();
	}
}

}
