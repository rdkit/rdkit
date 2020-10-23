/*
 * Logger.h
 *
 *  Created on: Apr 4, 2013
 *      Author: gjones
 */

#ifndef LOGGER_H_
#define LOGGER_H_

#include <stdlib.h>
#include <iostream>
#include <string>
#include "Util.h"

namespace GarethUtil {

using namespace std;

/*
 * Our logger class for generic reporting.
 *
 * See http://www.drdobbs.com/cpp/logging-in-c/201804215 for ideas
 */
class Reporter {
public:
	enum ReportingLevel {
		TRACE, DEBUG, DETAIL, NORMAL, INFO, WARN, FATAL
	};
private:
	static ostream * fileStream;
	static ReportingLevel minReportingLevel;
	ReportingLevel reportingLevel;
	Reporter(const Reporter &);
	Reporter & operator=(const Reporter &);
	std::ostringstream os;
public:
	Reporter() : reportingLevel(NORMAL){};
	~Reporter();
	std::ostringstream& get(int level);

	static void output(const std::string & message);

	static ostream & getFileStream() {
		return *fileStream;
	}


	static void setFileStream(ostream & fileStream_) {
		fileStream = &fileStream_;
	}

	static  const ReportingLevel  getMinReportingLevel() {
		return minReportingLevel;
	}

	static void setMinReportingLevel(const ReportingLevel reportingLevel_) {
		minReportingLevel = reportingLevel_;
	}

	std::ostringstream & get(ReportingLevel level = NORMAL);

	const static std::string levelToString(ReportingLevel level) {
		static const std::string labels[] = { "TRACE", "DEBUG", "DETAIL",
				"NORMAL", "INFO", "WARNING", "FATAL" };
		return labels[level];
	}
};

}

// turn off logging for TRACE and DEBUG completely for optimized code,
// and TRACE for debug code
// (or whatever MIN_REPORTING_LEVEL is set to)

#ifndef MIN_REPORTING_LEVEL
#ifdef NDEBUG
#define MIN_REPORTING_LEVEL GarethUtil::Reporter::ReportingLevel::DETAIL
#else
#define MIN_REPORTING_LEVEL GarethUtil::Reporter::ReportingLevel::DEBUG
#endif
#endif

#define REPORT(level) \
if (level < MIN_REPORTING_LEVEL) ;\
else if (level < GarethUtil::Reporter::getMinReportingLevel()) ; \
else GarethUtil::Reporter().get(level) << __FILE__ << ":" << __LINE__ << " "


#endif /* LOGGER_H_ */
