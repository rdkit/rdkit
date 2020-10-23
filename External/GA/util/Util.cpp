/*
 * Util.cpp
 *
 *  Created on: Apr 5, 2013
 *      Author: gjones
 */

#include <unistd.h>
#include <sstream>
#include <string>
#include <cstdlib>
#include <glob.h>

#include <boost/algorithm/string.hpp>

#include "Util.h"

namespace GarethUtil {

using namespace std;

string currentTime() {
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer, 80, "%m-%d-%Y %H:%M:%S", timeinfo);

	return buffer;
}

bool startsWith(string str, string prefix) {
	if (prefix.length() > str.length())
		return false;
	return str.compare(0, prefix.length(), prefix) == 0;
}

string getUserName() {
	const int bufsize = 100;
	char buffer[bufsize];

	if (!getlogin_r(buffer, bufsize))
		return string(buffer);
	else
		return string("");
}

string & removeTrailingLF(string & line) {
	if (!line.empty() && line[line.length() - 1] == '\r')
		line.erase(line.length() - 1);
	return line;
}

string & trim(string & str) {
	boost::trim(str);
	return str;
}

string & toUpperCase(string & str) {
	boost::to_upper(str);
	return str;
}

string & toLowerCase(string & str) {
	boost::to_lower(str);
	return str;
}

bool equals(const string & str1, const string & str2) {
	return str1.compare(str2) == 0;
}

bool equalsIgnoreCase(const string & str1, const string & str2) {
	return boost::iequals(str1, str2);
}

bool endsWith(const string & str, const string & suffix) {
	if (suffix.length() > str.length())
		return false;
	return str.compare(str.length() - suffix.length(), string::npos, suffix)
			== 0;
}

bool equals(const double d1, const double d2, const int ulp) {
	// the machine epsilon has to be scaled to the magnitude of the values used
	// and multiplied by the desired precision in ULPs (units in the last place)
	return abs(d1 - d2)
			< numeric_limits<double>::epsilon() * std::abs(d1 + d1) * ulp
	// unless the result is subnormal
			|| abs(d1 - d2) < numeric_limits<double>::min();
}

bool equals(const double d1, const double d2, const double epsilon) {
	return equals(d1, d2, 1) || abs(d1 - d2) < epsilon;
}

bool equals(const double d1, const double d2) {
	return equals(d1, d2, 1);
}

boost::optional<string> getEnv(const string &name) {
	const char *value = getenv(name.c_str());
	if (value == nullptr) {
		return boost::none;
	}
	return string(value);
}

std::vector<string> glob(const string& pattern) {
	glob_t glob_result;
	glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
	std::vector<string> ret;
	for (unsigned int i = 0; i < glob_result.gl_pathc; ++i) {
		ret.push_back(string(glob_result.gl_pathv[i]));
	}
	globfree(&glob_result);
	return ret;
}

void printMatrix(const Eigen::MatrixXd & m) {
	cout << "Matrix " << endl << m;
}

void printVector(const Eigen::VectorXd & v) {
	cout << "Vector " << endl << v;

}

}
