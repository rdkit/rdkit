// TEMPORARY PLACEHOLDER
#ifndef cs_STRINGUTILS_H
#define cs_STRINGUTILS_H

#include <cctype>    // std::tolower
#include <algorithm> // std::equal

namespace cs {
inline bool ichar_equals(char a, char b)
{
    return std::tolower(static_cast<unsigned char>(a)) ==
           std::tolower(static_cast<unsigned char>(b));
}

inline bool str_compareNoCase(const std::string& a, const std::string& b)
{
    return a.size() == b.size() &&
           std::equal(a.begin(), a.end(), b.begin(), ichar_equals);
}
 
}

#endif
