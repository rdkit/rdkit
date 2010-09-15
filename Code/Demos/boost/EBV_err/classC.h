#include <map>
#include <iostream>
#include <string>
#include <vector>

#ifndef CLASSC_H
#define CLASSC_H

#ifdef WIN32
#pragma warning (disable: 4786) // warning: long & complicated stl warning
#pragma warning (disable: 4788) // warning: long & complicated stl warning
#pragma warning (disable: 4660)
#pragma warning (disable: 4275) // warning: non dll-interface class used as...
#pragma warning (disable: 4305) // warning: truncation from 'const double' to 'const float' 
#endif

typedef std::pair<std::string, int> STR_INT;
typedef std::vector<STR_INT> PAIR_VECT;

class classC {
  
 public:
  classC() {
    setProp("useless", 10);
  };
  ~classC() {}
  
  void printC() const {
    if (hasProp("useless")) {
      std::cout << "has useless\n";
    }
  }

  bool hasProp(std::string key) const {
    PAIR_VECT::const_iterator pvi;
    for (pvi = dp_props.begin(); pvi != dp_props.end(); pvi++) {
      if (pvi->first == key) {
        return true;
      }
    }
    return false; 
  }

  void setProp(std::string key, int val) {
    std::pair<std::string, int> newp(key, val);
    dp_props.push_back(newp);
  }

 private:
  PAIR_VECT dp_props;
  
};

#endif
