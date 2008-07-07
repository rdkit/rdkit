#include <map>
#include <iostream>
#include <string>


#ifndef CLASSA_H
#define CLASSA_H

#ifdef WIN32
#pragma warning (disable: 4786) // warning: long & complicated stl warning
#pragma warning (disable: 4788) // warning: long & complicated stl warning
#pragma warning (disable: 4660)
#pragma warning (disable: 4275) // warning: non dll-interface class used as...
#pragma warning (disable: 4305) // warning: truncation from 'const double' to 'const float' 
#endif

class classA {
 public:
  classA() {
    setProp("useless", 10);
  };
  ~classA() {}
  
  void printA() const {
    if (hasProp("useless")) {
      std::cout << "has useless\n";
    }
  }

  bool hasProp(std::string key) const {
    return (dp_props.find(key) != dp_props.end());
  }

  void setProp(std::string key, int val) {
    dp_props[key] = val;
  }

 private:
  std::map<std::string, int> dp_props;
};

#endif
