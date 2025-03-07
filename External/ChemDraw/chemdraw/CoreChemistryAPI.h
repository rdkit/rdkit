#ifndef CORE_CHEISTRY_API_H
#define CORE_CHEISTRY_API_H
#include <RDGeneral/export.h>
#include <boost/lexical_cast.hpp>
#include <string>

#ifdef __linux
// conda-clang can set this, the linux version works just
//  fine now, so unset mac specific stuff
# ifdef TARGET_OS_MAC
#  define TARGET_OS_MAC 0
# endif
#endif

#ifdef RDKIT_CHEMDRAW_BUILD
#define CORE_CHEMISTRY_API RDKIT_EXPORT_API
#else
#define CORE_CHEMISTRY_API RDKIT_IMPORT_API
#endif

enum Dimensions2or3 {
  kIn2D = 2,
  kIn3D = 3
};

template<typename T>
void DeleteAndNull(T* ptr) {
  delete ptr;
  ptr = nullptr;
}


namespace cs {
  template<class T>
  std::string NumToStr(T value) {
    return boost::lexical_cast<std::string>(value);
  }
  
  inline size_t StrToNum(const std::string &value) {
    return boost::lexical_cast<size_t>(value);
  }

  inline double StrToDub(const std::string &value) {
    return boost::lexical_cast<double>(value);
  }
}

template<class T, class U> T FAST_dynamic_cast(U* obj) {
  return dynamic_cast<T>(obj);
}
#endif
