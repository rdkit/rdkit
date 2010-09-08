//
//  Copyright (C) 2003 Rational Discovery LLC
//

#include <boost/python.hpp>


class ClassA{
public:
  int get4() { return 4; };
};

class ClassB{
public:
  ClassA *returnOther() { return new ClassA; };
  int get3() { return 3; };
  int acceptOther(ClassA *other) { return other->get4();};
};
