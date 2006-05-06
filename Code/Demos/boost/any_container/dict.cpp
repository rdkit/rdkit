
#include <map>
#include <boost/any.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/smart_ptr.hpp>
#include <string>
#include <iostream>
#include <Invariant/Invariant.h>

int idx=0;

class Contained {
public:
  Contained() : _idx(idx) {
    std::cout << "Contained CTOR:" << _idx << std::endl;
    idx++;
  }
  ~Contained() {
    std::cout << "Contained DTOR:" << _idx << std::endl;
  }
  Contained( const Contained &other) : _idx(other._idx) {
    std::cout << "Contained Copy CTOR:" << _idx << std::endl;
  }
  Contained &operator=( const Contained &other) {
    _idx = other._idx;
    std::cout << "Contained Assign:" << _idx << std::endl;
    return *this;
  }
  int _idx;
};


  template <typename T>
  struct larger_of {
    T operator()(T arg1,T arg2) { return arg1>arg2 ? arg1 : arg2; };
  };


  struct charptr_functor {
    bool operator()(const char* s1, const char* s2) const
    {
      return strcmp(s1, s2) < 0;
    };
  };
  

class Dict {
  public:
    typedef std::map<const char *,boost::any,charptr_functor> DataType;
    Dict() {_data.clear();};
    Dict(const Dict &other) { _data = other._data; };
    Dict &operator=(const Dict &other) {
      _data = other._data;
      return *this;
    };
    
    bool hasVal(const char *what) const {
      //DataType::const_iterator i;
      //std::cerr << "\t\tpre: ";
      //for(i=_data.begin();i!=_data.end();i++){
      //std::cerr << (*i).first << " ";
      //}
      //std::cerr << std::endl;
      //std::cerr << "\thasVal: " << what << ": " << _data.count(what) << std::endl;
      //std::cerr << "\t\tpost: ";
      //for(i=_data.begin();i!=_data.end();i++){
      //std::cerr << (*i).first << " ";
      //}
      //std::cerr << std::endl;
      return _data.count(what);
    };
    bool hasVal(const std::string &what) const {
      return hasVal(what.c_str());
    };
    
    //
    //  We're going to try and be somewhat crafty about this getVal stuff to make these
    //  containers a little bit more generic.  The normal behavior here is that the
    //  value being queried must be directly castable to type T.  We'll robustify a
    //  little bit by trying that and, if the cast fails, attempting a couple of 
    //  other casts, which will then be lexically cast to type T.
    //
    template <typename T>
    void getVal(const char *what,T &res) const {
      PRECONDITION(hasVal(what),"getVal called on non-existent key");
      const boost::any &val = _data.find(what)->second;
      res = boost::any_cast<T>(val);
    };

    void getVal(const char *what,std::string &res) const {
      PRECONDITION(hasVal(what),"getVal called on non-existent key");
      const boost::any &val = _data.find(what)->second;
      try{
	res = boost::any_cast<std::string>(val);
      } catch (const boost::bad_any_cast &) {
	if(val.type()==typeid(int)){
	  res = boost::lexical_cast<std::string>(boost::any_cast<int>(val));
	} else if(val.type()==typeid(long)){
	  res = boost::lexical_cast<std::string>(boost::any_cast<long>(val));
	} else if(val.type()==typeid(float)){
	  res = boost::lexical_cast<std::string>(boost::any_cast<float>(val));
	} else if(val.type()==typeid(double)){
	  res = boost::lexical_cast<std::string>(boost::any_cast<double>(val));
	} else if(val.type()==typeid(const char *)){
	  res = std::string(boost::any_cast<const char *>(val));
	} else {
	  throw;
	}
      }
    };


    template <typename T>
    void getVal(const std::string &what,T &res) const { getVal(what.c_str(),res); };

    template <typename T>
    void setVal(const char *what,T &val){
      _data[what] = boost::any(val);
    };
    template <typename T>
    void setVal(const std::string &what,T &val) { setVal(what.c_str(),val); };

    void clearVal(const char *what) { _data.erase(what); };

    void reset() { _data.clear(); };

  private:
    DataType _data;

  };


void test1(){
  // basic containers
  std::cout << "----------- TEST1 -----------" << std::endl;

  Contained c;
  TEST_ASSERT(c._idx==0);
  Contained d(c);
  TEST_ASSERT(d._idx==0);
  Contained e;
  TEST_ASSERT(e._idx==1);
  std::cout << "assign" << std::endl;
  e=c;
  TEST_ASSERT(e._idx==0);


  std::cout <<" Done" << std::endl;
}

void test2(){
  // shared pointers with containers
  std::cout << "----------- TEST2 -----------" << std::endl;

  typedef boost::shared_ptr<Contained> CONT_SPTR;
  CONT_SPTR s_d(new Contained());
  TEST_ASSERT(s_d.use_count()==1);
  CONT_SPTR s_e=s_d;
  TEST_ASSERT(s_e.use_count()==2);
  s_d.reset();
  TEST_ASSERT(s_e.use_count()==1);
  s_e.reset();

  std::cout <<" Done" << std::endl;
}


void test3(){
  // basic dict stuff
  std::cout << "----------- TEST3 -----------" << std::endl;

  Dict dict;
  TEST_ASSERT(!dict.hasVal("foo"));

  int tmp=4;
  dict.setVal("foo",tmp);
  TEST_ASSERT(dict.hasVal("foo"));
  tmp=0;
  dict.getVal("foo",tmp);
  TEST_ASSERT(tmp==4);

  Dict d2(dict);
  TEST_ASSERT(d2.hasVal("foo"));
  d2.getVal("foo",tmp);
  TEST_ASSERT(tmp==4);

  Dict d3;
  d3=d2;
  TEST_ASSERT(d3.hasVal("foo"));
  d3.getVal("foo",tmp);
  TEST_ASSERT(tmp==4);
  
  std::cout <<" Done" << std::endl;
}

void test4(){
  // dict containing Contained objects
  std::cout << "----------- TEST4 -----------" << std::endl;

  Dict dict;
  TEST_ASSERT(!dict.hasVal("foo"));

  Contained tmp;
  int hold=tmp._idx;
  std::cout << " set" << std::endl;
  dict.setVal("foo",tmp);
  std::cout << " query" << std::endl;
  TEST_ASSERT(dict.hasVal("foo"));
  std::cout << " get" << std::endl;
  dict.getVal("foo",tmp);
  TEST_ASSERT(tmp._idx==hold);

  std::cout << " copy" << std::endl;
  Dict d2(dict);
  TEST_ASSERT(d2.hasVal("foo"));
  d2.getVal("foo",tmp);
  TEST_ASSERT(tmp._idx==hold);


  Dict d3;
  d3=d2;
  TEST_ASSERT(d3.hasVal("foo"));
  d3.getVal("foo",tmp);
  TEST_ASSERT(tmp._idx==hold);

  
  std::cout <<" Done" << std::endl;
}

void test5(){
  // dict containing sptrs to Contained objects
  std::cout << "----------- TEST5 -----------" << std::endl;
  typedef boost::shared_ptr<Contained> CONT_SPTR;

  Dict dict;
  TEST_ASSERT(!dict.hasVal("foo"));

  CONT_SPTR tmp(new Contained());
  int hold=tmp->_idx;
  std::cout << " set" << std::endl;
  dict.setVal("foo",tmp);
  std::cout << " query" << std::endl;
  TEST_ASSERT(dict.hasVal("foo"));
  std::cout << " get" << std::endl;
  dict.getVal("foo",tmp);
  TEST_ASSERT(tmp->_idx==hold);

  std::cout << " copy" << std::endl;
  Dict d2(dict);
  TEST_ASSERT(d2.hasVal("foo"));
  d2.getVal("foo",tmp);
  TEST_ASSERT(tmp->_idx==hold);


  Dict d3;
  d3=d2;
  TEST_ASSERT(d3.hasVal("foo"));
  d3.getVal("foo",tmp);
  TEST_ASSERT(tmp->_idx==hold);

  
  std::cout <<" Done" << std::endl;
}


int main(){
  test1();
  test2();
  test3();
  test4();
  test5();
  
}
