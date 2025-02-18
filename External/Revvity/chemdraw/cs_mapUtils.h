// TEMPORARY PLACE HOLDER

namespace cs {
  template<typename Key, typename Value>
  class csInitMap {
  public:
    std::map<Key,Value> storage;
    typedef std::pair<Key,Value> pair;
    csInitMap() : storage() {}

    typename std::map<Key,Value>::iterator find(const Key &key) {
      return storage.find(key);
    }
    Value & operator[](const Key&key) {
      return storage[key];
    }
    typename std::map<Key,Value>::iterator end() { return storage.end(); }
    
    csInitMap<Key,Value>& operator<<(const pair &) {
      return *this;
    }
  //private:
  //  friend csInitMap<Key,Value> & operator<<(csInitMap<Key,Value> &map, const csInitMap<Key,Value>::pair& p);
  };

//template<typename Key, typename Value>
//csInitMap<Key,Value>& operator<<(csInitMap<Key,Value>& os, const std::pair<Key, Value>& p) {
//    return os;
//}
}
