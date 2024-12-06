// TEMPORARY PLACE HOLDER
#ifndef CDMAP_H
#define CDMAP_H
template<class Key, class Value>
class CDMap {
public:
  std::map<Key,Value> storage;
  CDMap<Key,Value>() : storage() {}
  bool Contains(const Key &key) const {
    return storage.find(key) != storage.end();
  }
  Value & operator [](const Key &key)       { return storage[key]; }
  //  const Value & operator [](const Key &key) const {
  //    return storage[key]; }
  void insert(const Key &key, Value &value) {
    storage[key] = value;
  }
  
};
#endif
