// TEMP HOLDER
#include "cs_assert.h"
// stupid implementation
template<class T>
struct auto_buffer{
  std::vector<T> buffer;

  auto_buffer<T>(unsigned num_elements) : buffer(num_elements) {
    ASSERT(num_elements);
  }
 
   
  size_t length() const { return buffer.size(); }
        T* ptr() { return &buffer[0]; }
  const T* ptr() const { return &buffer[0]; }

        T  operator[](unsigned i) { return buffer[i]; }
  const T operator[](unsigned i) const { return buffer[i]; }
  
};
