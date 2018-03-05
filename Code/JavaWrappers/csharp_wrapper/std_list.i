/* -----------------------------------------------------------------------------
 * See the LICENSE file for information on copyright, usage and redistribution
 * of SWIG, and the README file for authors - http://www.swig.org/release.html.
 *
 * std_list.i
 * ----------------------------------------------------------------------------- */

%include <std_common.i>

%{
#include <list>
#include <stdexcept>
%}

namespace std {
    
    template<class T> class list {
      public:
        typedef size_t size_type;
        typedef T value_type;
        typedef const value_type& const_reference;
        list();
        size_type size() const;
        %rename(isEmpty) empty;
        bool empty() const;
        void clear();
        %rename(add) push_back;
        void push_back(const value_type& x);
        %extend {
            const_reference get(int i) throw (std::out_of_range) {
                int size = int(self->size());
                int j;
                if (i>=0 && i<size) {
                  std::list< T >::const_iterator p;  
                  p=self->begin(); 
                  for (j=0; j<i; j++) {p++;}
                  return (*p);   
                }
                else
                    throw std::out_of_range("list index out of range");
            }

            bool equals(const list<T> &o){
              if(self->size()==o.size()){
                std::list< T >::const_iterator sIt=self->begin();
                std::list< T >::const_iterator oIt=o.begin();
                while(sIt != self->end()){
                  if(*sIt != *oIt) return false;
                  ++sIt;
                  ++oIt;
                }
                return true;
              } else {
                return false;
              }
            }

        }
   };
}

%define specialize_std_list(T)
#warning "specialize_std_list - specialization for type T no longer needed"
%enddef

