//  (C) Copyright Greg Landrum 2003-2006.
//  Derived from vector_indexing_suite.hpp,
//  which is (C) Copyright Joel de Guzman 2003.
//  Permission to copy, use, modify, sell and distribute this software
//  is granted provided this copyright notice appears in all copies. This
//  software is provided "as is" without express or implied warranty, and
//  with no claim as to its suitability for any purpose.

#ifndef LIST_INDEXING_SUITE_GL20042_HPP
# define LIST_INDEXING_SUITE_GL20042_HPP

# include <boost/python/suite/indexing/vector_indexing_suite.hpp>
# include <boost/python/suite/indexing/indexing_suite.hpp>
# include <boost/python/suite/indexing/container_utils.hpp>
# include <boost/python/iterator.hpp>

#include "pyint_api.h"

namespace boost { namespace python {
            
    // Forward declaration
    template <class Container, bool NoProxy, class DerivedPolicies>
    class list_indexing_suite;
    
    namespace detail
    {
        template <class Container, bool NoProxy>
        class final_list_derived_policies 
            : public list_indexing_suite<Container, 
                NoProxy, final_list_derived_policies<Container, NoProxy> > {};
    }

    // The vector_indexing_suite class is a predefined indexing_suite derived 
    // class for wrapping std::vector (and std::vector like) classes. It provides
    // all the policies required by the indexing_suite (see indexing_suite).
    // Example usage:
    //
    //  class X {...};
    //
    //  ...
    //
    //      class_<std::vector<X> >("XVec")
    //          .def(vector_indexing_suite<std::vector<X> >())
    //      ;
    //
    // By default indexed elements are returned by proxy. This can be
    // disabled by supplying *true* in the NoProxy template parameter.
    //
    template <
        class Container, 
        bool NoProxy = false,
        class DerivedPolicies 
            = detail::final_list_derived_policies<Container, NoProxy> >
    class list_indexing_suite 
        : public indexing_suite<Container, DerivedPolicies, NoProxy>
    {
    public:
    
        typedef typename Container::value_type data_type;
        typedef typename Container::value_type key_type;
        typedef typename Container::size_type index_type;
        typedef typename Container::size_type size_type;
        typedef typename Container::iterator iterator_type;
        

      static 
        typename mpl::if_<
            is_class<data_type>
          , data_type&
          , data_type
        >::type
        get_item(Container& container, index_type i)
        {
	  iterator_type pos=moveToPos(container,i);
	  return *pos;
        }
        static object 
        get_slice(Container& container, index_type from, index_type to){
	  Container res;
	  iterator_type beg=moveToPos(container,from);
	  iterator_type end=moveToPos(container,to);
	  std::copy(beg,end,res.begin());
	  return object(res);
	};
        static void 
        set_item(Container& container, index_type i, data_type const& v){
	  iterator_type pos=moveToPos(container,i);
	  *pos = v;
	};

        static void 
        set_slice(Container& container, index_type from, 
		  index_type to, data_type const& v){
	  iterator_type beg=moveToPos(container,from);
	  iterator_type end=moveToPos(container,to);
	  container.erase(beg,end);
	  // FIX: didn't we just invalidate this iterator?
	  container.insert(beg,v);
	  
	};
        template <class Iter>
        static void 
        set_slice(Container& container, index_type from, 
		  index_type to, Iter first, Iter last){
	  iterator_type beg=moveToPos(container,from);
	  iterator_type end=moveToPos(container,to);
	  container.erase(beg,end);
	  // FIX: didn't we just invalidate this iterator?
	  container.insert(beg,first,last);
	};

        static void 
        delete_item(Container& container, index_type i)
        { 
	  iterator_type pos=moveToPos(container,i);
	  container.erase(pos);
        }
        
        static void 
        delete_slice(Container& container, index_type from, index_type to){
	  iterator_type beg=moveToPos(container,from);
	  iterator_type end=moveToPos(container,to);
	  container.erase(beg,end);
	}
        
        static size_t
        size(Container& container)
        {
            return container.size();
        }
        
        static bool
        contains(Container& container, key_type const& key)
        {
            return std::find(container.begin(), container.end(), key)
                != container.end();
        }
        
        static index_type
        get_min_index(Container& container)
        { 
            return 0;
        }

        static index_type
        get_max_index(Container& container)
        { 
            return container.size();
        }
      
        static bool 
        compare_index(Container& container, index_type a, index_type b)
        {
            return a < b;
        }
        
        static void 
        append(Container& container, data_type const& v)
        { 
            container.push_back(v);
        }
        template <class Iter>
        static void 
        extend(Container& container, Iter first, Iter last)
        { 
            container.insert(container.end(), first, last);
        }


        static index_type
        convert_index(Container& container, PyObject* i_)
        { 
            extract<long> i(i_);
            if (i.check())
            {
                long index = i();
                if (index < 0)
                    index += DerivedPolicies::size(container);
                if (index >= long(container.size()) || index < 0)
                {
                    PyErr_SetString(PyExc_IndexError, "Index out of range");
                    throw_error_already_set();
                }
                return index;
            }
            
            PyErr_SetString(PyExc_TypeError, "Invalid index type");
            throw_error_already_set();
            return index_type();
        }

    private:
#if 0    
        static void
        base_append(Container& container, object v)
        {
            extract<data_type&> elem(v);
            // try if elem is an exact Data
            if (elem.check())
            {
                DerivedPolicies::append(container, elem());
            }
            else
            {
                //  try to convert elem to data_type
                extract<data_type> elem(v);
                if (elem.check())
                {
                    DerivedPolicies::append(container, elem());
                }
                else
                {
                    PyErr_SetString(PyExc_TypeError, 
                        "Attempting to append an invalid type");
                    throw_error_already_set();
                }
            }
        }
        
        static void
        base_extend(Container& container, object v)
        {
            std::vector<data_type> temp;
            container_utils::extend_container(temp, v);
            DerivedPolicies::extend(container, temp.begin(), temp.end());
        }
#endif
      static iterator_type
      moveToPos(Container &container, index_type i){
	iterator_type pos;
	index_type idx=0;
	pos=container.begin();
	while(idx<i && pos != container.end()){
	  pos++;
	  idx++;
	}
	if(pos==container.end()){
	  PyErr_SetObject(PyExc_IndexError,PyInt_FromLong(i));
	  python::throw_error_already_set();
	}
	return pos;
      }

    };
       
}} // namespace boost::python 

#endif 
