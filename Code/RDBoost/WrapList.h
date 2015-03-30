#include <list>
#include <algorithm>
#include <boost/python.hpp>

// from
// http://stackoverflow.com/questions/6776888/wrapping-a-list-of-structs-with-boost-python
//  * modified to allow ptr classes and actually compile as intended
template<class T>
struct listwrap
{
    typedef typename std::list<T> L;
    typedef typename L::value_type value_type;
    typedef typename L::iterator iter_type;

    static void add(L & x, value_type const& v)
    {
        x.push_back(v);
    }

    static bool in(L const& x, value_type const& v)
    {
        return std::find(x.begin(), x.end(), v) != x.end();
    }

    static int index(L const& x, value_type const& v)
    {
        int i = 0;
        for(typename L::const_iterator it=x.begin(); it!=x.end(); ++it,++i)
            if( *it == v ) return i;

        PyErr_SetString(PyExc_ValueError, "Value not in the list");
        throw boost::python::error_already_set();
    }

    static void del(L& x, int i)
    {
        if( i<0 ) 
            i += x.size();

        iter_type it = x.begin();
        for (int pos = 0; pos < i; ++pos)
            ++it;

        if( i >= 0 && i < (int)x.size() ) {
            x.erase(it);
        } else {
            PyErr_SetString(PyExc_IndexError, "Index out of range");
            boost::python::throw_error_already_set();
        }
    }

    static value_type& get(L& x, int i)
    {
        if( i < 0 ) 
            i += x.size();

        if( i >= 0 && i < (int)x.size() ) {
            iter_type it = x.begin(); 
            for(int pos = 0; pos < i; ++pos)
                ++it;
            return *it;                             
        } else {
            PyErr_SetString(PyExc_IndexError, "Index out of range");
            throw boost::python::error_already_set();
        }
    }

    static void set(L& x, int i, value_type const& v)
    {
        if( i < 0 ) 
            i += x.size();

        if( i >= 0 && i < (int)x.size() ) {
            iter_type it = x.begin(); 
            for(int pos = 0; pos < i; ++pos)
                ++it;
            *it = v;
        } else {
            PyErr_SetString(PyExc_IndexError, "Index out of range");
            boost::python::throw_error_already_set();
        }
    }
};


template<class T>
void export_ConstSTLListOfPtrs(const char* typeName)
{
    using namespace boost::python;

    class_<std::list<T> >(typeName)
        .def("__len__", &std::list<T>::size)
        .def("__getitem__", &listwrap<T>::get,
             python::return_internal_reference<1,
             python::with_custodian_and_ward_postcall<0,1> >())
        .def("__contains__", &listwrap<T>::in)
        .def("__iter__", iterator<std::list<T> >())
        .def("index", &listwrap<T>::index);
}
