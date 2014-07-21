#ifndef RDKIT_IMPORT_ARRAY_INCLUDED
#define RDKIT_IMPORT_ARRAY_INCLUDED

#if PY_MAJOR_VERSION >= 3
  inline void * rdkit_import_array()
#else
  inline void rdkit_import_array()
#endif
  {
    // numpy's import_array is defined as a macro that expands into a block
    // statement that inlines a return. In python3 it returns a NULL value
    // (to comply with the Py_InitModule signature) so it can't be called 
    // directly from within the BOOST_PYTHON_MODULE init function (that
    // returns void)
    import_array();
  }

#endif
