// $Id: SDMolSupplier.cpp 585 2008-03-30 13:36:56Z glandrum $
//
//  Copyright (C) 2003-2008  Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#define NO_IMPORT_ARRAY
#include <boost/python.hpp>
#include <string>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/algorithm/string.hpp>

namespace io=boost::iostreams;

//ours
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/RDKitBase.h>
#include <RDBoost/Wrap.h>

#include "MolSupplier.h"

namespace python = boost::python;

namespace RDKit {
  // Note that this returns a pointer to the supplier itself, so be careful
  // that it doesn't get deleted by python!
  ForwardSDMolSupplier *MolSupplIter(ForwardSDMolSupplier *suppl){
    return suppl;
  }

  ROMol *MolSupplNext(ForwardSDMolSupplier *suppl){
    ROMol *res=0;
    if (!suppl->atEnd()) {
      try {
        res=suppl->next();
      } catch(...){
        res=0;
      }
    }
    if(!res && suppl->atEnd()) {
      PyErr_SetString(PyExc_StopIteration,"End of supplier hit");
      throw boost::python::error_already_set();
    }
    return res;
  }

  ForwardSDMolSupplier *createForwardSupplier(std::string filename,bool sanitize,
                                              bool removeHs){
    std::vector<std::string> splitName;
    boost::split(splitName,filename,boost::is_any_of("."));
    io::filtering_istream *strm=new io::filtering_istream();
    if(splitName.back()=="gz"){
      strm->push(io::gzip_decompressor());
    } else if(splitName.back()=="bz"){
      strm->push(io::bzip2_decompressor());
    } else if(splitName.back()=="sdf"){
    } else {
      std::string errorTxt="Unrecognized extension: "+splitName.back();
      throw_value_error(errorTxt);
    }

    strm->push(io::file_source(filename));
    
    ForwardSDMolSupplier *res=new ForwardSDMolSupplier(strm,true,sanitize,removeHs);
    return res;
  }
  
  std::string fsdMolSupplierClassDoc="A class which supplies molecules from an SD file.\n \
\n \
  Usage examples:\n \
\n \
    1) Lazy evaluation: the molecules are not constructed until we ask for them:\n \
       >>> suppl = SDMolSupplier('in.smi')\n \
       >>> for mol in suppl:\n \
       ...    mol.GetNumAtoms()\n \
\n \
  Properties in the SD file are used to set properties on each molecule.\n\
  The properties are accessible using the mol.GetProp(propName) method.\n\
\n";
  struct compressedsdmolsup_wrap {
    static void wrap() {
      python::class_<ForwardSDMolSupplier,boost::noncopyable>("ForwardSDMolSupplier",
                                                              fsdMolSupplierClassDoc.c_str(),
                                                              python::no_init)
	.def("__iter__", (ForwardSDMolSupplier *(*)(ForwardSDMolSupplier *))&MolSupplIter,
	     python::return_internal_reference<1>() )
	.def("next", (ROMol *(*)(ForwardSDMolSupplier *))&MolSupplNext,
	     "Returns the next molecule in the file.  Raises _StopIteration_ on EOF.\n",
	     python::return_value_policy<python::manage_new_object>())
	;
      python::def("CompressedSDMolSupplier",createForwardSupplier,
                  (python::arg("fileName"),python::arg("sanitize")=true,
                   python::arg("removeHs")=true),
                  python::return_value_policy<python::manage_new_object>());
    };
  };
}

void wrap_compressedsdsupplier() {
  RDKit::compressedsdmolsup_wrap::wrap();
}
