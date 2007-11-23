// $Id$
//
// Copyright (c) 2004-2006 greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <boost/python.hpp>
#define PY_ARRAY_UNIQUE_SYMBOL rdbase_array_API
#include <RDBoost/Wrap.h>
#include "Numeric/arrayobject.h"

#include <RDGeneral/RDLog.h>
#include <boost/log/functions.hpp>
#if defined(BOOST_HAS_THREADS)
#include <boost/log/extra/functions_ts.hpp>
#endif



namespace python = boost::python;
namespace logging = boost::logging;

std::string _version(){
  return "$Id$";
}

void EnableLog(std::string spec){
  logging::enable_logs(spec);
}
void DisableLog(std::string spec){
  logging::disable_logs(spec);
}
void AttachFileToLog(std::string spec,std::string filename,int delay=100){
#if defined(BOOST_HAS_THREADS)
  logging::manipulate_logs(spec)
    .add_appender(logging::ts_appender(logging::write_to_file(filename),
				       delay));
#else
  logging::manipulate_logs(spec)
    .add_appender(logging::write_to_file(filename));

#endif
}
void LogMessage(std::string spec,std::string msg){
  logging::logger theLog(spec);
  if(theLog.is_enabled(logging::level::default_)){
    *(theLog.stream().stream()) << msg;
  }
}



BOOST_PYTHON_MODULE(rdBase)
{
  python::scope().attr("__doc__") =
    "Module containing basic definitions for wrapped C++ code\n"
    "  No user-serviceable parts inside.\n"
    "\n"
    ;
  RDLog::InitLogs();
  RegisterVectorConverter<int>();
  RegisterVectorConverter<unsigned>();
  RegisterVectorConverter<double>();
  RegisterVectorConverter<std::string>(1);
  RegisterVectorConverter< std::vector<int> >();
  RegisterVectorConverter< std::vector<unsigned> >();
  RegisterVectorConverter< std::vector<double> >();

  RegisterListConverter< int >();
  RegisterListConverter< std::vector<int> >();

  
  import_array();

  python::register_exception_translator<IndexErrorException>(&translate_index_error);
  python::register_exception_translator<ValueErrorException>(&translate_value_error);

  python::def("_version",_version);

  python::def("EnableLog",EnableLog);
  python::def("DisableLog",DisableLog);
  python::def("AttachFileToLog",AttachFileToLog,"Causes the log to write to a file",
	      (python::arg("spec"),python::arg("filename"),python::arg("delay")=100));
  python::def("LogMessage",LogMessage);
}
