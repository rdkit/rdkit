
// Copyright (c) 2004-2019 greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/python.h>
#include <iostream>
#include <fstream>
#include <RDBoost/Wrap.h>
#include <RDBoost/python_streambuf.h>
#include <RDGeneral/versions.h>
#include <RDGeneral/Invariant.h>
#include <cstdlib>

#include <RDGeneral/RDLog.h>
#if 0
#include <boost/log/functions.hpp>
#if defined(BOOST_HAS_THREADS)
#include <boost/log/extra/functions_ts.hpp>
#endif
#endif

namespace python = boost::python;
namespace logging = boost::logging;

std::string _version() { return "$Id$"; }

// std::ostream wrapper around Python's stderr stream
struct PyErrStream : std::ostream, std::streambuf {
  static thread_local std::string buffer;

  PyErrStream() : std::ostream(this) {
    // All done!
  }

  int overflow(int c) override {
    write(c);
    return 0;
  }

  void write(char c) {
    if (c == '\n') {
      PyGILStateHolder h;
      PySys_WriteStderr("%s\n", buffer.c_str());
      buffer.clear();
    } else {
      buffer += c;
    }
  }
};

// std::ostream wrapper around Python's logging module
struct PyLogStream : std::ostream, std::streambuf {
  static thread_local std::string buffer;
  PyObject *logfn = nullptr;

  PyLogStream(std::string level) : std::ostream(this) {
    PyObject *module = PyImport_ImportModule("logging");
    PyObject *logger = nullptr;

    if (module != nullptr) {
      logger = PyObject_CallMethod(module, "getLogger", "s", "rdkit");
      Py_DECREF(module);
    }

    if (logger != nullptr) {
      logfn = PyObject_GetAttrString(logger, level.c_str());
      Py_DECREF(logger);
    }

    if (PyErr_Occurred()) {
      PyErr_Print();
    }
  }

  ~PyLogStream() {
    if (!_Py_IsFinalizing()) {
      Py_XDECREF(logfn);
    }
  }

  int overflow(int c) override {
    write(c);
    return 0;
  }

  void write(char c) {
    if (logfn == nullptr) {
      return;
    }

    if (c == '\n') {
      PyGILStateHolder h;
      PyObject *result = PyObject_CallFunction(logfn, "s", buffer.c_str());
      Py_XDECREF(result);
      buffer.clear();
    } else {
      buffer += c;
    }
  }
};

// per-thread buffers for the Python loggers
thread_local std::string PyErrStream::buffer;
thread_local std::string PyLogStream::buffer;

void LogToPythonLogger() {
  static PyLogStream debug("debug");
  static PyLogStream info("info");
  static PyLogStream warning("warning");
  static PyLogStream error("error");

  rdDebugLog = std::make_shared<logging::rdLogger>(&debug);
  rdInfoLog = std::make_shared<logging::rdLogger>(&info);
  rdWarningLog = std::make_shared<logging::rdLogger>(&warning);
  rdErrorLog = std::make_shared<logging::rdLogger>(&error);
}

void LogToPythonStderr() {
  static PyErrStream debug;
  static PyErrStream info;
  static PyErrStream warning;
  static PyErrStream error;

  rdDebugLog = std::make_shared<logging::rdLogger>(&debug);
  rdInfoLog = std::make_shared<logging::rdLogger>(&info);
  rdWarningLog = std::make_shared<logging::rdLogger>(&warning);
  rdErrorLog = std::make_shared<logging::rdLogger>(&error);
}

void WrapLogs() {
  static PyErrStream debug;    //("RDKit DEBUG: ");
  static PyErrStream error;    //("RDKit ERROR: ");
  static PyErrStream warning;  //("RDKit WARNING: ");
  static PyErrStream info;     //("RDKit INFO: ");

  if (!rdDebugLog || !rdInfoLog || !rdErrorLog || !rdWarningLog) {
    RDLog::InitLogs();
  }

  rdDebugLog->SetTee(debug);
  rdInfoLog->SetTee(info);
  rdWarningLog->SetTee(warning);
  rdErrorLog->SetTee(error);
}

void EnableLog(std::string spec) { logging::enable_logs(spec); }

void DisableLog(std::string spec) { logging::disable_logs(spec); }

std::string LogStatus() { return logging::log_status(); }

void AttachFileToLog(std::string spec, std::string filename, int delay = 100) {
  (void)spec;
  (void)filename;
  (void)delay;
#if 0
#if defined(BOOST_HAS_THREADS)
  logging::manipulate_logs(spec)
    .add_appender(logging::ts_appender(logging::write_to_file(filename),
				       delay));
#else
  logging::manipulate_logs(spec)
    .add_appender(logging::write_to_file(filename));

#endif
#endif
}

void LogDebugMsg(const std::string &msg) {
  // NOGIL nogil;
  BOOST_LOG(rdDebugLog) << msg << std::endl;
}

void LogInfoMsg(const std::string &msg) {
  // NOGIL nogil;
  BOOST_LOG(rdInfoLog) << msg << std::endl;
}

void LogWarningMsg(const std::string &msg) {
  // NOGIL nogil;
  BOOST_LOG(rdWarningLog) << msg << std::endl;
}

void LogErrorMsg(const std::string &msg) {
  // NOGIL nogil;
  BOOST_LOG(rdErrorLog) << msg << std::endl;
}

void LogMessage(std::string spec, std::string msg) {
  if (spec == "rdApp.error") {
    LogErrorMsg(msg);
  } else if (spec == "rdApp.warning") {
    LogWarningMsg(msg);
  } else if (spec == "rdApp.info") {
    LogInfoMsg(msg);
  } else if (spec == "rdApp.debug") {
    LogDebugMsg(msg);
  }
}

class BlockLogs : public boost::noncopyable {
 public:
  BlockLogs() : m_log_setter{new RDLog::LogStateSetter} {}
  ~BlockLogs() = default;

  BlockLogs *enter() { return this; }

  void exit(python::object exc_type, python::object exc_val,
            python::object traceback) {
    RDUNUSED_PARAM(exc_type);
    RDUNUSED_PARAM(exc_val);
    RDUNUSED_PARAM(traceback);
    m_log_setter.reset();
  }

 private:
  std::unique_ptr<RDLog::LogStateSetter> m_log_setter;
};

namespace {
struct python_streambuf_wrapper {
  typedef boost_adaptbx::python::streambuf wt;

  static void wrap() {
    using namespace boost::python;
    class_<wt, boost::noncopyable>("streambuf", no_init)
        .def(init<object &, std::size_t>(
            (arg("python_file_obj"), arg("buffer_size") = 0),
            "documentation")[with_custodian_and_ward_postcall<0, 2>()]);
  }
};

struct python_ostream_wrapper {
  typedef boost_adaptbx::python::ostream wt;

  static void wrap() {
    using namespace boost::python;
    class_<std::ostream, boost::noncopyable>("std_ostream", no_init);
    class_<wt, boost::noncopyable, bases<std::ostream>>("ostream", no_init)
        .def(init<object &, std::size_t>(
            (arg("python_file_obj"), arg("buffer_size") = 0)));
  }
};

void seedRNG(unsigned int seed) { std::srand(seed); }
}  // namespace

BOOST_PYTHON_MODULE(rdBase) {
  python::scope().attr("__doc__") =
      "Module containing basic definitions for wrapped C++ code\n"
      "\n";
  RDLog::InitLogs();
  RegisterVectorConverter<int>();
  RegisterVectorConverter<unsigned>();
  RegisterVectorConverter<double>();
  RegisterVectorConverter<std::string>(1);
  RegisterVectorConverter<std::vector<int>>();
  RegisterVectorConverter<std::vector<unsigned>>();
  RegisterVectorConverter<std::vector<double>>();

  RegisterListConverter<int>();
  RegisterListConverter<std::vector<int>>();
  RegisterListConverter<std::vector<unsigned int>>();

  python::register_exception_translator<IndexErrorException>(
      &translate_index_error);
  python::register_exception_translator<ValueErrorException>(
      &translate_value_error);
  python::register_exception_translator<KeyErrorException>(
      &translate_key_error);

#if INVARIANT_EXCEPTION_METHOD
  python::register_exception_translator<Invar::Invariant>(
      &translate_invariant_error);
#endif

  python::def("_version", _version,
              "Deprecated, use the constant rdkitVersion instead");

  python::scope().attr("rdkitVersion") = RDKit::rdkitVersion;
  python::scope().attr("boostVersion") = RDKit::boostVersion;
  python::scope().attr("rdkitBuild") = RDKit::rdkitBuild;

  python::def("LogToCppStreams", RDLog::InitLogs,
              "Initialize RDKit logs with C++ streams");
  python::def("LogToPythonLogger", LogToPythonLogger,
              "Initialize RDKit logs with Python's logging module");
  python::def("LogToPythonStderr", LogToPythonStderr,
              "Initialize RDKit logs with Python's stderr stream");
  python::def("WrapLogs", WrapLogs,
              "Tee the RDKit logs to Python's stderr stream");

  python::def("EnableLog", EnableLog);
  python::def("DisableLog", DisableLog);
  python::def("LogStatus", LogStatus);

  python::def("LogDebugMsg", LogDebugMsg,
              "Log a message to the RDKit debug logs");
  python::def("LogInfoMsg", LogInfoMsg, "Log a message to the RDKit info logs");
  python::def("LogWarningMsg", LogWarningMsg,
              "Log a message to the RDKit warning logs");
  python::def("LogErrorMsg", LogErrorMsg,
              "Log a message to the RDKit error logs");
  python::def("LogMessage", LogMessage, "Log a message to any rdApp.* log");

  python::def("AttachFileToLog", AttachFileToLog,
              "Causes the log to write to a file",
              (python::arg("spec"), python::arg("filename"),
               python::arg("delay") = 100));

  python::def("SeedRandomNumberGenerator", seedRNG,
              "Provides a seed to the standard C random number generator\n"
              "This does not affect pure Python code, but is relevant to some "
              "of the RDKit C++ components.",
              (python::arg("seed")));

  python_streambuf_wrapper::wrap();
  python_ostream_wrapper::wrap();

  python::class_<BlockLogs, boost::noncopyable>(
      "BlockLogs",
      "Temporarily block logs from outputting while this instance is in scope.",
      python::init<>())
      .def("__enter__", &BlockLogs::enter,
           python::return_internal_reference<>())
      .def("__exit__", &BlockLogs::exit);
}
