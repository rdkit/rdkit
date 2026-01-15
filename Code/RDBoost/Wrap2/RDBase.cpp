
// Copyright (c) 2004-2019 greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <fstream>
#include <RDBoost/Wrap_nb.h>
#include <RDBoost/python_streambuf_nb.h>
#include <RDGeneral/versions.h>
#include <RDGeneral/Invariant.h>
#include <cstdlib>

#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/tuple.h>

#include <RDGeneral/RDLog.h>

namespace nb = nanobind;
using namespace nb::literals;
namespace logging = boost::logging;

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

  ~PyLogStream() override {
#if PY_VERSION_HEX < 0x30d0000
    if (!_Py_IsFinalizing()) {
#else
    if (!Py_IsFinalizing()) {
#endif
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

  void exit(nb::object exc_type, nb::object exc_val, nb::object traceback) {
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

  static void wrap(nb::module_ &m) {
    nb::class_<wt>(m, "streambuf")
        .def(nb::init<nb::object &, std::size_t>(), "python_file_obj"_a,
             "buffer_size"_a = 0, nb::keep_alive<1, 2>());
  }
};

struct python_ostream_wrapper {
  typedef boost_adaptbx::python::ostream wt;

  static void wrap(nb::module_ &m) {
    nb::class_<std::ostream>(m, "std_ostream");
    nb::class_<wt, std::ostream>(m, "ostream")
        .def(nb::init<nb::object &, std::size_t>(), "python_file_obj"_a,
             "buffer_size"_a = 0);
  }
};

void seedRNG(unsigned int seed) { std::srand(seed); }
#if 0
/// Simple Boost.Python custom converter from pathlib.Path to std::string
template <typename T = std::string>
struct path_converter {
  path_converter() {
    python::converter::registry::push_back(&path_converter::convertible,
                                           &path_converter::construct,
                                           boost::python::type_id<T>());
  }

  /// Check PyObject is a pathlib.Path
  static void *convertible(PyObject *object) {
    // paranoia
    if (object == nullptr) {
      return nullptr;
    }
    python::object boost_object(python::handle<>(python::borrowed(object)));

    std::string object_classname = boost::python::extract<std::string>(
        boost_object.attr("__class__").attr("__name__"));
    // pathlib.Path is always specialized to the below derived classes
    if (object_classname == "WindowsPath" || object_classname == "PosixPath") {
      return object;
    }

    return nullptr;
  }

  /// Construct a std::string from pathlib.Path using its own __str__ attribute
  static void construct(
      PyObject *object,
      boost::python::converter::rvalue_from_python_stage1_data *data) {
    void *storage =
        ((boost::python::converter::rvalue_from_python_storage<T> *)data)
            ->storage.bytes;
    python::object boost_object{python::handle<>{python::borrowed(object)}};
    new (storage)
        T{boost::python::extract<std::string>{boost_object.attr("__str__")()}};
    data->convertible = storage;
  }
};

/// Convert a Python str to a std::string_view
template <typename T = std::string_view>
struct string_view_from_python_converter {
  string_view_from_python_converter() {
    python::converter::registry::push_back(
        &string_view_from_python_converter::convertible,
        &string_view_from_python_converter::construct, python::type_id<T>());
  }

  /// Check PyObject is a str
  static void *convertible(PyObject *object) {
    if (object != nullptr && PyUnicode_Check(object)) {
      return object;
    }
    return nullptr;
  }

  /// Construct a std::string_view from the internal string of the PyObject.
  /// This shouldn´t fail, as we block the Python thread until the C++ code
  /// returns, so the PyObject and the internal string should remain valid.
  static void construct(
      PyObject *object,
      python::converter::rvalue_from_python_stage1_data *data) {
    const char *tmp = PyUnicode_AsUTF8(object);

    void *storage = ((python::converter::rvalue_from_python_storage<T> *)data)
                        ->storage.bytes;
    new (storage) T{tmp};
    data->convertible = storage;
  }
};

struct string_view_to_python_converter {
  static PyObject *convert(const std::string_view &s) {
    return python::incref(python::str(s.data(), s.size()).ptr());
  }
};
#endif

}  // namespace

NB_MODULE(rdBase, m) {
  m.doc() = "Module containing basic definitions for wrapped C++ code";
  RDLog::InitLogs();
  // RegisterVectorConverter<int>();
  // RegisterVectorConverter<unsigned>();
  // RegisterVectorConverter<size_t>("UnsignedLong_Vect");
  // RegisterVectorConverter<boost::uint64_t>("VectSizeT");

  // RegisterVectorConverter<double>();
  // RegisterVectorConverter<std::string>(1);
  // RegisterVectorConverter<std::vector<int>>();
  // RegisterVectorConverter<std::vector<unsigned>>();
  // RegisterVectorConverter<std::vector<double>>();
  // RegisterVectorConverter<std::vector<std::string>>("VectorOfStringVectors");

  // RegisterVectorConverter<std::pair<int, int>>("MatchTypeVect");

  // path_converter();
  // string_view_from_python_converter();

  // RegisterListConverter<int>();
  // RegisterListConverter<std::vector<int>>();
  // RegisterListConverter<std::vector<unsigned int>>();

  nb::exception<IndexErrorException>(m, "IndexErrorException",
                                     PyExc_IndexError);
  nb::exception<ValueErrorException>(m, "ValueErrorException",
                                     PyExc_ValueError);
  nb::exception<KeyErrorException>(m, "KeyErrorException", PyExc_KeyError);
#if INVARIANT_EXCEPTION_METHOD
  nb::exception<Invar::Invariant>(m, "InvariantException", PyExc_RuntimeError);
#endif

  // boost::python::to_python_converter<std::string_view,
  //                                    string_view_to_python_converter>();

  // python::def("_version", _version,
  //             "Deprecated, use the constant rdkitVersion instead");

  m.attr("rdkitVersion") = RDKit::rdkitVersion;
  m.attr("boostVersion") = RDKit::boostVersion;
  m.attr("rdkitBuild") = RDKit::rdkitBuild;

  m.attr("_serializationEnabled") =
#ifdef RDK_USE_BOOST_SERIALIZATION
      true;
#else
      false;
#endif
  m.attr("_iostreamsEnabled") =
#ifdef RDK_USE_BOOST_IOSTREAMS
      true;
#else
      false;
#endif
  m.attr("_multithreadedEnabled") =
#ifdef RDK_BUILD_THREADSAFE_SSS
      true;
#else
      false;
#endif

  m.def("LogToCppStreams", RDLog::InitLogs,
        "Initialize RDKit logs with C++ streams");
  m.def("LogToPythonLogger", LogToPythonLogger,
        "Initialize RDKit logs with Python's logging module");
  m.def("LogToPythonStderr", LogToPythonStderr,
        "Initialize RDKit logs with Python's stderr stream");
  m.def("WrapLogs", WrapLogs, "Tee the RDKit logs to Python's stderr stream");

  m.def("EnableLog", EnableLog, "spec"_a);
  m.def("DisableLog", DisableLog, "spec"_a);
  m.def("LogStatus", LogStatus);

  m.def("LogDebugMsg", LogDebugMsg, "msg"_a,
        "Log a message to the RDKit debug logs");
  m.def("LogInfoMsg", LogInfoMsg, "msg"_a,
        "Log a message to the RDKit info logs");
  m.def("LogWarningMsg", LogWarningMsg, "msg"_a,
        "Log a message to the RDKit warning logs");
  m.def("LogErrorMsg", LogErrorMsg, "msg"_a,
        "Log a message to the RDKit error logs");
  m.def("LogMessage", LogMessage, "spec"_a, "msg"_a,
        "Log a message to any rdApp.* log");

  m.def("AttachFileToLog", AttachFileToLog, "spec"_a, "filename"_a,
        "delay"_a = 100, "Causes the log to write to a file");

  m.def("SeedRandomNumberGenerator", seedRNG, "seed"_a,
        "Provides a seed to the standard C random number generator\n"
        "This does not affect pure Python code, but is relevant to some "
        "of the RDKit C++ components.");

  python_streambuf_wrapper::wrap(m);
  python_ostream_wrapper::wrap(m);

  nb::class_<BlockLogs>(
      m, "BlockLogs",
      "Temporarily block logs from outputting while this instance is in scope.")
      .def(nb::init<>())
      .def("__enter__", &BlockLogs::enter, nb::rv_policy::reference)
      .def("__exit__", &BlockLogs::exit);
}
