
// Copyright (c) 2004-2026 greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <array>
#include <cstdlib>
#include <fstream>

#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/string.h>

#include <RDBoost/Wrap_nb.h>
#include <RDBoost/python_capsule_nb.h>
#include <RDBoost/python_streambuf_nb.h>
#include <RDGeneral/versions.h>
#include <RDGeneral/Invariant.h>
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

  PyLogStream(const std::string &level) : std::ostream(this) {
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
    if (!nb::is_alive()) {
      return;
    }
    Py_XDECREF(logfn);
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
  using logger_array_t = std::array<std::shared_ptr<logging::rdLogger>, 4>;

  constexpr const char *rdkLoggerCapsule = "__rdkit_loggers__";

  auto capsuleData = getCapsulePtr(rdkLoggerCapsule);
  if (capsuleData) {
    // We already created the loggers, so just reenable them
    auto loggers = static_cast<logger_array_t *>(capsuleData);

    rdDebugLog = (*loggers)[0];
    rdInfoLog = (*loggers)[1];
    rdWarningLog = (*loggers)[2];
    rdErrorLog = (*loggers)[3];

  } else {
    // Create the loggers and store them in a capsule that will
    // destroy them when Python is shutting down, and also allows
    // us to retrieve them later if we need to reenable them.

    auto makeLogger = [](const char *level) {
      return std::make_shared<logging::rdLogger>(new PyLogStream(level), true);
    };

    auto cleanUpFn = [](void *ptr) noexcept {
      // Reset log endpoints to stdout/stderr before destroying
      // the current loggers
      RDLog::InitLogs();

      auto loggers = static_cast<logger_array_t *>(ptr);
      if (!loggers) {
        return;
      }

      // By this point, this should be the only instance
      // of the loggers shared_ptrs, but reset them to make
      // sure the streams are destroyed.
      for (auto logger : *loggers) {
        logger.reset();
      }

      delete loggers;
    };

    rdDebugLog = makeLogger("debug");
    rdInfoLog = makeLogger("info");
    rdWarningLog = makeLogger("warning");
    rdErrorLog = makeLogger("error");

    auto logStreamV =
        new logger_array_t{rdDebugLog, rdInfoLog, rdWarningLog, rdErrorLog};

    // Store the loggers in a capsule so we can retrieve them later
    nb::capsule loggerCapsule(logStreamV, cleanUpFn);
    installCapsule(loggerCapsule, rdkLoggerCapsule);
  }
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

struct PyCaptureErrorLog : boost::noncopyable {
  PyCaptureErrorLog *enter() { return this; }
  void exit(nb::object /*exc_type*/, nb::object /*exc_val*/,
            nb::object /*traceback*/) {
    if (m_capturer) {
      m_messages = m_capturer->messages();
      m_capturer.reset();
    }
  }
  std::string messages() const {
    return m_capturer ? m_capturer->messages() : m_messages;
  }

 private:
  std::unique_ptr<RDLog::CaptureErrorLog> m_capturer{
      new RDLog::CaptureErrorLog};
  std::string m_messages;
};

namespace {

void seedRNG(unsigned int seed) { std::srand(seed); }

}  // namespace

NB_MODULE(rdBase, m) {
  m.doc() = "Module containing basic definitions for wrapped C++ code";
  RDLog::InitLogs();
  nb::exception<IndexErrorException>(m, "IndexErrorException",
                                     PyExc_IndexError);
  nb::exception<ValueErrorException>(m, "ValueErrorException",
                                     PyExc_ValueError);
  nb::register_exception_translator(
      [](const std::exception_ptr &p, void * /* unused */) {
        try {
          std::rethrow_exception(p);
        } catch (const KeyErrorException &e) {
          PyErr_SetString(PyExc_KeyError, e.key().c_str());
        }
      });
#if INVARIANT_EXCEPTION_METHOD
  nb::register_exception_translator(
      [](const std::exception_ptr &p, void * /* unused */) {
        try {
          std::rethrow_exception(p);
        } catch (const Invar::Invariant &e) {
          PyErr_SetString(PyExc_RuntimeError, e.toUserString().c_str());
        }
      });
#endif

  m.attr("rdkitVersion") = RDKit::rdkitVersion;
  m.attr("boostVersion") = RDKit::boostVersion;
  m.attr("rdkitBuild") = RDKit::rdkitBuild;
  m.attr("_wrapperType") = "nanobind";

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

  nb::class_<BlockLogs>(
      m, "BlockLogs",
      "Temporarily block logs from outputting while this instance is in scope.")
      .def(nb::init<>())
      .def("__enter__", &BlockLogs::enter, nb::rv_policy::reference_internal)
      .def("__exit__", &BlockLogs::exit, "excType"_a = nb::none(),
           "excValue"_a = nb::none(), "traceback"_a = nb::none());

  nb::class_<PyCaptureErrorLog>(
      m, "CaptureErrorLog",
      R"DOC(Captures messages from rdErrorLog while this instance is in scope.
      Can be used as a context manager. The ``messages`` property is
      accessible both inside the context and after it exits.
      Nesting is supported: inner captures shadow outer ones.

      Example::

        with rdBase.CaptureErrorLog() as capture:
            rdkit_function_that_may_fail()
        print(capture.messages)
  )DOC")
      .def(nb::init<>())
      .def("__enter__", &PyCaptureErrorLog::enter,
           nb::rv_policy::reference_internal)
      .def("__exit__", &PyCaptureErrorLog::exit, "excType"_a = nb::none(),
           "excValue"_a = nb::none(), "traceback"_a = nb::none())
      .def_prop_ro("messages", &PyCaptureErrorLog::messages,
                   "Messages captured from rdErrorLog.");
}
