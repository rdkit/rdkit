//
//  This file is part of the CCTBX distribution:
//     http://cctbx.sourceforge.net/
//  Downloaded from here: 
//    http://cctbx.svn.sourceforge.net/viewvc/cctbx/trunk/boost_adaptbx/python_streambuf.h?revision=13619
//
//  Copyright (c) 2006, The Regents of the University of
//  California, through Lawrence Berkeley National Laboratory (subject to
//  receipt of any required approvals from the U.S. Dept. of Energy).  All
//  rights reserved.
//
//  The license is here:
//    http://cctbx.svn.sourceforge.net/viewvc/cctbx/trunk/boost_adaptbx/LICENSE_2_0.txt?revision=5148
//
#ifndef BOOST_ADAPTBX_PYTHON_STREAMBUF_H
#define BOOST_ADAPTBX_PYTHON_STREAMBUF_H

#include <boost/python/object.hpp>
#include <boost/python/str.hpp>
#include <boost/python/extract.hpp>

#include <boost/optional.hpp>
#include <boost/utility/typed_in_place_factory.hpp>

//#include <tbxx/error_utils.hpp>
#include <RDGeneral/Invariant.h>

#include <streambuf>
#include <iostream>

namespace boost_adaptbx { namespace python {

namespace bp = boost::python;

/// A stream buffer getting data from and putting data into a Python file object
/** The aims are as follow:

    - Given a C++ function acting on a standard stream, e.g.

      \code
      void read_inputs(std::istream& input) {
        ...
        input >> something >> something_else;
      }
      \endcode

      and given a piece of Python code which creates a file-like object,
      to be able to pass this file object to that C++ function, e.g.

      \code
      import gzip
      gzip_file_obj = gzip.GzipFile(...)
      read_inputs(gzip_file_obj)
      \endcode

      and have the standard stream pull data from and put data into the Python
      file object.

    - When Python \c read_inputs() returns, the Python object is able to
      continue reading or writing where the C++ code left off.

    - Operations in C++ on mere files should be competitively fast compared
      to the direct use of \c std::fstream.


    \b Motivation

      - the standard Python library offer of file-like objects (files,
        compressed files and archives, network, ...) is far superior to the
        offer of streams in the C++ standard library and Boost C++ libraries.

      - i/o code involves a fair amount of text processing which is more
        efficiently prototyped in Python but then one may need to rewrite
        a time-critical part in C++, in as seamless a manner as possible.

    \b Usage

    This is 2-step:

      - a trivial wrapper function

        \code
          using boost_adaptbx::python::streambuf;
          void read_inputs_wrapper(streambuf& input)
          {
            streambuf::istream is(input);
            read_inputs(is);
          }

          def("read_inputs", read_inputs_wrapper);
        \endcode

        which has to be written every time one wants a Python binding for
        such a C++ function.

      - the Python side

        \code
          from boost.python import streambuf
          read_inputs(streambuf(python_file_obj=obj, buffer_size=1024))
        \endcode

        \c buffer_size is optional. See also: \c default_buffer_size

  Note: references are to the C++ standard (the numbers between parentheses
  at the end of references are margin markers).
*/
class streambuf : public std::basic_streambuf<char>
{
  private:
    typedef std::basic_streambuf<char> base_t;

  public:
    /* The syntax
        using base_t::char_type;
       would be nicer but Visual Studio C++ 8 chokes on it
    */
    typedef base_t::char_type   char_type;
    typedef base_t::int_type    int_type;
    typedef base_t::pos_type    pos_type;
    typedef base_t::off_type    off_type;
    typedef base_t::traits_type traits_type;

    // work around Visual C++ 7.1 problem
    inline static int
    traits_type_eof() { return traits_type::eof(); }

    /// The default size of the read and write buffer.
    /** They are respectively used to buffer data read from and data written to
        the Python file object. It can be modified from Python.
    */
    const static std::size_t default_buffer_size=1024;

    /// Construct from a Python file object
    /** if buffer_size is 0 the current default_buffer_size is used.
    */
    streambuf(
      bp::object& python_file_obj,
      std::size_t buffer_size_=0)
    :
      py_read (getattr(python_file_obj, "read",  bp::object())),
      py_write(getattr(python_file_obj, "write", bp::object())),
      py_seek (getattr(python_file_obj, "seek",  bp::object())),
      py_tell (getattr(python_file_obj, "tell",  bp::object())),
      buffer_size(buffer_size_ != 0 ? buffer_size_ : default_buffer_size),
      write_buffer(0),
      pos_of_read_buffer_end_in_py_file(0),
      pos_of_write_buffer_end_in_py_file(buffer_size),
      farthest_pptr(0)
    {
      TEST_ASSERT(buffer_size != 0);
      /* Some Python file objects (e.g. sys.stdout and sys.stdin)
         have non-functional seek and tell. If so, assign None to
         py_tell and py_seek.
       */
      if (py_tell != bp::object()) {
        try {
          off_type py_pos = bp::extract<off_type>(py_tell());
          if(py_seek != bp::object()){
            /* Make sure we can actually seek. 
               bzip2 readers from python have a seek method, but it fails
               when they are in write mode.
             */
            py_seek(py_pos);
          }
        }
        catch (bp::error_already_set&) {
          py_tell = bp::object();
          py_seek = bp::object();
          /* Boost.Python does not do any Python exception handling whatsoever
             So we need to catch it by hand like so.
           */
          PyErr_Clear();
        }
      }

      if (py_write != bp::object()) {
        // C-like string to make debugging easier
        write_buffer = new char[buffer_size + 1];
        write_buffer[buffer_size] = '\0';
        setp(write_buffer, write_buffer + buffer_size);  // 27.5.2.4.5 (5)
        farthest_pptr = pptr();
      }
      else {
        // The first attempt at output will result in a call to overflow
        setp(0, 0);
      }

      if (py_tell != bp::object()) {
        off_type py_pos = bp::extract<off_type>(py_tell());
        pos_of_read_buffer_end_in_py_file = py_pos;
        pos_of_write_buffer_end_in_py_file = py_pos;
      }
    }

    /// Mundane destructor freeing the allocated resources
    virtual ~streambuf() {
      if (write_buffer) delete[] write_buffer;
    }

    /// C.f. C++ standard section 27.5.2.4.3
    /** It is essential to override this virtual function for the stream
        member function readsome to work correctly (c.f. 27.6.1.3, alinea 30)
     */
    virtual std::streamsize showmanyc() {
      int_type const failure = traits_type::eof();
      int_type status = underflow();
      if (status == failure) return -1;
      return egptr() - gptr();
    }

    /// C.f. C++ standard section 27.5.2.4.3
    virtual int_type underflow() {
      int_type const failure = traits_type::eof();
      if (py_read == bp::object()) {
        throw std::invalid_argument(
          "That Python file object has no 'read' attribute");
      }
      read_buffer = py_read(buffer_size);
      char *read_buffer_data;
      bp::ssize_t py_n_read;
      if (PyBytes_AsStringAndSize(read_buffer.ptr(),
                                   &read_buffer_data, &py_n_read) == -1) {
        setg(0, 0, 0);
        throw std::invalid_argument(
          "The method 'read' of the Python file object "
          "did not return a string.");
      }
      off_type n_read = (off_type)py_n_read;
      pos_of_read_buffer_end_in_py_file += n_read;
      setg(read_buffer_data, read_buffer_data, read_buffer_data + n_read);
      // ^^^27.5.2.3.1 (4)
      if (n_read == 0) return failure;
      return traits_type::to_int_type(read_buffer_data[0]);
    }

    /// C.f. C++ standard section 27.5.2.4.5
    virtual int_type overflow(int_type c=traits_type_eof()) {
      if (py_write == bp::object()) {
        throw std::invalid_argument(
          "That Python file object has no 'write' attribute");
      }
      farthest_pptr = std::max(farthest_pptr, pptr());
      off_type n_written = (off_type)(farthest_pptr - pbase());
      bp::str chunk(pbase(), farthest_pptr);
      py_write(chunk);
      if (!traits_type::eq_int_type(c, traits_type::eof())) {
        py_write(traits_type::to_char_type(c));
        n_written++;
      }
      if (n_written) {
        pos_of_write_buffer_end_in_py_file += n_written;
        setp(pbase(), epptr());
        // ^^^ 27.5.2.4.5 (5)
        farthest_pptr = pptr();
      }
      return traits_type::eq_int_type(
        c, traits_type::eof()) ? traits_type::not_eof(c) : c;
    }

    /// Update the python file to reflect the state of this stream buffer
    /** Empty the write buffer into the Python file object and set the seek
        position of the latter accordingly (C++ standard section 27.5.2.4.2).
        If there is no write buffer or it is empty, but there is a non-empty
        read buffer, set the Python file object seek position to the
        seek position in that read buffer.
    */
    virtual int sync() {
      int result = 0;
      farthest_pptr = std::max(farthest_pptr, pptr());
      if (farthest_pptr && farthest_pptr > pbase()) {
        off_type delta = pptr() - farthest_pptr;
        int_type status = overflow();
        if (traits_type::eq_int_type(status, traits_type::eof())) result = -1;
        if (py_seek != bp::object()) py_seek(delta, 1);
      }
      else if (gptr() && gptr() < egptr()) {
        if (py_seek != bp::object()) py_seek(gptr() - egptr(), 1);
      }
      return result;
    }

    /// C.f. C++ standard section 27.5.2.4.2
    /** This implementation is optimised to look whether the position is within
        the buffers, so as to avoid calling Python seek or tell. It is
        important for many applications that the overhead of calling into Python
        is avoided as much as possible (e.g. parsers which may do a lot of
        backtracking)
    */
    virtual
    pos_type seekoff(off_type off, std::ios_base::seekdir way,
                     std::ios_base::openmode which=  std::ios_base::in
                                                   | std::ios_base::out)
    {
      /* In practice, "which" is either std::ios_base::in or out
         since we end up here because either seekp or seekg was called
         on the stream using this buffer. That simplifies the code
         in a few places.
      */
      int const failure = off_type(-1);

      if (py_seek == bp::object()) {
        throw std::invalid_argument(
          "That Python file object has no 'seek' attribute");
      }

      // we need the read buffer to contain something!
      if (which == std::ios_base::in && !gptr()) {
        if (traits_type::eq_int_type(underflow(), traits_type::eof())) {
          return failure;
        }
      }

      // compute the whence parameter for Python seek
      int whence;
      switch (way) {
        case std::ios_base::beg:
          whence = 0;
          break;
        case std::ios_base::cur:
          whence = 1;
          break;
        case std::ios_base::end:
          whence = 2;
          break;
        default:
          return failure;
      }

      // Let's have a go
      boost::optional<off_type> result = seekoff_without_calling_python(
        off, way, which);
      if (!result) {
        // we need to call Python
        if (which == std::ios_base::out) overflow();
        if (way == std::ios_base::cur) {
          if      (which == std::ios_base::in)  off -= egptr() - gptr();
          else if (which == std::ios_base::out) off += pptr() - pbase();
        }
        py_seek(off, whence);
        result = off_type(bp::extract<off_type>(py_tell()));
        if (which == std::ios_base::in) underflow();
      }
      return *result;
    }

    /// C.f. C++ standard section 27.5.2.4.2
    virtual
    pos_type seekpos(pos_type sp,
                     std::ios_base::openmode which=  std::ios_base::in
                                                   | std::ios_base::out)
    {
      return streambuf::seekoff(sp, std::ios_base::beg, which);
    }

  private:
    bp::object py_read, py_write, py_seek, py_tell;

    std::size_t buffer_size;

    /* This is actually a Python string and the actual read buffer is
       its internal data, i.e. an array of characters. We use a Boost.Python
       object so as to hold on it: as a result, the actual buffer can't
       go away.
    */
    bp::object read_buffer;

    /* A mere array of char's allocated on the heap at construction time and
       de-allocated only at destruction time.
    */
    char *write_buffer;

    off_type pos_of_read_buffer_end_in_py_file,
             pos_of_write_buffer_end_in_py_file;

    // the farthest place the buffer has been written into
    char *farthest_pptr;


    boost::optional<off_type> seekoff_without_calling_python(
      off_type off,
      std::ios_base::seekdir way,
      std::ios_base::openmode which)
    {
      boost::optional<off_type> const failure;

      // Buffer range and current position
      off_type buf_begin, buf_end, buf_cur, upper_bound;
      off_type pos_of_buffer_end_in_py_file;
      if (which == std::ios_base::in) {
        pos_of_buffer_end_in_py_file = pos_of_read_buffer_end_in_py_file;
        buf_begin = reinterpret_cast<std::streamsize>(eback());
        buf_cur = reinterpret_cast<std::streamsize>(gptr());
        buf_end = reinterpret_cast<std::streamsize>(egptr());
        upper_bound = buf_end;
      }
      else if (which == std::ios_base::out) {
        pos_of_buffer_end_in_py_file = pos_of_write_buffer_end_in_py_file;
        buf_begin = reinterpret_cast<std::streamsize>(pbase());
        buf_cur = reinterpret_cast<std::streamsize>(pptr());
        buf_end = reinterpret_cast<std::streamsize>(epptr());
        farthest_pptr = std::max(farthest_pptr, pptr());
        upper_bound = reinterpret_cast<std::streamsize>(farthest_pptr) + 1;
      }
      else {
        CHECK_INVARIANT(0,"unreachable code");
      }

      // Sought position in "buffer coordinate"
      off_type buf_sought;
      if (way == std::ios_base::cur) {
        buf_sought = buf_cur + off;
      }
      else if (way == std::ios_base::beg) {
        buf_sought = buf_end + (off - pos_of_buffer_end_in_py_file);
      }
      else if (way == std::ios_base::end) {
        return failure;
      }
      else {
        CHECK_INVARIANT(0,"unreachable code");
      }

      // if the sought position is not in the buffer, give up
      if (buf_sought < buf_begin || buf_sought >= upper_bound) return failure;

      // we are in wonderland
      if      (which == std::ios_base::in)  gbump(buf_sought - buf_cur);
      else if (which == std::ios_base::out) pbump(buf_sought - buf_cur);
      return pos_of_buffer_end_in_py_file + (buf_sought - buf_end);
    }

  public:

    class istream : public std::istream
    {
      public:
        istream(streambuf& buf) : std::istream(&buf)
        {
          exceptions(std::ios_base::badbit);
        }

        ~istream() { if (this->good()) this->sync(); }
    };

    class ostream : public std::ostream
    {
      public:
        ostream(streambuf& buf) : std::ostream(&buf)
        {
          exceptions(std::ios_base::badbit);
        }

        ~ostream() { if (this->good()) this->flush(); }
    };
};

    //std::size_t streambuf::default_buffer_size = 1024;

struct streambuf_capsule
{
  streambuf python_streambuf;

  streambuf_capsule(
    bp::object& python_file_obj,
    std::size_t buffer_size=0)
  :
    python_streambuf(python_file_obj, buffer_size)
  {}
};

struct ostream : private streambuf_capsule, streambuf::ostream
{
  ostream(
    bp::object& python_file_obj,
    std::size_t buffer_size=0)
  :
    streambuf_capsule(python_file_obj, buffer_size),
    streambuf::ostream(python_streambuf)
  {}

  ~ostream()
  {
    try {
      if (this->good()) this->flush();
    }
    catch (bp::error_already_set&) {
      PyErr_Clear();
      throw std::runtime_error(
        "Problem closing python ostream.\n"
        "  Known limitation: the error is unrecoverable. Sorry.\n"
        "  Suggestion for programmer: add ostream.flush() before"
        " returning.");
    }
  }
};

}} // boost_adaptbx::python

#endif // GUARD
