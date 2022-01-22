
/*
 *  Copyright (c) 2019, Novartis Institutes for BioMedical Research Inc.
 *  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *     * Neither the name of Novartis Institutes for BioMedical Research Inc.
 *       nor the names of its contributors may be used to endorse or promote
 *       products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


%include "std_string.i"


%{
#include <../RDStreams/streams.h>
#include <sstream> 
%}

// We just need to pass the pointers around here
%nodefaultctor istream;
%nodefaultdtor istream;
%nodefaultctor filtering_istream;
%nodefaultdtor filtering_istream;

namespace std {
class istream;
}

#ifdef RDK_USE_BOOST_IOSTREAMS
%extend RDKit::gzstream {
    std::istream* _GetStream() { return (std::istream*)$self; }
    std::string Dump() {
      std::ostringstream stream;
      std::copy(std::istreambuf_iterator<char>(*$self),
                std::istreambuf_iterator<char>(),
                std::ostreambuf_iterator<char>(stream));

      return stream.str();
    }
}
%typemap(javacode) RDKit::gzstream %{
  private SWIGTYPE_p_std__istream streamRef;
  public SWIGTYPE_p_std__istream GetStream() {
     if (streamRef == null)
       streamRef = _GetStream();
     return streamRef;
  }
%}
#endif

%include <../RDStreams/streams.h>


