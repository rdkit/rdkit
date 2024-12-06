//
//  Copyright (c) 2024 Glysade Inc and other RDkit contributors
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <CDXMLParser.h>
#include <CDXStdObjects.h>

// #define DEBUG 1
namespace {
void Recurse(CDXObject &object, int depth) {
  for(auto child : object.ContainedObjects()) {
    for(int i=0; i<depth; ++i)
      std::cerr << ".";
    std::cerr << "id:" << child.second->GetObjectID() << " tag: " << (CDXDatumID) child.second->GetTag() << std::endl;
    
    if(child.second->GetAttributes()) {
      for(auto attr : *child.second->GetAttributes()) {
        for(int i=0; i<depth; ++i)
          std::cerr << ">";
        std::cerr << (CDXDatumID) attr.second->GetTag() << std::endl;
      }
    }
    Recurse(*child.second, depth+1);
    }
    
  }
}
namespace RDKit {
void ChemDrawToMols(const std::string &filename) {
  CDXMLParser parser;
  const bool HaveAllXml = true;
  std::fstream chemdrawfile(filename);
  if (!chemdrawfile.is_open() ) {
    // error
    return;
  }
  std::stringstream stream;
  stream << chemdrawfile.rdbuf();
  std::string data(stream.str());
  parser.XML_Parse(data.c_str(), static_cast<int>(data.size()), HaveAllXml);
  std::unique_ptr<CDXDocument> document = parser.ReleaseDocument();
  if(!document) {
    // error
    return;
  }
  Recurse(*document, 0);
}
}
