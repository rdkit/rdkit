//
//
//  Copyright (C) 2023 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MolWriters.h"

#include <fstream>
#include <memory>
#include <string>

#include <maeparser/Writer.hpp>

#include <RDGeneral/BadFileException.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/RDLog.h>

using namespace schrodinger;

namespace RDKit {

MaeWriter::MaeWriter(const std::string &fileName) {
  auto *tmpStream = new std::ofstream(fileName.c_str());
  if (!(*tmpStream) || (tmpStream->bad())) {
    delete tmpStream;
    std::ostringstream errout;
    errout << "Bad output file " << fileName;
    throw BadFileException(errout.str());
  }
  dp_ostream.reset(static_cast<std::ostream*>(tmpStream));
}

MaeWriter::MaeWriter(std::ostream* outStream) : dp_ostream{outStream} {
  PRECONDITION(outStream, "null stream");
  if (outStream->bad()) {
    throw FileParseException("Bad output stream");
  }
}

MaeWriter::MaeWriter(std::shared_ptr<std::ostream> outStream)
    : dp_ostream{std::move(outStream)} {
  PRECONDITION(outStream, "null stream");
  if (outStream->bad()) {
    throw FileParseException("Bad output stream");
  }
}

void MaeWriter::open() { dp_writer.reset(new mae::Writer(dp_ostream)); }

void MaeWriter::setProps(const STR_VECT& propNames) {
  if (d_molid > 0) {
    BOOST_LOG(rdWarningLog) << "WARNING: Setting property list after a few "
                               "molecules have been written\n";
  }
  d_props = propNames;
}

void MaeWriter::flush() {
  PRECONDITION(dp_ostream, "no output stream");
  try {
    dp_ostream->flush();
  } catch (...) {
    try {
      if (dp_ostream->good()) {
        dp_ostream->setstate(std::ios::badbit);
      }
    } catch (const std::runtime_error&) {
    }
  }
}

void MaeWriter::close() {
  if (dp_writer) {
    dp_writer.reset();
  }
  if (dp_ostream) {
    flush();
  }
  dp_ostream.reset();
}

void MaeWriter::write(const ROMol &mol, int confId) {
  (void)mol;
  (void)confId;
}

}  // namespace RDKit