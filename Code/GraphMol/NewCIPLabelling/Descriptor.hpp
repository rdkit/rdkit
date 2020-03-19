//
//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once

namespace RDKit {
namespace NewCIPLabelling {

static const std::string CIP_LABEL_KEY = "cip.label";
static const std::string CONF_INDEX = "conf.index";

/**
 * Defines a descriptor which can be assigned to an atom to indicate the type of
 * chirality (if there is any). Each descriptor defines it's general @{link
 * Type} which can be useful when comparing centres of different geometry.
 *
 */
enum class Descriptor {
  NONE,  // Unspecified
  UNKNOWN,
  ns,  // Other
  /**
   * Tetrahedral
   */
  R,
  S,
  r,
  s,
  /**
   * Cis/Trans
   */
  seqTrans,
  seqCis,
  E,
  Z,
  /* Axial */
  M,
  P,
  m,
  p,

  SP_4,
  TBPY_5,
  OC_6
};

static std::string to_string(const Descriptor& desc) {
  switch (desc) {
    case Descriptor::NONE:
      return "NONE";
    case Descriptor::UNKNOWN:
      return "UNKNOWN";
    case Descriptor::ns:
      return "ns";
    case Descriptor::R:
      return "R";
    case Descriptor::S:
      return "S";
    case Descriptor::r:
      return "r";
    case Descriptor::s:
      return "s";
    case Descriptor::seqTrans:
      return "seqTrans";
    case Descriptor::seqCis:
      return "seqCis";
    case Descriptor::E:
      return "E";
    case Descriptor::Z:
      return "Z";
    case Descriptor::M:
      return "M";
    case Descriptor::P:
      return "P";
    case Descriptor::m:
      return "m";
    case Descriptor::p:
      return "p";
    case Descriptor::SP_4:
      return "SP_4";
    case Descriptor::TBPY_5:
      return "TBPY_5";
    case Descriptor::OC_6:
      return "OC_6";
  }
  throw std::runtime_error("Unknown descriptor");
}

}  // namespace NewCIPLabelling
}  // namespace RDKit