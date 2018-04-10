//
//  Copyright (C) 2015,2016 Greg Landrum and NextMove Software
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_SEQUENCEPARSE_H_
#define _RD_SEQUENCEPARSE_H_
#include <string>

namespace RDKit {
class RWMol;

// \brief construct a molecule from a sequence string (currently only supports
// peptides)
/*!
 *   \param seq      - the string to be processed
 *   \param sanitize - toggles sanitization and stereochemistry perception of
 *the molecule
 *   \param lowerD   - if set, lower case letters will be parsed as the d form
 *of the corresponding amino acid
 *
 */
RWMol *SequenceToMol(const char *seq, bool sanitize, bool lowerD);
//! \overload
RWMol *SequenceToMol(const std::string &seq, bool sanitize, bool lowerD);

// \brief construct a protein, RNA or DNA molecule from a sequence string
/*!
 *   \param seq      - the string to be processed
 *   \param sanitize - toggles sanitization and stereochemistry perception of
 *the molecule
 *   \param flavor   -
 *      0 Protein, L amino acids (default)
 *      1 Protein, D amino acids
 *      2 RNA, no cap
 *      3 RNA, 5' cap
 *      4 RNA, 3' cap
 *      5 RNA, both caps
 *      6 DNA, no cap
 *      7 DNA, 5' cap
 *      8 DNA, 3' cap
 *      9 DNA, both caps
 *
 */
RWMol *SequenceToMol(const char *seq, bool sanitize = true, int flavor = 0);
//! \overload
RWMol *SequenceToMol(const std::string &seq, bool sanitize = true,
                     int flavor = 0);

// \brief construct a molecule from a FASTA string (currently only supports
// peptides)
/*!
 *   \param seq      - the string to be processed
 *   \param sanitize - toggles sanitization and stereochemistry perception of
 *the molecule
 *   \param lowerD   - if set, lower case letters will be parsed as the d form
 *of the corresponding amino acid
 *
 */
RWMol *FASTAToMol(const char *seq, bool sanitize, bool lowerD);
//! \overload
RWMol *FASTAToMol(const std::string &seq, bool sanitize, bool lowerD);

// \brief construct a protein, DNA or RNA molecule from a FASTA string
/*!
 *   \param seq      - the string to be processed
 *   \param sanitize - toggles sanitization and stereochemistry perception of
 *the molecule
 *   \param flavor   -
 *      0 Protein, L amino acids (default)
 *      1 Protein, D amino acids
 *      2 RNA, no cap
 *      3 RNA, 5' cap
 *      4 RNA, 3' cap
 *      5 RNA, both caps
 *      6 DNA, no cap
 *      7 DNA, 5' cap
 *      8 DNA, 3' cap
 *      9 DNA, both caps
 *
 */
RWMol *FASTAToMol(const char *seq, bool sanitize = true, int flavor = 0);
//! \overload
RWMol *FASTAToMol(const std::string &seq, bool sanitize = true, int flavor = 0);

// \brief construct a molecule from a HELM string (currently only supports
// peptides)
/*!
 *   \param seq      - the string to be processed
 *   \param sanitize - toggles sanitization and stereochemistry perception of
 *the molecule
 *
 */
RWMol *HELMToMol(const char *helm, bool sanitize = true);
//! \overload
RWMol *HELMToMol(const std::string &helm, bool sanitize = true);
}

#endif
