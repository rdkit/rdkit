//
// Copyright (C) David Cosgrove 2023
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Options for Rascal Clustering.  In general, the option names and defaults
// are taken from the paper:
// 'A Line Graph Algorithm for Clustering Chemical Structures Based
// on Common Substructural Cores', JW Raymond, PW Willett.
// https://match.pmf.kg.ac.rs/electronic_versions/Match48/match48_197-207.pdf
// https://eprints.whiterose.ac.uk/77598/

#include <RDGeneral/export.h>

#ifndef RASCALCLUSTEROPTIONS_H
#define RASCALCLUSTEROPTIONS_H

namespace RDKit {
namespace RascalMCES {

struct RDKIT_RASCALMCES_EXPORT RascalClusterOptions {
  double similarityCutoff = 0.7; /* Similarity cutoff for clustering.  Initial
                                    clusters will have molecule pairs of at
                                    least this similarity. */
  double a = 0.05; /* penalty score for each unconnected component in MCES */
  double b = 2.0;  /* weight of matched bonds over matched atoms */
  unsigned int minFragSize =
      3; /* minimum number of atoms in a fragment for it to
            be included in the MCES.  Also p in the paper. */
  double minIntraClusterSim = 0.9; /* two pairs of molecules are included in the
                                      same cluster if the similarity between
                                      their MCESs is greater than this. S_a
                                      in the paper */
  double clusterMergeSim = 0.6;    /* two clusters are merged if fraction of
                                      molecules they have in common is greater than
                                      this. S_b in the paper */
  unsigned int maxNumFrags = 2; /* The maximum number of fragments in any MCES.
                                   Otherwise the MCES can be a lot of small
                                   fragments scattered across the molecule - it
                                   tries too hard to find a match, sometimes */
  int numThreads = -1; /* The number of threads to use.  If > 0, will use that
                          number.  If <= 0, will use the number of hardware
                          threads plus this number.  So if the number of
                          hardware threads is 8, and numThreads is -1, it will
                          use 7 threads. */
};
}  // namespace RascalMCES
}  // namespace RDKit
#endif  // RASCALCLUSTEROPTIONS_H
