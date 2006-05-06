# $Id$
#
# Copyright (C) 2003-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" Python functions for path-based (like Daylight) molecular fingerprinting


Schematic Algorithm:

  1) find all its paths within a particular range of distances

  2) generate a hash for each path

  3) use the hash to seed a random number generator (RNG)

  4) generate a series of values from the RNG and set the
     corresponding bits in the fingerprint.

During the actual computation, we're not really so stupid as to store
all the paths.


 **Notes**

 - The algorithm for setting bits is deterministic (one would
   certainly hope so), but irreversible.  There is no way of backing
   out information about structure from which bits are set by a
   molecule.  Which leads to the next point:

 - Like the Daylight algorithm, a path will always set the same set of
   bits, so we have the guarantee that if molecule M1 is a subgraph of
   molecule M2, all the bits set by M1 will also be set by M2
   (provided the same parameters are used to generate the
   fingerprints, of course).

 - This code is here only as a prototype, it really needs to all be in
   C++

"""
raise  ImportError,"This functionality has been moved into Chem.rdchem"
