# Release_2019.03.1
(Changes relative to Release_2018.09.1)

## REALLY IMPORTANT ANNOUNCEMENT
- As of this realease (2019.03.1) the RDKit no longer supports Python 2. Please read
this rdkit-discuss post to learn what your options are if you need to keep
using Python 2: https://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg08354.html

## Backwards incompatible changes
- The fix for github #2245 means that the default behavior of the MaxMinPicker is now truly random. 
  If you would like to reproduce the previous behavior, provide a seed value of 42.


# Release_2018.09.1
(Changes relative to Release_2018.03.1)

## Deprecations
- As part of the changes and cleaning up done for #1836 many of the `#defines`
  used in the RDKit have been renamed.
    - `USE_BUILTIN_POPCOUNT` -> `RDK_OPTIMIZE_NATIVE`
    - `RDK_THREADSAFE_SSS` -> `RDK_BUILD_THREADSAFE_SSS`
    - `BUILD_COORDGEN_SUPPORT` -> `RDK_BUILD_COORDGEN_SUPPORT`
    - `BUILD_AVALON_SUPPORT` -> `RDK_BUILD_AVALON_SUPPORT`
    - `BUILD_INCHI_SUPPORT` -> `RDK_BUILD_INCHI_SUPPORT`
    - `BUILD_SLN_SUPPORT` -> `RDK_BUILD_SLN_SUPPORT`
    - `RDK_CAIRO_BUILD` -> `RDK_BUILD_CAIRO_SUPPORT`

## Documentation updates
We have moved to using Sphinx's autodoc to create the Python API documentation.
`epydoc`, the tool we used to use, is no longer actively developed and only supports
Python2. There will undoubtedly be problems associated with the change; if you notice
anything missing in the documetation or something that's really badly formatted,
please either let us know or submit a PR with a fix.

## Backwards incompatible changes
This release includes a set of changes to make the default arguments to common
functions less error prone (github #1679).
- GetAtomSmiles() now generates isomeric SMILES by default.
- The ringMatchesRingOnly option to the FindMCS() function now applies to
  atom-atom matches as well as bond-bond matches.
- The Python functions EmbedMolecule() and EmbedMultipleConfs() now use the
  ETKDG algorithm by default instead of standard distance geometry.

## Highlights:
- This release includes two contributions from the Google Summer of Code:
   - A new generalized fingerprint generator.
   - An integration/port of MolVS to the core RDKit.
  The API on both of these may change a bit with future releases.
- The rdkit.Chem.Draw module now includes functions for depicting fingerprint
  bits. Morgan and RDKit fingerprints are currently supported.

## Acknowledgements:
Boran Adas, Francois Berenger, Thomas Blaschke,  Brian Cole, Andrew Dalke, Guillaume Godin,
Brice Hoff, Brian Kelley, Karl Leswing, Susan Leung, Pat Lorton, Josh Meyers, Hirotomo Moriwaki,
Dan Nealschneider, Noel O'Boyle, Pavel Raiskup, Sereina Riniker, Ricardo Rodriguez Schmidt,
Stephen Roughley, Roger Sayle, Takayuki Serizawa, Rim Shayakhmetov, Gregory Simm, Jon Sorenson,
Matt Swain, Kiran Telukunta, Paulo Tosco, Alain Vaucher, Maciej Wójcikowski, '0xDECAFC0FFEE',
'jaechanglim', 'paconius'

## Contrib updates:
- The FastCluster code has been updated.

## New Features and Enhancements:
  - expose MolChemicalFeature.{get,set}ActiveConformer()  to python
 (github issue #1636 from greglandrum)
  - added Scripts/create_deb_packages.sh
 (github pull #1655 from UnixJunkie)
  - Start to use Catch2 for new tests
 (github pull #1732 from greglandrum)
  - Switch DbCLI scripts from optparse to argparse
 (github issue #1778 from greglandrum)
  - Add EEM partial charges
 (github pull #1828 from greglandrum)
  - Add header file providing access to RDKit compile time flags
 (github issue #1836 from greglandrum)
  - add control over the coordinate precision to coordgen
 (github pull #1847 from greglandrum)
  - Add Properties interface to ChemicalReactions
 (github pull #1848 from greglandrum)
  - Switch Python API documentation generation over to Sphinx autodoc
 (github pull #1849 from greglandrum)
  - expose MolOps::assignStereochemistryFrom3D() to Python
 (github issue #1850 from greglandrum)
  - bivariate_normal deprecation in mlab.py of matplotlib
 (github issue #1851 from telukir)
  - Expose minAtomRingSize() and minBondRingSize() to Python wrappers
 (github pull #1859 from mwojcikowski)
  - enable building DLLs on Windows
 (github pull #1861 from ptosco)
  - Fix compatibility with Boost 1.67+
 (github pull #1864 from mcs07)
  - Don't manually set RDConfig paths in conda env
 (github pull #1865 from mcs07)
  - Make svg xmlns prefix use more consistent
 (github pull #1866 from mcs07)
  - Add custom 3D Descriptors
 (github pull #1867 from greglandrum)
  - Add initial Maestro format Supplier using maeparser library
 (github pull #1872 from lorton)
  - add queryAtomNonHydrogenDegree() query operation
 (github issue #1873 from greglandrum)
  - Add an auto-populated file with cmake config options
 (github pull #1874 from greglandrum)
  - Custom property VSA
 (github pull #1884 from sriniker)
  - Swap maeparser and coordgen dependencies to use releases
 (github issue #1887 from greglandrum)
  - expose MolDraw2DSVG.tagAtoms() to python
 (github pull #1897 from greglandrum)
  - allow the cleanup step of Mol2 parsing to be disabled
 (github pull #1898 from greglandrum)
  - Allow Atom.GetAtomSmarts() to return isomeric SMILES
 (github pull #1902 from greglandrum)
   - Allow coordgen and maeparser to be built as static libraries
 (github pull #1909 from ptosco)
  - Support reaction_to_svg() in cartridge
 (github issue #1916 from greglandrum)
  - Addresses several minor warning messages during the build
 (github pull #1935 from d-b-w)
  - Some optimization of the queries constructed from SMARTS
 (github pull #1940 from greglandrum)
  - Add ring/chain match constraints options to AdjustQueryParameters()
 (github issue #1943 from greglandrum)
  - roc calculation naming problem
 (github pull #1975 from 0xDECAFC0FFEE)
  - Fingerprinting functions should call assignStereochemistry() when necessary
 (github issue #1993 from bricehoff)
  - Dev/GSOC2018_MolVS_Integration
 (github pull #2002 from susanhleung)
  - GSoC 2018 - Fingerprints
 (github pull #2005 from Boranadas)
  - port fingerprint bit rendering code from CheTo
 (github pull #2008 from greglandrum)
  - PgSQL: add support for PostgreSQL_CONFIG cmake var
 (github pull #2014 from praiskup)
  - Add missing boost header for v1.64
 (github pull #2016 from gncs)
  - Enhanced stereo read/write support in SDF files.
 (github pull #2022 from d-b-w)
  - IPythonConsole hooks should copy the original docstring
 (github issue #2025 from adalke)
  - Allow dumping interchange information into SVG files
 (github pull #2030 from greglandrum)
  - MCS: add test for ring--non-ring matches at the atom level
 (github issue #2034 from greglandrum)
  - Ability to generate a list of possible smiles representation for a given molecule
 (github issue #2042 from thegodone)
  - add scoring test (relevant to #1975)
 (github pull #2064 from greglandrum)
  - actually run the SmilesWriterNoNames() test
 (github pull #2067 from greglandrum)
  - Add a default for maximum products generated by a reaction (maxProduc…
 (github pull #2069 from bp-kelley)
  - Add user-defined literals for parsing SMILES and SMARTS
 (github pull #2070 from greglandrum)
  - move rdInfoLog to stderr
 (github pull #2073 from greglandrum)
  - add confId argument to MolChemicalFeatureFactor::getFeaturesForMol()
 (github issue #2077 from greglandrum)
  - Expose a CMake flag to build RDKit with -rpath
 (github pull #2084 from coleb)
  - Dev/expose setquery to python
 (github pull #2088 from bp-kelley)
  - Updated .gitignore with files generated outside of build directory.
 (github pull #2095 from ricrogz)
  - Address compile warnings & trivial improvements
 (github pull #2097 from ricrogz)
  - Coordgen: stop printing Templates location
 (github pull #2102 from greglandrum)
  - Update Docs For CalcBEDROC
 (github pull #2103 from lilleswing)

## Bug Fixes:
  - Cannot find rings for ChEBI 50252
 (github issue #299 from greglandrum)
  - Feature request: #defines to test RDKit version
 (github issue #1454 from baoilleach)
  - Atoms in residue read from pdb lose their AtomPDBResidueInfo after reaction (Bug)
 (github issue #1632 from hjuinj)
  - SMARTS parse failure for some queries involving Hs
 (github issue #1719 from greglandrum)
  - Conformer indexing bug in RDFreeSASA.cpp?
 (github issue #1786 from paconius)
  - allow libs to not be installed
 (github pull #1832 from greglandrum)
  - RWMol::addAtom(Atom,bool) missing from Java wrappers
 (github issue #1837 from greglandrum)
  - RWMol::clear now calls ROMol::initMol
 (github pull #1844 from bp-kelley)
  - Allow types.h to be included in applications that use /D_USE_MATH_DEFINES
 (github pull #1846 from d-b-w)
  - Fixes failing Python tests on Windows
 (github pull #1855 from ptosco)
  - Allow building on cygwin using -std=c++11
 (github pull #1856 from greglandrum)
  - Stop using the cmake Boost:: targets
 (github pull #1858 from greglandrum)
  - R-group Decomposition: allow H replacements when matchOnlyAtRgroups is set
 (github pull #1871 from bp-kelley)
  - Mark cartridge functions as being parallel safe
 (github issue #1886 from greglandrum)
  - Fixes locale handling on Windows
 (github pull #1892 from ptosco)
  - get the pandas tests working with pandas 0.23
 (github pull #1896 from greglandrum)
  - "make install" appears to miss RDBoost/export.h
 (github issue #1903 from baoilleach)
  - Fix curl fallback for downloading files
 (github pull #1904 from d-b-w)
  - Bond stereo information not output to SMARTS
 (github issue #1906 from greglandrum)
  - Library .so names missing RDKit?
 (github issue #1913 from baoilleach)
  - Negated atom number queries in SMARTS should not set atomic number of query atom
 (github issue #1920 from greglandrum)
  - memory leak in Get3DDistanceMatrix
 (github issue #1924 from jaechanglim)
  - Atom with bond to itself is accepted by the SMILES parser.
 (github issue #1925 from tblaschke)
  - Possibly incorrect aromatic SMILES generated for structure
 (github issue #1928 from baoilleach)
  - Using the coordgen library seems to cause a seg fault
 (github issue #1929 from JoshuaMeyers)
  - Aromaticity failure in 7-ring with charged radical carbon
 (github issue #1936 from bp-kelley)
  - Fix embarassing bug, check the counter each iteration
 (github pull #1939 from bp-kelley)
  - RuntimeError when importing rdkit.Chem.Descriptors with Python 3.7
 (github issue #1948 from drkeoni)
  - Query features in products of rxn files not properly handled
 (github issue #1950 from greglandrum)
  - ReactionToSmarts broken for multi-component templates
 (github issue #1955 from bp-kelley)
  - update knime urls in overview.md
 (github pull #1966 from greglandrum)
  - CXN extended SMILES labels are not applied to the correct atom in SMILES with explicit H
 (github issue #1968 from baoilleach)
  - MolFromSmarts MolToSmarts fails to roundtrip on patterns with chirality
 (github issue #1985 from bp-kelley)
  - QueryAtoms with atom list queries should not have the atomic number set
 (github issue #1988 from greglandrum)
  - RemoveHs not properly updating double bond stereoatoms
 (github issue #1990 from shayakhmetov)
  - Error while parsing empty atom list in Mol files.
 (github issue #2000 from drkeoni)
  - Cleanup step of sanitization sometimes sets undesired formal charges
 (github issue #2020 from avaucher)
  - GetBondSmiles() returns nothing for wedged bonds when allBondsExplicit is true
 (github issue #2029 from greglandrum)
  - PMIs and NPRs return same value between different conformer
 (github issue #2037 from philopon)
  - Failure to parse V3K mol file with bonds to multi-center linkage points
 (github issue #2040 from greglandrum)
  - patch a memory allocation problem in the maeparser v1.0.0
 (github pull #2044 from greglandrum)
  - CIPRank values from JSONDataToMols are not unsigned
 (github issue #2046 from greglandrum)
  - switch to v1.0.1 of the maeparser
 (github pull #2048 from greglandrum)
  - Update fastcluster code
 (github pull #2051 from greglandrum)
  - Fix memory leak in Dict operator=
 (github pull #2061 from bp-kelley)
  - Atom mapping lost after call to MergeQueryHs()
 (github issue #2062 from greglandrum)
  - Drawing racemic bond stereo as crossed bonds should be the default
 (github issue #2063 from coleb)
  - Moved test.h from RDBoost to RDGeneral for consistency with export.h
 (github pull #2074 from ptosco)
  - CTABs behave differently than before regarding stereo
 (github issue #2082 from bp-kelley)
  - BitInfo not complete for RDKFingerprint
 (github issue #2115 from greglandrum)

## Removed code:
- Remove the deprecated MolDrawing code
 (github pull #2111 from greglandrum)


# Release_2018.03.1
(Changes relative to Release_2017.09.1)

## C++11 notes

Starting with this release, the RDKit core C++ code is written in modern C++;
for this release that means C++11. This means that the compilers used to build
it cannot be completely ancient. Here are the minimum tested versions:
- g++ v4.8: though note that the SLN parser code cannot be built with v4.8. It
  will automatically be disabled when this older compiler is used.
- clang v3.9: it may be that older versions of the compiler also work, but we
  haven't tested them.
- Visual Studio 2015: it may be that older versions of the compiler also work,
  but we haven't tested them.

## Backwards incompatible changes

This release includes a set of changes to make the default arguments to common
functions less error prone (github #1679).
- MolToSmiles() now generates isomeric SMILES by default.
- The embedding code now uses the ETKDG method by default.
- MolToMolBlock() will now by default generate a set of 2D coordinates for
  molecules when the includeStereo option is set to True. The changes are made
  to a copy of the molecule; the molecule itself will not be modified.
- The Mol file (and SDF) parser now determines atomic stereochemisty based on
  the 3D coordinates provided (if 3D coordinates are provided).
- The SMILES parser now supports CXSMILES by default (assuming that additional
  text that looks like CXSMILES extensions is there).

In every case the old behavior can be obtained by providing an optional argument
to the function(s) mentioned.

## Acknowledgements:
Boran Adas, José Emilio Sánchez Aparicio, Patrick Avery, Jason Biggs, Brian
Cole, Andrew Dalke, JW Feng, Peter Gedeck, Guillaume Godin, Richard Hall, Thomas
Heavy, Gareth Jones, Brian Kelley, Karl Leswing, Susan Leung, Chris Morris, Dan
Nealschneider, Noel O'Boyle, Axel Pahl, Pavel Polishchuk, Sereina Riniker, Jeff
van Santen, Roger Sayle, Martin Šícho, Matt Swain, Paolo Tosco, Sam Webb, Maciej
Wójcikowski, Nicola Zonta, 'clinntt', 'hjuinj', 'iwatobipen',

## Highlights:
  - An initial version of an integration with Schrodinger's coordgen library is
    included. This produces much better 2D coordinates for complex molecules.
  - Thanks to the move to modern C++ the RDKit is now faster and uses less
    memory
  - A number of improvements were made to the PDB reader
  - v2 of the ETKDG torsions and potentials is now available

## Contrib updates:
  - Implementation of Peter Ertl's IFG method for identifying functional groups
    from Guillaume Godin and Richard Hall
  (github pull #1813 from thegodone)

## New Features and Enhancements:
  - Support InChi 1.05
 (github issue #1533 from greglandrum)
  - Update CPack to create .deb files correctly
 (github pull #1580 from psavery)
  - Initial commit of EnumerateHeterocycles
 (github pull #1588 from coleb)
  - Version 2 of ETKDG
 (github pull #1597 from sriniker)
  - GetMolFrags now optionally returns atom indices along with mols
 (github pull #1602 from ptosco)
  - NP Likeness with confidence value
 (github pull #1608 from apahl)
  - Adding an option to EnumerateStereoisomers to only return unique isomers
 (github pull #1612 from coleb)
  - Add function wedgeBond()
  (github issue #1615 from greglandrum)  
  - Dev/substructlibrary docs
 (github pull #1620 from bp-kelley)
  - Turns off exception throwing for certain classes Rlabel sanitization.
 (github pull #1621 from bp-kelley)
  - Add an "MDL" aromaticity model
 (github issue #1622 from hjuinj)
  - Add support for %(NNN) notation for ring closures
 (github pull #1624 from baoilleach)
  - Enable windows build that uses cairo
 (github pull #1628 from greglandrum)
  - [MRG] Fix PDB reader + add argument to toggle proximity bonding
 (github pull #1629 from mwojcikowski)
  - Improve AddHs for molecules read from PDB
 (github pull #1647 from mwojcikowski)
  - Improved regression test for ETKDG version 2
 (github pull #1640 from sriniker)
  - RDKit interpretation of atom stereo SMILES is different from 4 other toolkits
 (github issue #1652 from coleb)
  - Treat bonds in PDB CONECT records explicitly, but make blacklisted ones zero-order.
 (github pull #1658 from mwojcikowski)
  - There is no need to enforce that (i, j) and (k, l) be bonded when setting a i, j, k, l dihedral
 (github pull #1673 from ptosco)
  - Make default arguments to common functions less error prone
 (github issue #1679 from greglandrum)
  - Add Fast cluster script
 (github pull #1683 from iwatobipen)
  - Update embedded InChI to v1.05
 (github pull #1684 from mcs07)
  - Add `AllChem.MMFFGetMoleculeForceField().CalcGradient()` to Python wrappers
 (github issue #1688 from theavey)
  - Play nice with naughty MOL blocks
 (github issue #1689 from jw-feng)
  - Make the defaults for some functions less error prone.
 (github pull #1690 from greglandrum)
  - implemented Python wrappers for computing PMI axes and moments
 (github pull #1700 from ptosco)
  - Enable range-based for loops for molecules
 (github pull #1701 from bp-kelley)
  - Support some cactvs extensions to SMARTS
 (github pull #1704 from greglandrum)
  - Integrate Coordgen
 (github pull #1708 from greglandrum)
  - Removes ATOM/BOND_SPTR in boost::graph in favor of raw pointers
 (github pull #1713 from greglandrum)
  - Set atomic properties from SMARTS
 (github pull #1716 from greglandrum)
  - Allow installation of Python tests to facilitate testing installations
 (github pull #1724 from greglandrum)
  - setAromaticity() should work even if there are aromatic atoms present
 (github issue #1730 from greglandrum)
  - Use uint32 atom and bond indices
 (github pull #1742 from greglandrum)
  - Switch from boost::thread to std::thread
 (github pull #1745 from greglandrum)
  - switch to using std::regex in the SLN parser
 (github pull #1746 from greglandrum)
  - replace the usage of rdk_auto_ptr with std::unique_ptr
 (github pull #1752 from greglandrum)
  - getMolBoundsMatrix() should do triangle bound smoothing by default
 (github issue #1763 from greglandrum)
  - Added Morgan feature fingerprints to Java API
 (github pull #1764 from jones-gareth)
  - Reaction fingerprints not exposed in Java wrapper
 (github issue #1776 from webbres)
  - add Tversky index calculation for shapes
 (github pull #1777 from susanhleung)
  - Add MolToInchiKey function()
 (github pull #1784 from greglandrum)
  - speedup the NumBitsInCommon operation
 (github pull #1785 from greglandrum)
  - Stop putting brackets around * atoms in SMILES
 (github pull #1788 from greglandrum)
  - Support for a JSON-based molecule interchange format
 (github pull #1798 from greglandrum)

## Bug Fixes:
  - Fixes Java wrapper build error with Boost 1.64
 (github pull #1613 from ptosco)
  - AssignStereochemistry cleanIt=True incorrectly removing new CIS/TRANS bond stereo
 (github issue #1614 from coleb)
  - switch to using a specific freesasa version
 (github pull #1619 from greglandrum)
  - Add support for %(NNN) notation for ring closures
 (github pull #1624 from baoilleach)
  - Draw._moltoSVG() raises an exception
 (github issue #1625 from greglandrum)
  - MolDrawCairo2D does not build on windows
 (github issue #1627 from greglandrum)
  - Enable windows build that uses cairo
 (github pull #1628 from greglandrum)
  - don't always download the FreeSASA source
 (github issue #1630 from greglandrum)
  - Make sure EmbedMultipleConfs is deterministic for very large seeds and a seed of 0
 (github pull #1635 from coleb)
  - from rdkit.Chem import AllChem has grown a six dependency
 (github issue #1637 from bp-kelley)
  - Fixing bug in IPythonConsole SVG rendering introduced in 1027d4469545653180fff9a38dc8224bd50e8b0d
 (github pull #1641 from coleb)
  - changes required to allow replacing the obsolete __conda_version__ in conda-rdkit
 (github pull #1644 from ptosco)
  - GetConformerRMSMatrix does not work if some conformers were removed
 (github issue #1650 from DrrDom)
  - EnumerateLibrary with initFromString called twice doesn't clear the reaction
 (github issue #1657 from bp-kelley)
  - Missed symmetrization in R-Group decomposition
 (github issue #1659 from greglandrum)
  - Use numpy not numeric for boost 1.65+ - fixes #1581
 (github pull #1664 from mcs07)
  - Support valence 7 for As, Sb, and Bi
 (github issue #1668 from greglandrum)
  - Fix: GetDonor2FeatVects heavy atoms confusion
 (github pull #1676 from josan82)
  - Acetylenic hydrogens not given appropriate 2D coordinates
 (github issue #1691 from jasondbiggs)
  - Warning on import of rgroup decomposition package
 (github issue #1695 from greglandrum)
  - AUTOCORR2D.h not installed unless RDK_BUILD_DESCRIPTORS3D but is required
 (github issue #1702 from baoilleach)
  - Dative bonds interfere with kekulization and the perception of aromaticity
 (github issue #1703 from greglandrum)
  - Fix/rgroup prefer matching nonhs over hs
 (github pull #1707 from bp-kelley)
  - bonds that are STEREOCIS or STEREOTRANS cannot be depickled
 (github issue #1710 from greglandrum)
  - Get queries from the new cactvs SMARTS extensions to pickle correctly
 (github pull #1712 from greglandrum)
  - fix an irritating cmake problem
 (github pull #1715 from greglandrum)
  - Added dependency from Boost headers to PgSQL CMakeLists.txt
 (github pull #1717 from ptosco)
  - Updates python test runner to always use sys.executable
 (github pull #1721 from bp-kelley)
  - - make bond stereo detection in rings consistent
 (github pull #1727 from ptosco)
  - xlocale.h not needed to compile with clang
 (github issue #1728 from adalke)
 - BreakBRICSBonds() not preserving stereochemistry
 (github issue #1734 from greglandrum)
  - rdmolfiles.CanonicalRankAtoms segfaults on 0 atom molecules
 (github issue #1735 from lilleswing)
  - deprecated apply() function causes GetRDKFingerprint to fail in Python 3
 (github issue #1747 from clinntt)
  - Stop dereferencing end() iterators
 (github pull #1748 from greglandrum)
  - out of range fromAtom causes GetMorganFingerprintAsBitVect to segfault
 (github issue #1749 from adalke)
  - Generated SMARTS does not contain atomic chiral tags
 (github issue #1756 from greglandrum)
  - make the build work even if boost::serialization is disabled
 (github pull #1767 from greglandrum)
  - Fix typo in GetBoolProp documentation
 (github pull #1770 from jvansan)
  - Fingerprint segfaults with branchedPaths=False and useHs=False
 (github issue #1793 from chrishmorris)
  - Fix python linkage (primarily for conda builds)
 (github pull #1808 from greglandrum)
  - removeHs() should not remove H atoms that are contributing to the definition of a stereo bond
 (github pull #1810 from d-b-w)
  - global EmbedParameters objects should not be writeable in SWIG wrappers
 (github issue #1826 from greglandrum)
  - RDKit crashes when MolsToGridImage function is called with an empty iterable.
 (github issue #1829 from martin-sicho)


# Release_2017.09.1
(Changes relative to Release_2017.03.1)

## Important
- The fix for bug #1567 changes the way fragment SMILES are canonicalized.
  MolFragmentToSmiles() and canonicalizeFragment() will now often return
  different results
- The fix for bug #1604 changes the behavior of QueryAtom::setQuery(), which
  now deletes the current query before setting the new value. If you are using
  QueryAtom::setQuery() from C++ (or possibly Java), be sure that you are not
  also deleting that memory.

## Acknowledgements:
Brian Cole, Peter Gedeck, Guillaume Godin, Jan Halborg Jensen, Malitha Kabir,
Tuomo Kalliokoski, Brian Kelley, Noel O'Boyle, Matthew O'Meara, Pavel
Polishchuk, Cameron Pye, Christian Ribeaud, Stephen Roughley, Patrick Savery,
Roger Sayle, Nadine Schneider, Gregor Simm, Matt Swain, Paolo Tosco, Alain
Vaucher, Sam Webb, 'phenethyl', 'xiaotaw'

## Highlights:
- The new R-Group decomposition code provides a flexible and powerful tool for
  building R-group tables or datasets look in $RDBASE/Docs/Notebooks for
  example notebooks showing how to use this.
- Drawing of chemical reactions has been greatly improved and is now done using
  the C++ rendering code.
- The MaxMinPicker is dramatically faster.
- New descriptors: the QED descriptor has been added as have a large collection
  of new 3D descriptors and implementations of the USR and USRCAT fingerprints.

## New Features and Enhancements:
  - Bring back USR and USRCAT descriptors
 (github pull #1417 from greglandrum)
  - Generate a warning for conflicting bond directions
 (github issue #1423 from greglandrum)
  - expose and test GetDrawCoords()
 (github pull #1427 from greglandrum)
  - Improvement suggestions for SaltRemover
 (github issue #1431 from ribeaud)
  - Remove obsolete scripts from Scripts dir
 (github pull #1440 from greglandrum)
  - Support drawing reactions from C++
 (github pull #1444 from greglandrum)
  - QED code with unit test file
 (github pull #1445 from gedeck)
  - Add support for other datatypes to  ConvertToNumpyArray
 (github issue #1447 from pyeguy)
  - - updated FindCairo.cmake
 (github pull #1455 from ptosco)
  - - fixes PgSQL CMakeLists.txt to enable conda build on Windows
 (github pull #1457 from ptosco)
  - Some cleanups to make Travis builds faster
 (github pull #1464 from greglandrum)
  - ExplainPairScore does not support includeChirality=True
 (github issue #1466 from xiaotaw)
  - Add a collection of new 3D descriptors
 (github pull #1467 from greglandrum)
  - Update cartridge documentation to use ChEMBL 23
 (github issue #1491 from greglandrum)
  - First entry of the SubstructLibrary module
 (github pull #1493 from bp-kelley)
  - assorted fixes to get the current master branch to build on Windows
 (github pull #1495 from ptosco)
  - Support assignment of stereochemistry tags to bonds from 3D structure  
 (github issue #1497 from gncs)
  - Support black and white molecule drawing
 (github issue #1510 from greglandrum)
  - Missing def_readwrite for backgroundColour in rdMolDraw2D.cpp
 (github issue #1519 from goraj)
  - Adds canonicalization of atom maps
 (github pull #1521 from bp-kelley)
  - Implement stereoisomer enumeration
 (github pull #1531 from greglandrum)
  - Add a MolBundle class
 (github pull #1537 from greglandrum)
  - Provide support for color palettes in MolDraw2D
 (github pull #1546 from greglandrum)
  - A few reaction drawing tweaks
 (github pull #1549 from greglandrum)
  - R group improvements
 (github pull #1552 from greglandrum)
  - Add a canned Atom query for heavy atom degree
 (github issue #1563 from greglandrum)
  - Adds FreeSASA adapter
 (github pull #1565 from bp-kelley)
  - Added C++ version of getBestRMS()
 (github pull #1568 from psavery)
  - SMILES lexer optimization/enhancement
 (github pull #1575 from greglandrum)
  - Update IPythonConsole and PandasTools to use new drawing code
 (github pull #1577 from greglandrum)
  - Squashes warnings on cygwin
 (github pull #1578 from bp-kelley)
  - Support continuous highlighting in drawMolecules().
 (github pull #1579 from greglandrum)
  - Enhanced Similarity Maps depiction
 (github pull #1594 from gerebtzoff)

## Bug Fixes:
  - RDKit gets stuck on PubChem CID 102128817
 (github issue #1281 from TuomoKalliokoski)
  - MMP code not including molecules with no cuts
 (github issue #1406 from greglandrum)
  - Fixes PandasTools to also work with pandas 0.20
 (github pull #1410 from bp-kelley)
  - csharp input files out of date
 (github issue #1413 from greglandrum)
  - Fix cxsmiles parse on VS2008
 (github pull #1415 from mcs07)
  - MaxMinPicker picking non-existent element
 (github issue #1421 from greglandrum)
  - _isCallable clashes with Celery
 (github issue #1434 from momeara)
  - Impossible to build the RDKit from source without Python installed
 (github issue #1435 from greglandrum)
  - RemoveHs() removes H atom attached to dummy if it came from AddHs()
 (github issue #1439 from DrrDom)
  - fix a couple failing windows tests related to temp file removal
 (github pull #1446 from greglandrum)
  - SanitizeRxn fails with a runtime exception when unused Rlabels are in product
 (github issue #1448 from bp-kelley)
  - String module conversion bug
 (github pull #1452 from coleb)
  - GetConformerRMS() documentation is misleading
 (github pull #1459 from greglandrum)
  - URANGE_CHECK not doing its job in RWMol::addBond
 (github issue #1461 from baoilleach)
  - ExplainPairScore does not support includeChirality=True
 (github issue #1466 from xiaotaw)
  - MolToSmarts does not include atom-map or isotope info for molecules built from SMILES
 (github issue #1472 from greglandrum)
  - AdjustQueryProperties() removing properties from dummy atoms
 (github issue #1474 from greglandrum)
  - Fixes lookup for HELM Monomer 'D'
 (github pull #1477 from bp-kelley)
  - Aromatic rings composed solely of dummy atoms should not be kekulized
 (github issue #1478 from bp-kelley)
  - Directly specify rotor model used in QED.
 (github pull #1483 from bp-kelley)
  - Unicode problem with pidPS tests on Mac
 (github issue #1490 from greglandrum)
  - Pattern fingerprint setting bad bits with degree zero atoms
 (github issue #1496 from greglandrum)
  - Remove xlocale header
 (github pull #1501 from greglandrum)
  - Fixes atom documentation
 (github pull #1505 from bp-kelley)
  - TypeError from PandasTools.SaveXlsxFromFrame
 (github issue #1507 from pyeguy)
  - Removes trailing spaces after \ to fix windows compilation errors
 (github pull #1516 from bp-kelley)
  - prepareMolForDrawing() not in SWIG wrappers
 (github issue #1522 from greglandrum)
  - Bond is missing IsInRing methods in Java wrapper
 (github issue #1535 from sroughley)
  - Fixes blanking of non-query atom data when QueryAtomData was being pi…
 (github pull #1541 from bp-kelley)
  - ChemicalReaction code not calling setNoImplicit() when H counts are set.
 (github issue #1544 from greglandrum)
  -  Fixes failing build with MSVC
 (github pull #1547 from ptosco)
  - Kekulization error with cores from R-Group Decomposition
 (github issue #1550 from greglandrum)
  - Fixes double free for Dict::update
 (github pull #1571 from bp-kelley)
  - QueryAtom::setQuery() should delete the old query first
 (github pull #1604 from greglandrum)


# Release_2017.03.1
(Changes relative to Release_2016.09.1)

## Important
- The fix for bug #879 changes the definition of the layered fingerprint.
  This means that all database columns using layered fingerprints as well as
  all substructure search indices should be rebuilt.
- All C++ library names now start with RDKit (see #1349).

## Acknowledgements:
Brian Cole, David Cosgrove, JW Feng, Berend Huisman, Peter Gedeck, 'i-tub',
Jan Holst Jensen, Brian Kelley, Rich Lewis, Brian Mack, Eloy Felix Manzanares,
Stephen Roughley, Roger Sayle, Nadine Schneider, Gregor Simm, Matt Swain,
Paolo Tosco, Riccardo Vianello, Hsiao Yi

## Highlights:
  - It's now possible (though not the default) to pickle molecule properties
  with the molecule
  - There's a new, and still in development, "Getting started in C++" document.
  - A lot of the Python code has been cleaned up

## New Features and Enhancements:
  - Add removeHs option to MolFromSmiles()
 (github issue #554 from greglandrum)
  - support a fixed bond length in the MolDraw2D code
 (github issue #565 from greglandrum)
  - Pattern fingerprint should set bits for single-atom fragments.
 (github issue #879 from greglandrum)
  - Reviewed unit tests of rdkit.ML - coverage now 63.1%
 (github pull #1148 from gedeck)
  - Reviewed unit tests of rdkit.VLib - coverage now 67.1%
 (github pull #1149 from gedeck)
  - Removes exponetial numBonds behavior
 (github pull #1154 from bp-kelley)
  - Exposes normalize option to GetFlattenedFunctionalGroupHierarchy
 (github pull #1165 from bp-kelley)
  - Expose RWMol.ReplaceBond to Python
 (github pull #1174 from coleb)
  - Review of rdkit.Chem.Fraggle code
 (github pull #1184 from gedeck)
  - Add support for dative bonds.
 (github pull #1190 from janholstjensen)
  - Python 3 compatibility (issue #398)
 (github pull #1192 from gedeck)
  - 1194: Review assignments of range in Python code
 (github pull #1195 from gedeck)
  - Moved GenerateDepictionMatching[23]DStructure from Allchem.py to C++
 (github pull #1197 from DavidACosgrove)
  - Review rdkit.Chem.pharm#D modules
 (github pull #1201 from gedeck)
  - Find potential stereo bonds should return any
 (github pull #1202 from coleb)
  - Gedeck coverage sim div filters
 (github pull #1208 from gedeck)
  - Gedeck review unit test inchi
 (github pull #1209 from gedeck)
  - Coverage rdkit.Dbase
 (github pull #1210 from gedeck)
  - Coverage rdkit.DataStructs
 (github pull #1211 from gedeck)
  - UnitTestPandas works on Python3
 (github pull #1213 from gedeck)
  - Cleanup and improvement to test coverage of PandasTools
 (github pull #1215 from gedeck)
  - Cleanup of rdkit.Chem.Fingerprints
 (github pull #1217 from gedeck)
  - Optimization of UFF and MMFF forcefields
 (github pull #1218 from ptosco)
  - Support for ChemAxon Extended SMILES/SMARTS
 (github issue #1226 from greglandrum)
  - Improved test coverage for rdkit.Chem.Fingerprints
 (github pull #1243 from gedeck)
  - Adding a few tests for coverage utils
 (github pull #1244 from gedeck)
  - Make Pandastools modifications to generic RDkit functionality more obvious
 (github pull #1245 from gedeck)
  - Rename test file and cleanup
 (github pull #1246 from gedeck)
  - Review of rdkit.Chem.MolKey
 (github pull #1247 from gedeck)
  - Review tests in rdkit.Chem.SimpleEnum
 (github pull #1248 from gedeck)
  - Move execution of DocTests in rdkit.Chem into a UnitTest file
 (github pull #1256 from gedeck)
  - Review code in rdkit.Chem.Suppliers
 (github pull #1258 from gedeck)
  - Add python wraps
 (github pull #1259 from eloyfelix)
  - Rename file UnitTestDocTests in rdkitChem
 (github pull #1263 from gedeck)
  - Gedeck rdkit chem unit test surf
 (github pull #1267 from gedeck)
  - cleanup rdkit.Chem.Lipinski and rdkit.Chem.GraphDescriptors
 (github pull #1268 from gedeck)
  - Address Issue #1214
 (github pull #1275 from gedeck)
  - Dev/pickle properties
 (github pull #1277 from bp-kelley)
  - Remove unused test boilerplate
 (github pull #1288 from gedeck)
  - Refactored the script SDFToCSV
 (github pull #1289 from gedeck)
  - Dev/rdmmpa api update
 (github pull #1291 from bp-kelley)
  - Fix/rogers fixes
 (github pull #1293 from bp-kelley)
  - Remove expected (error) output during unit tests
 (github pull #1298 from gedeck)
  - Refactor FeatFinderCLI and add unittests
 (github pull #1299 from gedeck)
  - Refactor BuildFragmentCatalog - 1
 (github pull #1300 from gedeck)
  - Review of rdkit.Chem code - 1
 (github pull #1301 from gedeck)
  - Minor cleanup in rdkit.Chem
 (github pull #1304 from gedeck)
  - Start using py3Dmol in the notebook
 (github pull #1308 from greglandrum)
  - Add the option to match formal charges to FMCS
 (github pull #1311 from greglandrum)
  - Review of rdkit.Chem.Subshape
 (github pull #1313 from gedeck)
  - Review rdkit.Chem.UnitTestSuppliers
 (github pull #1315 from gedeck)
  - Add cis/trans tags to double bonds
 (github pull #1316 from greglandrum)
  - MolDraw2D: make custom atom labels easier
 (github issue #1322 from greglandrum)
  - MolDraw2D: allow DrawMolecules() to put all molecules in one pane
 (github issue #1325 from greglandrum)
  - Refactoring rdkit.Chem.SATIS
 (github pull #1329 from gedeck)
  - Minor cleanup of rdkit.Chem.SaltRemover
 (github pull #1330 from gedeck)
  - Review rdkit.chem.FunctionalGroups and rdkit.Chem.UnitTestSuppliers
 (github pull #1331 from gedeck)
  - Get the tests working with python 3.6
 (github pull #1332 from greglandrum)
  - add "RDKit" to the beginning of all library names
 (github pull #1349 from greglandrum)
  - Fix/sanitizerxn merge hs
 (github pull #1367 from bp-kelley)
  - Update AllChem.py
 (github pull #1378 from BerendHuisman)

## New Java Wrapper Features:

## Bug Fixes:
  - python2 code in python3 install
 (github issue #1042 from kcamnairb)
  - Fixes #1162 (resMolSupplierTest failing with boost 1.62)
 (github pull #1166 from ptosco)
  - add missing $RDKLIBS to cartridge build
 (github pull #1167 from rvianello)
  - Include <boost/cstdint.hpp> for uint64_t
 (github pull #1168 from mcs07)
  - replace std::map::at with std::map::find
 (github pull #1169 from mcs07)
  - Fix Trajectory GetSnapshot behaviour after Clear
 (github pull #1172 from mcs07)
  - Add Contrib dir to RDPaths
 (github pull #1176 from mcs07)
  - RDThreads.h: No such file or directory
 (github issue #1177 from gncs)
  - this now builds with vs2008
 (github pull #1178 from greglandrum)
  - Add information on building RDkit on macOS using conda
 (github pull #1180 from gedeck)
  - new sequence capabilities not available from either Python or Java
 (github issue #1181 from greglandrum)
  - Gets the reaction sanitization code working correctly on 32bit systems
 (github pull #1187 from greglandrum)
  - Adds RDProps to c# wrapper
 (github pull #1188 from bp-kelley)
  - fix compatibility with PostgreSQL 9.2
 (github pull #1189 from greglandrum)
  - Fixes memory leak in closeCheckMolFiles, fixes valgrind read issue in…
 (github pull #1200 from bp-kelley)
  - Support valences of 4 and 6 for Te
 (github issue #1204 from hsiaoyi0504)
  - Stereochemistry not output to SMILES when allHsExplicit=True
 (github issue #1219 from greglandrum)
  - Remove deprecated string module functions
 (github pull #1223 from gedeck)
  - Turns on -fpermissive for gcc >= 6 and boost < 1.62
 (github pull #1225 from bp-kelley)
  - all-atom RMSD used to prune conformers in embedding code, docs say heavy-atom RMSD is used
 (github issue #1227 from greglandrum)
   - FindPotentialStereoBonds() failure
 (github issue #1230 from greglandrum)
  - make the Pandas version checking more robust
 (github pull #1239 from greglandrum)
  - Failure to embed larger aromatic rings
 (github issue #1240 from greglandrum)
   - fixed build failure on Windows due to missing link to library
 (github pull #1241 from ptosco)
  - fixed a test failure on Windows due to CR+LF encoding
 (github pull #1242 from ptosco)
  - MolFromMolBlock sanitizing when it should not be
 (github issue #1251 from greglandrum)
  - PMI descriptors incorrect
 (github issue #1262 from greglandrum)
  - Reactions don't modify isotope unless chemical element is specified for the product
 (github issue #1266 from i-tub)
  - Do not include the 3D descriptors in rdkit.Chem.Descriptors.descList
 (github issue #1287 from greglandrum)
  - ring stereochemistry perception failing for spiro centers
 (github issue #1294 from greglandrum)
  - Property pickling test failing on windows
 (github issue #1348 from greglandrum)
  - Fixes overflow error in boost when compiler chooses int for enum type
 (github pull #1351 from bp-kelley)
  - Hybridization type of group 1 metals
 (github issue #1352 from richlewis42)
  - bad python docs for some distance geometry functions
 (github issue #1385 from greglandrum)
  - Bond from reactant not added to product
 (github issue #1387 from greglandrum)
  - int32_t with no namespace in MolPickler.h
 (github issue #1388 from greglandrum)

## Contrib updates:
  - Chemical reaction role assignment code from Nadine Schneider
 (github pull #1185 from NadineSchneider)

## Deprecated code (to be removed in a future release):
- rdkit.Chem.MCS: please use rdkit.Chem.rdFMCS instead

# Release_2016.09.1
(Changes relative to Release_2016.03.1)

## Important
- The adjustQueryParameters structure has been rewritten. If you were using
  this in your code, you'll need to update your code.

## Acknowledgements:
Brian Cole, Piotr Dabrowski, Jan Domanski, Peter Gedeck, Richard Hall, Brian
Kelley, Joos Kiener, 'maddogcz', John Mayfield, 'michalsta', Michal Nowotka,
'philopon', Nico Pulver, Sereina Riniker, Stephen Roughley, Roger Sayle, Nadine
Schneider, Gianluca Sforna, Peter Shenkin, Paolo Tosco, David Turbert, Riccardo
Vianello, Maciek Wojcikowski  

## Highlights:
 - New AdjustQueryProperties() (adjustQueryProperties() in C++) for fine-tuning substructure queries.
 - ReactionEnum functionality adds significant new capabilities for doing library enumeration
 - Similarity searches with the PostgreSQL cartridge are substantially faster
 - New 3D descriptors (requires Eigen3 if you are building them yourself)

## New Features and Enhancements:
  - Trajectory/Snapshot objects
 (github pull #863 from ptosco)
  - Adds Avalon fingerprints to default set
 (github pull #871 from bp-kelley)
  - Adds the default index to the building block templates
 (github pull #873 from bp-kelley)
  - Pandas: Allow reading just the properties from SDF file
 (github pull #883 from mwojcikowski)
  - Dev/filtercatalog functional groups
 (github pull #885 from bp-kelley)
  - Dev/preprocessrxn cpp
 (github pull #892 from bp-kelley)
  - Rollup of warning squashing (with some tests diagnostics thrown in)
 (github pull #895 from bp-kelley)
  - Adds RDAny (smaller generic holder) Updates all used dictionaries
 (github pull #896 from bp-kelley)
  - expose FPS functions to SWIG
 (github pull #897 from greglandrum)
  - Add SaveFile method to the PyMol wrapper
 (github pull #898 from greglandrum)
  - Add a MultiFPBReader class
 (github pull #909 from greglandrum)
  - Improved Python doc strings for Trajectory/Snapshot objects
 (github pull #912 from ptosco)
  - Added support for building the gzip'd stream test
 (github pull #914 from ptosco)
  - Improved Trajectory Python doc strings
 (github pull #915 from ptosco)
  - improve error reporting for kekulization failures
 (github pull #919 from greglandrum)
  - Feat/github934
 (github pull #939 from greglandrum)
  - Add support for a simplified aromaticity model.
 (github pull #942 from greglandrum)
  - Dev/moldescriptors callables
 (github pull #944 from bp-kelley)
  - Dev/cleanup warnings
 (github pull #948 from greglandrum)
  - Modifications to enable building with MinGW compilers
 (github pull #960 from ptosco)
  - Made DistGeomHelpers test robust against small 3D coordinate variations
 (github pull #961 from ptosco)
  - Adds aromatization and reaction options to AdjustQuery
 (github pull #965 from bp-kelley)
  - Improved planarity for ETKDG
 (github pull #967 from sriniker)
  - Fixes built-in popcount in PgSQL cartridge on Windows
 (github pull #978 from ptosco)
  - A variety of drawing-related changes
 (github pull #986 from greglandrum)
  - Get pango 2D depiction to work with cairocffi
 (github pull #998 from ptosco)
  - Adds Atom atom map and rlabel apis
 (github pull #1004 from bp-kelley)
  - Dev/chemtransforms chirality
 (github pull #1006 from bp-kelley)
  - Added the option to label deuterium and tritium as D and T
 (github pull #1011 from ptosco)
  - Adds replaceCore function that takes a matchVect
 (github pull #1013 from bp-kelley)
  - Add an initial version of wavy bonds
 (github pull #1014 from greglandrum)
  - remove a compiler warning
 (github pull #1019 from greglandrum)
  - Make the Contrib directory available in RDConfig
 (github pull #1024 from NadineSchneider)
  - Adds some additional canned atom and bond query definitions
 (github pull #1047 from greglandrum)
  - Draw crossed bonds
 (github pull #1052 from greglandrum)
  - Alex/struct checker apr15
 (github pull #1054 from bp-kelley)
  - MolDraw2D: allow the amount of padding around atom labels to be adjusted.
 (github issue #1056 from greglandrum)
  - Add multiple molecule drawing to the C++ interface
 (github pull #1059 from greglandrum)
  - add pickle support to FilterCatalog
 (github pull #1063 from greglandrum)
  - Issue #1066: Improved .gitignore file
 (github pull #1068 from gedeck)
  - Cleanup of Scaffolds Python code
 (github pull #1069 from gedeck)
  - Consistent formatting of Python code
 (github issue #1071 from gedeck)
  - Improved test coverage of Python code
 (github pull #1081 from gedeck)
  - Improved test coverage of rdkit.DataStructs
 (github pull #1083 from gedeck)
  - Add some 3D molecular descriptors (requires Eigen3)
 (github pull #1084 from greglandrum)
  - Conformer GetPos returns a numpy array rather than a tuple of tuples
 (github pull #1087 from jandom)
  - make the 3D descriptors available in the Descriptors module
 (github pull #1097 from greglandrum)
  - Documentation update.
 (github pull #1100 from greglandrum)
  - Provide SVG output from the cartridge
 (github pull #1109 from greglandrum)
  - Allow the output of ROMol::debugMol() to show up in jupyter
 (github pull #1110 from greglandrum)
  - Dev/reaction enumeration
 (github pull #1111 from bp-kelley)
  - yapf formatting of recent changes to Python code in rdkit and Code
 (github pull #1120 from gedeck)
  - Add a parameters structure for controlling the embedding options.
 (github pull #1121 from greglandrum)
  - add more detailed error reporting when python tests fail in TestRunner.py
 (github pull #1122 from greglandrum)
  - add support for a default constructor to the python-exposed RWMol class
 (github pull #1129 from greglandrum)
  - The RunStruchk function is not exposed in pyAvalonTools
 (github issue #1130 from pulveni1)
  - SSSR performance improvements to support larger systems
 (github pull #1131 from coleb)
  - Script PythonFormat.py will test the RDkit python code for conformance with the agreed format using yapf
 (github pull #1133 from gedeck)
  - support additional trans-uranic elements
 (github pull #1134 from greglandrum)
  - Expanded sequence support
 (github pull #1140 from greglandrum)
  - add UGM link to documentation
 (github pull #1142 from David-Turbert)
  - Remove iPythonConsole configuration for normal Draw tests
 (github pull #1146 from gedeck)
  - Wrap DetectBondStereoChemistry in python
 (github pull #1156 from coleb)

## New Database Cartridge Features:
  - Provide SVG output from the cartridge
 (github pull #1109 from greglandrum)
  - Add cartridge support for adjustQueryProperties()
 (github pull #949 from greglandrum)
 - refactoring of the postgresql cartridge
(github pull #992 from rvianello)

## New Java Wrapper Features:
  - Expose filtermatch to swig
 (github pull #1117 from bp-kelley)
  - adjustQueryProperties()
  - Java wrappers for Trajectory/Snapshot objects
  (github pull #977 from ptosco)
  - Added getAlignmentTransform to ROMol.i to expose in Java SWIG wrapper
 (github pull #1155 from sroughley)

## Bug Fixes:
  - initialization of the PeriodicTable object should be made thread-safe
 (github issue #381 from greglandrum)
  - AssignAtomChiralTagsFromStructure() not recognizing chiral S
 (github issue #607 from greglandrum)
  - Fixed a few typos in Code/PgSQL/rdkit/CMakeLists.txt
 (github pull #867 from ptosco)
  - MergeQueryHs explicit H warning when no explicit Hs were actually used
 (github issue #868 from bp-kelley)
  - Fixes regression in python api CalcNumRotatableBonds
 (github pull #870 from bp-kelley)
  - Single atoms setting radius 1 bits in Morgan fingerprints
 (github issue #874 from greglandrum)
  - Providing subImgSize argument to MolsToGridImage() causes drawing failure
 (github issue #876 from greglandrum)
  - javadoc failure on CentOS 7
 (github pull #878 from ptosco)
  - adjust cartridge tests after the fix for #874
 (github pull #884 from greglandrum)
  - bugreport: invalid handling of negation of aromaticity when parsing SMARTS
 (github issue #893 from michalsta)
  - Fixes depictor problem with empty fragments
 (github pull #894 from greglandrum)
  - Fix building with G++ on Mac OS X
 (github pull #900 from johnmay)
  - linked additional libs to fix a build failure on Windows
 (github pull #901 from ptosco)
  - Rdkit 2016_03_1 generate SVG typo in Python bindings
 (github issue #903 from maddogcz)
  - PAINS filters update fails when old Python is installed
 (github issue #904 from greglandrum)
  - rdMolDraw2D.PrepareMolForDrawing() should not default to forceCoords=True
 (github issue #906 from greglandrum)
  - AddHs() using 3D coordinates with 2D conformations
 (github issue #908 from greglandrum)
  - ugly coordinates generated for peptide chains
 (github issue #910 from greglandrum)
  - Cartridge: makefile not using -O2 for C code.
 (github issue #920 from greglandrum)
  - Removes incorrect setting of hasNonPodData
 (github pull #923 from bp-kelley)
  - cleanups of RDLog's tee behavior
 (github pull #926 from greglandrum)
  - initialize boost::once_flag properly
 (github pull #927 from greglandrum)
  - sys not imported in IPythonConsole.py
 (github issue #928 from greglandrum)
  - AddTee is now SetTee
 (github pull #930 from bp-kelley)
  - mistake in SVG generated for wedged bonds
 (github issue #932 from greglandrum)
  - PandasTools AttributeError with pandas-0.18.1
 (github issue #933 from philopon)
  - Jupyter Notebooks: Issue with PyMol.MolViewer on Windows
 (github issue #936 from kienerj)
  - Subshape module: Not Python 3 compatible
 (github issue #937 from kienerj)
  - property dictionaries leaking memory
 (github issue #940 from greglandrum)
  - Bug when removing stereo info?
 (github pull #946 from mnowotka)
  - Distorted aromatic rings from ETKDG
 (github issue #952 from greglandrum)
  - MolDraw2D: default color should not be cyan
 (github issue #953 from greglandrum)
  - GetPropNames() no longer working on Atoms or Bonds
 (github issue #955 from greglandrum)
  - Kekulization issues post successful smiles parsing
 (github issue #962 from bp-kelley)
  - Fixes includes for older boost/gcc
 (github pull #966 from bp-kelley)
  - ugly conformations can be generated for highly constrained ring systems
 (github issue #971 from greglandrum)
  - Cleanup bad conformations
 (github pull #973 from greglandrum)
  - Unnecessary warnings in rxn.validate()
 (github issue #975 from greglandrum)
  - Minor fix to Code/GraphMol/Wrap/testTrajectory.py
 (github pull #979 from ptosco)
  - prepareMolForDrawing(): Do not add Hs to some three-coordinate Cs
 (github issue #982 from greglandrum)
  - MolDraw2D: wedged bonds between chiral centers drawn improperly
 (github issue #983 from greglandrum)
  - Fix format-security GCC warning
 (github pull #984 from giallu)
  - MolDraw2D scaling problems
 (github issue #985 from greglandrum)
  - RIght-justified elements in RCSB SDF files can now be parsed
 (github pull #994 from ptosco)
  - Right-justified elements in RCSB SDF files raise an exception
 (github issue #995 from ptosco)
  - ChemReactions: Bugfix in copy constructor
 (github pull #996 from NadineSchneider)
  - PgSQL README typos
 (github pull #997 from ptosco)
  - Fixes rounding errors in test
 (github pull #1001 from bp-kelley)
  - Fixes middle-justified symbols in sd files, adds M_CHG tests
 (github pull #1002 from bp-kelley)
  - fix compatibility issues with postgres < 9.5 (#1000)
 (github pull #1005 from rvianello)
  - Fixes MMFF94 aromaticity perception and ChemicalForceFields.MMFFHasAllMoleculeParams()
 (github pull #1007 from ptosco)
  - fixes typo which breaks the PostgreSQL cartridge build on Windows
 (github pull #1008 from ptosco)
  - Fix Inchi being hardcoded into PostgreSQL
 (github pull #1009 from ptosco)
  - Support ETKDG from within the SWIG wrappers
 (github pull #1010 from greglandrum)
  - move definition of computedPropName to namespace RDKit::detail
 (github issue #1017 from greglandrum)
  - fix non-inchi build
 (github pull #1018 from greglandrum)
  - Fixes #1018
 (github pull #1020 from ptosco)
  - GetSSSR interrupts by segmentation fault
 (github issue #1023 from PiotrDabr)
  - FMCS fix for Windows DLLs
 (github pull #1030 from ptosco)
  - Cause ImportError from failed dlopen of the rdBase.so shared library to propagate.
 (github pull #1032 from coleb)
  - typos in MMPA hash code
 (github issue #1044 from greglandrum)
  - MolOps::cleanUp() being called by CTAB parser even when sanitization isn't on
 (github issue #1049 from greglandrum)
  - Bond::BondDir::EITHERDOUBLE not exposed to python
 (github issue #1051 from greglandrum)
  - add python3 compatibility
 (github pull #1057 from greglandrum)
  - doc updates from Dave Cosgrove
 (github pull #1060 from greglandrum)
  - Fix leak with renumberAtoms() in the SWIG wrappers
 (github pull #1064 from greglandrum)
  - Timings on Windows with Python 3
 (github pull #1067 from ptosco)
  - computeInitialCoords() should call the SSSR code before it calls assignStereochemistry()
 (github issue #1073 from greglandrum)
  - Remove duplicates doesn't work on first column in rdkit.Dbase.DbUtils.GetData
 (github issue #1082 from gedeck)
  - clear up a bunch of windows warnings
 (github pull #1086 from greglandrum)
  - MolsToGridImage barfs on '&' in labels, at least with useSVG=True
 (github issue #1090 from shenkin)
  - Fixes csharp build for 64 bit systems
 (github pull #1098 from bp-kelley)
  - Cartridge: some C++ functions returning pointers to local storage
 (github issue #1106 from greglandrum)
  - Check for doubles after other integer types when reporting properties
 (github pull #1115 from bp-kelley)
  - Replace has_key use in Python (#issue1042)
 (github pull #1132 from gedeck)
  - fix moldraw2d headers installation path
 (github pull #1143 from giallu)
  - Remove iPythonConsole configuration for normal Draw tests
 (github pull #1146 from gedeck)
  - Adds SWIGWIN definition in WIN32 if not 32bit
 (github pull #1158 from bp-kelley)
  - Fix/java win64 memoryleak
 (github pull #1159 from bp-kelley)

## Deprecated code (to be removed in a future release):
  - rdkit.VLib python module
  - SanitizeRxn parameters "ChemDrawRxnAdjustParams" has been renamed to
    "MatchOnlyAtRgroupsAdjustParams".  These settings did not reflect
    how integrations with SciQuest or the Perkin Elmer ELN behaved and
    were confusing to users (especially since they were not explicit)

## Removed code:

## Contrib updates:
  - added an implementation of the Gobbi Atom-Atom-Path (AAP) similarity
 (github pull #1015 from Richard-Hall)

## Other:

# Release_2016.03.1
(Changes relative to Release_2015.09.2)

## Important
In order to build the RDKit, it is now necessary to have at least v1.7 of numpy installed.

## Acknowledgements:
Note: The RDKit has the wonderful "problem" that there are a lot of
contributors and it's tough for me to capture them all to put together release
notes. I don't even know many of the contributors (which is *awesome!*)
The names here come largely from what I pull in an automated way from github.
In cases where there's no real name listed in github, I either guessed
or used just the github alias in quotes. If I got it wrong, please let me know!

Josep Arus, Nik Bates-Haus, Andrew Dalke, 'DoliathGavid', 'elcaceres', Peter
Gedeck, James Jeffryes, Brian Kelley, Juuso Lehtivarjo, Rich Lewis, Daniel Lowe,
'maddogcz', Kozo Nishida, Michal Nowotka, Axel Pahl, Steven Roughley, Alexander
Savelyev, Nadine Schneider, Gianluca Sforna, Teague Sterling, Nik Stiefl, Matt
Swain, Eric Ting, Paolo Tosco, Samo Turk, Riccardo Vianello

## Highlights:
- Improvements to the build system: it's now much easier to build with InChI
  and/or Avalon support since cmake now knows how to fetch the appropriate
  source code for you. Building the PostgreSQL cartridge is now integrated into
  normal build process.
- Some improvements to molecule rendering and Jupyter notebook integration: The
  new `Draw.PrepareMolForDrawing()` function takes care of standard tasks like
  wedging bonds, kekulization, and adding chiral Hs. `Draw.MolsToGridImage()`
  can generate SVGs and uses the new molecular drawing code for PNGs when
  possible. The Jupyter notebook integration uses the new drawing code when
  possible.
- Error and warning messages from the C++ core can now be displayed in the
  Jupyter notebook

## Bug Fixes:
  - Sanitizer rejects higher valency halides
 (github issue #115 from dan2097)
  - Bad E/Z assignment from ctab
 (github issue #188 from greglandrum)
  - bad E/Z assignment from ctab
 (github issue #192 from greglandrum)
  - Documentation is still python2 specific.
 (github issue #374 from greglandrum)
  - SVG export - Python 3 support
 (github issue #398 from maddogcz)
  - FragmentOnBonds() producing incorrect chirality
 (github issue #511 from greglandrum)
  - Rings containing all dummy atoms with single bonds are flagged as aromatic
 (github issue #518 from greglandrum)
  - IPython integration broken with latest Jupyter
 (github issue #666 from samoturk)
  - Added missing include/forward declarations
 (github pull #668 from ptosco)
  - Fixes a memory leak in fragmentMol
 (github pull #669 from bp-kelley)
  - resetVect option being ignored by reaccsToFingerprint()
 (github issue #671 from greglandrum)
  - failure in AddHs when addCoords is true and coords are all zero
 (github issue #678 from greglandrum)
  - 404 error for the link to Installation instructions
 (github issue #679 from EricTing)
  - Fix java8 build
 (github pull #681 from greglandrum)
  - Smiles containing "[as]" do not parse.
 (github issue #682 from greglandrum)
  - SMARTS reaction triggers invariant violation on chiral compounds
 (github issue #685 from JamesJeffryes)
  - partially specified chiral substructure queries don't work properly
 (github issue #688 from bp-kelley)
  - ExactMolWt ignoring the mass of the electron
 (github issue #694 from greglandrum)
  - Bad 1-4 bounds matrix elements in highly constrained system
 (github issue #696 from greglandrum)
  - More ChEMBL molecules that fail bounds smoothing
 (github issue #697 from greglandrum)
  - Molecule serialization doesn't read/write atomic numbers above 128
 (github issue #713 from greglandrum)
  - AddHs cip rank is declared <int> should be unsigned int?
 (github issue #717 from bp-kelley)
  - ensure line endings are handled consistently for all users
 (github pull #729 from rvianello)
  - Fixes return type of operator[] (fails on later clangs)
 (github pull #733 from bp-kelley)
  - Fix/thread safe localeswitcher line endings
 (github pull #743 from bp-kelley)
  - Fixes Boost 1.46 issues with type traits
 (github pull #748 from bp-kelley)
  - PR #749 causes seg faults on windows
 (github issue #750 from greglandrum)
  - Fixes notebook problems with newer jupyter installs
 (github pull #753 from bp-kelley)
  - Double bond geometry loss on calling removeHs
 (github issue #754 from sroughley)
  - Bug fix to getShortestPath
 (github pull #757 from JLVarjo)
  - reversed stereochemistry with sulfoxides and ring closures
 (github issue #760 from greglandrum)
  - libRDBoost.so.1: undefined symbol
 (github issue #762 from kozo2)
  - Removed -Xdoclint:none flag when packing org.RDKitDoc.jar
 (github pull #763 from undeadpixel)
  - AnyBond specification treated as single when joining rings in SMARTS
 (github issue #766 from teaguesterling)
  - CanonicalRankAtomsInFragment() leaks when called from Python
 (github issue #769 from greglandrum)
  - MolCanvas2D drawing upside down
 (github issue #774 from greglandrum)
  - Drawing single-atom molecules hangs.
 (github issue #781 from greglandrum)
  - chiral lexical order for ring closure after branch
 (github issue #786 from adalke)
  - surface -> self.surface
 (github pull #787 from mnowotka)
  - Chem.MolToSmarts param misnomer
 (github issue #792 from elcaceres)
  - Fixes MolToSmarts python docs
 (github pull #793 from bp-kelley)
  - npscorer.py: Py3 compat and importable from other locations
 (github #801 from apahl)
  - Pre-condition Violation: bad bond type
 (github issue #805 from nbateshaus)
  - rooted atom fingerprint non identical for the same molecules
 (github issue #811 from nisti74)
  - test60RunSingleReactant() not being run
 (github issue #825 from greglandrum)
  - PostgreSQL bug fixes
 (github pull #835 from ptosco)
  - Crash while running testGithub497() on Windows
 (github pull #842 from ptosco)
  - Return value of NumRadicalElectrons and NumValenceElectrons should be integer
 (github issue #846 from gedeck)
  - Fixed a bug in getUFFAngleBendParams()
 (github pull #850 from ptosco)
  - Lines used to wedge bonds are too thick
 (github issue #852 from greglandrum)
  - Fix out of range dereference in MCS code.
 (github pull #857 from DoliathGavid)
  - Atom symbols in wrong order if bond comes from right
 (github issue #860 from greglandrum)

## New Features and Enhancements:
  - switch to using new version of avalon toolkit
 (github issue #382 from greglandrum)
  - MolDraw2D: Expand basic drawing api
 (github issue #417 from greglandrum)
  - MolDraw2D: add options
 (github issue #424 from greglandrum)
  - fixed FutureWarning in PeriodicTable.py
 (github pull #665 from richlewis42)
  - first pass, using google style
 (github pull #672 from greglandrum)
  - Use sets instead of and map. Minor comments cleanup.
 (github pull #675 from DoliathGavid)
  - Dev/squash msvc14 warnings
 (github pull #684 from bp-kelley)
  - Fix/stop unnecessary filtercatalog updates
 (github pull #690 from bp-kelley)
  - Add RDK_USE_BOOST_SERIALIZATION configure option (On by default)
 (github pull #691 from bp-kelley)
  - Minor optimizations of the force field minimization code, fix for issue 696
 (github pull #693 from greglandrum)
  - Include cis/trans stereochemistry when useChirality=true with the morgan fingerprints
 (github issue #695 from greglandrum)
  - Fixed a couple of compilation warnings in Resonance.cpp/Resonance.h
 (github pull #701 from ptosco)
  - Dev/run single reactant
 (github pull #705 from bp-kelley)
  - Updates CMAKE_SOURCE_DIR to CMAKE_CURRENT_SOURCE_DIR
 (github pull #707 from bp-kelley)
  - Make LocaleSwitcher threadsafe
 (github issue #710 from greglandrum)
  - Exposes Get/Set Double, Int, Uint and bool props to molecules
 (github pull #711 from bp-kelley)
  - Speed up molblock generation
 (github pull #712 from greglandrum)
  - Expose generateOneProductSet?
 (github issue #721 from DoliathGavid)
  - Add a reader for FPB files (still experimental)
 (github pull #724 from greglandrum)
  - replace std::map::at with std::map::operator[]
 (github pull #730 from rvianello)
  - Fix/get double prop get props asdict
 (github pull #734 from bp-kelley)
  - Add support for Tversky similarity to the FPB reader
 (github pull #735 from greglandrum)
  - Fix ConformerParser to use const std::string &
 (github pull #737 from mcs07)
  - Fix/expose invariant exception
 (github pull #740 from bp-kelley)
  - Support CTABs where the second letter in atom symbols is capitalized
 (github issue #741 from greglandrum)
  - Adds support for capturing RDLogs in Python StdErr streams
 (github pull #749 from bp-kelley)
  - Allow adding Hs only to atoms matching a query operator
 (github issue #758 from greglandrum)
  - Add argument to addHs allowing only certain Hs to be considered
 (github pull #759 from greglandrum)
  - avoid the multiple definition of rdkitVersion/boostVersion
 (github pull #761 from rvianello)
  - cleanup possible pythonObjectToVect leaks in python wrappers
 (github issue #764 from greglandrum)
  - Stop possible core leaks in pythonObjectToVect()
 (github pull #770 from greglandrum)
  - Add C++ function to prepare mol for rendering
 (github issue #771 from greglandrum)
  - Prefer wedging bonds to Hs
 (github issue #772 from greglandrum)
  - Add prepareMolForDrawing() function to C++
 (github pull #775 from greglandrum)
  - Support blanks in MolsToGridImage()
 (github issue #776 from greglandrum)
  - A number of small additions and features to the drawing code
 (github pull #802 from greglandrum)
  - Support larger isotope deltas in the chirality assignment
 (github issue #803 from greglandrum)
  - Adds option RDK_USE_COMPLEX_ROTOR_DEFINITION
 (github pull #810 from bp-kelley)
  - add Draw.MolsToSVGGrid()
 (github pull #817 from greglandrum)
  - make Hs black instead of gray
 (github pull #819 from greglandrum)
  - Fix alignMols so that it takes into account of QueryAtoms and QueryBonds
 (github pull #821 from DoliathGavid)
  - feat/github831: Add getText() static method.
 (github pull #832 from greglandrum)
  - Add an unfolded count-based version of the RDKFingerprint
 (github pull #838 from NadineSchneider)
  - Add some utils functions to ChemReactions
 (github pull #840 from NadineSchneider)
  - Autodetect boost c++ library and compile with matching one
 (github pull #845 from bp-kelley)
  - Add automatic downloads of junit.jar
 (github pull #859 from greglandrum)

## New Database Cartridge Features:
  - support providing InChI (or InChI key) generation options in cartridge
 (github pull #755 from greglandrum)
  - building the cartridge is now integrated with the cmake build system
 (github pull #785 from ptosco)

## New Java Wrapper Features:
  - Add a bit more control over the lazy MaxMin picker to the java layer
 (github pull #791 from greglandrum)
  - Ensure reduceProductToSideChains exposed in Java/Swig
 (github issue #744 from bp-kelley)

## Deprecated code (to be removed in next release):

## Removed code:

## Contrib updates:

## Other:

# Release_2015.09.2
(Changes relative to Release_2015.09.1)

## Acknowledgements:
Brian Kelley, Paolo Tosco, Riccardo Vianello

## Bug Fixes:
  - Fixed a post-decrement which causes a crash when compiling under Windows with MSVC 9
  (from ptosco)
  - Fixes a memory leak in fragmentMol
  (github #669 from bp-kelley)
  - MMPA compile error with Microsoft Visual C++ Compiler for Python 2.7
  (github #655 from rvianello)

# Release_2015.09.1
(Changes relative to Release_2015.03.1)

## Acknowledgements:

Pierre Bhoorasingh, Gungor Budak, Andrew Dalke, JP Ebejer, Peter Ertl,
Jan Holst Jensen, Brian Kelley, Joos Kiener, Noel O'Boyle, Axel Pahl,
Sereina Riniker, Alexander Savelyev, Roger Sayle, Nadine Schneider,
Andrew Stretton, Paolo Tosco, Samo Turk, JL Varjo, Riccardo Vanello,
Maciek Wojcikowski

## Highlights:

  - Addition of parsers/writers for sequence notation, FASTA, and basic HELM
  - Improved conformation generation based on experimental torsional parameters
  - Much better filtering of generated conformations to ensure they
    match the chirality of the input structure
  - New method for enumerating molecular resonance structures
  - Addition of a molecular FilterCatalog data structure

## Bug Fixes:
  - Draw.MolToImage(mol) does not work for Python 3, because cairo for Python 3 has not yet implemented Surface.create_for_data
 (github issue #460 from apahl)
  - SDWriter.close() fails when the underlying stream has an error
 (github issue #498 from greglandrum)
  - Hexafluorophosphate cannot be handled
 (github issue #510 from greglandrum)
  - Labels of highlighted atoms are invisible
 (github issue #519 from NadineSchneider)
  - MolDraw2D: Fix in highlighting atoms
 (github pull #521 from NadineSchneider)
  - Cartridge: index failing for molecular equality in some circumstances
 (github issue #525 from greglandrum)
  - Bad ring finding in a complex fused ring
 (github issue #526 from greglandrum)
  - Fixed crash upon closing a gzip/bzip2 stream opened in binary mode for use with SDWriter under Python3
 (github pull #531 from ptosco)
  - Regression: _smilesAtomOutputOrder incorrect for dot disconnected molecules
 (github issue #532 from baoilleach)
  - Fix #532 - smilesAtomOutputOrder incorrect
 (github pull #535 from baoilleach)
  - Fix Python3 encoding for FilterCatalog/Entry serialization
 (github pull #537 from bp-kelley)
  - Catalog::setCatalogParams needs to be virtual now
 (github pull #538 from bp-kelley)
  - Bonds in allyl cation are not conjugated
 (github issue #539 from greglandrum)
  - Fixes GitHub issue 539
 (github pull #540 from ptosco)
  - SaltBuilder._initPatterns incorrectly handles SMARTS errors
 (github issue #541 from adalke)
  - merging query Hs failing on recursive SMARTS
 (github issue #544 from greglandrum)
  - Crash in O3A alignment when running multi-threaded
 (github issue #546 from greglandrum)
  - PyMol.py: changed xmlrpc import to have Python2/Python3 compatibility
 (github pull #547 from apahl)
  - fix MMPA bugs for some tests
 (github pull #548 from AlexanderSavelyev)
  - TFD fix for single bonds adjacent to triple bonds
 (github pull #550 from sriniker)
  - Canonicalization paper Aug2015
 (github pull #552 from NadineSchneider)
  - Chirality not affected by atom-map index
 (github issue #553 from adalke)
  - Implicit Hs should not appear in depictions for query atoms.
 (github issue #556 from greglandrum)
  - Fix issue with invalid reactions throwing NameError
 (github pull #559 from bp-kelley)
  - InChI radicals not properly converted
 (github issue #562 from pierrelb)
  - MMPA code not python3 compatible
 (github issue #564 from greglandrum)
  - mmpa: fix bug with num_of_cuts > 2 case
 (github pull #566 from AlexanderSavelyev)
  - Incorrect stereochemistry after embedding
 (github issue #568 from greglandrum)
  - changed FrameToGridImage() so that the dataframe index can be used as legend
 (github pull #570 from apahl)
  - MolDraw2DCairo.GetDrawingText() doesn't work with Python3
 (github issue #571 from greglandrum)
  - addBond return value
 (github issue #572 from JLVarjo)
  - Process aborts during ForwardSDMolSupplier gc when the file object is closed
 (github issue #579 from adalke)
  - Fix/update pains filter catalog
 (github pull #581 from bp-kelley)
  - Importing PandasTools on Windows fails due to Salts.txt
 (github issue #583 from baoilleach)
  - renumberAtoms() not setting conformer dimensionality properly
 (github issue #584 from greglandrum)
  - stereo atoms property on double bonds not being updated properly with insertMol
 (github issue #608 from greglandrum)
  - UFF Atom type not properly assigned to lanthanides
 (github issue #613 from greglandrum)
  - segfault from MolToInchi when bad bond stereochem info is present
 (github issue #614 from greglandrum)
  - MQN12 (heavy atom count) seems to be always 0.
 (github issue #623 from kienerj)
  - fmcs: fix issue with initial seed for PR 580
 (github pull #624 from AlexanderSavelyev)
  - fmcs: fixes #631 with chiralirty
 (github pull #634 from AlexanderSavelyev)
  - Fixes sprintf not found on some gcc compiles
 (github pull #635 from bp-kelley)
  - Fix AppVeyor and Travis UnitTests
 (github pull #636 from bp-kelley)
  - Fixes #629 - python GetSubstructureMatch thread safety
 (github pull #637 from bp-kelley)
  - Fix regressions occurring when building with msvc9
 (github pull #638 from rvianello)
  - Fix/python gil release on rdkit threadsafe sss only
 (github pull #639 from bp-kelley)
  - ctest not running some of the python tests.
 (github issue #643 from greglandrum)
  - Issue643
 (github pull #646 from greglandrum)
  - Fix/various bug fixes filtercatalog and bad operator
 (github pull #648 from bp-kelley)
  - unlocking or destroying a locked mutex in the substructure matcher
 (github issue #653 from greglandrum)
  - MMPA compile error with Microsoft Visual C++ Compiler for Python 2.7
 (github issue #655 from rvianello)
  - new canon: fix in special symmetry invariant
 (github pull #663 from NadineSchneider)

## New Features and Enhancements:
  - enable popcount by default for cartridge
 (github issue #428 from greglandrum)
  - Added the RSMD calculation over N molecules in the cookbook.
 (github pull #495 from malteseunderdog)
  - Modified force field constraint tests to be more robust
 (github pull #503 from ptosco)
  - Fix mol drawing on Python3 (issue #460)
 (github pull #504 from apahl)
  - mmpa first version was added
 (github pull #505 from AlexanderSavelyev)
  - Forcefield tests now use RDKit::feq() instead of RDKit::round()
 (github pull #506 from ptosco)
  - SDMolSupplier(), setData() and strictParsing
 (github pull #507 from ptosco)
  - Improvements to LoadSDF and WriteSDF
 (github pull #513 from samoturk)
  - updates to mmpa. reduce number of smiles parsing
 (github pull #515 from AlexanderSavelyev)
  - Some enhancements for the new canonicalization
 (github pull #520 from NadineSchneider)
  - mmpa: remove boost_regex dependency at all. Add new test references
 (github pull #527 from AlexanderSavelyev)
  - Support getting atoms involved in Pharm2D bits
 (github issue #530 from greglandrum)
  - Optimized MMFF::MMFFOptimizeMoleculeConfs()
 (github pull #534 from ptosco)
  - RDKit learns how to filter PAINS/BRENK/ZINC/NIH via FilterCatalog
 (github pull #536 from bp-kelley)
  - Expose Conformer's copy constructor to Python
 (github issue #545 from greglandrum)
  - PyMol.py: changed xmlrpc import to have Python2/Python3 compatibility
 (github pull #547 from apahl)
  - Update PAINS smarts to validated set.  Always use mergeQueryHs when reading Smarts
 (github pull #549 from bp-kelley)
  - Add a curated set of PAINS filters to the RDKit
 (github issue #555 from greglandrum)
  - Update pains from wehi_pains.csv
 (github pull #560 from bp-kelley)
  - Change MMPA to use CanonicalRankAtoms instead of _CIPRank
 (github issue #561 from adalke)
  - Add adjustQuery() function to make substructure queries more specific.
 (github issue #567 from greglandrum)
  - changed sascorer.py to enable import from different location
 (github pull #569 from apahl)
  - add export_values() to enums in the python wrapper where it's sensible to do so
 (github issue #573 from greglandrum)
  - Support Sheridan's BT and BP fingerprints
 (github issue #574 from greglandrum)
  - Locale dependent float for _GasteigerCharge
 (github issue #577 from adalke)
  - fmcs: implement adding an initial seed structure
 (github pull #580 from AlexanderSavelyev)
  - Where possible convert docs from RST to MD
 (github issue #585 from greglandrum)
  - [UGM2015] Autodetection of numThreads
 (github issue #586 from mwojcikowski)
  - Generating SVG Images with Transparent Background
 (github issue #587 from gungorbudak)
  - updates to PandasTools.LoadSDF
 (github pull #599 from adalke)
  - Can control mol depiction size with PandasTools.molSize = (200,200).
 (github pull #600 from samoturk)
  - pandas performance and functionality improvements
 (github pull #601 from adalke)
  - Adding documentation for installation with conda.
 (github pull #602 from strets123)
  - Automatic atom reordering in TorsionFingerprints
 (github pull #603 from sriniker)
  - remove bare excepts
 (github pull #606 from adalke)
  - Added option to use SVG rendering in pandas data frames
 (github pull #609 from samoturk)
  - Handle inserting molecules with conformers into molecules without conformers
 (github issue #610 from greglandrum)
  - If wedged bonds are already present, write them to mol blocks
 (github issue #611 from greglandrum)
  - Dev/filter catalog java wrapper
 (github pull #612 from bp-kelley)
  - Support extended reduced graphs
 (github issue #615 from greglandrum)
  - port NumSpiroCenters and NumBridgeheadAtoms descriptors to C++
 (github issue #617 from greglandrum)
  - Add parser/writer for peptide sequences
 (github issue #620 from greglandrum)
  - Add parser/writer for FASTA
 (github issue #621 from greglandrum)
  - Add parser/writer for HELM
 (github issue #622 from greglandrum)
  - Migrate std::string APIs to const std::string &
 (github pull #627 from bp-kelley)
  - Chirality tests
 (github pull #628 from sriniker)
  - Improvements of TFD
 (github pull #630 from sriniker)
  - Added ResonanceMolSupplier
 (github pull #632 from ptosco)
  - Image tostring/fromstring methods replaced by tobytes/frombytes
 (github pull #644 from rvianello)
  - ET(K)DG implementation
 (github pull #647 from sriniker)
  - Build multithreading support by default when boost::threads is present
 (github issue #649 from greglandrum)
  - Minimal SWIG Java wrappers for ResonanceMolSupplier
 (github pull #657 from ptosco)

## New Database Cartridge Features:
  - Support for PostgreSQL v8.x has been removed
  - NumSpiroCenters and NumBridgeheadAtoms added

## New Java Wrapper Features:
  - Support for FilterCatalogs
  - Support for ResonanceMolSuppliers

## Deprecated code (to be removed in next release):

## Removed code:
  - rdkit/DataStructs/BitVect.py : a C++ version is used, this was only present for historical reasons
  - rdkit/DataStructs/SparseIntVect.py : a C++ version is used, this was only present for historical reasons

## Contrib updates:
  - Addition of Peter Ertl's Natural Product Likeness score.

## Other:
  - Much of the documentation has been translated from RST to MD

# Release_2015.03.1
(Changes relative to Release_2014.09.2)

## IMPORTANT

 - This release has a new algorithm for canonical atom ordering. This
   means that canonical SMILES generated with the new version will be
   different from those generated with previous versions.

## Acknowledgements:

David Cosgrove, Andrew Dalke, JP Ebejer, Niko Fechner, Igor Filippov,
Patrick Fuller, David Hall, J Bach Hardie, Jan Holst Jensen, Brian
Kelley, Rich Lewis, John May, Michael Reutlinger, Sereina Riniker,
Alexander Savelyev, Roger Sayle, Nadine Schneider, Gianluca Sforna,
Paolo Tosco, Samo Turk, JL Varjo, Riccardo Vianello

## Highlights:

  - A new algorithm for generating canonical atom orders. The new
    algorithm is faster and more robust than the old one.
  - C++-based molecule drawing code, allows consistent molecule
    renderings from C++, Java, and Python. This will become the
    default renderer in the next release.
  - General performance improvements and reduction of memory usage.
  - Torsion Fingerprints for comparing 3D conformations to each other
  - MCS code now available from the PostgreSQL cartridge

## Bug Fixes:

  - fmcs: use the zero bond type instead of other
    (github issue #368 from AlexanderSavelyev)
  - bug fixed in TorsionFingerprints.GetTFDMatrix
    (github issue #376 from sriniker)
  - rdkit.Chem.MolDb.Loader_sa does not work with python3.4
    (github issue #390)
  - ChemReactions: Bugfix/Memleak-fix in runReactants
    (github issue #394 from NadineSchneider)
  - TorsionConstraint bug fix
    (github issue #395 from ptosco)
  - Fixed LoadSDF to keep 3D info
    (github pull #401 from samoturk)
  - Incorrect expected absolute stereochemistries in a test
    (github issue #407)
  - GetSubstructMatches() consumes all available memory
    (github issue #409)
  - Bugfix in the chirality handling of chemical reactions
    (github pull #410 from NadineSchneider)
  - fixed two minor bugs in MMFF code
    (github issue #416 from ptosco)
  - Varied sanitise rejection with hydrogen representation
    (github issue #418)
  - Fixed two trivial bugs
    (github issue #419 from ptosco)
  - RDKit learns how to roundtrip MOLFile values.
    (github issue #420 from bp-kelley)
  - cairoCanvas.py python3 compatibility
    (github issue #426 from bach-hardie)
  - Chem.FragmentOnSomeBonds() should update implicit H count on aromatic heteroatoms when  addDummies is False
    (github issue #429)
  - Chem.FragmentOnSomeBonds() crashes if the bond list is empty
    (github issue #430)
  - Increase limit for smallest ring size
    (github issue #431 from cowsandmilk)
  - Problems caused by heteroatoms with radicals in aromatic rings.
    (github issue #432)
  - Conversion from mol to InChI getting ring stereo wrong
    (github issue #437)
  - Metals in mol2 blocks handled incorrectly
    (github issue #438)
  - Fixed a few tests failing when the Windows git client is used
    (github pull #439 from ptosco)
  - Fixed tests failing on Windows when retrieving sources through the  Windows git client
    (github pull #440 from ptosco)
  - Conformers not copied correctly in renumberAtoms
    (github issue #441)
  - Radicals are not correctly assigned when reading from SMILES
    (github issue #447)
  - moving RDBoost/Exceptions.h to RDGeneral
    (github pull #458 from rvianello)
  - Add code to accept blank SMILES input.
    (github issue #468 from janholstjensen)
  - Accept empty SMILES and SMARTS
    (github issue #470)
  - Fix MolFile Atom Line List Writer
    (github issue #471 from bp-kelley)
  - Moved 3D ipython renderer to dependency
    (github pull #472 from patrickfuller)
  - Windows build/test failure fixes
    (github issue #473 from ptosco)
  - install missing FMCS/Graph.h header file
    (github pull #478 from rvianello)
  - added const qualifiers to some methods of Atom/Bond iterator classes
    (github pull #479 from rvianello)
  - Bugfix in SmilesWrite and some additional tests for getMolFrags function
    (github issue #482 from NadineSchneider)
  - Removed a while(1) {} in BFGSOpt.h which might result in an infinite loop
    (github issue #484 from ptosco)
  - Gasteiger charge calculation fails with hexavalent sulfur #485
    (github issue #485)
  - SmilesWriter not creating automatic name values for molecules read from CTABs
    (github issue #488)
  - Fix StringType access, remove unused imports
    (github issue #494 from bp-kelley)

## New Features:

  - Isomeric fix in PandasTools
    (github issue #369 from samoturk)
  - added SaveXlsxFromFrame - first version
    (github issue #371 from samoturk)
  - New feature: Torsion fingerprints
    (github issue #372 from sriniker)
  - Add function to extract a molecule with a single (particular) conformer from a multiconf mol
    (github issue #384)
  - Added C++ and Python helpers to retrieve force-field parameters
    (github issue #387 from ptosco)
  - Support providing all options to SWIG-wrapped FMCS
    (github issue #397)
  - Add function to request if calling UpdatePropertyCache is necessary
    (github issue #399)
  - Add function to test if UpdatePropertyCache is necessary
    (github issue #400 from NadineSchneider)
  - Substructure highlighting in pandas dataframe (fixes #362)
    (github issue #403 from nhfechner)
  - Added SDF Export to PandasTools and adjusted SDF Import
    (github issue #404 from nhfechner)
  - support count-based avalon fingerprint
    (github issue #408)
  - New C++ drawing code
    (github issue #412 from greglandrum)
  - RDKit learns how to compute code coverage for tests
    (github issue #413 from bp-kelley)
  - Dictionary access is saniztized and optimized.
    (github issue #414 from bp-kelley)
  - expose MolOps::rankAtoms() and MolOps::rankAtomsInFragment() to python
    (github issue #421)
  - Dev/expose rank atoms to python
    (github issue #422 from bp-kelley)
  - rdkit learns how to wrap a proper RWMol in python
    (github issue #423 from bp-kelley)
  - Docs: document "magic" property values that are used throughout the code
    (github issue #425)
  - MolDraw2D: expose class to Python
    (github issue #433)
  - RDKit learns how to query properties on Atoms
    (github issue #442 from bp-kelley)
  - Issue445: provide helper functions for multithreaded evaluation of some expensive functions
    (github issue #448 from greglandrum)
  - RDKit learns how to release the GIL in python
    (github pull #449 from bp-kelley)
  - Dev/property queries
    (github pull #450 from bp-kelley)
  - Support a confId argument to the atom pair fingerprinting code
    (github issue #453)
  - Save the atom bookmarks so we can deconvolute reaction products.
    (github pull #454 from bp-kelley)
  - Cartridge: Support converting no structs to InChI
    (github issue #455)
  - RDKit learns how to expose ChemicalReaction copy constructor to python
    (github pull #456 from bp-kelley)
  - chirality flag was implemented for fmcs() function
    (github pull #466 from AlexanderSavelyev)
  - Support copy and deepcopy properly for at least Mol and RWMol
    (github issue #467)
  - Cartridge: add qmol_from_smiles() and qmol_from_ctab()
    (github issue #469)
  - restore java and python wrappers. New parameter (matchChiralTag)
    (github issue #477 from AlexanderSavelyev)
  - Added a Python wrapper for getShortestPath()
    (github issue #487 from ptosco)
  - Dev/merge query hs no unmapped atoms
    (github issue #490 from bp-kelley)

## New Database Cartridge Features:

  - The MCS code is now available within the cartridge
  - The functions qmol_from_smiles() and qmol_from_ctab()

## New Java Wrapper Features:

  - The new molecule rendering code is accessible from the SWIG wrappers.

## Deprecated code (to be removed in next release):

  - C++: The functionality in $RDBASE/Code/GraphMol/MolDrawing has been
    superseded by the new drawing code in $RDBASE/Code/GraphMol and will
    be removed in the next release.
  - Python:
     - rdkit/Dbase/Pubmed
     - rdkit/Chem/fmcs (this has been superseded by the C++ implementation)
  - Cartridge: support for v8.x of PostgreSQL (v8.4 is no longer
    supported by the PostgreSQL project team as of July 2014)

## Removed code:

  - the method Atom::setMass() has been removed. Please use setIsotope()
    instead.
  - the methods involving an atom's dativeFlag have been removed.

## Contrib updates:

## Other:
  - Python 2.6 support is deprecated. Starting with the next RDKit
release, we will only support python 2.7 and python 3.4 and
higher. Python 2.6 has not been supported since October 2013. If you
believe that you are stuck with python 2.6 because of the version of
RHEL you are running, please read this post to learn about your
options:
http://www.curiousefficiency.org/posts/2015/04/stop-supporting-python26.html
  - The RDKit molecule objects (ROMol and RWMol) now require about 15%
less memory to store


# Release_2014.09.2
(Changes relative to Release_2014.09.1)

## Acknowledgements:
Sereina Riniker, Nadine Schneider, Paolo Tosco

## Bug Fixes:
- SMILES parser doing the wrong thing for odd dot-disconnected construct
 (github issue #378)
- O3A code generating incorrect results for multiconformer molecules
 (github issue #385)
- Suppliers probably not safe to read really large files
 (github issue #392)
- Torsion constraints not worrking properly for negative torsions
 (github issue #393)
- Core leak in reactions
 (github issue #396)


# Release_2014.09.1
(Changes relative to Release_2014.03.1)

## Acknowledgements:
Andrew Dalke, James Davidson, Jan Domanski, Patrick Fuller, Seiji
Matsuoka, Noel O'Boyle, Sereina Riniker, Alexander Savelyev, Roger
Sayle, Nadine Schneider, Matt Swain, Paolo Tosco, Riccardo Vianello,
Richard West

## Bug Fixes:
- Bond query information not written to CTAB
 (github issue #266)
- Bond topology queries not written to CTABs
 (github issue #268)
- Combined bond query + topology query not correctly parsed from CTAB
 (github issue #269)
- SWIG wrapped suppliers leak memory on .next()
 (github issue #270)
- SWIG wrappers don't build with SWIG 3.0.x
 (github issue #277)
- core leak from DataStructs.ConvertToNumpyArray
 (github issue #281)
- MolTransforms not exposed to Java wrapper
 (github issue #285)
- Seg fault in ReactionFromRxnBlock
 (github issue #290)
- BitInfo from GetHashedMorganFingerprint() has non-folded values
 (github issue #295)
- bad inchi for chiral S when the molecule is sanitized
 (github issue #296)
- Cannot generate smiles for ChEBI 50252
 (github issue #298)
- Either molecule-molecule substruct matching is wrong *OR* the docs for Atom::Match incorrect
 (github issue #304)
- fluorine F-F  gives segmentation fault with MMFF forcefield
 (github issue #308)
- cartridge: MACCS similarity wrong when using the builtin popcount and the index
 (github issue #311)
- Substructure Search via SMARTS implicit hydrogens
 (github issue #313)
- SMARTS output for [x] is wrong
 (github issue #314)
- Bonds not being set up properly in renumberAtoms
 (github issue #317)
- Python 2 code in python 3 branch
 (github issue #326)
- Linking error with ICC 15.0 on Linux
 (github issue #327)
- Using explicit hydrogens in the SMILES lead to the same AP FP for two different molecules
 (github issue #334)
- memory leaks when smiles/smarts parsers fail
 (github issue #335)
- No double bond stereo perception from CTABs when sanitization is turned off
 (github issue #337)
- missing MACCS key 44 might be found
 (github issue #352)
- Hydrogens in mol blocks have a valence value set
 (github issue #357)
- Computed props on non-sanitized molecule interfering with substructure matching
 (github issue #360)
- Fixed a weakness in the angular restraint code
 (github pull #261 from ptosco)
- A few fixes to improve MMFF/UFF robustness
 (github pull #274 from ptosco)
- Static webGL rendering fix
 (github pull #287 from patrickfuller)
- Revert #include ordering in SmilesMolSupplier.cpp
 (github pull #297 from mcs07)
- Add missing include for RDDepict::compute2DCoords
 (github pull #301 from baoilleach)
- Herschbach-Laurie fallback implemented to fix GitHub 308
 (github pull #312 from ptosco)
- Issue #320 Making GetBestRMS more idiot-proof
 (github pull #322 from jandom)
- Update URLs to InChI API after inchi-trust.org website redesign.
 (github pull #341 from rwest)

## New Features:
- Should be able to do intramolecular bond breaking in reactions.
 (github issue #58)
- Support reactions in cartridge
 (github issue #223)
- Documentation of Inchi methods
 (github issue #240)
- add DescribeQuery() to Bond python wrapper
 (github issue #267)
- support avalon fingerprint in cartridge
 (github issue #286)
- support partial fragmentation with fragmentOnSomeBonds
 (github issue #288)
- Add calcNumHeterocycles() descriptor
 (github issue #351)
- C++ implementation of FMCS algorithm
- Reordering feature for Butina clustering
 (github pull #302 from sriniker)
- Changes and new functions for the calculation of RMS values between conformers of a molecule
 (github pull #306 from sriniker)
- Extended chemical reaction functionality and add chemical reactions to cartridge
 (github pull #315 from NadineSchneider)
- Custom color to highlight atoms in Mol2Image
 (github pull #316 from jandom)
- Several different fingerprint algorithms for chemical reactions are now available
- add Chem.Draw.MolToQPixmap
 (github pull #355 from mojaie)


## New Database Cartridge Features:
- *NOTE:* the configuration variable rdkit.ss_fp_size has been renamed to rdkit.sss_fp_size
- Chemical reactions and several operations on them are now supported
- Avalon fingerprints now supported (when support has been compiled in)


## New Java Wrapper Features:
- FMCS implementation exposed
- Fingerprints for chemical reactions
- Possible core leak in some of the MolSuppliers was fixed

Deprecated modules (to be removed in next release):
- Projects/SDView4
- rdkit/utils/
  - GUIDs.py
  - GenLicense.py
  - Licensing.py
  - PiddleTools.py
  - PilTools.py
  - REFile.py
  - SliceTools.py
- rdkit/Logger

Removed modules:

## Contrib updates:

## Other:
- The RDKit now supports both python3 and python2.
- There is now conda integration for the RDKit.
- SMILES generation is substantially faster


# Release_2014.03.1
(Changes relative to Release_2013.09.2)

## IMPORTANT
 - Due to a bug fix in the rotatable bonds definition, the default
   rotatable bond calculation returns different values than before.
   This also affects MQN descriptor #18.

## Acknowledgements:
Paul Czodrowski, James Davidson, Markus Elfring, Nikolas Fechner, Jan Holst Jensen, Christos Kannas, Sereina Riniker, Roger Sayle, Paolo Tosco, Samo Turk, Riccardo Vianello, Maciej Wójcikowski, Toby Wright

## Bug Fixes:
- Dict::DataType declaration causing problems with C++11 std::lib
 (github issue 144)
- Pre-condition Violation in AllChem.Compute2DCoords
 (github issue 146)
- GetSimilarityMapForFingerprint() fails when similarity is zero
 (github issue 148)
- PatternFingerprint failure for substructure matching
 (github issue 151)
- query atoms don't match ring queries
 (github issue 153)
- Incorrect SMILES generated after calling MMFF parameterization
 (github issue 162)
- Problems with Xe from SDF
 (github issue 164)
- Radicals not being used in atom--atom matches
 (github issue 165)
- Cannot skip sanitization while reading PDB
 (github issue 166)
- Distance Geometry embedding not thread safe
 (github issue 167)
- O3A::align() and O3A::trans() now return "true" RMSD value
 (github pull 173)
- RangeError when pre-incrementing or decrementing AtomIterators
 (github issue 180)
- ctabs do not contain wedged bonds for chiral s
 (github issue 186)
- ctabs do not contain "A" when appropriate
 (github issue 187)
- Problems round-tripping Al2Cl6 via CTAB
 (github issue 189)
- Don't merge Hs onto dummies
 (github issue 190)
- Wavy bonds to Hs in CTABs should affect the stereochemistry of attached double bonds
 (github issue 191)
- Rendering binary compounds as ClH, FH, BrH or IH rather than putting H first.
 (github issue 199)
- Fixed data race condition in Code/GraphMol/MolAlign/testMolAlign.cpp
 (github pull 202)
- Re-prepared SDF/SMILES files of the MMFF validation suite + a fix
 (github pull 205)
- Problems round-tripping P with non-default valence.
 (github issue 206)
- Added a stricter definition of rotatable bonds as a new function in the ...
 (github pull 211)
- Code/GraphMol/AddHs patch proposal
 (github pull 212)
- Fix: Changed getNumReactantTemplates to GetNumReactantTemplates.
 (github pull 219)
- aromatic B ("b") causes errors from SMARTS parser
 (github issue 220)
- Segmentation fault for MMFF optimization with dummy atoms
 (github issue 224)
- isMoleculeReactantOfReaction() and isMoleculeProductOfReaction() not useable from SWIG wrappers
 (github issue 228)
- cartridge: mol_from_ctab() ignores argument about keeping conformers
 (github issue 229)
- Reaction not correctly preserving chirality on unmapped atoms.
 (github issue 233)
- AllChem.AssignBondOrdersFromTemplate() fails with nitro groups
 (github issue 235)
- Fix molecule dataframe rendering in pandas 0.13.x
 (github pull 236)
- Dummy labels copied improperly into products
 (github issue 243)
- Two bugfixes in MMFF code
 (github pull 248)
- seg fault when trying to construct pharmacophore with no conformation
 (github issue 252)
- EmbedMolecule() should not create a conformer for molecules that have zero atoms
 (github issue 256)
- cartridge: dice similarity calculation does not use USE_BUILTIN_POPCOUNT flag
 (github issue 259)
- cartridge: similarity calculations wrong for maccs fps when USE_BUILTIN_POPCOUNT flag is set
 (github issue 260)

## New Features:
- Expose gasteiger charge calculation to SWIG
 (github issue 152)
- Added additional functionality to PandasTools
 (github pull 155)
- Add MMFFHasAllMoleculeParams() to SWIG interface
 (github issue 157)
- O3A code should throw an exception if the parameterization is not complete.
 (github issue 158)
- Support zero order bonds
 (github issue 168)
- Add attachmentPoint argument to ReplaceSubstructs
 (github issue 171)
- Forcefield constraints (distances, angles, torsions, positions)
 (github pull 172)
- Add kekulize flag to SDWriter
 (github issue 174)
- Support operator= for RWMol
 (github issue 175)
- Get GetAromaticAtoms() and GetQueryAtoms() working from python
 (github issue 181)
- Support richer QueryAtom options in Python
 (github issue 183)
- Support writing V3000 mol blocks
 (github issue 184)
- Allow disabling the building of tests
 (github issue 185)
- Expand DbCLI to allow updating databases
 (github issue 197)
- Code refactoring and enhancement to allow for O3A alignment according to atom-based Crippen logP contribs
 (github pull 201)
- call MolOps::assignStereochemistry() with flagPossibleStereoCenters true from within the molecule parsers.
 (github issue 210)
- Support passing of file-like PandasTools.LoadSDF
 (github pull 221)
- Reaction SMARTS parser should support agents
 (github issue 222)
- Add function to MolOps to allow a molecule to be split into fragments based on a query function
  This is useable from python as Chem.SplitMolByPDBResidues() and Chem.SplitMolByPDBChainId()
 (github issue 234)
- Adding option useCounts for Morgan fingerprints
 (github pull 238)
- support SimpleEnum functionality for adding recursive queries to reactions
 (github issue 242)
- Additional functions for bit vectors
 (github pull 244)
- Support of RDK fingerprints added to SimilarityMaps
 (github pull 246)
- add get3DDistance
- support 3D distances in the atom pair fingerprints
 (github issue 251)
- added MolOps::get3DDistanceMat() (Chem.Get3DDistanceMatrix() from python)


## New Database Cartridge Features:
- Support configuration of fingerprint sizes in cartridge.
 (github issue 216)
- Add mol_to_ctab(mol, bool default true) to Postgres cartridge.
 (github pull 230)
- Adds sum formula function to PG cartridge.
 (github pull 232)

## New Java Wrapper Features:

Deprecated modules (to be removed in next release):

Removed modules:
- The CMIM integration (previously available to python in the rdkit.ML.FeatureSelect package)
  has been removed due to license incompatibility.

## Contrib updates:
- Added Contribution to train ChEMBL-based models
 (github pull 213)
- ConformerParser functionality
 (github pull 245)

## Other:
- The Open3DAlign code is considerably faster.
- The SMILES parsing code is faster.
- Fix Bison 3.x incompabtility
 (github pull 226)
- Add Travis support
 (github pull 227)
- port of rdkit.ML  bindings from Python/C API to boost::python
 (github pull 237)
- The code now builds more easily using the Anaconda python distribution's
  conda package manager
 (github pull 247)

# Release_2013.09.2
(Changes relative to Release_2013.09.1)

## Acknowledgements:
Andrew Dalke, JP Ebejer, Daniel Moser, Sereina Riniker, Roger Sayle, Manuel Schwarze, Julia Weber

## Bug Fixes:
- cannot pickle unsanitized molecules
  (github issue 149)
- Problems reading PDB files when locale is DE
  (github issue 170)
- calling RWMol::clear() leaves property dict empty
  (github issue 176)
- zero atom molecule generates exception in MolToSmiles when
  rootedAtAtom is provided
  (github issue 182)
- bond orders not being set when PDB files are read
  (github issue 194)
- GenMACCSKeys() raises an exception with an empty molecule
  (github issue 195)
- left-justified SDF bond topology of "0" raises an unexpected
  exception
  (github issue 196)

# Release_2013.09.1
(Changes relative to Release_2013.06.1)

## Acknowledgements:
James Davidson, JP Ebejer, Nikolas Fechner, Grégori Gerebtzoff, Michal Nowotka, Sereina Riniker, Roger
Sayle, Gianluca Sforna, Matthew Szymkiewicz, Paolo Tosco, Dan Warner,

## IMPORTANT
 - Due to a bug fix in the parameter set, the MolLogP and MolMR
   descriptor calculators now return different values for molecules
   with pyrrole (or pyrrole-like) Ns.

## Bug Fixes:
 - The pymol ShowMol method can now handle molecules with more than
   999 atoms (they are sent via PDB)
 - Various stability improvements to the Pandas integration.
   (github issues 129 and 51)
 - Some RDKit methods require python lists and don't allow passing
   numpy arrays or pandas series directly
   (github issue 119)
 - mol2 parser not setting E/Z flags on double bonds
   (github issue 114)
 - Incorrect angle terms in UFF
   (github issue 105)
 - Problems with stereochemistry flags and PathToSubmol()
   (github issue 103)
 - Bad Crippen atom type for pyrrole H
   (github issue 92)
 - PandasTools tests fail with Pandas v0.12
   (github issue 91)
 - Isotope information not affecting chirality
   (github issue 90)
 - properties are not read from SDF files with V3000 mol blocks.
   (github issue 88)
 - assignStereochemistry does not remove bond wedging that shouldn't be there.
   (github issue 87)
 - Drawing code modifies wedge bonds in reactions
   (github issue 86)
 - Stereochemistry not perceived when parsing CTABs unless sanitization is done.
   (github issue 82)
 - 2D rendering issue for epoxide ( CAS  70951-83-6)
   (github issue 78)
 - PandasTools doctests should be failing, but are not
   (github issue 75)
 - Better handling of radicals to/from mol files
   (github issue 73)
 - Benzothiazolium structure can be parsed from ctab, but the SMILES generated cannot be processed.
   (github issue 72)
 - Chem.MolFromInch hangs on CID 23691477 and CID 23691480
   (github issue 68)
 - Chem.MolFromInchi on CHEMBL104337 leads to segmentation fault
   (github issue 67)
 - "Could not embed molecule." (The Anthony Conundrum)
   (github issue 55)

## New Features:
 - Add fragmentOnBonds() to python wrapper
   (github issue 142)
 - Allow renumbering atoms in a molecule.
   (github issue 140)
 - MMFF94 and MMFF94S support
 - implementation of the Open3DAlign rigid alignment algorithm
 - Support for reading and writing PDB files
 - The python function AllChem.AssignBondOrdersFromTemplate() can be
   used to assign bond orders from a reference molecule to the bonds
   in another molecule. This is helpful for getting bond orders
   correct for PDB ligands.
   (github issue 135)
 - Bond lengths, angles, and torsions can now be queries and adjusted.
   (github issue 132)
 - Implementation of similarity maps
   (github issue 94)
 - Python implementation of the Fraggle similarity algorithm.
   See Jameed Hussain's presentation from the 2013 UGM for details:
   https://github.com/rdkit/UGM_2013/blob/master/Presentations/Hussain.Fraggle.pdf?raw=true
 - SparseIntVects now support -=, +=, /=, and *= with ints from C++
   and Python
 - support \\ in SMILES
   (github issue 136)
 - Support a similarity threshold in DbCLI
   (github issue 134)
 - Support construction molecules from other molecules in the python wrapper
   (github issue 133)
 - support tversky similarity in DbCLI
   (github issue 130)
 - support tversky similarity in cartridge
   (github issue 121)
 - support reading and writing reactionComponentType and reactionComponentNumber from ctabs
   (github issue 118)
 - Add in-place forms of addHs(), removeHs(), and mergeQueryHs()
   (github issue 117)
 - modify MolOps::cleanUp() to support this azide formulation: C-N=N#N
   (github issue 116)
 - Dihedral rotation exposed in python
   (github issue 113)
 - Support for cairocffi (cairo drop-in replacement that plays nicely with virtualenv)
   (github issue 80)
 - Grey color for Hydrogens
   (github issue 97)
 - Improvements to the Dict interface in C++
   (github issue 74)
 - customizable drawing options
   (github issue 71)
 - Add method for setting the chiral flag in mol blocks
   (github issue 64)
 - New descriptors added (Python only for now):
   MaxPartialCharge(),MinPartialCharge(),MaxAbsPartialCharge(),MinAbsPartialCharge(),
   MaxEStateIndex(),MinEStateIndex(),MaxAbsEStateIndex(),MinAbsEStateIndex()

## New Database Cartridge Features:

## New Java Wrapper Features:
 - MMFF support
 - PDB reading and writing
 - Open3DAlign support


Deprecated modules (to be removed in next release):

Removed modules:

## Contrib updates:
 - The MMPA implementation has been updated
   See Jameed Hussain's tutorial from the 2013 UGM for details:
   https://github.com/rdkit/UGM_2013/tree/master/Tutorials/mmpa_tutorial
   [Jameed Hussain]
 - An implementation of Ertl and Schuffenhauer's Synthetic
   Accessibility score is available in Contrib/SA_Score
   [Peter Ertl, Greg Landrum]
 - Command line scripts for the Fraggle similarity algorithm
   See Jameed Hussain's presentation from the 2013 UGM for details:
   https://github.com/rdkit/UGM_2013/blob/master/Presentations/Hussain.Fraggle.pdf?raw=true
   [Jameed Hussain]

## Other:
 - Some of the changes to UFF deviate from the published force
   field. Specifics of the changes, and the reasoning behind them, are
   in Paolo Tosco's 2013 RDKit UGM presentation:
   https://github.com/rdkit/UGM_2013/blob/master/Presentations/Tosco.RDKit_UGM2013.pdf?raw=true
 - Reaction drawing has been improved. Support for reaction drawing
   has been added to the IPython notebook.


# Release_2013.06.1
(Changes relative to Release_2013.03.2)

Administrivia note:
In the course of this release cycle, development was moved over
entirely to github. The sourceforge svn repository no longer contains
an up-to-date version of the code.

## Acknowledgements:
Andrew Dalke, JP Ebejer, Nikolas Fechner, Roger Sayle, Riccardo Vianello,
Yingfeng Wang, Dan Warner

## Bug Fixes:
 - The docs for Descriptors.MolWt are now correct (GitHub #38)
 - Molecules coming from InChi now have the correct molecular
   weight. (GitHub #40)
 - RemoveAtoms() no longer leads to problems in canonical SMILES
   generation when chiral ring atoms are present. (GitHub #42)
 - Atom invariants higher than the number of atoms in the molecule can
   now be provided to the atom pairs and topological torsions
   fingerprinters. (GitHub #43)
 - A typo with the handling of log levels was fixed in the python
   wrapper code for InChI generation. (GitHub #44)
 - Stereochemistry no longer affects canonical SMILES generation if
   non-stereo SMILES is being generated. (GitHub #45)
 - The ExactMolWt of [H+] is no longer zero. (GitHub #56)
 - The MPL canvas now has an addCanvasDashedWedge() method. (GitHub
   #57)
 - RWMol::insertMol() now copies atom coordinates (if
   present). (GitHub #59)
 - The "h" primitive in SMARTS strings now uses the method
   getTotalNumHs(false) instead of getImplicitValence().
   (GitHub #60)
 - bzip2 files now work better with the SDWriter class. (GitHub #63)
 - a crashing bug in InChI generation was fixed. (GitHub #67)

## New Features:
 - Sanitization can now be disabled when calling GetMolFrags() from
   Python (GitHub #39)
 - Bond.GetBondTypeAsDouble() has been added to the python
   wrapper. (GitHub #48)
 - The fmcs code now includes a threshold argument allowing the MCS
   that hits a certain fraction of the input molecules (instead of all
   of them) to be found. The code has also been synced with the most
   recent version of Andrew Dalke's version.
 - Atoms now have a getTotalValence() (GetTotalValence() from Python)
   method. (GitHub #61)
 - R labels from Mol files now can go from 0-99
 - chiral flags in CTABs are now handled on both reading and writing.
   The property "_MolFileChiralFlag" is used.


## New Database Cartridge Features:

## New Java Wrapper Features:
 - {Get,Set}Prop() methods are now available for both Atoms and
   Bonds. (GitHub #32)


Deprecated modules (to be removed in next release):

Removed modules:
 - rdkit.utils.pydoc_local

## Other:
 - the handling of flex/bison output files as dependencies has been
   improved (GitHub #33)
 - the molecule drawing code should now also work with pillow (a fork of
   PIL)
 - the PANDAS integration has been improved.  


# Release_2013.03.2
(Changes relative to Release_2013.03.1)

## Acknowledgements:
Manuel Schwarze

## Bug Fixes:
 - The hashed topological torsion fingerprints generated are now the
   same as in previous rdkit versions. (GitHub issue 25)


# Release_2013.03.1
(Changes relative to Release_2012.12.1)

## IMPORTANT

 - The algorithm for hashing subgraphs used in the RDKit fingerprinter
   has changed. The new default behavior will return different
   fingerprints than previous RDKit versions. This affects usage from
   c++, python, and within the postgresql cartridge. See the "## Other"
   section below for more details.

## Acknowledgements:
Paul Czodrowski, Andrew Dalke, Jan Domanski, Jean-Paul Ebejer, Nikolas
Fechner, Jameed Hussain, Stephan Reiling, Sereina Riniker, Roger
Sayle, Riccardo Vianello

## Bug Fixes:
 - removeBond now updates bond indices (sf.net issue 284)
 - dummy labels are no longer lost when atoms are copied (sf.net issue
   285)
 - more specific BRICS queries now match before less specific ones
   (sf.net issue 287, github issue 1)
 - molAtomMapNumber can now be set from Python (sf.net issue 288)
 - the legend centering for molecular image grids has been improved
   (sf.net issue 289)
 - make install now includes all headers (github issue 2)
 - InChIs generaged after clearing computed properties are now correct
   (github issue 3)
 - Reacting atoms that don't change connectivity no longer lose
   stereochemistry (github issue 4)  
 - Aromatic Si is now accepted (github issue 5)  
 - removeAtom (and deleteSubstructs) now correctly updates stereoAtoms
   (github issue 8)  
 - [cartridge] pg_dump no longer fails when molecules cannot be
    converted to SMILES (github issue 9)  
 - a canonicalization bug in MolFragmentToSmiles was fixed (github issue 12)  
 - atom labels at the edge of the drawing are no longer cut off (github issue 13)  
 - a bug in query-atom -- query-atom matching was fixed (github issue 15)  
 - calling ChemicalReaction.RunReactants from Python with None
   molecules no longer leads to a seg fault. (github issue 16)  
 - AllChem.ReactionFromSmarts now generates an error message when called
   with an empty string.
 - Writing CTABs now includes information about atom aliases.
 - An error in the example fdef file
   $RDBASE/Contrib/M_Kossner/BaseFeatures_DIP2_NoMicrospecies.fdef
   has been fixed. (github issue 17)
 - Quantize.FindVarMultQuantBounds() no longer generates a seg fault
   when called with bad arguments. (github issue 18)
 - The length of SDMolSuppliers constructed from empty files is no
   longer reported as 1. (github issue 19)
 - Partial charge calculations now work for B, Si, Be, Mg, and Al.
   (github issue 20)
 - Two logging problems were fixed (github issues 21 and 24)
 - Molecules that have had kappa descriptors generated can now be
   written to SD files (github issue 23)

## New Features:
 - The handling of chirality in reactions has been reworked and
   improved. Please see the RDKit Book for an explanation.
 - Atom-pair and topological-torsion fingerprints now support the
   inclusion of chirality in the atom invariants.
 - A number of new compositional descriptors have been added:
   calcFractionCSP3, calcNum{Aromatic,Aliphatic,Saturated}Rings,
   calcNum{Aromatic,Aliphatic,Saturated}Heterocycles,
   calcNum{Aromatic,Aliphatic,Saturated}Carbocycles
 - An implementation of the molecular quantum number (MQN) descriptors
   has been added.
 - RDKFingerprintMol now takes an optional atomBits argument which is
   used to return information about which bits atoms are involved in.
 - LayeredFingerprintMol no longer takes the arguments tgtDensity and
   minSize. They were not being used.
 - LayeredFingerprintMol2 has been renamed to PatternFingerprintMol
 - The substructure matcher can now properly take stereochemistry into
   account if the useChirality flag is provided.
 - The module rdkit.Chem.Draw.mplCanvas has been added back to svn.
 - A new module integrating the RDKit with Pandas (rdkit.Chem.PandasTools)
   has been added.

## New Database Cartridge Features:
 - The new compositional descriptors are available:
   calcFractionCSP3, calcNum{Aromatic,Aliphatic,Saturated}Rings,
   calcNum{Aromatic,Aliphatic,Saturated}Heterocycles,
   calcNum{Aromatic,Aliphatic,Saturated}Carbocycles
 - MACCS fingerprints are available
 - the substruct_count function is now available
 - substructure indexing has improved. NOTE: indexes on molecule
    columns will need to be rebuilt.

## New Java Wrapper Features:
 - The new compositional descriptors are available:
   calcFractionCSP3, calcNum{Aromatic,Aliphatic,Saturated}Rings,
   calcNum{Aromatic,Aliphatic,Saturated}Heterocycles,
   calcNum{Aromatic,Aliphatic,Saturated}Carbocycles
 - The molecular quantum number (MQN) descriptors are available
 - MACCS fingerprints are available
 - BRICS decomposition is available.

Deprecated modules (to be removed in next release):

Removed modules:

## Other:
 - RDKit fingerprint generation is now faster. The hashing algorithm
   used in the RDKit fingerprinter has changed.
 - Force-field calculations are substantially faster (sf.net issue 290)
 - The core of the BRICS implementation has been moved into C++.
 - The MACCS fingerprint implementation has been moved into
   C++. (contribution from Roger Sayle)
 - New documentation has been added: Cartridge.rst, Overview.rst,
   Install.rst

# Release_2012.12.1
(Changes relative to Release_2012.09.1)

## IMPORTANT

## Acknowledgements:
Andrew Dalke, James Davidson, Robert Feinstein, Nikolas Fechner,
Nicholas Firth, Markus Hartenfeller, Jameed Hussain, Thorsten Meinl,
Sereina Riniker, Roger Sayle, Gianluca Sforna, Pat Walters, Bernd
Wiswedel

## Bug Fixes:
 - Using parentheses for zero-level grouping now works in reaction
   SMARTS. This allows intramolecular reactions to be expressed.
 - SMILES generated for molecules with ring stereochemistry
   (e.g. N[C@H]1CC[C@H](CC1)C(O)=O) are now canonical. (issue 40)
 - SKP lines in a CTAB property block are now properly handled. (issue
   255)
 - The molecular drawing code now shows dotted lines for Any bonds.
   (issue 260)
 - ROMol::debugMol() (ROMol.DebugMol() in Python) now reports isotope
   information. (issue 261)
 - The molecular drawing code now correctly highlights wedged bonds.
   (issue 262)
 - RWMol::addAtom() now adds atoms to conformers.
   (issue 264)
 - TDT files with atomic coordinates now have those coordinates in the
   correct order. (issue 265)
 - A ring-finding error/crash has been fixed. (issue 266)
 - Dummy atoms now have a default valence of 0 and no maximim
   valence. (issue 267)
 - The Python code no longer throws string exceptions. (issue 268)
 - Invalid/unrecognized atom symbols in CTABs are no longer
   accepted. (issue 269)
 - Chem.RDKFingerprint now accepts atom invariants with values larger
   than the number of atoms. (issue 270)
 - The code should now all work when the locale (LANG) is set to
   values other than "C" or one of the English locales. (issue 271)
 - Two-coordinate Hs are no longer removed by
   MolOps::removeHs(). (issue 272)
 - R groups read from CTABs are now marked using setIsotope() instead
   of setMass(). (issue 273)
 - Hs present in the molecule graph no longer incorrectly impact
   substructure matches. (issue 274)
 - Murcko decomposition of molecules with chiral ring atoms now
   works. (issue 275)  
 - Methane now shows up in molecular drawings. (issue 276)
 - '&' in SLN properties is now correctly handled. (issue 277)
 - Molecules with string-valued molAtomMapNumber atomic properties can
   now be serialized. (issue 280)
 - SMARTS strings containing a dot in a recursive piece are now
   properly parsed. (issue 281)
 - The SMILES and SLN parsers no longer leak memory when sanitization
   of the result molecule fails. (issue 282)
 - The cairo canvas drawing code now works with PIL v1.1.6 as well as
   more recent versions.

## New Features:
 - RDKit ExplicitBitVects and DiscreteValueVects can now be directly
   converted into numpy arrays.
 - Rogot-Goldberg similarity has been added.
 - C++: BitVects and SparseIntVects now support a size() method.
 - C++: DiscreteValueVects now support operator[].
 - An initial version of a SWIG wrapper for C# has been added.
 - Support for easily adding recursive queries to molecules and
   reactions has been added. More documentation is required for this
   feature.
 - To allow more control over the reaction, it is possible to flag reactant
   atoms as being protected by setting the "_protected" property on those
   atoms. Flagged atoms will not be altered in the reaction.
 - Atoms and Bonds now support a ClearProp() method from python.
 - The new Python module rdkit.ML.Scoring.Scoring includes a number of
   standard tools for evaluating virtual screening experiments: ROC
   curve generation, AUC, RIE, BEDROC, and Enrichment.
 - The function RDKit::Descriptors::getCrippenAtomContribs()
   (rdkit.Chem.rdMolDescriptors._CalcCrippenContribs() from Python)
   can now optionally return atom-type information as ints or text.


## New Database Cartridge Features:
- The Chi and Kappa descriptors are now available

## New Java Wrapper Features:
- The Chi and Kappa descriptors are now available

Deprecated modules (to be removed in next release):

Removed modules:
- The old SWIG wrapper code in $RDBASE/Code/Demos/SWIG has been
  removed. The SWIG wrappers are now in $RDBASE/Code/JavaWrappers

## Other:
- The C++ code for drawing molecules previously found in
  $RDBASE/Code/Demos/RDKit/Draw has been moved to
  $RDBASE/Code/GraphMol/MolDrawing
- Calculation of the Chi and Kappa descriptors has been moved into
  C++.
- To make builds easier, the thread-safety of the recursive-smarts
  matcher has been made optional. The build option is
  RDK_BUILD_THREADSAFE_SSS.
- There are two new entries in the Contrib directory:
  * Contrib/PBF : An implementation of the Plane of Best Fit
    contributed by Nicholas Firth.
  * Contrib/mmpa : An implementation of GSK's matched molecular pairs
    algorithm contributed by Jameed Hussain
- A new "Cookbook" has been added to the documentation to provide
  a collection of recipes for how to do useful tasks.


# Release_2012.09.1
(Changes relative to Release_2012.06.1)

## IMPORTANT
 - Some of the bug fixes affect the generation of SMILES. Canonical
   SMILES generated with this version of the RDKit will be different
   from previous versions.
 - The fix to Issue 252 (see below) will lead to changes in calculated
   logP and MR values for some compounds.
 - The fix to Issue 254 (see below) will lead to changes in some
   descriptors and geometries for sulfur-containing compounds.
 - The fix to Issue 256 (see below) has changed the name of the
   optional argument to mol.GetNumAtoms from onlyHeavy to
   onlyExplicit. For compatibility reasons, Python code that uses
   explicitly uses onlyHeavy will still work, but it will generate
   warnings. This compatibility will be removed in a future release.

## Acknowledgements:
Gianpaolo Bravi, David Cosgrove, Andrew Dalke, Fabian Dey, James
Davidson, JP Ebejer, Gabriele Menna, Stephan Reiling, Roger Sayle,
James Swetnam

## Bug Fixes:
- The molecules that come from mergeQueryHs() now reset the RingInfo
  structure. (issue 245)
- The output from MurckoScaffold.MakeScaffoldGeneric no longer
  includes stereochemistry or explicit Hs. (issue 246)
- D and T atoms in CTABs now have their isotope information
  set. (issue 247)
- Some problems with ring finding in large, complex molecules have
  been fixed. (issue 249)
- The "rootedAtAtom" argument for FindAllSubgraphsOfLengthN is now
  handled properly. (issue 250)
- Bonds now have a SetProp() method available in Python. (issue 251)
- A number of problems with the Crippen atom parameters have been
  fixed. (issue 252)
- Ring closure digits are no longer repeated on the same atom in
  SMILES generated by the RDKit. (issue 253)
- Non-ring sulfur atoms adjacent to aromatic atoms are no longer set
  to be SP2 hybridized. This allows them to be stereogenic. (issue
  254)
- The combineMols() function now clears computed properties on the
  result molecule.
- A couple of problems with the pickling functions on big endian
  hardware were fixed.
- The molecule drawing code now uses isotope information
- Superscript/Subscript handling in the agg canvas has been improved.
- SKP lines in CTABS are now propertly handled. (Issue 255)
- The name of the optional argument to mol.GetNumAtoms has been
  changed from onlyHeavy to onlyExplicit. The method counts the number
  of atoms in the molecular graph, not the number of heavy
  atoms. These numbers happen to usually be the same (which is why
  this has taken so long to show up), but there are exceptions if Hs
  or dummy atoms are in the graph. (Issue 256)
- Unknown bonds in SMILES are now output using '~' instead of '?'. The
  SMILES parser now recognizes '~' as an "any bond" query. (Issue 257)
- Lines containing only white space in SDF property blocks are no
  longer treated as field separators.
- Transition metals and lanthanides no longer have default valences
  assigned.

## New Features:
- The RDKit now has a maximum common substructure (MCS) implementation
  contributed by Andrew Dalke. This is currently implemented in Python
  and is available as: from rdkit.Chem import MCS Documentation is
  available as a docstring for the function MCS.FindMCS and in the
  GettingStarted document.
- A few new functions have been added to rdkit.Chem.Draw:
  MolsToImage(), MolsToGridImage(), ReactionToImage()
- CalcMolFormula() now provides the option to include isotope
  information.
- The RDKit and Layered fingerprinters both now accept "fromAtoms"
  arguments that can be used to limit which atoms contribute to the
  fingerprint.
- Version information is now available in the Java wrapper.
- The descriptor NumRadicalElectrons is now available.
- The PyMol interface now supports a GetPNG() method which returns the
  current contents of the viewer window as an PIL Image object.
- Molecules (ROMol in C++, rdkit.Chem.Mol in Python) now have a
  getNumHeavyAtoms() method.
- Component-level grouping (parens) can be used in reaction SMARTS.


## New Database Cartridge Features:
- support for molecule <-> pickle conversion via the functions
  mol_to_pkl, mol_from_pkl, and is_valid_mol_pkl.
- support for bit vector <-> binary text conversion via the functions
  bfp_to_binary_text, bfp_from_binary_text

## New Java Wrapper Features:

Deprecated modules (to be removed in next release):

Removed modules:

## Other:
- During this release cycle, the sourceforge project was updated to
  their new hosting system. This explains the change in bug/issue
  ids.
- the SMILES parser is now substantially faster.
- The molecular drawings generated by Code/Demo/RDKit/Draw/MolDrawing.h
  have been improved.
- There is now demo code availble for using the C++ drawing code
  within Qt applications. (contributed by David Cosgrove)
- The directory $RDBASE/Regress now contains sample data and
  scripts for benchmarking the database cartridge.
- Fused-ring aromaticity is now only considered in rings of up to size
  24.
- It is no longer necessary to have flex and bison installed in order
  to build the RDKit.


# Release_2012.06.1
(Changes relative to Release_2012.03.1)

## IMPORTANT
 - Some of the bug fixes affect the generation of SMILES. Canonical
   SMILES generated with this version of the RDKit will be different
   from previous versions.

## Acknowledgements:
Andrew Dalke, JP Ebejer, Igor Filippov, Peter Gedeck, Jan Holst
Jensen, Adrian Jasiński, George Papadatos, Andrey Paramonov, Adrian
Schreyer, James Swetnam

## Bug Fixes:
 - Radicals are now indicated in molecular depictions. (Issue 3516995)
 - Calling .next() on an SDMolSupplier at eof no longer results in an
   infinite loop. (Issue 3524949)
 - Chirality perception no longer fails in large molecules.
   (Issue 3524984)
 - problem creating molblock for atom with four chiral nbrs
   (Issue 3525000)
 - A second sanitization leads to a different molecule.
   (Issue 3525076)
 - can't parse Rf atom in SMILES
   (Issue 3525668)
 - generates [HH2-] but can't parse it
   (Issue 3525669)
 - improper (re)perception of 1H-phosphole
   (Issue 3525671)
 - ForwardSDMolSupplier not skipping forward on some errors
   (Issue 3525673)
 - SMILES/SMARTS parsers don't recognize 0 atom maps
   (Issue 3525776)
 - R group handling in SMILES
   (Issue 3525799)
 - Canonical smiles failure in symmetric heterocycles
   (Issue 3526810)
 - Canonical smiles failure with "extreme" isotopes
   (Issue 3526814)
 - Canonical smiles failure with many symmetric fragments
   (Issue 3526815)
 - Canonical smiles failure with dependent double bonds
   (Issue 3526831)
 - Build Fails Due to Missing include in Code/RDBoost/Wrap.h
   (Issue 3527061)
 - Incorrect template parameter use in std::make_pair
   (Issue 3528136)
 - Canonicalization failure in cycle
   (Issue 3528556)
 - incorrect values reported in ML analysis
   (Issue 3528817)
 - Cartridge does not work on 32bit ubuntu 12.04
   (Issue 3531232)
 - Murcko Decomposition generates unuseable molecule.
   (Issue 3537675)
 - A few memory leaks were fixed in the Java Wrappers
 - The exact mass of molecules with non-standard isotopes is now
   calculated correctly.
 - The default (Euclidean) distance metric should now work with Butina
   clustering.
 - Some bugs in the depictor were fixed.
 - AvalonTools bug with coordinate generation for mols with no
   conformers fixed.

## New Features:
 - ChemicalFeatures now support an optional id
 - Isotope handling has been greatly improved. Atoms now have a
   getIsotope() (GetIsotope() in Python) method that returns zero if
   no isotope has been set, the isotope number otherwise.
 - The function MolFragmentToSmiles can be used to generate canonical
   SMILES for pieces of molecules.
 - The function getHashedMorganFingerprint (GetHashedMorganFingerprint
   in Python) has been added.

## New Database Cartridge Features:
 - The functions mol_from_smiles(), mol_from_smarts(), and
   mol_from_ctab() now return a null value instead of generating an
   error when the molecule processing fails. This allows molecule
   tables to be constructed faster.
 - The functions mol_to_smiles() and mol_to_smarts() have been added.
 - Creating gist indices on bit-vector fingerprint columns is faster.
 - The indexing fingerprint for molecular substructures has been changed.
   The new fingerprint is a bit slower to generate, but is
   considerably better at screening. More information here:
   http://code.google.com/p/rdkit/wiki/ImprovingTheSubstructureFingerprint

## New Java Wrapper Features:

Deprecated modules (to be removed in next release):
 - Support for older (pre9.1) postgresql versions.

Removed modules:
 - rdkit.Excel
 - the code in $RDBASE/Code/PgSQL/RDLib
 - rdkit.Chem.AvailDescriptors : the same functionality is now available
   in a more useable manner from rdkit.Chem.Descriptors

## Other:
 - Similarity calculations on ExplicitBitVectors should now be much faster
 - Use of [Xa], [Xb], etc. for dummy atoms in SMILES is no longer
   possible. Use the "*" notation and either isotopes (i.e. [1*],
   [2*]) or atom maps (i.e. [*:1], [*:2]) instead.
 - Initial work was done towards make the RDKit work on big endian
   hardware (mainly changes to the way pickles are handled)
 - Canonical SMILES generation is now substantially faster.

# Release_2012.03.1
(Changes relative to Release_2011.12.1)

## IMPORTANT
 - The atom-atom match behavior for non-query atoms has been changed.
    This affects the results of doing substructure matches using
    query molecules that are not constructed from SMARTS.

## Acknowledgements:
JP Ebejer, Paul Emsley, Roger Sayle, Adrian Schreyer, Gianluca Sforna,
Riccardo Vianello

## Bug Fixes:
- the older form of group evaluations in Mol blocks is now correctly
  parsed. (Issue 3477283)
- some problems with handling aromatic boron were fixed. (Issue 3480481)
- the SD writer no longer adds an extra $$$$ when molecule parsing
  fails (Issue 3480790)
- molecules in SD files that don't contain atoms are now parsed
  without warnings and their properties are read in. (Issue 3482695)
- it's now possible to embed molecules despite failures in the triangle
  bounds checking (Issue 3483968)
- Isotope information in Mol blocks is now written to M ISO lines
  instead of going in the atom block. (Issue 3494552)
- Better 2D coordinates are now generated for neighbors of atoms with
  unspecified hybridization. (Issue 3487469)
- Dummy atoms and query atoms are now assigned UNSPECIFIED hybridization
  instead of SP. (Issue 3487473)
- Error reporting for SMARTS involving recursion has been improved.
  (Issue 3495296)
- Some problems of queries and generating SMARTS for queries were resolved.
  (Issues 3496759, 3496799, 3496800)
- It's now possible to do database queries with SMARTS that use the index.
  (Issue 3493156).
- A series of problems related to thread safety were fixed.
- Tracking the lifetime of owning molecules across the C++/Python
  border is now being handled better (Issue 3510149)
- A bug with ring-finding in some complex fused ring systems was fixed.
  (Issue 3514824)
- The AllChem module now imports successfully even if the SLN parser
hasn't been built.

## New Features:
- The molecular sanitization is now configurable using an optional
  command-line argument.
- It's now possible to get information from the sanitization routine
  about which operation failed.
- Suppliers support GetLastItemText()
- ComputeDihedralAngle() and ComputeSignedDihedralAngle() were added
  to the rdkit.Geometry module.
- computeSignedDihedralAngle() was added to the C++ API
- ChemicalReactions now support a GetReactingAtoms() method
- the Mol file and Mol block parsers, as well as the SD suppliers,
  now support an optional "strictParsing" argument.
  When this is set to False, problems in the structure of the
  input file are ignored when possible
- EditableMols return the index of the atom/bond added by AddAtom/AddBond
- rdkit.Chem.Draw.MolToImage() now supports an optional "legend" argument
- The MolToSmiles function now supports an optional "allBondsExplicit" argument.

## New Database Cartridge Features:
- the functions mol_from_smiles() and mol_from_smarts() were added

## New Java Wrapper Features:
- the diversity picker now supports an optional random-number seed

Deprecated modules (to be removed in next release):
- rdkit.Excel

Removed modules:
- rdkit.ML.Descriptors.DescriptorsCOM
- rdkit.ML.Composite.CompositeCOM

## Other:
- Assigning/cleaning up stereochemistry is now considerably
faster. This makes standard molecule construction faster.


# Release_2011.12.1
(Changes relative to Release_2011.09.1)

## IMPORTANT
 - The functions for creating bit vector fingerprints using atom pairs
   and topological torsions have been changed. The new default
   behavior will return different fingerprints than previous RDKit
   versions. This affects usage from c++, python, and within the
   postgresql cartridge. See the "## Other" section below for more
   details.
 - Due to a bug fix in the parameter set, the MolLogP and MolMR
   descriptor calculators now return different values for some
   molecules. See the "## Bug Fixes" section below for more details.
 - To make storage more efficient, the size of the fingerprint
   used to store morgan fingerprints in the database cartridge
   has been changed from 1024 bits to 512 bits. If you update
   the cartridge version all morgan and featmorgan fingerprints
   and indices will need to be re-generated.

## Acknowledgements:
Andrew Dalke, JP Ebejer, Roger Sayle, Adrian Schreyer, Gianluca
Sforna, Riccardo Vianello, Toby Wright

## Bug Fixes:
- molecules with polymeric S group information are now rejected by the
  Mol file parser. (Issue 3432136)
- A bad atom type definition and a bad smarts definition were fixed in
  $RDBASE/Data/Crippen.txt. This affects the values returned by the
  logp and MR calculators. (Issue 3433771)
- Unused atom-map numbers in reaction products now produce warnings
  instead of errors. (Issue 3434271)
- rdMolDescriptors.GetHashedAtomPairFingerprint() now works. (Issue
  3441641)
- ReplaceSubstructs() now copies input molecule conformations to the
  output molecule. (Issue 3453144)
- three-coordinate S and Se are now stereogenic (i.e. the
  stereochemistry of O=[S@](C)F is no longer ignored). (Issue 3453172)

## New Features:
- Integration with the new IPython graphical canvas has been
  added. For details see this wiki page:
http://code.google.com/p/rdkit/wiki/IPythonIntegration
- Input and output from Andrew Dalke's FPS format
  (http://code.google.com/p/chem-fingerprints/wiki/FPS) for
  fingerprints.
- The descriptor CalcNumAmideBonds() was added.

## New Database Cartridge Features:
- Support for PostgreSQL v9.1
- Integration with PostgreSQL's KNN-GIST functionality. (Thanks to
  Adrian Schreyer)
- the functions all_values_gt(sfp,N) and all_values_lt(sfp,N) were
  added.

## New Java Wrapper Features:
- A function for doing diversity picking using fingerprint similarity.
- support for the Avalon Toolkit (see below)

Deprecated modules (to be removed in next release):
- rdkit.Excel
- rdkit.ML.Descriptors.DescriptorsCOM
- rdkit.ML.Composite.CompositeCOM

Removed modules:
- rdkit.WebUtils
- rdkit.Reports
- rdkit.mixins

## Other:
- Improvements to the SMARTS parser (Roger Sayle)
- The atom-pair and topological-torsion fingerprinting functions that
  return bit vectors now simulate counts by setting multiple bits in
  the fingerprint per atom-pair/torsion. The number of bits used is
  controlled by the nBitsPerEntry argument, which now defaults to 4.
  The new default behavior does a much better job of reproducing the
  similarities calculated using count-based fingerprints: 95% of
  calculated similarities are within 0.09 of the count-based value
  compared with 0.22 or 0.17 for torsions and atom-pairs previously.
  To get the old behavior, set nBitsPerEntry to 1.
- Optional support has been added for the Avalon Toolkit
(https://sourceforge.net/projects/avalontoolkit/) to provide an
alternate smiles canonicalization, fingerprint, and 2D coordination
generation algorithm.
- The SLN support can now be switched off using the cmake variable
RDK_BUILD_SLN_SUPPORT.
- There are now instructions for building the RDKit and the SWIG
wrappers in 64bit mode on windows.

# Release_2011.09.1
(Changes relative to Release_2011.06.1)

## IMPORTANT
 - A bug in the definition of the Lipinski HBD descriptor was fixed in
   this release. The descriptor Lipinski.NHOHCount will return
   different values for molecules containing Ns or Os with more than
   one attached H.

## Acknowledgements:
Eddie Cao, Richard Cooper, Paul Czodrowski, James Davidson, George
Papadatos, Riccardo Vianello  

## Bug Fixes:
 - A problem with interpretation of stereochemistry from mol files was
   fixed (Issue 3374639)
 - Sterochemistry information for exocyclic double bonds in mol blocks
   is no longer lost. (Issue 3375647)
 - linear double bonds from mol files now have their stereochemistry
   set correctly(Issue 3375684)
 - Chirality for phosphates and sulfates is not longer automatically
   removed. (Issue 3376319)
 - A bug with the reading of query information from mol files was
   fixed. (Issue 3392107)
 - Sterochemistry is now cleaned up after processing mol2
   files. (Issue 3399798)
 - mergeQueryHs now correctly handles atoms with multiple Hs (Issue
   3415204)
 - mergeQueryHs now correctly handles atoms without initial query
   information (Issue 3415206)
 - the calcLipinskiHBD() (equivalent to Lipinski.NHOHCount) descriptor
   now correctly handles Ns or Os with multiple Hs. (Issue 3415534)
 - Morgan fingerprints generated using the fromAtoms argument now have
   all bits from those atoms set.(Issue 3415636)
 - A problem with the way MolSuppliers handle the EOF condition when
   built with the most recent versions of g++ was fixed.
 - Translation of RDKit stereochemistry information into InChI
   stereochemistry information is improved.

## New Features:

## New Database Cartridge Features:
 - molecules can now be built from mol blocks using the function
   mol_from_ctab(). The corresponding is_valid_ctab() function was
   also added.
 - the radius argument is now optional for the functions morganbv_fp,
   morgan_fp, featmorganbv_fp, and featmorgan_fp. The default radius
   for all four functions is 2.

Deprecated modules (to be removed in next release):

Removed modules:

## Other:
 - The documentation in $RDBASE/Docs/Book has been migrated to use
   Sphinx instead of OpenOffice.
 - The optional InChI support can now be built using a system
   installation of the InChI library.



# Release_2011.06.1
(Changes relative to Release_2011.03.2)

## Acknowledgements:
 - Eddie Cao, Andrew Dalke, James Davidson, JP Ebejer, Gianluca
   Sforna, Riccardo Vianello, Bernd Wiswedel

## Bug Fixes:
 - A problem with similarity values between SparseIntVects that
   contain negative values was fixed. (Issue 3295215)
 - An edge case in SmilesMolSupplier.GetItemText() was fixed. (Issue
   3299878)
 - The drawing code now uses dashed lines for aromatic bonds without
   kekulization. (Issue 3304375)
 - AllChem.ConstrainedEmbed works again. (Issue 3305420)  
 - atomic RGP values from mol files are accessible from python (Issue
   3313539)
 - M RGP blocks are now written to mol files. (Issue 3313540)
 - Atom.GetSymbol() for R atoms read from mol files is now
   correct. (Issue 3316600)
 - The handling of isotope specifications is more robust.
 - A thread-safety problem in SmilesWrite::GetAtomSmiles() was fixed.
 - some of the MACCS keys definitions have been corrected
 - Atoms with radical counts > 2 are no longer always written to CTABs
   with a RAD value of 3. (Issue 3359739)

## New Features:
 - The smiles, smarts, and reaction smarts parsers all now take an additional
   argument, "replacements", that carries out string substitutions pre-parsing.
 - There is now optional support for generating InChI codes and keys
   for molecules.
 - the atom pair and topological torsion fingerprint generators now
   take an optional "ignoreAtoms" argument
 - a function to calculate exact molecular weight was added.
 - new java wrappers are now available in $RDBASE/Code/JavaWrappers
 - the methods getMostCommonIsotope() and getMostCommonIsotopeMass()
   have been added to the PeriodicTable class.

## New Database Cartridge Features:
 - Support for generating InChIs and InChI keys
   (if the RDKit InChI support is enabled).

Deprecated modules (to be removed in next release):
 - The original SWIG wrappers in $RDBASE/Code/Demos/SWIG are deprecated

Removed modules:

## Other:
 - The quality of the drawings produced by both the python molecule drawing
   code and $RDBASE/Code/Demos/RDKit/Draw is better.
 - the python molecule drawing code will now use superscripts and
   subscripts appropriately when using the aggdraw or cairo canvases
   (cairo canvas requires pango for this to work).
 - $RDBASE/Code/Demos/RDKit/Draw now includes an example using cairo
 - A lot of compiler warnings were cleaned up.
 - The error reporting in the SMILES, SMARTS, and SLN parsers was improved.
 - the code for calculating molecular formula is now in C++
   (Descriptors::calcMolFormula())


# Release_2011.03.2
(Changes relative to Release_2011.03.1)

## Bug Fixes:
 - A problem in the refactored drawing code that caused the
   rdkit.Chem.Draw functionality to not work at all was fixed.


# Release_2011.03.1
(Changes relative to Release_2010.12.1)

## Acknowledgements:
 - Eddie Cao, James Davidson, Kirk DeLisle, Peter Gedeck, George
   Magoon, TJ O'Donnell, Gianluca Sforna, Nik Stiefl, Bernd Wiswedel

## Bug Fixes:
 - The performance of SSSR finding for molecules with multiple highly-fused
   ring systems has been improved. (Issue 3185548)
 - Isotope information is now correctly saved when molecules are
   serialized (pickled). (Issue 3205280)
 - Generating SMILES for a molecule no longer changes the
   molecule. This fixes a round-trip bug with certain highly symmetric
   molecules read from SD files. (Issue 3228150)
 - Another bounds-matrix generation bug for highly (con)strained
   systems was fixed. (Issue 3238580)
 - Conformation information is now better handled by deleteSubstructs(),
   replaceSubstructs(), and replaceCore().

## New Features:
 - the rdkit.Chem.Draw package has been significantly refactored.
 - Code for doing Murcko decomposition of molecules has been
   added. From Python this is in the module:
   rdkit.Chem.Scaffolds.MurckoScaffold
   It's available in C++ in the GraphMol/ChemTransforms area.
 - rdkit.Chem.AllChem.TransformMol() now takes optional arguments
   allowing the conformation to be transformed to be specified and
   other existing conformations to be preserved.
 - Calculations for most of the descriptors in rdkit.Chem.Lipinski and
   rdkit.Chem.MolSurf have been moved into C++. The python API is the
   same, but the calculations should be somewhat faster.
 - Extensive feature additions to the SWIG-based java wrapper.
 - The Chem.ReplaceCore() function is now better suited for use
   in R-group decomposition.
 - The Morgan fingerprinting code can now return information about
   which atoms set particular bits.
 - The function pathToSubmol() now copies coordinate information
   from conformations (if present). The function is also now available
   from Python
 - The path and subgraph finding code now takes an optional
   rootedAtAtom argument to allow only paths/subgraphs that start at a
   particular atom to be generated.
 - The function findAtomEnvironmentOfRadiusN has been added to allow
   circular atom environments to be located in molecules.
 - MolOps::assignStereochemistry now can also flag potential
   stereocenters that are not specified.

## New Database Cartridge Features:
 - the descriptor-calculation functions mol_numrotatablebonds(),
   mol_numheteroatoms(), mol_numrings(), and mol_tpsa() have been
   added.

Deprecated modules (to be removed in next release):

Removed modules:

## Other:
 - In C++, the functions CalcCrippenDescriptors and CalcAMW have been
   renamed calcCrippenDescriptors and calcAMW to make them consistent
   with the other descriptor calculators.
 - The molecule serialization (pickling) format has been changed. The
   new format is more compact.



# Release_2010.12.1
(Changes relative to Release_2010.09.1)

## IMPORTANT
 - Due to changes made to the fingerprinting code, RDKit and layered
   fingerprints generated with this release are not compatible with
   those from previous releases. For users of the database cartridge:
   you will need to re-generate RDKit fingerprint columns and any
   indices on molecule tables.

## Acknowledgements:
 - Eddie Cao, Andrew Dalke, James Davidson, Kirk DeLisle, Peter Gedeck,
   TJ O'Donnell, Gianluca Sforna, Nik Stiefl, Riccardo Vianello

## Bug Fixes:
 - The depiction code no longer crashes with single-atom templates
   (issue 3122141)
 - Aromatic bonds in the beginning of a SMILES branch are now
   correctly parsed (issue 3127883)
 - A crash when generating 2d constrained coordinates was fixed (issue
   3135833)
 - Stereochemistry no longer removed from double bonds in large
   rings. (issue 3139534)
 - Atom mapping information no longer in reaction products (issue
   3140490)  
 - Smiles parse failure with repeated ring labels and dot disconnects
   fixed (issue 3145697)
 - a bug causing the molecule drawing code to not use the cairo canvas
   when it's installed was fixed
 - the SMILES generated for charged, aromatic Se or Te has been fixed
   (issue 3152751)
 - PropertyMols constructed from pickles and then written to SD files
   will now include the properties in the SD file.
 - SMILES can now be generated correctly for very large molecules
   where more than 50 rings are open at once. (issue 3154028)

## New Features:
 - All molecular descriptor calculators are now pulled in by the
   rdkit.Chem.Descriptors module. So you can do things like:
   Descriptors.MolLogP(mol) or Descriptors.fr_amide(mol)
 - Atom-map numbers in SMILES are now supported. They can be accessed
   as the atomic "molAtomMapNumber" property. (issue 3140494)
 - It's now possible to tell the RDKit to generate non-canonical
   SMILES via an optional argument to MolToSmiles. This is faster than
   generating canonical SMILES, but is primarity intended for
   debugging/testing. (issue 3140495)
 - The function GenerateDepictionMatching2DStructure() has been added
   to the rdkit.Chem.AllChem module to make generating
   template-aligned depictions easier.
 - Generating FCFP-like fingerprints is now more straightforward via
   the useFeatures optional argument to GetMorganFingerprint()
 - Extensive changes were made to the layered fingerprinting code to
   allow better coverage of queries.
 - Functionality for stripping common salts from molecules has been
   added in rdkit.Chem.SaltRemover. The salts themselves are defined
   in $RDBASE/Data/Salts.txt
 - Functionality for recognizing common functional groups has been
   added in rdkit.Chem.FunctionalGroups. The functional groups
   themselves are defined in
   $RDBASE/Data/Functional_Group_Hierarchy.txt

## New Database Cartridge Features:
 - The cartridge now supports SMARTS queries.
 - The functions is_valid_{smiles,smarts}() are now available
   (issue 3097359).
 - The operator @= is now supported for testing molecule equality.
   (issue 3120707)
 - The functions featmorgan_fp() and featmorganbv_fp() are now
   available for generating FCFP-like fingerprints.

Deprecated modules (to be removed in next release):
 - rdkit.Chem.AvailDescriptors : the same functionality is now available
   in a more useable manner from rdkit.Chem.Descriptors (see above).

Removed modules:

## Other:
 - RDKit support has been added to the Knime data mining and reporting
   tool. More information is available from the knime.org community
   site: http://tech.knime.org/community/rdkit
   Thanks to Thorsten, Bernd, Michael, and the rest of the crew at
   knime.com for making this possible.
 - RPMs to allow easy installation of the RDKit on Fedora/CentOS/RHEL
   and similar systems are now available. Thanks to Gianluca Sforna
   for doing this work.
 - The database cartridge now statically links the RDKit libraries.
   This should make installation easier.
 - The RDKit fingerprinter now by default sets 2 bits per hashed
   subgraph instead of 4. The old behavior can be regained by setting
   nBitsPerHash to 4.

# Release_2010.09.1
(Changes relative to Release_Q22010_1)

## IMPORTANT
 - Due to changes made to the layered fingerprinting code,
   fingerprints generated with this release are not compatible with
   fingerprints from earlier releases.
 - The default arguments to the Morgan fingerprinting code will yield
   fingerprints that are not backwards compatible.

## Acknowledgements:
 - Andrew Dalke, James Davidson, Paul Emsley, Peter Gedeck,
   Uwe Hoffmann, Christian Kramer, Markus Kossner, TJ O'Donnell,
   Gianluca Sforna, Nik Stiefl, Riccardo Vianello

## Bug Fixes:
 - A typo in the parameters for the Crippen clogp calculator was
   fixed. (issue 3057201)
 - some problems in the layered fingerprinting code were fixed. (issue
   3030388)
 - a bug in the ring-finding code that could lead to incorrect results
   or crashes in large molecules was fixed.  
 - the Murtagh clustering code should now execute correctly on recent
   versions of the MacOS.
 - some problems with the cairo canvas were fixed
 - a problem with matching non-default isotope SSS queries for molecules
   read in from CTABs was fixed (issue 3073163).
 - a problem with calculating AMW for molecules with non-default isotopes
   was fixed.

## New Features:
 - a PostgreSQL cartridge for similarity and substructure searching
   has been added to the RDKit distribution.
 - The Morgan fingerprinting code accepts additional arguments that
   control whether or not bond order and chirality are taken into
   account. By default chirality is ignored and the bond order is
   used. Another change with the MorganFPs is that ring information is
   now included by default.
 - 2D coordinates can now be generated for chemical reactions.
 - The functions IsMoleculeReactantOfReaction and
   IsMoleculeProductOfReaction have been added to the C++
   interface. From python these are methods of the ChemicalReaction
   class:
   rxn.IsMoleculeReactant and rxn.IsMoleculeProduct
 - The default bond length for depiction can now be changed.
 - FCFP-like fingerprints can now be generated with the Morgan
   fingerprinting code by starting with feature invariants.
 - The close() method has been added to MolWriters.
 - Morgan, atom-pair, and topological-torsion fingerprints can now
   also be calculated as bit vectors.
 - RDKit and layered fingerprints can now be generated using only
   linear paths.
 - the function findAllPathsOfLengthMtoN() was added

Deprecated modules (to be removed in next release):

Removed modules:
 - rdkit/qtGui
 - rdkit/RDToDo
 - Projects/SDView

## Other:
 - As of this release a new version numbering scheme is being used:
   YYYY.MM.minor. An example, this release was done in Sept. of 2010
   so it's v2010.09.1.
 - the RDBASE environment variable is no longer required. It will be
   used if set, but the code should work without it
 - The directory Contrib/M_Kossner contains two new contributions from
   Markus Kossner.
 - A change was made to the subgraph matching code that speeds up
   substructure searches involving repeated recursive queries.
 - the deprecated registerQuery argument has been removed from the
   substructure matching functions.
 - the empty header files AtomProps.h and BondProps.h have been
   removed.
 - in order to simplify the build process the test databases are now
   in svn
 - some python functions to calculate descriptors (i.e. pyMolWt,
   pyMolLogP, etc.) that have C++ equivalents have been removed to
   clean up the interface
 - the PIL canvas should no longer generate warnings
 - Thanks to the help of Gianluca Sforna and Riccardo Vianello, it is
   now much easier to package and distribute the RDKit.
 - the bjam-based build system has been removed.

# Release_Q22010_1
(Changes relative to Release_Q12010_1)

## IMPORTANT
 - There are a couple of refactoring changes that affect people using
   the RDKit from C++. Please look in the ## Other section below for a list.
 - If you are building the RDKit yourself, changes made in this
   release require that you use a reasonably up-to-date version of
   flex to build it. Please look in the ## Other section below for more
   information.

## Acknowledgements:
 - Andrew Dalke, James Davidson, Kirk DeLisle, Thomas Heller, Peter Gedeck,
   Greg Magoon, Noel O'Boyle, Nik Stiefl,  

## Bug Fixes:
 - The depictor no longer generates NaNs for some molecules on
   windows (issue 2995724)
 - [X] query features work correctly with chiral atoms. (issue
   3000399)
 - mols will no longer be deleted by python when atoms/bonds returned
   from mol.Get{Atom,Bond}WithIdx() are still active. (issue 3007178)
 - a problem with force-field construction for five-coordinate atoms
   was fixed. (issue 3009337)
 - double bonds to terminal atoms are no longer marked as "any" bonds
   when writing mol blocks. (issue 3009756)
 - a problem with stereochemistry of double bonds linking rings was
   fixed. (issue 3009836)
 - a problem with R/S assignment was fixed. (issue 3009911)
 - error and warning messages are now properly displayed when cmake
   builds are used on windows.
 - a canonicalization problem with double bonds incident onto aromatic
   rings was fixed. (issue 3018558)
 - a problem with embedding fused small ring systems was fixed.
   (issue 3019283)

## New Features:
 - RXN files can now be written. (issue 3011399)
 - reaction smarts can now be written.
 - v3000 RXN files can now be read. (issue 3009807)
 - better support for query information in mol blocks is present.
   (issue 2942501)
 - Depictions of reactions can now be generated.
 - Morgan fingerprints can now be calculated as bit vectors (as
   opposed to count vectors.
 - the method GetFeatureDefs() has been added to
   MolChemicalFeatureFactory
 - repeated recursive SMARTS queries in a single SMARTS will now be
   recognized and matched much faster.
 - the SMILES and SMARTS parsers can now be run safely in
   multi-threaded code.

Deprecated modules (to be removed in next release):
 - rdkit/qtGui
 - Projects/SDView

Removed modules:
 - SVD code: External/svdlibc External/svdpackc rdkit/PySVD
 - rdkit/Chem/CDXMLWriter.py

## Other:
 - The large scale changes in the handling of stereochemistry were
   made for this release. These should make the code more robust.
 - If you are building the RDKit yourself, changes made in this
   release require that you use a reasonably up-to-date version of
   flex to build it. This is likely to be a problem on Redhat, and
   redhat-derived systems. Specifically: if your version of flex is
   something like 2.5.4 (as opposed to something like 2.5.33, 2.5.34,
   etc.), you will need to get a newer version from
   http://flex.sourceforge.net in order to build the RDKit.

 - Changes only affecting C++ programmers:
   - The code for calculating topological-torsion and atom-pair
     fingerprints has been moved from $RDBASE/Code/GraphMol/Descriptors
     to $RDBASE/Code/GraphMol/Fingerprints.
   - The naming convention for methods of ExplicitBitVect and
     SparseBitVect have been changed to make it more consistent with
     the rest of the RDKit.
   - the bjam-based build system should be considered
     deprecated. This is the last release it will be actively
     maintained.


# Release_Q12010_1
(Changes relative to Release_Q42009_1)

## Acknowledgements:
 - Andrew Dalke, Jean-Marc Nuzillard, Noel O'Boyle, Gianluca Sforna,
   Nik Stiefl, Anna Vulpetti

## Bug Fixes
 - Substantial improvements were made to the SLN parser
 - A bad depiction case was fixed. (issue 2948402)
 - Hs added to planar carbons are no longer in the same plane as the
   other atoms. (issue 2951221)
 - Elements early in the periodic table (e.g. Mg, Na, etc.) no longer
   have their radical counts incorrectly assigned. (issue 2952255)
 - Some improvements were made to the v3k mol file parser. (issue
   2952272)
 - Double bonds with unspecified stereochemistry are now correctly
   flagged when output to mol files. (issue 2963522)
 - A segmentation fault that occured when kekulizing modified
   molecules has been fixed. (issue 2983794)

## New Features
 - The MaxMin diversity picker can now be given a seed for the random
   number generator to ensure reproducible results.

## Other
 - the vflib source, which is no longer used, was removed from the
   External source tree. It's still available in svn at rev1323 or via
   this tarball:
   http://rdkit.svn.sourceforge.net/viewvc/rdkit/trunk/External/vflib-2.0.tar.gz?view=tar&pathrev=1323
 - the directory Contrib has been added to the RDKit distribution to
   house contributions that don't necessarily fit anywhere else. The
   first contribution here is a collection of scripts required to
   implement local-environment fingerprints contributed by Anna
   Vulpetti.
 - Some optimization work was done on the molecule initialization code:
   reading in molecules is now somewhat faster.
 - Some optimization work was done on the RDK and Layered fingerprinting code.

# Release_Q42009_1
(Changes relative to Release_Q32009_1)

## IMPORTANT
  - A bug fix in the SMARTS parser has changed the way atom-map
    numbers in Reaction SMARTS are parsed.
      Earlier versions of the RDKit required that atom maps be
      specified at the beginning of a complex atom query:
        [CH3:1,NH2]>>[*:1]O
      The corrected version only accepts this form:
        [CH3,NH2:1]>>[*:1]O
    This change may break existing SMARTS patterns.
  - A switch to using cmake as the build system instead of bjam has
    made the RDKit much easier to build.

## Acknowledgements
  - Andrew Dalke, Kirk DeLisle, David Hall, Markus Kossner, Adrian
    Schreyer, Nikolaus Stiefl, Jeremy Yang

## Bug Fixes
  - the SMARTS parser now correctly requires tha atom-map numbers be
    at the end of a complex atom query.
    (issue 1804420)
  - a bug in the way SMARTS matches are uniquified has been fixed
    (issue 2884178)

## New Features
  - The new SMARTS atomic query feature "x" (number of ring bonds) is
    now supported.
  - The proof-of-concept for a SWIG-based wrapper around the RDKit has
    been expanded a bit in functionality. Samples are now included for
    Java, C#, and Python.
  - Information about the current RDKit and boost versions is now
    available from C++ (file RDGeneral/versions.h) and Python
    (rdBase.rdkitVersion and rdBase.boostVersion)
  - The KNN code now supports weighted nearest-neighbors calculations
    with a radius cutoff.

## Other
  - The lapack dependency has been completely removed from the RDKit.
  - The supported build system for the RDKit is now cmake
    (http://www.cmake.org) instead of bjam. See the file INSTALL for
    the new installation instructions. Files for bjam are still
    included in the distribution but are deprecated and will be
    removed in a future version.


# Release_Q32009_1
(Changes relative to Release_Q22009_1)

## IMPORTANT
  - Due to bug fixes in the boost random-number generator, RDK
    fingerprints generated with boost 1.40 are not backwards
    compatible with those from earlier versions.

## Acknowledgements
  - Uwe Hoffmann, Nik Stiefl, Greg Magoon, Ari Gold-Parker,
    Akihiro Yokota, Kei Taneishi, Riccardo Vianello, Markus Kossner

## Bug Fixes
  - the canonOrient argument to the depiction code now works
    (issue 2821647)
  - typo in the depictor 2D embedding code fixed  
    (issue 2822883)
  - single aromatic atoms in chains now (correctly) fail sanitization
    (issue 2830244)
  - problem with embedding and fused rings fixed
    (issue 2835784)
  - crash when reading some large molecules fixed
    (issue 2840217)
  - trailing newline in TemplateExpand.py fixed
    (issue 2867325)
  - fingerprint incompatibility on 64bit machines fixed
    (issue 2875658)
  - PropertyMol properties are now written to SD files
    (issue 2880943)

## New Features
  - to the extent possible, reactions now transfer coordinates from
    reactant molecules to product molecules (issue 2832951)
## Other
  - the function DaylightFingerprintMol() has been removed
  - the outdated support for Interbase has been removed
  - the Compute2DCoords() function in Python now canonicalizes the
    orientation of the molecule by default.
  - the distance-geometry code should now generate less bad amide
    conformations. (issue 2819563)
  - the quality of distance-geometry embeddings for substituted- and
    fused-ring systems should be better.  

# Release_Q22009_1
(Changes relative to Release_Q12009_2)

## Acknowledgements
  - Uwe Hoffmann, Marshall Levesque, Armin Widmer

## Bug Fixes
  - handling of crossed bonds in mol files fixed (issue 2804599)
  - serialization bug fixed (issue 2788233)
  - pi systems with 2 electrons now flagged as aromatic (issue 2787221)
  - Chirality swap on AddHs (issue 2762917)
  - core leak in UFFOptimizeMolecule fixed (issue 2757824)

## New Features
  - cairo support in the mol drawing code (from Uwe Hoffmann) (issue 2720611)
  - Tversky and Tanimoto similarities now supported for SparseIntVects
  - AllProbeBitsMatch supported for BitVect-BitVect comparisons
  - ChemicalReactions support serialization (pickling) (issue 2799770)
  - GetAtomPairFingerprint() supports minLength and maxLength arguments
  - GetHashedTopologicalTorsionFingerprint() added
  - preliminary support added for v3K mol files
  - ForwardSDMolSupplier added
  - CompressedSDMolSupplier added (not supported on windows)
  - UFFHasAllMoleculeParams() added
  - substructure searching code now uses an RDKit implementation of
    the vf2 algorithm. It's much faster.
  - Atom.GetPropNames() and Bond.GetPropNames() now available from
    python
  - BRICS code now supports FindBRICSBonds() and BreakBRICSBonds()
  - atom labels Q, A, and * in CTABs are more correctly supported
    (issue 2797708)
  - rdkit.Chem.PropertyMol added (issue 2742959)
  - support has been added for enabling and disabling logs
    (issue 2738020)

## Other
  - A demo has been added for using the MPI with the RDKit
    ($RDBASE/Code/Demos/RDKit/MPI).
  - Embedding code is now better at handling chiral structures and
    should produce results for molecules with atoms that don't have
    UFF parameters.
  - the UFF code is more robust w.r.t. missing parameters
  - GetHashedAtomPairFingerprint() returns SparseIntVect instead of
    ExplicitBitVect
  - the CTAB parser (used for mol files and SD files) is faster
  - extensive changes to the layered fingerprinting code;
    fingerprinting queries is now possible
  - molecule discriminator code moved into $RDBASE/Code/GraphMol/Subgraphs
  - the SDView4 prototype has been expanded
  - $RDBASE/Regress has been added to contain regression and
    benchmarking data and scripts.
  - support for sqlalchemy has been added to $RDBASE/rdkit/Chem/MolDb
  - $RDBASE/Projects/DbCLI/SDSearch.py has been removed; use the
    CreateDb.py and SearchDb.py scripts in the same directory instead.
  - the BRICS code has been refactored  

# Release_Q12009_2
(Changes relative to Release_Q42008_1)

## IMPORTANT

 - The directory structure of the distribution has been changed in
   order to make installation of the RDKit python modules more
   straightforward. Specifically the directory $RDBASE/Python has been
   renamed to $RDBASE/rdkit and the Python code now expects that
   $RDBASE is in your PYTHONPATH. When importing RDKit Python modules,
   one should now do: "from rdkit import Chem" instead of "import
   Chem". Old code will continue to work if you also add $RDBASE/rdkit
   to your PYTHONPATH, but it is strongly suggested that you update
   your scripts to reflect the new organization.
 - For C++ programmers: There is a non-backwards compatible change in
   the way atoms and bonds are stored on molecules. See the *## Other*
   section for details.

## Acknowledgements
 - Kirk DeLisle, Noel O'Boyle, Andrew Dalke, Peter Gedeck, Armin Widmer

## Bug Fixes
 - Incorrect coordinates from mol2 files (issue 2727976)
 - Incorrect handling of 0s as ring closure digits (issues 2525792,
 and 2690982)
 - Incorrect handling of atoms with explicit Hs in reactions (issue 2540021)
 - SmilesMolSupplier.GetItemText() crashes (issue 2632960)
 - Incorrect handling of dot separations in reaction SMARTS (issue 2690530)
 - Bad charge lines in mol blocks for large molecules (issue 2692246)
 - Order dependence in AssignAtomChiralTagsFromStructure (issue 2705543)
 - Order dependence in the 2D pharmacophore code
 - the LayeredFingerprints now handle non-aromatic single ring bonds
   between aromatic atoms correctly.


## New Features
 - BRICS implementation
 - Morgan/circular fingerprints implementation
 - The 2D pharmacophore code now uses standard RDKit fdef files.
 - Atom parity information in CTABs now written and read. If present
   on reading, atom parity flags are stored in the atomic property
   "molParity".
 - An optional "fromAtoms" argument has been added to the atom pairs
   and topological torsion fingerprints. If this is provided, only atom
   pairs including the specified atoms, or torsions that either start
   or end at the specified atoms, will be included in the fingerprint.
 - Kekulization is now optional when generating CTABs. Since the MDL
   spec suggests that aromatic bonds not be used, this is primarily
   intended for debugging purposes.
 - the removeStereochemistry() (RemoveStereoChemistry() from Python)
   function has been added to remove all stereochemical information
   from a molecule.

## Other
 - The Qt3-based GUI functionality in $RDBASE/rdkit/qtGui and
   $RDBASE/Projects/SDView is deprecated. It should still work, but it
   will be removed in a future release. Please do not build anything
   new on this (very old and creaky) framework.
 - The function DaylightFingerprintMol() is now deprecated, use
   RDKFingerprintMol() instead.
 - For C++ programmers: The ROMol methods getAtomPMap() and
   getBondPMap() have been removed. The molecules themselves now support
   an operator[]() method that can be used to convert graph iterators
   (e.g. ROMol:edge_iterator, ROMol::vertex_iterator,
   ROMol::adjacency_iterator) to the corresponding Atoms and Bonds.
   New API for looping over an atom's bonds:
        ... molPtr is a const ROMol * ...
        ... atomPtr is a const Atom * ...
        ROMol::OEDGE_ITER beg,end;
        boost::tie(beg,end) = molPtr->getAtomBonds(atomPtr);
        while(beg!=end){
          const BOND_SPTR bond=(*molPtr)[*beg];
          ... do something with the Bond ...
          ++beg;
        }
  New API for looping over a molecule's atoms:
        ... mol is an ROMol ...
        ROMol::VERTEX_ITER atBegin,atEnd;
        boost::tie(atBegin,atEnd) = mol.getVertices();  
        while(atBegin!=atEnd){
          ATOM_SPTR at2=mol[*atBegin];
          ... do something with the Atom ...
          ++atBegin;
        }

# Release_Q42008_1
(Changes relative to Release_Q32008_1)

## IMPORTANT
 - A fix in the handling of stereochemistry in rings means that the
   SMILES generated with this release are different from those in
   previous releases. Note that the canonicalization algorithm does
   not work in cases of pure ring stereochemistry : the SMILES should
   be correct, but it is not canonical. Rings containing chiral
   centers should be fine.

## Acknowledgements:
 - Kirk DeLisle, Markus Kossner, Greg Magoon, Nik Stiefl

## Bug Fixes
 - core leaks in learning code (issue 2152622)
 - H-bond acceptor definitions (issue 2183240)
 - handling of aromatic dummies (issue 2196817)
 - errors in variable quantization (issue 2202974)
 - errors in information theory functions on 64 bit machines (issue 2202977)
 - kekulization problems (issue 2202977)
 - infinite loop in getShortestPaths() for disconnected structures (issue 2219400)
 - error in depictor for double bonds with stereochemistry connected
   to rings (issue 2303566)
 - aromaticity flags not copied to null atoms in reaction products
   (issue 2308128)
 - aromaticity perception in large molecule hangs (issue 2313979)
 - invariant error in canonicalization (issue 2316677)
 - mol file parser handling of bogus bond orders (issue 2337369)
 - UFF optimization not terminating when atoms are on top of each
   other (issue 2378119)
 - incorrect valence errors with 4 coordinate B- (issue 2381580)
 - incorrect parsing of atom-list queries with high-numbered atoms
   (issue 2413431)
 - MolOps::mergeQueryHs() crashing with non-query molecules. (issue
   2414779)

## New Features
 - SLN parser (request 2136703).
 - Mol2 parser : Corina atom types (request 2136705).
 - Building under mingw (request 2292153).
 - Null bonds in reaction products are replaced with the corresponding
   bond from the reactants (request 2308123).

## Other
 - a bunch of deprecation warnings from numpy have been cleaned up
   (issue 2318431)
 - updated documentation
 - some optimization work on the fingerprinter

# Release_Q32008_1
(Changes relative to Release_May2008_1)

## Acknowledgements:
 - Noel O'Boyle, Igor Filippov, Evgueni Kolossov, Greg Magoon

## Bug Fixes
 - A memory leak in the ToBase64 and FromBase64 wrapper functions was
   fixed.
 - The UFF atom typer has been made more permissive: it now will pick
   "close" atom types for things it does not recognize. (issue
   2094445)
 - The handling of molecules containing radicals has been greatly
   improved (issues 2091839, 2091890, 2093420)
 - Iterative (or secondary, or dependent) chirality is now supported,
   see this page for more information:
   http://code.google.com/p/rdkit/wiki/IterativeChirality
   (issue 1931470)
 - Isotope handling has been changed, this allows correct matching of
   SMARTS with specified isotopes. (issue 1968930)
 - Some problems with the MACCS key definitions have been
   fixed. (issue 2027446)
 - Molecules with multiple fragments can now be correctly
   embedded. (issue 1989539)
 - Adding multiple bonds between the same atoms in a molecule now
   produces an error. (issue 1993296)
 - The chemical reaction code now handles chiral atoms correctly in
   when applying reactions with no stereochem information
   provided. (issue 2050085)
 - A problem with single-atom cores in TemplateExpand.py has been
   fixed. (issue 2091304)  
 - A problem causing bicyclobutane containing molecules to not be
   embeddable has been fixed. (issue 2091864)
 - The default parameters for embedding are now molecule-size
   dependent. This should help with the embedding of large, and
   crowded molecules. (issue 2091974)
 - The codebase can now be built with boost 1.36. (issue 2071168)
 - A problem with serialization of bond directions was fixed.
   (issue 2113433)

## New Features
 - The RDKit can now be built under Darwin (Mac OS/X).
 - Tversky similarity can now be calculated. (request 2015633)
 - Many of the core datastructures now support equality comparison
   (operator==). (request 1997439)
 - Chirality information can now be assigned based on the 3D
   coordinates of a molecule using
   MolOps::assignChiralTypesFrom3D(). (request 1973062)
 - MolOps::getMolFrags() can now return a list of split molecules
   instead of just a list of atom ids. (request 1992648)
 - ROMol::getPropNames() now supports the includePrivate and
   includeComputed options. (request 2047386)


## Other
 - the pointers returned from Base64Encode/Decode are now allocated
   using new instead of malloc or calloc. the memory should be
   released with delete[].
 - the generation of invariants for chirality testing is now quite a
   bit faster; this results in faster parsing of molecules.
 - The use of C include files instead of their C++ replacements has
   been dramatically reduced.
 - The new (as of May2008) hashing algorithm for fingerprints is now
   the default in the python fingerprinting code
   (Chem.Fingerprints.FingerprintMols).
 - The functions MolOps::assignAtomChiralCodes() and
   MolOps::assignBondStereoCodes() are deprecated. Use
   MolOps::assignStereochemistry() instead.
 - The RDKit no longer uses the old numeric python library. It now
   uses numpy, which is actively supported.
 - By default Lapack++ is no longer used. The replacement is the boost
   numeric bindings: http://mathema.tician.de/software/boost-bindings.


# Release_May2008_1
(Changes relative to Release_Jan2008_1)

## IMPORTANT
 - A fix to the values of the parameters for the Crippen LogP
   calculator means that the values calculated with this version are
   not backwards compatible. Old values should be recalculated.
 - topological fingerprints generated with this version *may* not be
   compatible with those from earlier versions. Please read the note
   below in the "## Other" section.
 - Please read the point about dummy atoms in the "## New Features"
   section. It explains a change that affects backwards compatibility
   when dealing with dummy atoms.


## Acknowledgements:
 - Some of the bugs fixed in this release were found and reported by
   Adrian Schreyer, Noel O'Boyle, and Markus Kossner.

## Bug Fixes
 - A core leak in MolAlign::getAlignmentTransform was fixed (issue
   1899787)
 - Mol suppliers now reset the EOF flag on their stream after they run
   off the end (issue 1904170)
 - A problem causing the string "Sc" to not parse correctly in
   recursive SMARTS was fixed (issue 1912895)
 - Combined recursive smarts queries are now output correctly.
   (issue 1914154)
 - A bug in the handling of chirality in reactions was fixed (issue
   1920627)
 - Looping directly over a supplier no longer causes a crash (issue
   1928819)
 - a core leak in the smiles parser was fixed (issue 1929199)
 - Se and Te are now potential aromatic atoms (per the proposed
   OpenSmiles standard). (issue 1932365)
 - isotope information (and other atomic modifiers) are now correctly
   propagated by chemical reactions (issue 1934052)
 - triple bonds no longer contribute 2 electrons to the count for
   aromaticity (issue 1940646)
 - Two bugs connected with square brackets in SMILES were fixed
   (issues 1942220 and 1942657)
 - atoms with coordination numbers higher than 4 now have tetrahedral
   stereochemistry removed (issue 1942656)
 - Bond.SetStereo() is no longer exposed to Python (issue 1944575)
 - A few typos in the parameter data for the Crippen logp calculator
   were fixed. Values calculated with this version should be assumed
   to not be backwards compatible with older versions (issue 1950302)
 - Isotope queries are now added correctly (if perhaps not optimally)
   to SMARTS.
 - some drawing-related bugs have been cleared up.
 - A bug in Chem.WedgeMolBonds (used in the drawing code) that was
   causing incorrect stereochemistry in drawn structures was
   fixed. (issue 1965035)
 - A bug causing errors or crashes on Windows with [r<n>] queries was
   fixed. (issue 1968930)
 - A bug in the calculation of TPSA values in molecules that have Hs
   in the graph was fixed. (issue 1969745)

## New Features
 - Support for supplying dummy atoms as "[Du]", "[X]", "[Xa]", etc. is
   now considered deprecated. In this release a warning will be
   generated for these forms and in the next release the old form will
   generate errors. Note that the output of dummy atoms has also
   changed: the default output format is now "*", this means that the
   canonical SMILES for molecules containing dummies are no longer
   compatible with the canonical SMILES from previous releases.
   (feature request 186217)
 - Atom and bond query information is now serializable; i.e. query
   molecules can now be pickled and not lose the query
   information. (feature request 1756596)
 - Query features from mol files are now fully supported. (feature
   request 1756962)
 - Conformations now support a dimensionality flag. Dimensionality
   information is now read from mol blocks and TDT files. (feature request
   1906758)
 - Bulk Dice similarity functions have been added for IntSparseIntVect
   and LongSparseIntVect (feature request 1936450)
 - Exceptions are no longer thrown during molecule parsing. Failure in
   molecule parsing is indicated by returning None. Failure to *open* a
   file when reading a molecule throws BadFileExceptions (feature
   requests 1932875 and 1938303)
 - The various similarity functions for BitVects and SparseIntVects
   now take an optional returnDistance argument. If this is provided,
   the functions return the corresponding distance instead of
   similarity.
 - Some additional query information from Mol files is now translated
   when generating SMARTS. Additional queries now translated:
     - number of ring bonds
     - unsaturation queries
     - atom lists are handled better as well
   (feature request 1902466)
 - A new algorithm for generating the bits for topological
   fingerprints has been added. The new approach is a bit quicker and
   more robust than the old, but is not backwards compatible.
   Similarity trends are more or less conserved.
 - The molecule drawing code in Chem.Draw.MolDrawing has been modified
   so that it creates better drawings. A new option for drawing that
   uses the aggdraw graphics library has been added.
 - The RingInfo class supports two new methods: AtomRings() and
   BondRings() that return tuples of tuples with indices of the atoms
   or bonds that make up the molecule's rings.

## Other
 - Changes in the underlying boost random-number generator in version
   1.35 of the boost library may have broken backwards compatibility
   of 2D fingerprints generated using the old fingerprinter. It is
   strongly suggested that you regenerate any stored fingerprints (and
   switch to the new fingerprinter if possible). There is an explicit
   test for this in $RDBASE/Code/GraphMol/Fingerprints/test1.cpp
 - The unofficial and very obsolete version of John Torjo's v1
   boost::logging library that was included with the RDKit
   distribution is no longer used. The logging library has been
   replaced with the much less powerful and flexible approach of just
   sending things to stdout or stderr. If and when the logging library
   is accepted into Boost, it will be integrated.
 - The DbCLI tools (in $RDBASE/Projects/DbCLI) generate topological
   fingerprints using both the old and new algorithms (unless the
   --noOldFingerprints option is provided). The default search
   uses the newer fingerprint.
 - The directory $RDBASE/Data/SmartsLib contains a library of sample
   SMARTS contributed by Richard Lewis.


# Release_Jan2008_1
(Changes relative to Release_Aug2007_1)

## IMPORTANT
 - Bug fixes in the canonicalization algorithm have made it so that
   the canonical SMILES from this version are not compatible with
   those from older versions of the RDKit.
 - Please read the point about dummy atoms in the "## New Features"
   section. It explains a forthcoming change that will affect
   backwards compatibility when dealing with dummy atoms.
 - The build system has been completely changed. Makefiles and Visual
   Studio project files have been removed. See the "## Other" section for
   more info.

## Acknowledgements:
 - Adrian Schreyer uncovered and reported a number of the bugs fixed
   in this release.

## Bug Fixes
 - the Recap code no longer generates illegal fragments for
   highly-branched atoms. (issue 1801871)
 - the Recap code no longer breaks cyclic bonds to N
   (issue 1804418)
 - A bug in the kekulization of aromatic nitrogens has been fixed
   (issue 1811276)
 - Bugs in the Atom Type definitions for polar carbons and positive
   nitrogens in BaseFeatures.fdef have been fixed. (issue 1836242)
 - A crash in the sanitization of molecules that only have degree 4
   atoms has been fixed; it now generates an exception. The underlying
   problem with ring-finding in these systems is still present. (issue
   1836576)
 - Mol files for molecules that have more than 99 atoms or bonds are
   no longer incorrectly generated. (issue 1836615)
 - Problems with the sping PIL and PDF canvases have been cleared
   up. The PIL canvas still generates a lot of warnings, but the
   output is correct.
 - The query "rN" is now properly interpreted to be "atom whose
   smallest ring is of size N" in SMARTS queries. It was previously
   interpreted as "atom is in a ring of size N". (issue 1811276)
   This change required that the default feature definitions for
   aromaticity and lumped hydrophobes be updated.
 - The MolSuppliers (SDMolSupplier, TDTMolSupplier, SmilesMolSupplier)
   no longer fail when reading the last element. (issue 1874882)
 - A memory leak in the constructor of RWMols was fixed.
 - A problem causing rapid memory growth with Recap analysis was fixed.
   (issue 1880161)
 - The Recap reactions are no longer applied to charged Ns or Os
   (issue 1881803)
 - Charges, H counts, and isotope information can now be set in
   reactions. (issue 1882749)
 - The stereo codes from double bonds (used for tracking cis/trans)
   are now corrected when MolOps::removeHs is called. (issue 1894348)
 - Various small code cleanups and edge case fixes were done as a
   result of things discovered while getting the VC8 build working.

## New Features
 - The SparseIntVect class (used by the atom pairs and topological
   torsions) is now implemented in C++.
 - The package $RDKit/Python/Chem/MolDb has been added to help deal
   with molecular databases. (this was intended for the August release
   and overlooked)
 - The module $RDKit/Python/Chem/FastSDMolSupplier has been added to
   provide a fast (at the expense of memory consumption) class for
   working with SD files. (this was intended for the August release
   and overlooked)
 - A new directory $RDKit/Projects has been created to hold things
   that don't really fit in the existing directory structure.
 - The new project $RDKit/Projects/DbCLI has been added. This contains
   command-line scripts for populating molecular database and
   searching them using substructure or similarity.
 - The code for calculating some descriptors has been moved into C++
   in the new module Chem.rdMolDescriptors. The C++ implementation is
   considerably faster than the Python one and should be 100%
   backwards compatible.
 - The MaxMinPicker (in Code/SimDivPickers) supports two new options:
   1) the user can provide a set of initial picks and the algorithm
      will pick new items that are diverse w.r.t. to those
   2) the user can provide a function to calculate the distance matrix
      instead of calculating it in advance. This saves the N^2 step of
      calculating the distance matrix.
 - A new piece of code demo'ing the use of the RDKit to add chemical
   functionality to SQLite is in Code/Demos/sqlite. This will
   eventually move from Demos into Code/sqlite once some more
   functionality has been added and more testing is done.
 - The distance geometry embedding code now supports using random
   initial coordinates for atoms instead of using the eigenvalues of
   the distance matrix. The default behavior is still to use the
   eigenvalues of the distance matrix.
 - The function Recap.RecapDecompose now takes an optional argument
   where the user can specify the minimum size (in number of atoms)
   of a legal fragment. (feature request 180196)
 - Dummy atoms can be expressed using asterixes, per the Daylight spec.
   Dummy atoms are also now legal members of aromatic systems (e.g.
   c1cccc*1 is a legal molecule). Support for supplying dummy atoms
   as "[Du]", "[X]", "[Xa]", etc. is now considered deprecated. In
   the next release a warning will be generated for these forms and
   in the release after that the old form will generate errors. Note
   that the output of dummy atoms will also change: in the next release
   the default output format will be "*".
   (feature request 186217)
 - A proof of concept for doing a SWIG wrapper of RDKit functionality
   has been added in: $RDBASE/Code/Demos/SWIG/java_example. This isn't
   even remotely production-quality; it's intended to demonstrate that
   the wrapping works and isn't overly difficult.

## Other
 - The full set of tests is now easier to setup and run on new
   machines. (issue 1757265)
 - A new build system, using Boost.Build, has been put into place on
   both the windows and linux sides. The new system does much better
   dependency checking and handles machine-specific stuff a lot
   better. The new system has been tested using Visual Studio 2003,
   Visual Studio Express 2005, Ubuntu 7.10, and RHEL5.
 - The "Getting Started in Python" document has been expanded.
 - There's now an epydoc config file for building the python
   documentation ($RDBASE/Python/epydoc.config).

# Release_Aug2007_1
(Changes relative to Release_April2007_1)

## Bug Fixes
 - operators and SparseIntVects. (issue 1716736)
 - the Mol file parser now calculates cis/trans labels for double
   bonds where the two ends had the same substituents. (issue 1718794)
 - iterator interface to DiscreteValueVects and UniformGrid3D. (issue
   1719831)
 - improper removal of stereochemistry from ring atoms (issue
   1719053)
 - stereochemistry specifications and ring bonds.  (issue 1725068)
 - handling of aromatic bonds in template molecules for chemical
   reactions. (issue 1748846)
 - handling of unrecognized atom types in the descriptor calculation
   code. (issue 1749494)
 - ChemicalReactionException now exposed to Python. (issue 1749513)
 - some small problems in topological torsions and atom pairs

## New Features
 - The Atom Pairs and Topological Torsions code can now provide
   "explanations" of the codes. See $RDBASE/Python/Chem/AtomPairs for
   details.
 - The PointND class has been exposed to Python
 - The "Butina" clustering algorithm [JCICS 39:747-50 (1999)] is now
   available in $RDBase/Python/Ml/Cluster/Butina.py
 - A preliminary implementation of the subshape alignment algorithm is
   available.
 - The free version of MS Visual C++ is now supported.   
 - There is better support for queries in MDL mol files. (issue 1756962)
   Specifically: ring and chain bond queries; the not modifier for
   atom lists; R group labels.
 - An EditableMol class is now exposed to Python to allow molecules to
   be easily edited. (issue 1764162)
 - The RingInfo class is now exposed to Python.
 - The replaceSidechains and and replaceCore functions have been added
   in the ChemTransforms library and are exposed to Python as
   Chem.ReplaceSidechains and Chem.ReplaceCore.
 - pickle support added to classes: PointND
 - atoms and bonds now support the HasQuery() and GetSmarts() methods
   from Python.   

## Other
 - Similarity scores can now be calculated from Python in bulk
   (i.e. calculating the similarity between one vector and a list of
   others). This can be substantially faster than calling the
   individual routines multiple times. The relevant functions are
   BulkTanimotoSimilarity, BulkDiceSimilarity, etc.
 - The calculation of AtomPairs and TopologicalTorsions fingerprints
   is now a lot more efficient.
 - Optimization of the Dice metric implementation for SparseIntVects
 - The Visual Studio build files have been moved to the directories
   $RDBASE/Code/Build.VC71 and $RDBASE/Code/Build.VC80. This allows
   simultaneous support of both versions of the system and cleans up
   the source trees a bit.
 - Boost version 1.34 is now supported (testing has been done on 1.34 and 1.34.1).
 - Updates to the "Getting Started" documentation.

# Release_April2007_1
(Changes relative to Release_Jan2007_1)

## Bug Fixes
 - handing of isotope information in SMILES has been fixed
 - "implicit" hydrogens are now added to charged atoms explicitly when
   writing SMILES. (issue 1670149)
 - the 2D->3D code no longer generates non-planar conjugated 4-rings
   (e.g. C1=CC=C1). (issue 1653802)
 - removing explicit hydrogens no longer produces incorrect smiles
   (issue 1694023)
 - bit indices and signature lengths in the AtomPairs code no longer
   being calculated incorrectly. *NOTE* this changes the bits that are
   set, so if you have existing signatures, they will need to be
   regenerated.
 - Fixed a bug causing MolSuppliers to generate incorrect length
   information when a combination of random access and iterator
   interfaces are used. (issue 1702647)
 - Fixed a bug leading to incorrect line numbers in error messages
   from the SDMolSuppler. (issue 1695221)

## New Features
 - chemical reactions are now supported
 - there is a new entry point into the 2D depictor code,
   compute2DCoordsMimicDistMat(), that attempts to generate 2D
   depictions that are similar to the structure described by the
   distance matrix. There's also a shorthand approach for calling this
   to mimic a 3D structure available as:
   AllChem.GenerateDepictionMatching3DStructure()
 - DiscreteValueVect and UniformGrid3D now support the binary
   operators |, &, +, and -.
 - a reader/writer for TPL files has been added.
 - support has been added for MolCatalogs: hierarchical catalogs that
   can store complete molecules.
 - the protrude distance metric for shapes has been added
 - pickle support added to classes: UniformGrid, DiscreteValueVect,
   Point
 - added the class DataStructs/SparseIntVect to improve performance
   and clarity of the AtomPairs code

## Other
 - the non-GUI code now supports python2.5; the GUI code may work with
   python2.5, but that has not been tested
 - the Mol and SD file parsers have been sped up quite a bit.
 - the "Crippen" descriptors are now calculated somewhat faster.
 - in-code documentation updates
 - new documentation for beginners in $RDBASE/Docs/Book

# Release_Jan2007_1
(Changes relative to Release_Oct2006_1)

## Bug Fixes
  - zero-atom molecules now trigger an exception
  - dummy atoms are no longer labelled 'Xe'
  - core leak in the mol file writer fixed
  - mol files with multiple charge lines are now correctly parsed
  - a workaround was installed to prevent crashes in the regression
    tests on Windows when using the newest VC++ v7 series compiler.
    (http://sourceforge.net/tracker/index.php?func=detail&aid=1607290&group_id=160139&atid=814650)
  - chirality perception (which requires partial sanitization) is no
    longer done by the MolFileParser when sanitization is switched
    off.
  - Two potential memory corruption problems were fixed (rev's 150 and
    151).

## New Features
  - optional use of chirality in substructure searches
  - MolWriters can now all take a stream as an argument
  - Chiral terms can now be included in the DistanceGeometry
    embedding.

## Other
  - $RDBASE/Code/Demos/RDKit/BinaryIO is a demonstration of using
    boost IOStreams and the ROMol pickling mechanism to generate
    highly compressed, random-access files of molecules.
  - the Point code has been refactored
