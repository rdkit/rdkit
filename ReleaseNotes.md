# Release_2024.09.1
(Changes relative to Release_2024.03.1)

## Acknowledgements
(Note: I'm no longer attempting to manually curate names. If you would like to
see your contribution acknowledged with your name, please set your name in
GitHub)


## Highlights

## Backwards incompatible changes

## New Features and Enhancements:

## Bug Fixes:
  
## Cleanup work:
 
## Code removed in this release:

## Deprecated code (to be removed in a future release):
- AtomPairs.Utils.NumPiElectrons is deprecated in favor of Chem.GetNumPiElectrons.
AtomPairs.Utils.NumPiElectrons failed if the atom had outgoing dative bonds (see Issue #7318).

# Release_2024.03.1
(Changes relative to Release_2023.09.1)

## Acknowledgements
(Note: I'm no longer attempting to manually curate names. If you would like to
see your contribution acknowledged with your name, please set your name in
GitHub)

Mark Archibald, Armin Ariamajd, Chris Von Bargen, Jason Biggs, Jonathan Bisson,
Jan C. Brammer, Jessica Braun, Benoît Claveau, David Cosgrove, James Davidson,
Hussein Faara, Théophile Gaudin, Gareth Jones, Christoph Hillisch, Tad Hurst,
Kevin Keating, Brian Kelley, Joos Kiener, David Lounsbrough, Jeremy Monat, Dan
Nealschneider, Yoshinobu Ogura, Marta Pasquini, Yakov Pechersky, Patrick Penner,
Rachael Pirie, Ricardo Rodriguez-Schmidt, Nate Russell, Ivan Tubert-Brohman,
Matthew Seddon, Leonid Stolbov, Paolo Tosco, Riccardo Vianello, Franz Waibl,
Rachel Walker, sitanshubhunia, skystreet8, dehaenw, dhibbit, vslashg, nbehrnd,
MarioAndWario, levineds-meta

## Highlights
- An initial version of support for atropisomers has been added; this will be expanded in future releases.
- Support for using multiple threads has been added in a few more places: many operations in rdMolStandardize, the fingerprint generators, and GetBestRMS()/GetAllConformerBestRMS()
- The initial release of version 2 of the RDKit C++ API; we will continue to expand this in future releases. The new API makes it easier to write correct and memory safe code. The current API is still supported and will remain so for the forseeable future, but we encourage C++ developers to start using v2 of the API in their code.

## Backwards incompatible changes
- Two changes to improve the defaults for conformer generation: the functions EmbedMolecule() and EmbedMultipleConfis() now use ETKDGv3 by default (previously they were using ETKDGV1) and only consider heavy atoms when calculating RMSD for conformer pruning (previously Hs were alos considered).
- The way that the number of radical electrons is calculated for atoms coming from mol blocks has been changed. Systems like a `[CH]` marked as a `DOUBLET` will now have three radical electrons assigned. This is consistent with the value from SMILES.
- The validation classes in MolStandardize were refactored in order to offer a simpler and more consistent API. In the C++ implementation, the `MolVSValidations` base class was removed and consolidated into `ValidationMethod`. Consequently, the `validate` method replaced `run` in the subclasses related to MolVS (namely `NoAtomValidation`, `FragmentValidation`, `NeutralValidation`, and `IsotopeValidation`) and all subclasses of `ValidationMethod` are now required to implement a `copy` method. Moreover, `MolStandardize::ValidationErrorInfo` was redefined as an alias for `std::string`. The changes related to the MolVS validation methods were similarly implemented in the Python API.
- Metal atoms (really any atom which has a default valence of -1) now have their radical electron count set to zero if they form any bonds. Metal atoms/ions without bonds will continue to be assigned a radical count of either 1 or 0 if they do/do not have an odd number of valence electrons. It is not possible in a cheminformatics system to generally answer what the spin state of a metal atom should be, so we are taking a simple and easily explainable approach. If you know the spin state of your species, you can directly provide that information by calling SetNumRadicalElectrons().
- Chirality will now be perceived for three-coordinate atoms with a T-shaped coordination environment and the wedge in the stem of the T. If we are perceiving tetrahedral stereo, it's possible to interpret this unambiguously.
- Bug fixes in the v2 tautomer hash algorithm will change the output for some molecules. Look at PR #7200 for more details: https://github.com/rdkit/rdkit/pull/7200
- RMS pruning during conformer generation now symmetrizes conjugated terminal groups by default. This can be disabled with the parameter "symmetrizeConjugatedTerminalGroupsForPruning"

## New Features and Enhancements:
  - Support writing detailed SMARTS queries to CTABs using the SMARTSQ mechanism
 (github issue #5819 from greglandrum)
  - add more error checking to substance groups
 (github issue #5923 from greglandrum)
  - add maxRecursiveMatches to SubstructMatchParameters
 (github issue #6017 from greglandrum)
  - Removed some code duplication between Depictor.cpp and common.h
 (github pull #6799 from greglandrum)
  - Add support for writing chirality and stereo in MaeWriter
 (github pull #6810 from rachelnwalker)
  - Implement MinimalLib get_mcs() version that returns JSON
 (github pull #6812 from ptosco)
  - support generalized substructure search in the SubstructLibrary
 (github pull #6835 from greglandrum)
  - Support copying of GeneralizeQueryMolecules 
 (github issue #6851 from greglandrum)
  - Enable chemist-friendly depiction of R-groups
 (github pull #6866 from ptosco)
  - Allow building DetermineBonds without YAeHMOP support
 (github pull #6885 from greglandrum)
  - Add multithreading to getBestRMS and new getAllConformerBestRMS
 (github pull #6896 from greglandrum)
  - switch to catch2 v3
 (github pull #6898 from greglandrum)
  - minilib functions exposure: mmpa
 (github pull #6902 from StLeonidas)
  - atropisomer handling added
 (github pull #6903 from tadhurst-cdd)
  - Add multi-threaded versions of some MolStandardize operations
 (github pull #6909 from greglandrum)
  - Add (multithreaded) functions to the fingerprint generators for calculating multiple fingeprints in one call
 (github pull #6910 from greglandrum)
  - Add Python modules to generate stubs and automatically patch docstrings
 (github pull #6919 from ptosco)
  - Update molecular templates headers and drop bond-length tests
 (github pull #6960 from github-actions[bot])
  - Add in place and multithread support for more of the MolStandardize code
 (github pull #6970 from greglandrum)
  - Enable in-tree builds and improve overloaded constructor docstrings
 (github pull #6980 from ptosco)
  - Change the defaults for the conformer generation to be ETKDGv3
 (github pull #6985 from greglandrum)
  - Added fingerprints to GeneralizedSubstruct search and extended SWIG wrappers
 (github pull #6991 from jones-gareth)
  - Allow sanitization to be disabled in PandasTools.LoadSDF
 (github issue #7019 from christophhillisch)
  - Add Atom::hasValenceViolation (Take 2)
 (github pull #7030 from cdvonbargen)
  - Please consider exposing maxBondMatchPairs param in rdRascalMCES.RascalOptions()
 (github issue #7054 from nate-russell)
  - Copy stereo and substance groups during insertMol
 (github issue #7064 from cdvonbargen)
  - [v2 API] FileParsers 
 (github issue #7074 from greglandrum)
  - [v2 API] Reaction Parsers
 (github issue #7075 from greglandrum)
  - Rationalize attachment points
 (github issue #7078 from cdvonbargen)
  - refactoring of MolStandardize validation module
 (github pull #7085 from rvianello)
  - Add a 'force' option to MolStandardizer::Uncharger
 (github pull #7088 from rvianello)
  - support sanitization of reaction product templates
 (github pull #7095 from greglandrum)
  - Support atropisomers in the conformer generator
 (github pull #7098 from greglandrum)
  - Compatibility with pathlib.Path
 (github pull #7100 from PatrickPenner)
  - Add option to sanitize reaction components like molecules
 (github issue #7108 from MartaPasquini)
  - [v2 API] MRV parsers
 (github pull #7110 from greglandrum)
  - Add v2 API for the molecule CDXML parser
 (github pull #7113 from greglandrum)
  - Make addStereoAnnotation public
 (github issue #7140 from cdvonbargen)
  - optimize batch operations when editing molecules
 (github pull #7145 from bp-kelley)
  - V2 API for the MolSuppliers
 (github pull #7168 from greglandrum)
  - Improve output of debugMol
 (github pull #7172 from greglandrum)
  - update cookbook, draw molecule with atom indices
 (github pull #7173 from nbehrnd)
  - Colinear bonds in depiction cause stereo to be lost when converting to mol block 
 (github issue #7177 from mps-hlx)
  - Update MinimalLib Dockerfiles
 (github pull #7182 from ptosco)
  - allow perception of stereo from T-shaped structures
 (github pull #7183 from greglandrum)
  - switch the TFD code to use a fingerprint generator
 (github pull #7187 from greglandrum)
  - Don't reset computed properties if already empty
 (github pull #7188 from rachelnwalker)
   - Enhance molzip to properly handle RGroupDecompositions
 (github pull #7202 from bp-kelley)
  - Add some ExplicitBitVect operations to Swig
 (github pull #7204 from jones-gareth)
  - Some modernization of core GraphMol classes
 (github pull #7228 from greglandrum)
  - Custom decimal precision
 (github pull #7229 from PatrickPenner)
  - Add Double Cubic Lattice Volume (DCLV).
 (github pull #7234 from RPirie96)
  - feat(minilib): expose the options parameter in get_inchi
 (github pull #7240 from BenoitClaveau)
  - Postpone clearing computed properties until after all Hs removed
 (github pull #7241 from rachelnwalker)
  - Speed up cleanMolStereo
 (github pull #7244 from ricrogz)
  - add HetAtomProtomerv2
 (github pull #7253 from greglandrum)
  - Support zero order bonds in V3K CTABs
 (github pull #7269 from greglandrum)
  - add option to symmetrize conjugated terminal groups when RMS pruning conformers
 (github pull #7270 from greglandrum)

## Bug Fixes:
  - STEREOANY bonds lead to non-stable SMILES/SMARTS strings
 (github issue #5499 from ricrogz)
  - Chemical reactions with radicals cannot be pickled and unpickled.
 (github issue #5890 from sitanshubhunia)
  - Postgresql: exact search showing false with radicals from CXSMILES
 (github issue #6276 from sitanshubhunia)
  - CXSMILES: atom with labels should not also have `dummyLabel` property set
 (github issue #6309 from greglandrum)
  - Query Features: Different input format leads to a different molecule
 (github issue #6349 from kienerj)
  - non-physical radical counts being preserved
 (github issue #6370 from greglandrum)
  - MolEnumerator: use repeat counts for SRUs when present
 (github issue #6429 from greglandrum)
  - Unexpected non-matching ElementGraph hashes
 (github issue #6472 from jepdavidson)
  - Fixes for canonicalization, and stereochemistry
 (github pull #6743 from tadhurst-cdd)
  - MCS query incorrect when ringCompare=RingCompare.StrictRingFusion
 (github issue #6773 from d-b-w)
  - Fixes bug in get_sss_json()
 (github pull #6806 from ptosco)
  - SWIG builds failing on Windows
 (github pull #6808 from jones-gareth)
  - Double bonds should not be depicted as crossed bonds in the presence of wavy bonds
 (github issue #6816 from ptosco)
  - We should be able to run the tests without boost::iostreams
 (github issue #6818 from greglandrum)
  - Fix stereo bond corruption on RGD.
 (github pull #6832 from jones-gareth)
  - MurckoScaffold.MakeScaffoldGeneric() has issues with isotopes
 (github issue #6836 from dehaenw)
  - Fix unclosed resource in BuildFuncGroupHierarchy
 (github pull #6846 from ricrogz)
  - RGD: Fix doEnumeration true for cores that are not bundles
 (github pull #6857 from jones-gareth)
  - Fix build error when serialization is off.
 (github pull #6867 from vslashg)
  - Wavy bonds in mol blocks can't be stereo enumerated
 (github issue #6876 from bp-kelley)
  - CDXML read of AND1 group (specifying racemic center) gets associated into an OR1 group
 (github issue #6887 from pechersky)
  - Segfault in JSONToMols when "commonchem" is an int
 (github issue #6890 from i-tub)
  - reapplyMolBlockWedging() should retain ENDDOWNRIGHT, ENDUPRIGHT dirs
 (github issue #6893 from ptosco)
  - MMPA FragmentMol segfaults when new stereo perception is turned on
 (github issue #6900 from jasondbiggs)
  - PositionVariationOp::getVariationCounts() does unnecessary copies of vectors
 (github issue #6906 from whosayn)
  - Obtaining descriptors via Descriptors.descList results in duplication of SPS.
 (github issue #6928 from wsuzume)
  - Some Clang-specific build instructions skip some clang compilers on mac
 (github issue #6941 from whosayn)
  - With new stereo, removing H from an Imine double bond does not remove bond stereo
 (github issue #6944 from ricrogz)
  - FindMolChiralCenters should honor RDK_USE_LEGACY_STEREO_PERCEPTION
 (github issue #6945 from ricrogz)
  - generateDepictionMatching2DStructure does not optimally align when refPatt!=None, allowRGroups=False, alignOnly=True
 (github issue #6952 from ptosco)
  - SpacialScore ignores undefined bond stereo
 (github issue #6957 from jasondbiggs)
  - GetAtomPairFingerprint yields different rooted FP from generator 
 (github issue #6958 from ptosco)
  - DetermineBonds() for PH3 yields no bonding
 (github issue #6961 from dhibbit)
  - Highlights of triple bonds come out wrong
 (github issue #6968 from DavidACosgrove)
  - MaeMolSupplier cannot read dummy atoms from Maestro files
 (github issue #6973 from ricrogz)
  - Chem.FindMolChiralCenters function should not be sensitive to atom-map numbers
 (github issue #6975 from skystreet8)
  - Parsing a Mol leaks the "_needsDetectBondStereo" property
 (github issue #6981 from ricrogz)
  - SubstructMatch maxRecursiveMatches is not being honored
 (github issue #6983 from ricrogz)
  - HierarchicalClusterPicker::pick() randomly fails with Invariant Violation
 (github issue #7001 from ricrogz)
  - rdkit.Dbase doesn't work correctly with Python 3.12
 (github issue #7009 from rvianello)
  - "Inconsistent state" when manually sanitizing and assigning stereo when using the new stereo algorithm
 (github issue #7023 from ricrogz)
  - Spacing bug in compute2DCoordsForReaction
 (github issue #7028 from KevKeating)
  - Update distance bounds calculation for conjugated double bonds in macrocycles
 (github pull #7032 from fwaibl)
  - Middle line in triple bond drawn to incorrect point when a wedged bond is present
 (github issue #7036 from greglandrum)
  - fragmentation of mol loses any sgroups
 (github pull #7056 from tadhurst-cdd)
  - CSharp Wrapper ExtendedQueryMol  Read Access Violation
 (github issue #7069 from jones-gareth)
  - removing an atom should not remove all stereo groups involving that atom.
 (github issue #7071 from greglandrum)
  - Sanitizing and assigning stereo twice can change bond stereo with new stereo
 (github issue #7076 from ricrogz)
  - testConrec.cpp:130 fails on ARM64
 (github issue #7083 from bjonnh-work)
  - Wrong stereochemistry in embedded rings from stereospecific SMILES
 (github issue #7109 from brje01)
  - Quaternary nitrogens with hydrogens are not a candidate for stereo
 (github issue #7115 from bp-kelley)
  - Some metal centers get radical electrons
 (github issue #7122 from cdvonbargen)
  - AddHs sets "no implicit Hs" on the atoms were Hs are added
 (github issue #7123 from ricrogz)
  - ReplaceBond may cause valence issues in specific edge cases
 (github issue #7128 from ricrogz)
  - Adding Wedge/Dash bond neighboring a stereo double bond causes a Precondition Violation
 (github issue #7131 from ricrogz)
  - Stereo Annotation Appears Incorrect
 (github issue #7157 from lounsbrough)
  - Unexpected exact mass values are returned for radium and radon
 (github issue #7162 from markarchibald)
  - Adding missing headers in ReactionParser.h
 (github pull #7163 from tgaudin)
  - fix: add PandasTools support for pandas 2.2
 (github pull #7165 from AAriam)
  - Fix leaking Bonds on unmatched ring closures
 (github pull #7178 from ricrogz)
  - fix a problem with tautomeric systems being extended too far
 (github pull #7200 from greglandrum)
  - Fixes #7181
 (github pull #7206 from greglandrum)
  - Fix Uncharger applying to already neutralized perhalic groups
 (github pull #7211 from rvianello)
  - Fix `Chem.Randomize.py`
 (github pull #7232 from JanCBrammer)
  - SGroup fields without values cause weird properties
 (github issue #7246 from ricrogz)
  - Remove duplicate entry in fragment descriptors
 (github pull #7249 from levineds-meta)
  - RDKit fails to parse "M RAD" lines were radical is 0
 (github issue #7256 from ricrogz)
  - Writing StereoGroups to Mol files should break lines at 80 characters
 (github issue #7259 from ricrogz)
  - Update ring fusion cache when needed
 (github pull #7274 from ptosco)
  - Ring stereo in SMILES inverted after sanitization in molecule with fragments
 (github issue #7295 from greglandrum)
  
## Cleanup work:
  - Switch over to using pytest to run the python tests
 (github pull #5916 from greglandrum)
  - Redundant variable`hasCoreDummies` in R-group decomposition code
 (github issue #6779 from MarioAndWario)
  - cmake cleanup 
 (github pull #6814 from greglandrum)
  - Remove boost::regex support
 (github issue #6817 from greglandrum)
  - remove the deprecated python implementation of MolStandardize
 (github pull #6819 from greglandrum)
  - Update CI, remove some warnings
 (github pull #6882 from greglandrum)
  - Deprecate some of the ancient python-based ML code
 (github pull #6891 from greglandrum)
  - Remove boost::regex support #6817
 (github pull #6913 from whosayn)
  - Fix minimal build, allow building without boost::serialization
 (github pull #6932 from greglandrum)
  - Drop unrequired zlib include which may break the windows build
 (github pull #6966 from ricrogz)
  - Compile time and runtime deprecation warnings
 (github pull #7004 from greglandrum)
  - New tests for specical query atoms and atropisomers
 (github pull #7010 from tadhurst-cdd)
  - fix GCC 13.2 warnings about redundant move in return statement
 (github pull #7029 from rvianello)
  - fix check of python version when updating Filters.cpp
 (github pull #7035 from rvianello)
  - fix several warnings originating from the swig wrappers
 (github pull #7063 from rvianello)
  - lock the versions of a bunch of stuff used in the CI builds
 (github pull #7082 from greglandrum)
  - remove deprecated packages from rdkit.ML
 (github pull #7107 from greglandrum)
  - require SWIG 4.1+ at cmake config time
 (github pull #7139 from rvianello)
  - RGD code cleanup
 (github pull #7186 from ptosco)
  - remove the broken Dbase.DbReport module
 (github pull #7227 from greglandrum)
  - remove a bunch of std::endls
 (github pull #7233 from greglandrum)
  - Avoid rebuilding FreeSASA at every build for no good reason
 (github pull #7245 from ptosco)

## Code removed in this release:
- The python implementations of MolStandardize has been removed.
  Please use the implementation in `rdkit.Chem.MolStandardize.rdMolStandardize` instead.
- The rdkit.six module, a leftover from the days when we supported both python 2
  and python 3, has been removed
- The RDKit implementation of standard machine learning algorithms has been
  removed. The affected packages include: rdkit.ML.Composite, rdkit.ML.DecTree,
  rdkit.ML.KNN, rdkit.ML.ModelPackage, rdkit.ML.NaiveBayes, rdkit.ML.Neural
  rdkit.ML.{Analyze,Screen,Grow,Build}Composite, rdkit.ML.CompositeRun,
  rdkit.ML.EnrichPlot
- The Dbase.DbReport package was no longer working and has been removed.

## Deprecated code (to be removed in a future release):
- The PDBMolSupplier class has been deprecated and will be removed in the next release
- The legacy Python code for drawing molecules has been deprecated and will be removed in the next release. This includes the following modules in rdkit.Chem.Draw: aggCanvas, cairoCanvas, canvasbase, MolDrawing, mplCanvas, qtCanvas, spingCanvas; the functions Draw.MolToImageFile(), Draw.MolToMPL(), and Draw.MolToQPixmap(); the "canvas" argument to the function Draw.MolToImage(); and calling Draw.MolToFile() with imageTypes other than PNG or SVG, 

# Release_2023.09.1
(Changes relative to Release_2023.03.1)

## Acknowledgements
(Note: I'm no longer attempting to manually curate names. If you would like to
see your contribution acknowledged with your name, please set your name in
GitHub)

Jason Biggs, Jonathan Bisson, David Cosgrove, Andrew Dalke, Christian W.
Feldmann, Eloy Félix, Richard Gowers, Tadd Hurst, Gareth Jones, Eisuke
Kawashima, Brian Kelley, Joos Kiener, Juuso Lehtivarjo, John Mayfield, Vedran
Miletić, Jeremy Monat, Dan Nealschneider, Timothy Ngotiaoco, Axel Pahl, Rachael
Pirie, Ricardo Rodriguez-Schmidt, Ernst-Georg Schmid, Paolo Tosco, Ivan
Tubert-Brohman, Riccardo Vianello, Rachel Walker, Maciej Wójcikowski, pierred5,
lhyuen, paconius, BartlomiejF, thomp-j, wangyingxie, teltim, Meteor-han,
abefrandsen, 

## Highlights
- The new RascalMCES code adds a very fast maximum common substructure
  implementation for pairs of molecules.
- The RDKit core now supports "generalized substructure searching", making it
  easier to take link nodes, variable attachment points, and tautomer queries
  into account when doing substructure searches. This is now also supported in
  the PostgreSQL cartridge.
- The RDKit now has support for reading and writing MRV files.

## Backwards incompatible changes
- The CDXML parser now returns mols with reasonable coordinates and in
the same coordinate axes as the other RDKit file parsers. 
- All methods returning `JSMol` and `JSReaction` objects now return a
`nullptr` (`null` in JS) when faling to generate a valid object, while
previously they were returning objects whose `is_valid()` method would
return `false`. The new implementation avoids the overhead of having to
call `delete()` on invalid objects and was approved in a
[public discussion on the `rdkit-js` GitHub repository](
  https://github.com/rdkit/rdkit-js/discussions/336)
- In JS MinimalLib, `MolIterator` was renamed to `MolList`: since now it
includes `at()`, `append()`, `insert()` and `pop()` methods, `MolIterator`
felt inappropriate. This change should have minimal impact on existing
JS code since so far there was no constructor for this class.
The only place where JS code needs to be updated is when parsing the return
value of `JSMol::get_frags()`: the return value consists of an object with
two keys, `molIterator` and `mappings`. The `molIterator` key has now
been renamed to `molList`.
- The user-configurable `MCSParameters::FinalMatchChecker` function is now
called after the built-in `FinalMatchChecker` function, rather as
alternatively to the built-in `FinalMatchChecker` function. This was a
design flaw which is worth correcting.
- Setting `MCSParameters::Timeout` to 0 means no timeout, rather than 0s
timeout, which is rather pointless as it would cause MCS to be canceled
immediately.
- Result SMARTS strings generated by `FindMCS` when
`MCSParameters::MatchRingFusionStrict` is `true` now include ring membership
queries where appropriate in order to ensure more specific substructure
matches.
- In MCS Verbose statistics, `SingleBondExcluded` was renamed to
`IndividualBondExcluded` to avoid confusion, since single bond has a
different meaning in chemistry.
- The error messages from failed type conversions in calls to `GetProp()` now
differ slightly between compilers. Instead of always including "boost::bad_any
cast", they now need to be matched with the regex `[B,b]ad any[\ ,_]cast`
- The functions for determining connectivity in DetermineBonds now use a more
efficient method by default. To go back to the old behavior, set the useVdw argument
to True.
- The algorithm for perception of atomic stereochemistry from 2D structures has
been rewritten. The new algorithm is more accurate, which results in some
differences in perceived stereo between this release and the previous ones.
- Information about stereo groups is no longer used in the SMILES
canonicalization process if CXSMILES are not being generated.

## New Features and Enhancements:
  - Mols matrix to grid image
 (github pull #6080 from bertiewooster)
  - Reduce space overhead of enabling the inclusion of properties when serializing molecules
 (github issue #6312 from rvianello)
  - Add optional sortsupport methods to the PostgreSQL GiST indices
 (github pull #6313 from rvianello)
  - Add a new parameter to mol_adjust_query_properties for generic query parameters
 (github pull #6332 from bjonnh-work)
  - add DebugDraw() function
 (github pull #6340 from greglandrum)
  - Optimize GetPropsFromDict: use tags for conversion
 (github pull #6355 from bp-kelley)
  - Fix cleanupOrganometallics and reinstate to sanitisation
 (github pull #6357 from DavidACosgrove)
  - postgres cartridge: cleanup a few obsolete build options
 (github pull #6363 from rvianello)
  - Fixes some issues in the SubstructLibrary JS implementation
 (github pull #6385 from ptosco)
  - Support TautomerQuery and MolBundle queries in the cartridge
 (github pull #6393 from greglandrum)
  - JS: Implement in-place aromatisation/kekulisation and avoid undesired exception
 (github pull #6407 from ptosco)
  - Optionally expose MCS to JS and extend optional compilation to JSReaction and JSSubstructLibrary
 (github pull #6409 from ptosco)
  - Add a method "HasQuery()" to Mol class
 (github issue #6411 from kienerj)
  - Enable using JSSubstructLibrary without pattern fps
 (github pull #6431 from ptosco)
  - Add support for generalized substructure searching
 (github pull #6443 from greglandrum)
  - Add atom and bond property parameters to substruct matching
 (github pull #6453 from rachelnwalker)
  - Replace a try..catch block with an if clause
 (github pull #6488 from ptosco)
  - Add an in place version of most of the MolStandardize functionality
 (github pull #6491 from greglandrum)
  - Exposed partial sanitization options to MinimalLib (both JS and CFFI)
 (github pull #6519 from ptosco)
  - Optionally forward Enhanced Stereo Group ids
 (github pull #6560 from ricrogz)
  - Add support for dative bonds in MOL files
 (github pull #6566 from bjonnh-work)
  - RASCAL MCES
 (github pull #6568 from DavidACosgrove)
  - Support additional generic groups
 (github pull #6570 from bjonnh-work)
  - Add support for Marvin files
 (github pull #6575 from bjonnh-work)
  - Add a "rootedAtAtom" to MolToSmarts
 (github pull #6581 from ricrogz)
  - RGD to support tautomers of core
 (github issue #6609 from jones-gareth)
  - Fix some build warnings
 (github pull #6618 from ricrogz)
  - Exposed log capture functionality to MinimalLib
 (github pull #6628 from ptosco)
  - add option to use sequential random seeds in the conformer generator
 (github pull #6639 from greglandrum)
  - Major MCS refactoring: new features and bug fixes
 (github pull #6646 from ptosco)
  - Lasso highlights
 (github pull #6653 from DavidACosgrove)
  - extract continuous lines from the conrec code
 (github pull #6676 from greglandrum)
  - Allow some tolerance in flatness determination
 (github pull #6696 from ricrogz)
  - Do not trust the 2D/3D tag in ctab mol files
 (github pull #6697 from ricrogz)
  - expose the CDXML reaction parsers to python
 (github pull #6700 from greglandrum)
  - Add hasQueryHs
 (github pull #6702 from bp-kelley)
  - Exporting to mol marks imine bonds EITHERDOUBLE when imine H is implicit
 (github issue #6703 from ricrogz) 
  - Use the connect-the-dots algorithm by default in DetermineBonds
 (github pull #6740 from greglandrum)
  - Add function to calculate all 3D descriptors
 (github pull #6741 from RPirie96)
  - Add SpacialScore
 (github pull #6742 from apahl)
  - Extract the core matching logic into a separate function
 (github pull #6754 from ptosco)


## Bug Fixes:
  - rdFMCS.FindMCS uses huge amounts of memory for this pair of molecules when CompleteRingsOnly is True
 (github issue #3965 from i-tub)
  - PF6- still can not get Bad Conformer Id after the #510 issue fix 
 (github issue #5145 from wangyingxie)
  - Order dependence for rdFMCS.FindMCS with MatchFusedRingsStrict
 (github issue #5411 from pierred5)
  - Failed FMCS results with certain seed SMARTS and MatchFusedRings* parameters on
 (github issue #5440 from pierred5)
  - Seed SMARTS to FindMCS() produces incorrect MCS
 (github issue #5457 from pierred5)
  - FindMCS returns wrong result with monoatomic molecules
 (github issue #5510 from ptosco)
  - queryMol from FindMCS doesn't match mols used to generate MCS
 (github issue #6082 from paconius)
  - Crash when parsing InChI
 (github issue #6172 from eloyfelix)
  - Iteration over Python GetAtoms is 10-20x slower than need be
 (github issue #6208 from d-b-w)
  - refactor(python): replace deprecated unittest methods
 (github pull #6304 from e-kwsm)
  - generateDepictionMatching2DStructure: bonds to R groups should be generic when matching
 (github pull #6306 from ptosco)
  - MolToSmiles(canonical=False) creates the wrong _smilesBondOutputOrder property
 (github issue #6315 from adalke)
  - MolToMolBlock ignores unspecified information for double bonds in rings
 (github issue #6316 from mwojcikowski)
  - bump yaehmop version
 (github pull #6330 from greglandrum)
  - rdMolDraw2D.MolDraw2DCairo produces Pre-condition Violation: no draw context when SetColour, DrawRect or DrawLine was called.
 (github issue #6336 from lhyuen)
  - Added cstdint include
 (github pull #6338 from vedranmiletic)
  - remove the dependency from python distutils in the top CMakeLists.txt file
 (github pull #6339 from rvianello)
  - take drawOptions into account when exporting structure to xlsx format
 (github pull #6341 from ptosco)
  - Fix swig memory leak
 (github pull #6346 from jones-gareth)
  - Add inlines to ForceFieldHelpers header functions
 (github pull #6356 from JLVarjo)
  - Bug relating to this PF6- still can not get Bad Conformer Id
 (github issue #6365 from teltim)
  - straightenDepiction should not consider 0-degree rotations as multiples of 60
 (github pull #6367 from ptosco)
  - expose two missing EmbedFailureCauses tags to python
 (github pull #6372 from greglandrum)
  - Molfile Unsaturation Query Not Parsed Correctly
 (github issue #6395 from timothyngo)
  - MolDraw2D: chiral tag overlapping atom label
 (github issue #6397 from greglandrum)
  - MolDraw2D: increasing padding results in the legend not being displayed
 (github issue #6400 from greglandrum)
  - expose some missing CXSmiles flags to python
 (github pull #6415 from greglandrum)
  - V3000 structure segfaults when converting to SVG
 (github issue #6416 from ergo70)
  - WedgeMolBonds won't wedge/dash a 2nd bond when input already has a wedge/dash around the same chiral atom
 (github issue #6423 from ricrogz)
  - MolEnumerate should clear the reaction properties on its results
 (github issue #6432 from greglandrum)
  - RDKit hangs indefinitely when parsing not so big molblock
 (github issue #6434 from eloyfelix)
  - Removing Hs on a pyrrol-like structure throws kekulization error
 (github issue #6437 from ricrogz)
  - Molecules from CDXML Parser have inverted, unrealistic atomic coordinates
 (github issue #6461 from greglandrum)
  - CDXML Parser does not preserve information about bond wedging
 (github issue #6462 from greglandrum)
  - boost::bad_any_cast error when calling getProp<string> on properties set by applyMolListPropsToAtoms<int64_t>
 (github issue #6465 from rachelnwalker)
  - Allow systems like C/C=N/[H] to be stereogenic with the new chirality code
 (github pull #6473 from greglandrum)
  - Fix RWMol::addAtom docstring
 (github pull #6477 from d-b-w)
  - StereoGroup information should not impact canonicalization when CXSMILES isn't being generated
 (github issue #6479 from greglandrum)
  - Fix a few broken docstrings
 (github pull #6480 from ptosco)
  - pin numpy to 1.24.3
 (github pull #6483 from bp-kelley)
  - CMAKE_INSTALL_PREFIX not honored for Python files installation on Windows
 (github pull #6485 from ricrogz)
  - Fixed tests that weren't being run in testDepictor.py
 (github pull #6486 from rachelnwalker)
  - Get tests to work when building without exception support (i.e., legacy pure JS library)
 (github pull #6487 from ptosco)
  - Fixes rdkit-js/issues/347
 (github pull #6490 from ptosco)
  - Make sure that molecules are shown as images by PandasTools also when DataFrames are truncated horizontally
 (github pull #6496 from ptosco)
  - MolToMolBlock writes "either" stereo for double bonds which shouldn't be stereo
 (github issue #6502 from ricrogz)
  - Double bonds are not correctly drawn on sulfoximines
 (github issue #6504 from ptosco)
  - RegistrationHash.GetMolLayers() with v2 tautomer hash does not filter CX extensions
 (github issue #6505 from ricrogz)
  - Drop the s_m_color_rgb property from MaeWriter
 (github pull #6511 from ricrogz)
  - update avalontools version to incorporate bug fixes
 (github pull #6513 from ptosco)
  - update windows DLL CI build config
 (github pull #6535 from greglandrum)
  - Add MolEnumerator to C#
 (github pull #6542 from jones-gareth)
  - MolDraw2D: placement of bond or atom labels gets confused when atoms overlap
 (github issue #6569 from greglandrum)
  - partial MCS failure
 (github issue #6578 from greglandrum)
  - Fix vulnerabilities found by fuzzer.
 (github pull #6579 from thomp-j)
  - allow building the cartridge against PostgreSQL 16
 (github pull #6580 from ptosco)
  - SIGSEGV while calculating molecular descriptors after using salt remover.
 (github issue #6592 from BartlomiejF)
  - Add newline to ConstrainedEmbed docstring.
 (github pull #6596 from DavidACosgrove)
  - Avoid leaking memory in case exceptions are thrown while generating FPs
 (github pull #6630 from ptosco)
  - Pre-condition violation in canonicalization of dative bond adjacent to double bond
 (github issue #6633 from ricrogz)
  - Incorrect most abundant isotope for Vanadium
 (github issue #6638 from abefrandsen)
  - Simple imine-containing molecule causes an infinite loop in FindStereo.cpp
 (github issue #6640 from ptosco)
  - Mol file parser strips stereogenic H from imine bonds
 (github issue #6664 from ricrogz)
  - ROMol move constructor and assignment do not update SubstanceGroup ownership
 (github issue #6681 from ricrogz)
  - Flexicanvas cuts off bottom of reaction
 (github issue #6685 from DavidACosgrove)
  - make sure bond attachpt info is pickled
 (github pull #6698 from greglandrum)
  - Catch meanBondLen of 0.0
 (github pull #6699 from DavidACosgrove)
  - Add trans layout for double bonds in rings
 (github pull #6709 from d-b-w)
  - Segmentation fault in MMFF
 (github issue #6728 from Meteor-han)
  - Fix chirality handling when the chiral atom is the first one in a SMARTS
 (github pull #6730 from johnmay)
  - ConnectTheDots can segfault if all atoms do not have residue info
 (github issue #6756 from jasondbiggs)
  - _moltoimg() should honor drawOptions.prepareMolsBeforeDrawing
 (github issue #6792 from ptosco)

## Cleanup work:
 - adjustQueryProperties cleanup
 (github pull #6361 from ptosco)
  - Fix identical for-loop variable names
 (github pull #6391 from JLVarjo)
  - Deprecate JSMol::is_valid() and JSReaction::is_valid() and return nullptr instead
 (github pull #6392 from ptosco)
  - misc jswrapper.cpp cleanup
 (github pull #6449 from ptosco)
  - Deprecate the pure python MolStandardize implementations.
 (github pull #6548 from greglandrum)
  - clear up some compiler warnings
 (github pull #6627 from greglandrum)
  - switch from boost::any to std::any
 (github pull #6662 from greglandrum)
  - Fail CI builds on compiler warnings + some fixes
 (github pull #6675 from ricrogz)
  - fix some more leaks
 (github pull #6684 from ricrogz)
  - Some small cleanups from the UGM Hackathon
 (github pull #6744 from greglandrum)
  - some modernization of the memory handing in the canonicalization code
 (github pull #6763 from greglandrum)


## Code removed in this release:
- The `compare` and `callback` methods (deprecated since release 2021.01)
were removed from the `MCSProgress`, `MCSBondCompare` and `MCSAtomCompare`
Python classes of the `rdFMCS` module. Both `compare` and `callback` methods
were replaced by `__call__`.
- The `SetAtomTyper` and `SetBondTyper` methods (deprecated since release 2021.01)
were removed from the `MCSParameters` Python class of the `rdFMCS` module.
The methods were replace by read-write properties `AtomTyper` and `BondTyper`,
respectively.
## Deprecated code (to be removed in a future release):
- JSMol::is_valid() and JSReaction::is_valid() are now deprecated and always
return true, as invalid `JSMol` and `JSReaction` cannot exist anymore.
- The python implementations of MolStandardize will be removed in the next release.
Please use the implementation in `rdkit.Chem.MolStandardize.rdMolStandardize` instead.

# Release_2023.03.1
(Changes relative to Release_2022.09.1)

## Acknowledgements
(Note: I'm no longer attempting to manually curate names. If you would like to
see your contribution acknowledged with your name, please set your name in
GitHub)

Michael Banck, Christopher Von Bargen, Jason Biggs, Jonathan Bisson, Jacob
Bloom, shang chien, David Cosgrove, Iren Azra Azra Coskun, Andrew Dalke, Eloy
Félix, Peter Gedeck, Desmond Gilmour, Mosè Giordano, Emanuele Guidotti, Tad
Hurst, Gareth Jones, Calvin Josenhans, Maria Kadukova, Brian Kelley, Joos
Kiener, Chris Kuenneth, Martin Larralde, Krzysztof Maziarz, Jeremy Monat, Michel
Moreau, Rocco Moretti, Lucas Morin, Dan Nealschneider, Noel O'Boyle, Vladas
Oleinikovas, Rachael Pirie, Ricardo Rodriguez-Schmidt, Vincent F. Scalfani,
Gregor Simm, Marco Stenta, Georgi Stoychev, Paolo Tosco, Kazuya Ujihara,
Riccardo Vianello, Franz Waibl, Rachel Walker, Patrick Walters,
'dangthatsright', 'mihalyszabo88', 'Deltaus', 'radchenkods',
'josh-collaborationspharma', 'jkh', 'yamasakih'

## Highlights
- The 2D coordinate generation can now optionally use templates when working with complex ring systems. We will continue to improve this functionality in future releases.
- There's now a single function which allows you to calculate all available 2D descriptors for a molecule: Descriptors.CalcMolDescriptors() 
- Support for working with organometallic molecules has improved: drawings of these structures are now better and there's new code for switching back and forth between dative and multi-center views of the bonding in systems like ferrocene.
- The fingerprint generator code has been improved and expanded with the idea of allowing user to switch entirely to the new code for the supported fingerprint types: Morgan, RDKit, topological torsion, and atom pairs.

## Backwards incompatible changes

- The ring-finding functions will now run even if the molecule already has ring information. Older versions of the RDKit would return whatever ring information was present, even if it had been generated using a different algorithm.
- The ring-finding functions now no longer consider dative bonds as possible ring bonds by default. All of the ring-finding functions have a new optional argument `includeDativeBonds` which can be used to change this behavior
- Generating 2D coordinates no longer has the side effect of running ring finding on molecules.
- The canonical SMILES and CXSMILES generated for molecules with enhanced stereochemistry (stereo groups) is different than in previous releases. The enhanced stereochemistry information and the stereo groups themselves are now canonical. This does *not* affect molecules which do not have enhanced stereo and will not have any effect if you generate non-isomeric SMILES. This change also affects the output of the MolHash and RegistrationHash code when applied to molecules with enhanced stereo.
- The doIsomericSmiles parameter in Java and C# ROMol.MolToSmiles() now defaults to true (previously it was false), thus aligning to the C++ and Python behavior.
- Double bonds which are marked as crossed (i.e. `bond.GetBondDir() == Bond.BondDir.EITHERDOUBLE`) now have their BondStereo set to `Bond.BondStereo.STEREOANY` and the BondDir information removed by default when molecules are parsed or `AssignStereochemistry()` is called with the `cleanIt` argument set to True.
- The conformers generated for molecules with three-coordinate chiral centers will be somewhat different due to the fix for #5883.
- Molecules which come from Mol or SDF files will now always have the "_MolFileChiralFlag" property set to the value of the chiral flag in the CTAB. In previous versions the property was not set if the chiral flag was 0.


## Bug Fixes:
  - GetSubstructMatches uniquify and maxMatches don't work well together
 (github issue #888 from adalke)
  - DrawRDKBits raised RDKit error when it applied to the compounds that contains imidazole.
 (github issue #2164 from yamasakih) 
  - MolFromMol2File: O.co2 atom type correctness check ignores phosphate groups
 (github issue #3246 from chmnk)
  - Enhanced Stereo is lost when using GetMolFrags(m, asMols=True)
 (github issue #4845 from kienerj)
  - Segfault with coordgen v3.0.0
 (github issue #4845 from lucasmorin222)
  - Dative bond and alkali and alkaline earth metals
  (github issue #5120 from marcostenta)
  - RGD Stereochemistry in decomposed structure is not copied to the matching core
 (github issue #5613 from jones-gareth)
  - fp.ToList() fails for empty molecule
 (github issue #5677 from baoilleach)
  - SMILES and SMARTS parse bonds in a different order
 (github issue #5683 from ricrogz)
  - postgresql makefile needs to be updated to use c++17
 (github issue #5685 from mbanck)
  - Exception raised when reading very large SMILES file
 (github issue #5692 from DavidACosgrove)
  - Update warning message about aromaticity detection
 (github pull #5696 from d-b-w)
  - stop building catch_main when tests are disabled
 (github pull #5697 from greglandrum)
  - Make PandasTools.RGroupDecompositionToFrame re-entrant
 (github pull #5698 from greglandrum)
  - PandasTools.RGroupDecompositionToFrame() should call ChangeMoleculeRendering()
 (github issue #5702 from greglandrum)
  - MolDraw2D should automatically set bond highlight color when atom colors are changed
 (github issue #5704 from greglandrum)
  - Use correct `_WIN32` macro for checking Windows target
 (github pull #5710 from giordano)
  - Environment not set properly in chirality tests for MinGW builds
 (github pull #5711 from ptosco)
  - windows.h header should be lowercase
 (github pull #5712 from ptosco)
  - Fixes bond index parsing for w/c/t/ctu labels in CXSMILES/CXSMARTS
 (github pull #5722 from ricrogz)
  - Fix a deprecation warning in pythonTestDirRoot
 (github pull #5723 from ricrogz)
  - allowNontetrahedralChiralty should be honored when reading/writing SMILES
 (github pull #5728 from greglandrum)
  - CFFI/MinimalLib fixes
 (github pull #5729 from ptosco)
  - Allow setting custom FREETYPE_LIBRARY/FREETYPE_INCLUDE_DIRS through CMake
 (github pull #5730 from ptosco)
  - Missing update path for postgreSQL from 3.8 to 4.2
 (github issue #5734 from Deltaus)
  - Avoid passing a NULL pointer to CanSmiles()
 (github pull #5750 from ptosco)
  - CDXML reader incorrectly sets stereo on crossed double bonds
 (github issue #5752 from baoilleach)
  - `R` atom label information lost in molfile if not backed by a `M RGP` entry
 (github issue #5763 from eloyfelix)
  - Missing monomer labels when depicting `MON` SGroups
 (github issue #5767 from eloyfelix)
  - Wrongly oriented SGroup bracket
 (github issue #5768 from eloyfelix)
  - Adjust LocaleSwitcher on Windows when RDK_BUILD_THREADSAFE_SSS not set
 (github pull #5783 from roccomoretti)
  - KekulizationException in tautomer canonicalization
 (github issue #5784 from d-b-w)
  - ChemicalReactionToRxnBlock ignores separateAgents if forceV3000 is true
 (github issue #5785 from jacobbloom)
  - extend the allowed valences of the alkali earths
 (github pull #5786 from greglandrum)
  - Minimallib build (rdkit-js) not working for release 2022.09.2
 (github issue #5792 from MichelML)
  - Remove dependency on MSVC runtime DLL in MinGW builds
 (github pull #5800 from ptosco)
  - Update macOS target platform to 10.13
 (github pull #5802 from ptosco)
  - `R#` atom label information lost in molfile if not handled by the `RGP` spec
 (github issue #5810 from eloyfelix)
  - Stop using recommonmark in the documentation
 (github issue #5812 from greglandrum)
  - Properties with new lines can create invalid SDFiles
 (github issue #5827 from bp-kelley)
  - Allow building PgSQL RPM and DEB packages
 (github pull #5836 from ptosco)
  - Additional output is incorrect when FP count simulation is active
 (github issue #5838 from ptosco)
  - Explicit valence check SiPa13fails for certain SMILES
 (github issue #5849 from josh-collaborationspharma)
  - Set emsdk path for freetype in emscripten builds
 (github pull #5857 from ptosco)
  - DrawMorganBit fails by default 
 (github issue #5863 from eguidotti)
  - Fix #5810 in V2000 mol files.
 (github pull #5864 from eloyfelix)
  - Chemical drawings should be automatically enabled on Colab
 (github pull #5868 from kuelumbus)
  - use enhanced stereo when uniquifying in SimpleEnum
 (github pull #5874 from greglandrum)
  - Conformer Generation Fails for three-coordinate Ns with specified stereo
 (github issue #5883 from gncs)
  - Fix documentation example for KeyFromPropHolder
 (github pull #5886 from gedeck)
  - Allow unrecognized atom types when strictParsing=False
 (github pull #5891 from greglandrum)
  - DetermineBonds assigning methyl carbon as tetrahedral center
 (github issue #5894 from jasondbiggs)
  - numpy.float is no longer supported and causes exceptions 
 (github issue #5895 from PatWalters)
  - moldraw2DTest1 failure when building on aarch64
 (github issue #5899 from vfscalfani)
  - DetermineBondOrders running out of memory on medium-sized disconnected structure
 (github issue #5902 from jasondbiggs)
  - clear MDL Rgroup labels from core atoms when we aren't using them
 (github pull #5904 from greglandrum)
  - Conformer generator produces strange structure for substituted butadiene
 (github issue #5913 from gncs)
  - `MHFPEncoder::Distance` doesn't compute a (Jaccard) distance
 (github issue #5919 from althonos)
  - AvalonTools: Avoid that trailing garbage pollutes the fmemopen buffer
 (github pull #5928 from ptosco)
  - "not" queries in molfiles get inverted
 (github issue #5930 from d-b-w)
  - CalcTPSA() doesn't use options when caching
 (github issue #5941 from greglandrum)
  - Bad drawing of end points for dative bonds
 (github issue #5943 from DavidACosgrove)
  - Extremes of drawn ellipses not being calculated correctly.
 (github issue #5947 from DavidACosgrove)
  - Arrow heads of dative bonds are different sizes
 (github issue #5949 from DavidACosgrove)
  - stop caching ring-finding results
 (github pull #5955 from greglandrum)
  - Wrong bond endpoint when connecting to wedge bond in 2D image
 (github issue #5963 from stgeo)
  - Tiny change to get demo.html to load in legacy browsers
 (github pull #5964 from ptosco)
  - detect bad double bond stereo in conformer generation 
 (github pull #5967 from greglandrum)
  - drawing code should not generate kekulization errors
 (github issue #5974 from greglandrum)
  - Adjust expected test results for newer freetype versions
 (github pull #5979 from greglandrum)
  - CanonicalRankAtomsInFragment example in the documentation is not reproducible
 (github issue #5986 from chmnk)
  - Exception in RegistrationHash for molecules with bad bond directions
 (github pull #5987 from d-b-w)
  - Updated the GetMolHash docstring for accuracy
 (github pull #5988 from irenazra)
  - Fix a problem with pickling molecules with more than 255 rings
 (github pull #5992 from greglandrum)
  - Support Python 3.11
 (github pull #5994 from greglandrum)
  - Incorrect disconnection of CC(=O)O[Mg]OC(=O)C
 (github issue #5997 from DavidACosgrove)
  - PostgreSQL autovacuum stuck when molecules with query features are stored in mol columns
 (github issue #6002 from mihalyszabo88)
  - Remove `and` from C++ headers
 (github pull #6003 from d-b-w)
  - [PH3] incorrectly recognized as potential stereo center
 (github issue #6011 from greglandrum)
  - Potential nontetrahedral stereo is recognized when nontetrahedral stereo is disabled.
 (github issue #6012 from greglandrum)
  - MolEnumerator is not propagating query information to molecules
 (github issue #6014 from greglandrum)
  - Reactions do not propagate query information to products
 (github issue #6015 from greglandrum)
  - Error rendering to very small canvas
 (github issue #6025 from DavidACosgrove)
  - Bad double bond drawn for collinear atoms
 (github issue #6027 from DavidACosgrove)
  - Fix some minor leaks
 (github pull #6029 from ricrogz)
  - Cannot draw molecule which includes an atom with a `[!#X]` query (for any X)
 (github issue #6033 from ShangChien)
  - FragmentOnBonds may create unexpected radicals
 (github issue #6034 from ricrogz)
  - Calling MurckoScaffold on molecule causes bug in pickling
 (github issue #6036 from dangthatsright)
  - Bump maeparser and coordgen versions
 (github pull #6039 from ricrogz)
  - enhanced stereo is still included in CXSMILES if isomericSmiles=False
 (github issue #6040 from greglandrum)
  - Issues with ACS1996 drawing mode on a small canvas
 (github issue #6041 from DavidACosgrove)
  - Cyclobutyl group in a macrocycle triggers a stereo center
 (github issue #6049 from cdvonbargen)
  - stereogroups not combined when parsing CXSMILES
 (github issue #6050 from greglandrum)
  - Regression in depicting molecules with MDL query atoms
 (github issue #6054 from ptosco)
  - Do not include dative bonds in ring finding by default
 (github issue #6058 from DavidACosgrove)
  - Remove check for ring information from Atom::Match
 (github pull #6063 from fwaibl)
  - Correct docstring for minFontSize.
 (github pull #6066 from DavidACosgrove)
  - Minor code cleanup
 (github pull #6101 from ptosco)
  - Dummy atoms should not be considered to be metals for M and MH queries
 (github issue #6106 from greglandrum)
  - Drawing in ACS mode crops small images
 (github issue #6111 from DavidACosgrove)
  - Drawing in ACS1996 mode throws ValueError: Bad Conformer Id if molecule has no coords
 (github issue #6112 from DavidACosgrove)
  - Single atom or queries with hydrogens shouldn't trigger warning in mergeQueryHs
 (github issue #6119 from bp-kelley)
   - DetermineBonds fails for single H atom
 (github issue #6121 from gncs)
  - MinimalLib: avoid that assignStereochemistry() fails when removeHs=true
 (github pull #6134 from ptosco)
  - Round-tripping a reaction through pickle changes the outputs from RunReactants
 (github issue #6138 from kmaziarz)
  - RGD and enhanced stereochemistry
 (github issue #6146 from jones-gareth)
  - MaeMolSupplier requires bond block
 (github issue #6153 from cdvonbargen)
  - Incorrect doule bond drawing with MolDraw2DSVG
 (github issue #6160 from radchenkods)
  - BondDir not cleared from bonds that aren't stereoactive
 (github pull #6162 from greglandrum)
  - Crossed bond not correctly drawn
 (github issue #6170 from ptosco)
  - ReactionFromRxnBlock fails on bond with reacting center status set
 (github issue #6195 from jones-gareth)
  - Possible regression in the atom/bond highlighting code
 (github issue #6200 from ptosco)
  - Updated README to build the PostgreSQL cartridge + bug fix
 (github pull #6214 from ptosco)
  - Atoms may get flagged with non-tetrahedral stereo even when it is not allowed
 (github issue #6217 from ricrogz)
  - Fix TorsionFingerprints for 15 membered rings
 (github pull #6228 from kazuyaujihara)
  - Fix build warnings
 (github pull #6235 from ricrogz)
  - Tri-coordinate atom with implicit + neighbor H atom is found potentially chiral
 (github issue #6239 from ricrogz)
  - DativeBondsToHaptic doesn't set _MolFileBondEndPts correctly.
 (github issue #6252 from DavidACosgrove)
  - Round-tripping ferrocene through HapticBondsToDatives loses drawing highlights.
 (github issue #6253 from DavidACosgrove)
  - Using Chiral Tag instead of CIPCode to ensure preservation of chirality in addHs
 (github pull #6259 from HalflingHelper)
  - Update assignSubstructureFilters.py
 (github pull #6270 from OleinikovasV)
  - deal with deprecated DataFrame.append method
 (github pull #6272 from greglandrum)
  - compile-time error with GCC 12.2.1 on Fedora 36
 (github issue #6274 from rvianello)
  - Fix UnitTestPandasTools for running without pandas installed.
 (github pull #6299 from roccomoretti)
  - Aromatic dashes look bad
 (github pull #6303 from greglandrum)


## Cleanup work:
  - Do deprecations for the 2023.03 release
 (github pull #5675 from greglandrum)
  - run clang_format
 (github pull #5676 from greglandrum)
  - Cleanup work on documentation Makefile
 (github pull #5804 from greglandrum)
  - Refactor RGD moving function implementations from header to source files
 (github pull #5958 from jones-gareth)
  - Disable POPCNT on M1
 (github pull #6081 from bjonnh-work)
  - Remove spurious full stops from warnings.
 (github pull #6124 from DavidACosgrove)
  - Reformat Python code for 2023.03 release
 (github pull #6294 from ricrogz)
  - Reformat C/C++ code ahead of 2023.03 release
 (github pull #6295 from ricrogz)


## New Features and Enhancements:
  - mol V3000: multicenter dative bond
 (github issue #5121 from marcostenta)
  - add molecular filter examples
 (github pull #5647 from RPirie96)
   - Use templates in RDKit coordinate generation
 (github pull #5643 from rachelnwalker)
  - add MACCS fp to the MinimalLib
 (github pull #5707 from eloyfelix)
  - Enable additional parameters in prepareAndDrawMolecule() and expose them to CFFI/MinimalLib
 (github pull #5731 from ptosco)
  - add includeRedundantEnvironments support to GetMorganGenerator
 (github pull #5732 from greglandrum)
  - FingerprintGenerator refactoring
 (github pull #5748 from greglandrum)
  - Expose RDLog to SWIG wrappers
 (github pull #5749 from ptosco)
  - Add a timeout protection for CIP calculations
 (github pull #5772 from tadhurst-cdd)
  - Expose getMolFrags in CFFI and MinimalLib
 (github pull #5774 from ptosco)
  - Get the wrappers working with SWIG 4.0
 (github pull #5795 from greglandrum)
  - Update AvalonTools to version 2.0.4a
 (github pull #5796 from ptosco)
  - Add early example of drawing a molecule to Getting Started with the RDKit in Python
 (github pull #5803 from bertiewooster)
  - Enable get_molblock(details_json) from MinimalLib
 (github pull #5806 from ptosco)
  - Improvements to PandasTools.SaveXlsxFromFrame
 (github pull #5835 from ptosco)
  - swap boost::tuple to std::tuple
 (github pull #5851 from greglandrum)
  - Make it easy to calculate all 2D descriptors
 (github pull #5892 from greglandrum)
  - Introduces AvgIpc descriptor
 (github pull #5896 from greglandrum)
  - Add SMILES to each group abbreviation in Cookbook
 (github pull #5908 from bertiewooster)
  - Support SubstanceGroups and StereoGroups in JSON
 (github pull #5909 from greglandrum)
  - Add info about mergeHs to README.
 (github pull #5910 from DavidACosgrove)
  - Cookbook - update entry 1 and add entries 38 and 39
 (github pull #5918 from vfscalfani)
  - Allow the sources of conformer generation failures to be retrieved
 (github pull #5960 from greglandrum)
  - Create getExperimentalTorsions() function
 (github pull #5969 from greglandrum)
  - Molblock wedging improvements
 (github pull #5981 from ptosco)
  - MinimalLib JS functions to add/remove Hs in place
 (github pull #5984 from ptosco)
  - Adds Pat Walter's Chembl Filters extraction to the FilterCatalog
 (github pull #5991 from bp-kelley)
  - Add depiction coordinates to molzip
 (github pull #5993 from jones-gareth)
  - Enable using STL algorithms on ROMol atoms and bonds
 (github pull #5995 from ptosco)
  - Enable building MinimalLib as a plain JS file for usage in legacy/headless browsers
 (github pull #5999 from ptosco)
  - Allow WriteSDF to create v3000 SDF files
 (github pull #6004 from jkhales)
  - add maxRecursiveMatches to SubstructMatchParameters
 (github issue #6017 from greglandrum)
  - Expose fingerprint generator options to python
 (github pull #6024 from greglandrum)
  - Allow SMARTS of zero order bonds to match zero order bonds
 (github pull #6037 from d-b-w)
  - Change IUPAC metal->non-metal single bonds to dative
 (github pull #6038 from DavidACosgrove)
  - Add canonicalization of stereo groups (enhanced stereo)
 (github pull #6051 from greglandrum)
  - Improve MaeMolSupplier API
 (github pull #6053 from ricrogz)
  - Enable optional visualization of complex query atoms in a more compact form
 (github pull #6056 from ptosco)
  - Start a Maestro file (.mae) writer.
 (github pull #6069 from ricrogz)
  - Expose some stereochemistry-related functions to SWIG wrappers
 (github pull #6075 from ptosco)
  - Add option to only include shortest paths for topological torsion fingerprints
 (github pull #6090 from greglandrum)
  - Enable "smilesSymbol" substitution in SMARTS
 (github pull #6096 from ricrogz)
  - Add the option to wedge two bonds at chiral centers
 (github pull #6108 from greglandrum)
  - Another minor code cleanup
 (github pull #6109 from ptosco)
  - A few SWIG tweaks
 (github pull #6110 from ptosco)
  - Stereochemistry-related SWIG updates
 (github pull #6127 from ptosco)
  - SWIG pickling improvements (and other cleanup)
 (github pull #6133 from ptosco)
  - Improved handling of organometallics
 (github pull #6139 from DavidACosgrove)
  - expose two missing QueryAtom types to python
 (github pull #6158 from greglandrum)
  - Support Pseudoatoms like Pol and Mod in the RegistrationHash
 (github pull #6175 from irenazra)
  - Better name for areBondsLinear.
 (github pull #6196 from DavidACosgrove)
  - add features to allow drawing molecules in arbitrary positions on a large canvas
 (github pull #6210 from greglandrum)
  - Support chirality when determining if a molecule is a reaction reactant
 (github issue #6211 from jones-gareth)  
  - rdMolHash.MolHash function should allow customization of the CXSmiles via Chem.CXSmilesFields
 (github issue #6224 from irenazra)
  - Updated README for cartridge installation into conda PostgreSQL
 (github pull #6236 from ptosco)
  - Add a function to translate the MDL chiral flag into enhanced stereo groups
 (github issue #6241 from ricrogz)
  - Add support for generic matching in the PgSQL cartridge
 (github pull #6269 from bjonnh-work)
  - allowOptionalAttachments should also include terminal query atoms matching hydrogen
 (github pull #6280 from ptosco)
  - Exposed queryColour in MolDrawOptions
 (github pull #6282 from ptosco)
  - Add a new tautomer mol hash
 (github pull #6289 from glandrum)
  - has_coords() now reports whether coords are 2D or 3D if present
 (github pull #6297 from ptosco)
 - Improve the installation/testing instructions. 
 (github pull #6298 from roccomoretti)

## Code removed in this release:
- The `SmilesParserParams` option `useLegacyStereo` has been removed. Please use
  `SetUseLegacyStereoPerception()` instead. 
- The following JS methods:
  * generate_aligned_coords()
  * get_morgan_fp()
  * get_morgan_fp_as_uint8array()
  * get_pattern_fp()
  * get_pattern_fp_as_uint8array()
  which used to take several individual parameters have been removed. 
  Please use the versions which take a single JSON string parameter.
- The `PrintAsBase64PNGString` function in `PandasTools` has been removed.
  Please use `PrintAsImageString` instead.




## Deprecated code (to be removed in a future release):


# Release_2022.09.1
(Changes relative to Release_2022.03.1)

## Acknowledgements

Jonathan Bisson, Andy Cai, David Cosgrove, JP Ebejer, Aleš Erjavec, Peter
Gedeck, Mosè Giordano, Sreya Gogineni, Emanuele Guidotti, Hyeonki Hong, Gareth
Jones, Per Johnsson, Maria Kadukova, Eisuke Kawashima, Brian Kelley, Alan
Kerstjens, Michel Moreau, Dan Nealschneider, Santeri Puranen, Ricardo
Rodriguez-Schmidt, Guy Rosin Jeff van Santen, David Schaller, David W.H.
Swenson, Paolo Tosco, Antonio Trande, Ivan Tubert-Brohman, Alexandra Wahab,
Rachel Walker, balducci, GLPG-GT

## Highlights
- The new RegistrationHash module provides one of the last pieces required to
  build a registration system with the RDKit.
- This release includes an initial version of a C++ implementation of the
  xyz2mol algorithm for assigning bonds and bond orders based on atomic
  positions. This work was done as part of the 2022 Google Summer of Code.
- A collection of new functionality has been added to minimallib and is now
  accessible from JavaScript and other programming languages.

## Backwards incompatible changes
- Changes to the way atomic chirality is used in the canonicalization algorithm
  mean that canonical atom ranking and canonical SMILES generated with this
  RDKit version can be different from those generated with previous versions
- `GetBestRMS() and CalcRMS()` by default now treat terminal conjugated functional
  groups like carboxylate and nitro symmetrically. For example, the group
  `C(=[O:1])[O-:2]` can match in either orientation. The SMARTS pattern which is
  used to recognize affected groups is:
  `[{atomP};$([{atomP}]-[*]=[{atomP}]),$([{atomP}]=[*]-[{atomP}])]~[*]` where
  `{atomP}` is `O,N;D1`. The previous behavior can be restored using by setting
  the `symmetrizeConjugatedTerminalGroups` argument to false when calling
  `GetBestRMS() and CalcRMS()`
- The following `#defines` are no longer provided in/used by the C++ code or `RDConfig.h`:
  - `BUILD_COORDGEN_SUPPORT`: use `RDK_BUILD_COORDGEN_SUPPORT` instead
  - `RDK_THREADSAFE_SSS`: use `RDK_BUILD_THREADSAFE_SSS` instead
  - `USE_BUILTIN_POPCOUNT`: use `RDK_OPTIMIZE_POPCNT` instead
- The Python function `Chem.GetSSSR()` now returns the SSSR rings found instead
  of just returning the count of rings. This is consistent with
  `Chem.GetSymmSSSR()` and more useful.
- The SMILES parser will ignore the value of
  `SmilesParserParams.useLegacyStereo` unless it is set to `false`. See the
  deprecation note about `useLegacyStereo` below for more information.
- The CFFI function `set_2d_coords_aligned()` now takes an additional `char **match_json`
  parameter; if `match_json` is not not `NULL`, `*match_json` will point to a
  JSON string containing the atoms and bonds which are part of the match.
  It is up to the user to free this string.
- The aliphatic imine rule used in tautomer enumeration has been changed to more
  closely match the definition in the original paper.

## Bug Fixes:
  - H atoms in SGroups cause RDKit to clear SGroups after parsing
 (github issue #2716 from ricrogz)
  - RDKit misplaces stereochemistry/chirality information for small ring
 (github issue #2984 from d-b-w)
  - DrawMorganBit returns empty image for "isolated" fingerprints
 (github issue #4242 from eguidotti)
  - Cores with query atoms may fail to R-group-decompose molecules
 (github issue #4505 from ptosco)
  - Unable to serialize coordinates as floats in combination with *Props
 (github issue #4621 from santeripuranen)
  - Image Generation: Highlighting looks off when bondLineWidth is increased for PNG generation
 (github issue #5122 from rachelnwalker)
  - RDKit crashes on CIP label calculation
 (github issue #5142 from ricrogz)
  - Modified the JS tests to comply with older nodejs versions
 (github pull #5148 from ptosco)
  - Presence of exocyclic S/D, S/A or D/A query bonds prevents benzene from being recognized as aromatic
 (github issue #5152 from rachelnwalker)
  - Fix for RGD dummy atom bug in RDKit::replaceCore
 (github pull #5154 from jones-gareth)
  - KekulizeException of molecule from Smarts pattern with new RDKit release 
 (github issue #5156 from schallerdavid)
  - Very small fix to avoid an AttributeError
 (github pull #5163 from ptosco)
  - issue with V3000 SD files containing enhanced stereochemistry information
 (github issue #5165 from GLPG-GT)
  - Draw.MolToQPixmap raises a TypeError with Python 3.10
 (github issue #5166 from ales-erjavec)
  - Standardization via RDKit breaks molecules
 (github issue #5169 from malteseunderdog)
  - Multiple calls to BlockLogs() permanently disable logging
 (github issue #5172 from ricrogz)
  - Check architecture of the target system to optimise popcnt
 (github pull #5182 from giordano)
  - More consistently check for `WIN32` instead of `MSVC` in CMake files
 (github pull #5183 from giordano)
  - Atom indices inside double bond
 (github issue #5185 from DavidACosgrove)
  - Bug in IPython display after setting a molecule property
 (github issue #5192 from dwhswenson)
  - Zero & coordinate bonds are being taken into account for chirality
 (github issue #5196 from ricrogz)
  - FindPotentialStereo does not clean stereoflags from atoms which cannot be stereocenters
 (github issue #5200 from greglandrum)
  - Fixes array overflow in FMCS code
 (github pull #5205 from ricrogz)
  - PeriodicTable initialization is not thread safe
 (github issue #5207 from ricrogz)
  - Find and remove deprecated ifdefs
 (github issue #5210 from greglandrum)
  - Fix use of not thread safe function localtime() 
 (github pull #5211 from ricrogz)
  - Fix duplicate non thread safe check in VarianceDataForLabel
 (github pull #5212 from ricrogz)
  - RDKit::Utils::LocaleSwitcher is not thread safe
 (github issue #5214 from ricrogz)
  - Core with query atoms and no user definded attachment points may create poor decompostions
 (github issue #5222 from jones-gareth)
  - error: format not a string literal and no format arguments
 (github issue #5234 from sagitter)
  - Fix qt build under VS2019
 (github pull #5238 from ricrogz)
  - Precondition violation on chiral Atoms with zero order bonds
 (github issue #5239 from ricrogz)
  - pyForceFieldConstraints test failed
 (github issue #5252 from sagitter)
  - drawReaction() should not hit a PRECONDITION with prepareMolsBeforeDrawing=false
 (github issue #5259 from ptosco)
  - Atom annotations poorly placed on highlighted atoms
 (github issue #5269 from DavidACosgrove)
  - Make the aliphatic imine rule more closely match the definition in the paper
 (github pull #5270 from greglandrum)
  - rdMolDraw2D.PrepareMolForDrawing(None) causes segmentation fault
 (github issue #5298 from perjo)
  - MolStandardize: uncharger failing in molecules with zwitterionic sulfone
 (github issue #5317 from greglandrum)
  - MolStandardize: some operations throwing on non-standardized molecules
 (github issue #5318 from greglandrum)
  - MolStandardize: cleanup() function not correctly reassigning stereochemistry
 (github issue #5320 from greglandrum)
  - Runtime error when writing reaction with substance group to V3000 rxnblock
 (github issue #5324 from rachelnwalker)
  - MolFromMolBlock should correctly assign stereochemistry to 3D molecules
 (github issue #5327 from greglandrum)
  - assignChiralTypesFrom3D() ignores wiggly bonds
 (github issue #5328 from greglandrum)
  - Molzip segfaults instead of throwing an error when multiple bonds are formed to the same pairs of atoms
 (github issue #5334 from loluwot)
  - leftover debugging code makes build of 2022_03_3 tarball fail
 (github issue #5351 from balducci)
  - Prevent wedging ring bonds
 (github pull #5356 from greglandrum)
  - Class info info missing for for wavy bonds in SVGs.
 (github pull #5358 from DavidACosgrove)
  - MolToSmiles with doRandom=True raises errors with molecules containing fragments.
 (github issue #5372 from greglandrum)
  - Restore #5103 that was accidentally reverted in #5132
 (github pull #5381 from ptosco)
  - cairo error when using similarity maps
 (github issue #5383 from greglandrum)
  - MolDraw2DQt and MolDraw2DJS don't draw wavy bonds
 (github issue #5386 from greglandrum)
  - MolDraw2DQt draws brackets incorrectly
 (github issue #5389 from greglandrum)
  - PandasTools should not always bring in IPythonConsole
 (github pull #5392 from greglandrum)
  - Tautomer enumeration removes stereochemistry depending on how aromatic rings are defined in SMILES
 (github issue #5402 from greglandrum)
  - Incorrect perception of pseudoasymmetric centers on non-canonical molecules
 (github issue #5403 from ptosco)
  - Fix performance issue in RingUtils::checkFused
 (github pull #5410 from rachelnwalker)
  - Multi-coloured highlights go wrong with small fonts
 (github issue #5420 from DavidACosgrove)
  - Parsing a Mol block/file does not clear the "molTotValence" property from atoms
 (github issue #5423 from ricrogz)
  - Pre-condition Violation - Failed Expression: dir == Bond::ENDUPRIGHT || dir == Bond::ENDDOWNRIGHT
 (github issue #5433 from bjonnh-work)
  - MolZip doesn't preserve bond directions when zipping
 (github issue #5450 from bp-kelley)
  - Draw option noAtomLabels and explicit hydrogen works badly
 (github issue #5453 from DavidACosgrove)
  - Fix integer overflow in RGroupDecomposition strategy GA
 (github pull #5460 from bp-kelley)
  - Invalid number of radical electrons calculated for [Pr+4]
 (github issue #5462 from bjonnh-work)
  - CXSmiles isn't properly escaping floating point properties
 (github issue #5466 from bp-kelley)
  - Crossed trans bonds in rings
 (github issue #5486 from DavidACosgrove)
  - MolDraw2D::drawMolecules() should not crash on null molecules
 (github pull #5503 from ptosco)
  - Running kekulization on mols with query bonds will either fail or return incorrect results.
 (github issue #5505 from ricrogz)
  - Regression in the way aldehydes are drawn in current master
 (github issue #5511 from ptosco)
  - Drawing code gives segmentation fault.
 (github issue #5534 from DavidACosgrove)
  - RGD may yield poor depictions
 (github issue #5569 from jones-gareth)
  - Fix a problem when generating conformers for molecules with larger heteroatoms in conjugated 5-rings
 (github pull #5586 from greglandrum)
  - Avoid incurring into division by zero in normalizeDepiction
 (github pull #5619 from ptosco)
  - Allow properties to be displayed in jupyter when 3D rendering is enabled
 (github pull #5624 from greglandrum)
  - Checking whether double bond's controlling atoms are duplicate may cause an infinite loop.
 (github issue #5659 from ricrogz)

## Cleanup work:
  - update release notes, do deprecations
 (github pull #5161 from greglandrum)
  - Fix typo: quarternary --> quaternary
 (github pull #5243 from guyrosin)
  - fix doxygen comments
 (github pull #5254 from e-kwsm)
  - make the catch tests build faster
 (github pull #5284 from greglandrum)
  - make the logging tests more robust
 (github pull #5312 from greglandrum)
  - Update docs for ReplaceSubstructs
 (github pull #5343 from i-tub)
  - Minimallib NPM release fixes - wrong node bin command in Dockerfile and prepublish npm script command replacement
 (github pull #5349 from MichelML)
  - Update code formatting in GettingStartedInPython.rst
 (github pull #5350 from gedeck)
  - fix: rdkit.Chem.rdDistGeom.EmbedMultipleConfs docstring indentation
 (github pull #5404 from jvansan)
  - Remove obsolete files related to rdkit-js + write a rescoped README for MinimalLib
 (github pull #5432 from MichelML)
  - Fixes "DeprecationWarning: invalid escape sequence \C" when importing EnumerateStereoisomers
 (github pull #5478 from ricrogz)
  - Fix double to float downcasting warning
 (github pull #5479 from ricrogz)
  - Remove spurious doxygen tag
 (github pull #5488 from ptosco)
  - Add deprecation warning to docstring.
 (github pull #5498 from DavidACosgrove)
  - Remove unnecessary compiler flags that would be ignored anyway
 (github pull #5587 from ptosco)
  - Update AvalonTools to version 2.0.1
 (github pull #5591 from ptosco)

## New Features and Enhancements:
  - Initial support for non-tetrahedral stereochemistry
 (github pull #5084 from greglandrum)
  - Move to C++17
 (github pull #5155 from greglandrum)
  - Add multi-template isMoleculeXOfReaction overloads
 (github pull #5171 from AlanKerstjens)
  - Add support for Qt6
 (github pull #5203 from ricrogz)
  - Make atom, bond iterators usable in STL algorithms
 (github pull #5204 from ricrogz)
  - Make TautomerQuery serializable
 (github pull #5248 from greglandrum)
  - add boost::serialization support to ROMol
 (github pull #5249 from greglandrum)
  - Add support for Pol and Mod pseudoatoms
 (github pull #5264 from greglandrum)
  - Allow the compiler to generate default ChemicalReaction copy assignment operator
 (github pull #5265 from AlanKerstjens)
  - cdxml parser
 (github pull #5273 from bp-kelley)
  - Expose reaction drawing and additional FPs in MinimalLib
 (github pull #5277 from ptosco)
  - Expose mappings of atom/bond indices abbreviated mol->original mol
 (github pull #5300 from ptosco)
  - Start using string_views in the file parsers
 (github pull #5301 from greglandrum)
  - Add a global feature flag to enable the "new" stereo perception code
 (github pull #5309 from greglandrum)
  - Support conjugated terminal groups in GetBestRMS()
 (github pull #5322 from greglandrum)
  - Open sources Schrodinger's implementation of "molhash"
 (github pull #5360 from d-b-w)
  - Drop usage of CIP information from the canonicalization code
 (github pull #5385 from greglandrum)
  - GetSSSR should have same return type as GetSymmSSSR
 (github issue #5395 from chmnk)
  - Add prop related method into JSMol
 (github pull #5414 from hhk7734)
  - Add draw option to force use of wedge information in MolBlock if pres…
 (github pull #5417 from DavidACosgrove)
  - Add ACS1996 drawing style
 (github pull #5425 from DavidACosgrove)
  - Treat unspecified stereo as unknown
 (github issue #5436 from DavidACosgrove)
  - New version of AssignStereochemistry
 (github pull #5442 from greglandrum)
  - MolZip: add atom property labels for mol-zipping
 (github pull #5446 from bp-kelley)
  - Add atom property (int/unsigned int) to indicate which atoms to MolZip
 (github issue #5451 from bp-kelley)
  - Improved bond highlights
 (github pull #5484 from DavidACosgrove)
  - Include element name in atomic data
 (github pull #5524 from rachelnwalker)
  - Add Substance Groups and Stereo Groups to debugMol()
 (github pull #5526 from ricrogz)
  - make it easy to create DAT SGroups from Python
 (github pull #5544 from greglandrum)
  - Integrating xyz2mol Into The RDKit Core (GSoC 2022)
 (github pull #5557 from gosreya)
  - Switch yaehmop wrapper to basic cmake
 (github pull #5558 from greglandrum)
  - fix warnings ahead of 2022.09 release
 (github pull #5561 from ricrogz)
  - add mergeIsotopes option to mergeQueryHs
 (github pull #5563 from jones-gareth)
  - Add CXSMILES support for bond wedging and wiggly bonds
 (github pull #5575 from greglandrum)
  - Add GetBestTFDBetweenMolecules()
 (github pull #5577 from greglandrum)
  - Add ColorPalette_Vect to SWIG bindings
 (github pull #5580 from jones-gareth)
  - small changes to get the RDKit to build with C++20
 (github pull #5581 from greglandrum)
  - Improvements to 2D depiction and alignment/RMSD calculation
 (github pull #5598 from ptosco)
  - Expose additional functionality to SWIG wrappers
 (github pull #5614 from ptosco)
  - Add an RGroupDecomp aware  molzip to the FreeWilson Contribution
 (github pull #5615 from bp-kelley)
  - PandasTools and InteractiveRenderer improvements
 (github pull #5628 from ptosco)
  - Add updateMolDrawOptionsFromJSON()
 (github pull #5630 from ptosco)
  - InteractiveRenderer.setEnabled() improvements
 (github pull #5635 from ptosco)
  - Support stereo for double bonds in rings from CXSMILES
 (github pull #5636 from greglandrum)
  - add stop condition for arom calc of large ringysystems
 (github pull #5648 from alexwahab)
  - Speed up ring detection by reducing allocations
 (github pull #5654 from d-b-w)
  - Expose highlighAtomColors, highlighBondColors and highlightAtomRadii to CFFI and JS
 (github pull #5657 from ptosco)
  - Update obsolete SWIG definitions
 (github pull #5658 from ptosco)
  - Speed up ring detection by reducing count() calls
 (github pull #5663 from d-b-w)
  - Expose two SubstructUtils functions to SWIG wrappers
 (github pull #5666 from ptosco)
 

## Code removed in this release:
- The C++ class `RDLog::BlockLogs` has been removed. Please use the class `RDLog::LogStateSetter`. The Python class rdBase.BlockLogs() is still available and supported.
- Python function `rdkit.Chem.WrapLogs()` has been removed. Please use
  `rdkit.rdBase.LogToPythonStderr()`. `rdkit.rdBase.WrapLogs()` also exists, but
  unless you need the old teeing behavior, prefer the former.
- Python function `rdkit.Chem.LogWarningMsg()` has been removed. Please use
  `rdkit.rdBase.LogWarningMsg()`.
- Python function `rdkit.Chem.LogErrorMsg()` has been removed. Please use
  `rdkit.rdBase.LogErrorMsg()`.

## Deprecated code (to be removed in a future release):
- The `SmilesParserParams` option `useLegacyStereo` is deprecated and will be
  removed in the 2023.03 release. Please use `SetUseLegacyStereoPerception()`
  instead. In the meantime the SMILES parser will use only use the value of
  `SmilesParserParams.useLegacyStereo` if it is set to `false`, otherwise the
  value of the global `useLegacyStereoPerception` parameter will control the
  behavior of the SMILES parser.
- The following JS methods:
  * generate_aligned_coords()
  * get_morgan_fp()
  * get_morgan_fp_as_uint8array()
  * get_pattern_fp()
  * get_pattern_fp_as_uint8array()
  which used to take several individual parameters now take a single JSON string parameter.
  The overloads taking several individual parameters are now deprecated and will be
  removed in a future release.
- The `PrintAsBase64PNGString` function in `PandasTools` is deprecated and replaced
  by `PrintAsImageString`, which has a more appropriate name given it actually supports
  both PNG and SVG images.

# Release_2022.03.1
(Changes relative to Release_2021.09.1)

## Backwards incompatible changes
- When running in Jupyter Notebook, logs are now sent only to Python's
  standard error stream, and no longer include the `RDKit LEVEL` prefix.
- The Debug and Info logs are now disabled by default. If you would like to
  enable them within your code you can call `rdBase.EnableLog("rdApp.info")`
  and/or `rdBase.EnableLog("rdApp.debug")`.
- The MolHash functions now reassign stereochemistry after modifying the
  molecule and before calculating the hash. Previous versions would still
  include information about atom/bond stereochemistry in the output hash even if
  that no longer applies in the modified molecule.
- The rules for aromaticity in rings containing dummy atoms have been changed.
  The general intention of the new handling is that aromaticity will not be
  perceived for rings containing dummy atoms unless it's clear that the dummies
  should be aromatic. As an example: the SMILES `C1=C*2=CC=CC=*2C=C1` is
  perceived to be aromatic while the SMILES `C1=C*2C=CC=C*2C=C1` (which does not
  have any double bonds to the dummy atoms) is not; in previous RDKit releases
  both of these structures were aromatic. There's more information about this in
  the discussion of PR #4722 (https://github.com/rdkit/rdkit/pull/4722) and
  Issue #4721 (https://github.com/rdkit/rdkit/issues/4721).
- In the PostgreSQL cartridge the `mol_in()` function no longer performs full
  sanitization of the molecule. One consequence of this is that directly casting
  from strings to molecules also no longer does sanitization, so `select 'CN(=O)=O'::mol` 
  does not sanitize the molecule. If you want to convert a string to a molecule
  with full sanitization you can either cast to `text` first 
  (i.e. `select 'CN(=O)=O'::text::mol` or use the `mol_from_smiles()` function.
- The code to calculate bit vector topological torsion fingerprints for
  reactions no longer ignore the fingerprint size argument.
- The rules for tautomer enumeration in `MolStandardize` have been updated to
  more closely match the rules in the original publication. These changes
  primarily consist of making the rules more specific; the consequence is that
  less tautomers will be generated with this version. The previous rules can
  still be accessed via the function `GetV1TautomerEnumerator()` (Python) or
  `getV1TautomerEnumerator()` (C++)

## Highlights
- The RDKit can now integrate with the python logger: calling
  `rdBase.LogToPythonLogger()` enables this. All log messages are sent to a
  logger named "rdkit".
- The backend of the MolDraw2D code has been extensively refactored. This should
  be mostly invisible to RDKit users, but it makes supporting and extending that
  code much easier.
- Beilstein generics (i.e. things like "ARY", "ALK", or "CAL") are now supported
  when doing substructure queries. This is a first step towards enabling some
  really cool substructure search possibilities.

## Acknowledgements
Marcel Baltruschat, Jason Biggs, Kevin Burk, Cédric Bouysset, David Cosgrove,
Joel Duerksen, Jacob Gora, Gareth Jones, Toshiki Kataoka, Eisuke Kawashima,
Brian Kelley, Sonia Liggi, Niels Kristian Kjærgård Madsen, Hector
Martinez-Seara, Dan Nealschneider, Alex Rebert, Ricardo Rodriguez-Schmidt, Steve
Roughley, Roger Sayle, Nikolai Schapin, Ansgar Schuffenhauer, Kaushalesh Shukla,
Jon Sorenson, Ichiru Take, Paolo Tosco, Kazuya Ujihara, Fabio Urbina, Riccardo
Vianello Rachel Walker, Maciej Wójcikowski, SPKorhonen, yuri@FreeBSD,

## Code removed in this release:
- The `useCountSimulation` keyword argument for
  `rdFingerprintGenerator.GetMorganGenerator` and
  `rdFingerprintGenerator.GetAtomPairGenerator` has been removed. Please use the
  `countSimulation` keyword argument instead.
- The function `mol_from_smarts()` in the PostgreSQL cartridge has been removed.
  Please use the `qmol_from_smarts()` function instead.
- The `computeBalabanJ()` functions from the `MolOps` namespace were removed.
  These were not exposed to Python, so this will not affect any Python code.

## Bug Fixes:
  - Fix Flake8 erorrs
 (github pull #4252 from e-kwsm)
  - handle sqlalchemy deprecation
 (github pull #4625 from greglandrum)
  - Debug build of the postgres cartridge fails pg_regress tests for reaction.sql
 (github issue #4631 from rvianello)
  - fix parsing beyond the end of the input string in findMCSsmiles
 (github pull #4636 from rvianello)
  - Mem leak fixes
 (github pull #4637 from ricrogz)
  - Subsequent call to rdChemReactions.ChemicalReaction.RunReactants will block indefinitely.
 (github issue #4651 from goraj)
  - Fix docstring of ConstrainedEmbed
 (github pull #4666 from kazuyaujihara)
  - Postgres cartridge build fails under Ubuntu
 (github issue #4681 from SPKorhonen)
  - Molfile SDD records not properly displayed
 (github pull #4690 from jones-gareth)
  - RGD:  fix for cores with MOL  block atom lists
 (github pull #4695 from jones-gareth)
  - Wrong tautomers generated
 (github issue #4700 from sonial)
  - RGD align output core to input structure
 (github pull #4709 from jones-gareth)
  - TorsionFingerprints raises error with S(Cl)F4 group
 (github issue #4720 from kazuyaujihara)
  - Dummy atoms next to aromatic are always kekulized even when they should not
 (github issue #4721 from ptosco)
  - TautomerEnumerator will crash if copied with a callback set
 (github issue #4736 from ptosco)
  - Minor PandasTools cleanup
 (github pull #4744 from ptosco)
  - Reaction parser fails when CX extensions are present
 (github issue #4759 from greglandrum)
  - GetSimilarityMapFromWeights changes behavior of parameter "colorMap" depending on whether the parameter "draw2d" is provided or not
 (github issue #4763 from FabioUrbina)
  - Highlight bond width is different for different PNG image sizes
 (github issue #4764 from rachelnwalker)
  - Fixes crashing bug with finalSubstructChecks
 (github pull #4782 from greglandrum)
  - MDL query with aromatic bond sets aromatic flag on atoms even though they are not in an aromatic ring
 (github issue #4785 from ptosco)
  - [Cartridge]: qmol_from_ctab and qmol_from_smiles are sanitizing molecules
 (github issue #4787 from greglandrum)
  - AdjustQueryProperties() is inconsistent when adjustConjugatedFiveRings is set
 (github issue #4789 from greglandrum)
  - `Draw.MolToFile` to SVG raises "Can't kekulize" error
 (github issue #4792 from toslunar)
  - Ring double bonds written as crossed bonds after RGD
 (github issue #4809 from greglandrum)
  - Fix a number of crashing bugs in the python wrappers
 (github pull #4810 from greglandrum)
  - correctly tag unspecified branch-begin bonds in SMARTS
 (github pull #4811 from greglandrum)
  - AttributeError in PandasTools
 (github issue #4821 from greglandrum)
  - reloading PandasTools leads to infinite recursion
 (github issue #4823 from greglandrum)
  - ReplaceCore should set stereo on ring bonds when it breaks rings
 (github issue #4825 from greglandrum)
  - MolDraw2D::drawArc() starts at wrong angle
 (github issue #4836 from greglandrum)
  - MolDraw2D::drawArc() not exposed to Python
 (github issue #4837 from greglandrum)
  - align argument to MolDraw2D::DrawString() cannot be used from Python
 (github issue #4838 from greglandrum)
  - Fix RunFilterCatalog() thread counts.
 (github pull #4856 from xavierholt)
  - RGD: dummy atom in input structure is mishandled
 (github pull #4863 from jones-gareth)
  - Fix bug with wedges being drawn backwards
 (github pull #4868 from greglandrum)
  - Missing dependency on RDKit::RingDecomposerLib_static in RDKit::GraphMol_static
 (github issue #4875 from nielskm)
  - cannot parse coordinate bonds from CXSMARTS
 (github issue #4878 from greglandrum)
  - fix some leaks in the SWIG wrappers
 (github pull #4916 from greglandrum)
  - Fix some problems turned up by ossfuzz
 (github pull #4927 from greglandrum)
  - Fix i-files for RDK_USE_BOOST_IOSTREAMS=OFF
 (github pull #4933 from kazuyaujihara)
  - Fails on i386: non-constant-expression cannot be narrowed from type 'unsigned int' to 'npy_intp' (aka 'int') in initializer list
 (github issue #4934 from yurivict)
  - Fix SWIG wrappers for C#
 (github pull #4935 from kazuyaujihara)
  - Definition of eccentricity in documentation is wrong
 (github issue #4952 from drkeoni)
  - fix CMakeLists.txt extracting link line from python LDSHARED config var
 (github pull #4954 from rvianello)
  - add JavaGzStreamTests
 (github pull #4973 from kazuyaujihara)
  - Invalid SMARTS generated by MolToSmarts
 (github issue #4981 from nielskm)
  - Fix memory safety issues found by OSS-Fuzz
 (github pull #4983 from alpire)
  - AssignStereochemistry should remove nonchiral atoms from StereoGroups
 (github pull #4986 from greglandrum)
  - Transform3D#SetRotation() around arbitrary axis scales coordinates
 (github issue #4995 from sroughley)
  - Bad handling of dummy atoms in the CIP assignment code
 (github issue #4996 from greglandrum)
  - Bad handling of fragments in CIP code
 (github issue #4998 from greglandrum)
  - Drawing query atoms containing an AND query raises an exception
 (github issue #5006 from ptosco)
  - Bad tautomers produced for phosphorous compounds
 (github issue #5008 from NikSchap2107)
  - Fix warning when generating coordinates for ZOBs
 (github pull #5011 from d-b-w)
  - Code in the docstring for `FindMolChiralCenters()` doesn't work
 (github pull #5014 from greglandrum)
  - Mol images in DataFrames are drawn only once in Jupyter Lab
 (github issue #5017 from mrcblt)
  - Drop `gist_qmol_ops` in upgrade scripts in case it exists
 (github pull #5021 from mwojcikowski)
  - Remove extra newline from Kekulize error message.
 (github pull #5022 from xavierholt)
  - Neighboring Hs not taken into account in connectivity invariants
 (github issue #5036 from greglandrum)
  - smiles parsing error due to erroneous ring perception
 (github issue #5055 from AnsgarSchuffenhauer)
  - To INCHI conversion leaks on kekulization failure
 (github pull #5057 from ricrogz)
  - Add a CXSMILES option to the MolHash
 (github pull #5058 from greglandrum)
  - Make the RGD code work when `rgroupLabelling` is `Isotope`
 (github pull #5088 from greglandrum)
  - Compilation issue with catch.hpp
 (github issue #5089 from hseara)
  - pg_restore: error: COPY failed for table "mols": ERROR:  could not create molecule from SMILES
 (github issue #5095 from joelduerksen)
  - Removing H preserving only wedged ones strips all H
 (github issue #5099 from ricrogz)
 - fix a mistake in the enhanced stereochemistry substructure table
 (github issue #5101 from greglandrum)
  - NumRotatableBonds() incorrect for partially sanitized molecule
 (github issue #5104 from greglandrum)
  - Wiggly bonds don't override wedged bonds
 (github issue #5108 from greglandrum)

## Cleanup work:
  - Do the deprecations for the 2022.03  release
 (github pull #4626 from greglandrum)
  - clang-tidy: readability-simplify-boolean-expr
 (github pull #4639 from e-kwsm)
  - Clean-up Python 4815 - P1.1: Chem\AtomPairs
 (github pull #4859 from IchiruTake)
  - Clean-up Python #4815 - P1.2: Chem/ChemUtils
 (github pull #4860 from IchiruTake)
  - Clean-up Python #4815 - P1.3: Chem\Draw
 (github pull #4891 from IchiruTake)
  - Clean-up Python #4815 - P1.4: Chem\EState
 (github pull #4893 from IchiruTake)
  - Clean-up Python #4815 - P1.5: Chem\FeatMaps
 (github pull #4894 from IchiruTake)
  - Clean-up Python #4815 - P1.6: Chem\Features
 (github pull #4896 from IchiruTake)
  - Clean-up Python #4815 - P1.8: Chem\fmcs
 (github pull #4898 from IchiruTake)
  - Clean-up Python #4815 - P1.9: Chem\Fraggle
 (github pull #4906 from IchiruTake)
  - Clean-up Python #4815 - P1.10: Chem\MolDb
 (github pull #4907 from IchiruTake)
  - Clean-up Python #4815 - P1.11: Chem\MolKey
 (github pull #4910 from IchiruTake)
  - Clean-up Python #4815 - P1.12: Chem\MolStandardize
 (github pull #4911 from IchiruTake)
  - Clean-up Python #4815 - P1.13: Chem\Pharm2D & Chem\Pharm3D
 (github pull #4912 from IchiruTake)
  - Clean-up Python #4815 - P1.14: Chem\Scaffolds & Chem\SimpleEnum
 (github pull #4913 from IchiruTake)
  - Run clang-tidy (readability-braces-around-statements)
 (github pull #4977 from e-kwsm)
  - silence warnings in MSVC compliatons
 (github pull #5044 from bp-kelley)
  - Clean up the warning landscape
 (github pull #5048 from greglandrum)
  - Cleanup of python API documentation stubs
 (github pull #5105 from greglandrum)


## New Features and Enhancements:
  - Update coordgenlibs to v3.0.0
 (github pull #4638 from ricrogz)
  - Fix some compile-time warnings in the postgres cartridge code
 (github pull #4657 from rvianello)
  - Add some new color palettes to MolDraw2D
 (github pull #4668 from greglandrum)
  - Add support for Beilstein generics when doing substructure queries
 (github pull #4673 from greglandrum)
  - Remove unnecessary mutex in InChI wrapper
 (github pull #4680 from greglandrum)
  - Update mac CI builds
 (github pull #4738 from greglandrum)
  - Remove dead code
 (github pull #4739 from ptosco)
  - Refactor the memory management of the postgres cartridge cache module
 (github pull #4755 from rvianello)
  - Improve CMake integration of PgSQL build
 (github pull #4767 from ptosco)
  - Allow MolDraw2DCairo and MolDraw2DSVG to determine canvas size based on the molecule
 (github pull #4772 from greglandrum)
  - generate the sql update files in the binary (build) directory
 (github pull #4777 from rvianello)
  - Allow using heavy atoms only in the FragmentChooser
 (github pull #4791 from ptosco)
  - silence warnings in MSVC compilations
 (github pull #4796 from bp-kelley)
  - Support using the python logger
 (github issue #4840 from xavierholt)
  - Add support for quadruple bonds in SMILES
 (github issue #4842 from jasondbiggs)
  - Some refactoring of the depictor code
 (github pull #4865 from greglandrum)
  - Build documentation for 3 missing modules
 (github pull #4879 from ptosco)
  - Add access to query atom symbols from python or an option to render them in images
 (github issue #4880 from rachelnwalker)
   - Start adding move constructors and move-assignment operators
 (github pull #4909 from greglandrum)
  - Refactor mol draw2 d
 (github pull #4948 from DavidACosgrove)
  - [ENH]: Support greater use of `findAtomEnvironmentOfRadiusN()`
 (github pull #4970 from IchiruTake)
  - Move isEarlyAtom to a table to reduce lock contention in getPeriodicTable
 (github pull #4980 from bp-kelley)
  - Support writing V3000 reactions
 (github pull #4982 from greglandrum)
  - Make FMCS check bond stereo.
 (github pull #5009 from DavidACosgrove)
  - Improving atom colors for dark mode.
 (github pull #5038 from kaushaleshshukla)
  - enable the multithreaded LeaderPicker on linux
 (github pull #5043 from greglandrum)
  - Expose MolzipParams::atomSymbols to python
 (github pull #5054 from bp-kelley)
  - disable Info and Debug logs by default
 (github pull #5065 from greglandrum)
  - Add sanitize option to molzip
 (github pull #5069 from bp-kelley)
  - "Powered by RDKit" Badge
 (github pull #5085 from cbouy)
  - Add a couple of depiction helper functions and some JS bindings
 (github pull #5115 from ptosco)
  - Swig MolDraw2D cairo
 (github pull #5128 from jones-gareth)
  - Enables rdkit-structure-renderer.js in Jupyter Lab and Notebook
 (github pull #5132 from ptosco)

## Deprecated code (to be removed in a future release):
- Python function `rdkit.Chem.WrapLogs()` is deprecated in favor of
  `rdkit.rdBase.LogToPythonStderr()`.  `rdkit.rdBase.WrapLogs()` also exists,
  but unless you need the old teeing behavior, prefer the former.
- Python function `rdkit.Chem.LogWarning()` is deprecated in favor of
  `rdkit.rdBase.LogWarning()`.
- Python function `rdkit.Chem.LogError()` is deprecated in favor of
  `rdkit.rdBase.LogError()`.
- The C++ class `RDLog::BlockLogs` is deprecated in favor of the the class `RDLog::LogStateSetter`.

# Release_2021.09.1
(Changes relative to Release_2021.03.1)

## Backwards incompatible changes
- `RWMol.replaceAtom()` no longer removes `SubstanceGroups` which reference that atom.
- The `keepSGroups` argument to `RWMol.replaceBond()` now defaults to true.
- The SMARTS parser now by default accepts CXSMILES extensions and molecule
  names. SMARTS which previously failed to parse like `CCC fail` will now return
  valid molecules.
- Molecule names in SMILES and SMARTS are now parsed by default. Previously they
  were ignored.
- The `getParams()` function for retrieving UFF parameters now returns a const
  pointer instead of a standard pointer. This shouldn't affect the functionality
  of any code since the only method of the class is also const.

## Highlights
 - Single reactant/single product reactions can now be applied in-place. This
   really helps with the performance of these kinds of transformations.
 - The CFFI wrapper around MinimalLib greatly expands the number of possible
   places the RDKit can be used from.
 - A number of general enhancements and quality-of-life improvements were made
   to the PostgreSQL cartridge.

## Acknowledgements
Jason Biggs, Kit Choi, David Cosgrove, Eloy Félix, Harrison Green, Gareth Jones,
Eisuke Kawashima, Alan Kerstjens, Brian Kelley, John Konecny, Stephanie
Labouille, Rasmus Lundsgaard, Hadrien Mary, Michel Moreau, Dan Nealschneider,
Axel Pahl, Maximilian Peters, Alessio Ragno, Ricardo Rodriguez-Schmidt, Riccardo
Sabatini, Roger Sayle, Vincent F. Scalfani, Dan Skatov, David Slochower, Peter
St. John, Mihaly Szabo, Ichiru Take, Paolo Tosco, Ivan Tubert-Brohman, Kazuya
Ujihara, Alain Vaucher, Riccardo Vianello, Rachel Walker, Shuzhe Wang, Maciej
Wójcikowski, bzoracler, jungb-basf, charly828, yoer77,

## Code removed in this release:
- The minimizeOnly option for coordgen has been removed.

## Contrib updates:
 - Contribute FreeWilson analysis
 (github pull #4026 from bp-kelley)

## Bug Fixes:
  - cannot pass drawOptions to MolsToGridImage when using notebook
 (github issue #3101 from slochower)
  - Draw.MolToImage() cannot highlight with highlightMap (v > '2019.09.3' )
 (github issue #3616 from spideralessio)
  - EnumerateStereoisomers fail with STEREOANY bonds from molblock
 (github issue #3759 from TermeHansen)
  - Double bond stereo gets flipped by SMILES reader/writer
 (github issue #3967 from mwojcikowski)
  - SparseIntVect copy constructor and assignment operators not clearing existing data
 (github issue #3994 from AlanKerstjens)
  - MolFragmentToSmiles with kekuleSmiles=True raises AtomKekulizeException
 (github issue #3998 from kazuyaujihara)
  - update clang version for linux CI fuzzer builds
 (github pull #4012 from greglandrum)
  - Update coordgen to 2.0.3
 (github pull #4017 from d-b-w)
  - Get SWIG wrappers working with C# again
 (github pull #4020 from kazuyaujihara)
  - replaceSidechains creates aromatic dummy atoms
 (github pull #4022 from ptosco)
  - A set of fixes for problems caused by bad input
 (github pull #4033 from greglandrum)
  - Cleanup some problems found during an ASAN build
 (github pull #4054 from greglandrum)
  - Avoid that lone atoms which are part of a ring in one of the molecules become part of the MCS
 (github pull #4065 from ptosco)
  - StereoGroups not preserved by RenumberAtoms() 
 (github issue #4071 from greglandrum)
  - call to pyAvalonTools.Generate2DCoords results in an assert violation
 (github issue #4075 from rvianello)
  - Update boost download location in Dockerfile
 (github pull #4094 from greglandrum)
  - HCount field in v2000 Mol blocks ignored
 (github issue #4099 from riccardosabatini)
  - Reactions don't propagate bond properties
 (github issue #4114 from d-b-w)
  - RemoveStereochemistry should also remove stereogroups
 (github issue #4115 from greglandrum)
  - Avoid that MolStandardizer segfaults on empty mols
 (github pull #4119 from ptosco)
  - SEGV in RWMol::commitBatchEdit
 (github issue #4122 from hgarrereyn)
  - SEGV in ROMol::getAtomDegree if atom is not in graph
 (github issue #4127 from hgarrereyn)
  - SEGV from unsigned integer overflow in Conformer::setAtomPos
 (github issue #4128 from hgarrereyn)
  - HCOUNT from v3000 CTABS incorrectly interpreted
 (github issue #4131 from greglandrum)
  - Empty query produces empty match, but at the same time is considered non-matching
 (github issue #4138 from i-tub)
  - fixed AddBond documentation
 (github pull #4142 from Ashafix)
  - Possible bug with `EnumerateStereoisomers`
 (github issue #4144 from stephanielabouille)
  - Odd drawing behavior with radicals and MolsToGridImage
 (github issue #4156 from pstjohn)
  - pre-condition violation when sanitizing a de-pickled reaction
 (github issue #4162 from jasondbiggs)
  - Many of the PMI descriptors are not being recalculated for different conformers
 (github issue #4167 from greglandrum)
  - bug in MDLParser.cpp when reading a rxn file in v3000 format that contains agents
 (github issue #4183 from jungb-basf)
  - Potentially chiral bridgehead atoms not being identified.
 (github pull #4192 from greglandrum)
  - allow more recoverable V3000 parsing errors when strictParsing=false
 (github pull #4210 from greglandrum)
  - RGD: Fix memory leak with deleting array
 (github pull #4211 from bp-kelley)
  - UnfoldedRDKFingerprintCountBased returns a different fingerprint length for every molecule
 (github issue #4212 from greglandrum)
  - rdMolHash.MolHash fails on non-standard valences
 (github issue #4222 from ricrogz)
  - Fix a couple of problems with fingerprint count simulation
 (github pull #4228 from greglandrum)
  - Chem.MolFromSmiles using SmilesParserParams throws exceptions
 (github issue #4232 from greglandrum)
  - Parse failure for data groups in CXSMILES
 (github issue #4233 from greglandrum)
  - double bonds now have EITHER stereo if no coordinates are present
 (github pull #4239 from greglandrum)
  - Fix CMakeLists for FileParsers
 (github pull #4240 from kazuyaujihara)
  - Multiple ATTCHPT entries for one atom handled incorrectly
 (github issue #4256 from greglandrum)
  - Exception thrown by reionizer when dealing with Mg+2
 (github issue #4260 from greglandrum)
  - Fallback ring finding failing on molecules with multiple fragments
 (github issue #4266 from avaucher)
  - Make sure that ResonanceMolSupplier substructure matches are uniquified consistently
 (github pull #4274 from ptosco)
  - FindPotentialStereo() doesn't find *marked* ring stereo when flagPossible=False
 (github issue #4279 from greglandrum)
  - The normalization pattern for pyridine N-oxide is not specific enough
 (github issue #4281 from ptosco)
  - computeCanonicalTransform may generate non-canonical coords
 (github issue #4302 from ptosco)
  - Unreasonable calculation of implicit valence for atoms with query bonds
 (github issue #4311 from greglandrum)
  - call to AvalonTools::set2DCoords results in an assert violation
 (github issue #4330 from jasondbiggs)
  - MolBlock writer gives non-stereo double bonds "unspecified" parity
 (github issue #4345 from d-b-w)
  - Specified trans stereo being ignored during conformation generation in macrocycles
 (github issue #4346 from greglandrum)
  - Two MinGW build fixes and one MSVC build fix
 (github pull #4347 from ptosco)
  - Fixes RDK_BUILD_THREADSAFE_SSS=OFF build
 (github pull #4349 from ptosco)
  - clean up some leaks identified by an ASAN build
 (github pull #4354 from greglandrum)
  - Three more Windows build fixes
 (github pull #4356 from ptosco)
  - Specified grid spacing for ShapeTanimotoDistance is ignored.
 (github issue #4364 from greglandrum)
  - Need implicit H cleanup after rdMolEnumerator.Enumerate()
 (github issue #4381 from greglandrum)
  - rdMolEnumerator.Enumerate fails on variable attachment points with queries
 (github issue #4382 from greglandrum)
  - RDKit reaction produces wrong double bond stereochemistry
 (github issue #4410 from mwojcikowski)
  - "to-Python converter already registered; second conversion method ignored." warnings issued at import
 (github issue #4425 from ricrogz)
  - v2000 SGroups do not generate an "index" property
 (github issue #4434 from ricrogz)
  - pg_restore does not work with some mol type molecule
 (github issue #4442 from mihalyszabo88)
  - Building with static dependencies breaks CMake exports
 (github issue #4449 from ricrogz)
  - DataStruct vectors leak when iterating
 (github issue #4465 from bp-kelley)
  - SGroups: Additional SDT properties not decoded if FIELDNAME is empty
 (github issue #4476 from greglandrum)
  - Test failure in reaction.sql
 (github issue #4486 from yoer77)
  - Small rings can have STEREOANY/EITHERDOUBLE bonds
 (github issue #4494 from ricrogz)
  - OR queries involving aromatic atoms cannot be drawn
 (github issue #4496 from ptosco)
  - MolFromSmiles and MolFromSmarts incorrectly accepting input with spaces
 (github issue #4503 from greglandrum)
  - Native 2D layout engine may generate overlapping coordinates
 (github issue #4504 from ptosco)
  - SubstanceGroup labels sometimes overlap with atoms in image generation
 (github issue #4508 from rachelnwalker)
  - SGroups do not have a way of unsetting properties from Python
 (github issue #4514 from ricrogz)
  - operator<< is declared for AtomPDBResidueInfo but not defined
 (github issue #4535 from greglandrum)
  - Improve test coverage and some bug fixes
 (github pull #4536 from greglandrum)
  - Seg fault in MolDraw2D::drawMolecules()
 (github issue #4538 from greglandrum)
  - Salt removal forces sanitization
 (github issue #4550 from ricrogz)
  - fix a thread-safety bug in the UFF parameter loading
 (github pull #4553 from greglandrum)
  - GetSubstructMatches() loops at 43690 iterations.
 (github issue #4558 from ricrogz)
  - failure to parse CTAB with LINKNODE and SGROUP
 (github issue #4561 from greglandrum)
  - Requesting "density" fingerprint Hydrogen molecule fails with exception
 (github issue #4567 from ricrogz)
  - Incorrect double bond stereo in output SMILES around ring closures
 (github issue #4582 from greglandrum)
  - migrate the MHFP implementation to use boost::random
 (github pull #4603 from rvianello)
  - Fix EnumerateStereoisomers with tryEmbedding
 (github pull #4615 from kazuyaujihara)

## New Features and Enhancements:
  - Support Chemical Markup Language, CML, for writing
 (github pull #3024 from e-kwsm)
  - Add Eigen to ExternalProject and automatically download if RDK_BUILD_DESCRIPTORS3D
 (github pull #3075 from e-kwsm)
  - updates to postgreSQL cartridge
 (github pull #3976 from greglandrum)
  - Update Overview.md
 (github pull #3992 from charly828)
  - MinimalLib: add CFFI interface
 (github pull #4018 from greglandrum)
  - Contribute FreeWilson analysis
 (github pull #4026 from bp-kelley)
  - Allow partial deserialization of molecules
 (github pull #4040 from greglandrum)
  - Add datamol project
 (github pull #4046 from hadim)
  - Build BCUT when RDK_BUILD_DESCRIPTORS3D=OFF
 (github pull #4085 from kazuyaujihara)
  - Making RDKit minimallib (JS lib) available through the npm package manager
 (github pull #4086 from MichelML)
  - Normalize line endings in source code files
 (github pull #4104 from ptosco)
  - Allow MolToQPixmap to support PySide2
 (github pull #4110 from kitchoi)
  - update ChEMBL projects in Projects using RDKit
 (github pull #4116 from eloyfelix)
  - A collection of MolStandardize improvements
 (github pull #4118 from greglandrum)
  - Run clang-format against header files
 (github pull #4143 from e-kwsm)
  - Some miscellaneous MinimalLib enhancements
 (github pull #4169 from ptosco)
  - [MinimalLib] Add number of heavy atoms to descriptors
 (github issue #4184 from apahl)
  - Comments added to RGD core matching
 (github pull #4189 from jones-gareth)
  - Fix/rdprop integer conversions
 (github pull #4194 from bp-kelley)
  - RWMol cleanup
 (github pull #4198 from greglandrum)
  - Test comparing SVGs via hash code - ready for review
 (github pull #4199 from DavidACosgrove)
  - support getNumAtoms and getNumHeavyAtoms as Descriptors
 (github pull #4200 from greglandrum)
  - Canon.cpp canonicalizeDoubleBond function refactor
 (github pull #4204 from jfkonecn)
  - Add low level functions to bulk-update Substance Group atoms & bonds
 (github pull #4206 from ricrogz)
  - Run clang-tidy (modernize-pass-by-value)
 (github pull #4224 from e-kwsm)
  - Shift `Trajectory` and `Snapshot` constructors to methods on classes
 (github pull #4225 from bzoracler)
  - Allow depiction of "either" double bonds as "wiggly neighbors"
 (github issue #4238 from d-b-w)
  - Some cartridge enhancements
 (github pull #4271 from greglandrum)
  - [Enhancement]: Allow every bit-vect fixed-size fingerprints can be directly embedded || attached into the (Numpy) Array or return the Numpy Array
 (github issue #4273 from IchiruTake)
  - Turn MRV_COORDINATE_BOND_TYPE data Substance Groups into coordinate bonds
 (github pull #4299 from ricrogz)
  - Support using SubstructMatchParameters in RGD
 (github pull #4318 from greglandrum)
  - Add partial CXSMARTS support
 (github issue #4319 from greglandrum)
  - Support toggling components of CXSMILES output
 (github issue #4320 from greglandrum)
  - Switch to using InChI v1.06
 (github issue #4322 from greglandrum)
  - support using RGBA colors
 (github issue #4323 from greglandrum)
  - MinimalLib: return fingerprints as BitSets
 (github issue #4329 from dskatov)
  - Expose atomColourPalette as JSON drawOption
 (github pull #4337 from ptosco)
  - Improved sgroup output
 (github pull #4343 from greglandrum)
  - support get_json in the rdkitjs wrapper
 (github pull #4348 from greglandrum)
  - Enable using the URF library in Windows static builds
 (github pull #4357 from ptosco)
  - a few doxygen comment fixes
 (github pull #4368 from jasondbiggs)
  - Enable building Java wrappers with MinGW compilers
 (github pull #4384 from ptosco)
  - add exactmw to cartridge
 (github issue #4386 from greglandrum)
  - Some cartridge additions and fixes
 (github pull #4387 from greglandrum)
  - Remove SWIG kludge on Windows
 (github pull #4388 from ptosco)
  - CXSMILES improvements
 (github pull #4396 from greglandrum)
  - SubstructLibrary improvements
 (github pull #4403 from greglandrum)
  - Add 3 new examples to Cookbook.
 (github pull #4404 from vfscalfani)
  - cleanup the use of lambdas in the code
 (github pull #4432 from greglandrum)
  - Swap from RDUNUSED_PARAM to unnamed parameters
 (github pull #4433 from greglandrum)
  - Fix #4442 and other cartridge improvements
 (github pull #4448 from greglandrum)
  - RDKit Cartridge: qmol GiST support
 (github issue #4463 from mwojcikowski)
  - Add ToList method to Sparse/ExplicitBitVector
 (github pull #4467 from bp-kelley)
  - GiST support to qmol type
 (github pull #4470 from mwojcikowski)
  - be more tolerant when kekulizing
 (github pull #4492 from greglandrum)
  - Allow applying single-reactant/single-product reactions in place
 (github pull #4511 from greglandrum)
  - Add custom distance bounds parameter for ETKDG conformer generation
 (github pull #4516 from hjuinj)
  - cleanup some compiler warnings
 (github pull #4521 from greglandrum)
  - Improve test coverage and some bug fixes
 (github pull #4536 from greglandrum)
  - Additions to the jupyter integration
 (github pull #4541 from greglandrum)
  - another round of cartridge improvements
 (github pull #4543 from greglandrum)
  - Major speed-up of RGD scoring
 (github pull #4544 from ptosco)
  - Expose get_smarts to JS
 (github pull #4547 from ptosco)
  - Improve performance of removing substruct/tautomer duplicates
 (github pull #4560 from ricrogz)
  - Add support for SRUs to MolEnumerator
 (github pull #4563 from greglandrum)
  - Adds KeyFromPropHolder to hold user defined indices
 (github pull #4571 from bp-kelley)
  - add ROMol::atomNeighbors() and ROMol::atomBonds()
 (github pull #4573 from greglandrum)
  - Improvements to some constructors in python wrappers
 (github pull #4581 from greglandrum)
  - Add FreeSASA support to Windows builds
 (github pull #4584 from ptosco)
  - use V3K mol blocks in PNG metadata
 (github pull #4588 from greglandrum)
  - Prevent some loop variables from creating unnecessary copies
 (github pull #4610 from rvianello)
  - Rename the molLinkNode property to _molLinkNode
 (github pull #4614 from greglandrum)
  - Fix clang warning `-Wabsolute-value`
 (github pull #4616 from e-kwsm)

## Deprecated code (to be removed in a future release):
- The `useCountSimulation` keyword argument for
  `rdFingerprintGenerator.GetMorganGenerator` and
  `rdFingerprintGenerator.GetAtomPairGenerator` has been deprecated and will be
  removed in the next release. Please use the `countSimulation` keyword argument
  instead.
- The function `mol_from_smarts()` in the PostgreSQL cartridge has been
  deprecated and will be removed in the next release. Please use the
  `qmol_from_smarts()` function instead.
- The `computeBalabanJ()` functions from the `MolOps` namespace have been
  deprecated and will be removed in the next release. These have not been
  exposed to Python, so this will not affect any Python code.


# Release_2021.03.1
(Changes relative to Release_2020.09.1)

## Backwards incompatible changes
- The distance-geometry based conformer generation now by defaults generates
  trans(oid) conformations for amides, esters, and related structures. This can
  be toggled off with the `forceTransAmides` flag in EmbedParameters. Note that
  this change does not impact conformers created using one of the ET versions.
  (#3794)
- The conformer generator now uses symmetry by default when doing RMS pruning.
  This can be disabled using the `useSymmetryForPruning` flag in
  EmbedParameters. (#3813)
- Double bonds with unspecified stereochemistry in the products of chemical
  reactions now have their stereo set to STEREONONE instead of STEREOANY (#3078)
- The MolToSVG() function has been moved from rdkit.Chem to rdkit.Chem.Draw
  (#3696)
- There have been numerous changes to the RGroup Decomposition code which change
  the results. (#3767)
- In RGroup Decomposition, when onlyMatchAtRGroups is set to false, each molecule
  is now decomposed based on the first matching scaffold which adds/uses the
  least number of non-user-provided R labels, rather than simply the first
  matching scaffold.
  Among other things, this allows the code to provide the same results for both
  onlyMatchAtRGroups=true and onlyMatchAtRGroups=false when suitable scaffolds
  are provided without requiring the user to get overly concerned about the
  input ordering of the scaffolds. (#3969)
- There have been numerous changes to `GenerateDepictionMatching2DStructure()` (#3811)
- Setting the kekuleSmiles argument (doKekule in C++) to MolToSmiles will now
  cause the molecule to be kekulized before SMILES generation. Note that this
  can lead to an exception being thrown. Previously this argument would only
  write kekulized SMILES if the molecule had already been kekulized (#2788)
- Using the kekulize argument in the MHFP code will now cause the molecule to be
  kekulized before the fingerprint is generated. Note that becaues kekulization
  is not canonical, using this argument currently causes the results to depend
  on the input atom numbering. Note that this can lead to an exception being
  thrown. (#3942)
- Gradients for angle and torsional restraints in both UFF and MMFF were computed
  incorrectly, which could give rise to potential instability during minimization.
  As part of fixing this problem, force constants have been switched to using
  kcal/degree^2 units instead of kcal/rad^2 units, consistently with the fact that
  angle and dihedral restraints are specified in degrees. (#3975)

## Highlights
- MolDraw2D now does a much better job of handling query features like common
  query bond types, atom lists, variable attachment points, and link nodes. It
  also supports adding annotations at the molecule level, displaying brackets
  for Sgroups, rendering the ABS flag for stereochemistry, and a new "comic"
  mode.
- There are two new Contrib packages: NIBRStructureFilters and CalcLigRMSD
- There have been a number of improvements made to the R-Group Decomposition
  code which make it both more flexible and considerably faster

## Acknowledgements
Michael Banck, Christopher Von Bargen, Jason Biggs, Patrick Buder, Ivan
Chernyshov, Andrew Dalke, Xiaorui Dong, Carmen Esposito, Nicholas Firth, Enrico
Gandini, James Gayvert, Gareth Jones, Eisuke Kawashima, Steven Kearnes, Brian
Kelley, Mark Mackey, Niels Kristian Kjærgård Madsen, Luca Naef, Dan
Nealschneider, Jin Pan, Daniel Paoliello, António JM Ribeiro, Sereina Riniker,
Braxton Robbason, Jaime Rodríguez-Guerra, Ricardo Rodriguez-Schmidt, Steve
Roughley, Vincent F. Scalfani, Nadine Schneider, Philippe Schwaller, Dan Skatov,
Pascal Soveaux, Paolo Tosco, Kazuya Ujihara, Riccardo Vianello, Shuzhe Wang,
Piotr Wawrzyniak, Maciej Wójcikowski, Zhijiang Yang, Yutong Zhao
'driesvr', 'GintasKam', 'SPKorhonen', 'pkubaj', 'AnPallo', 'darintay',
'slchan3', 'Robins', 'sirbiscuit', 'amateurcat', 'noncomputable', 'yurivict',
'magattaca'

## Contrib updates:
  - Added NIBRStructureFilters: a set of substructure filters for hit-list triaging together with python code for applying them. The filters are described in the publication https://dx.doi.org/10.1021/acs.jmedchem.0c01332
   (github pull #3516 from NadineSchneider)
  - Added CalcLigRMSD: flexible python code for calculating RMSD between pre-aligned molecules
   (github pull #3812 from cespos)

## Bug Fixes:
  - Casting int to uint in MorganFingerprintHelper leads to unexpected behaviour.
 (github issue #1761 from SiPa13)
  - MolChemicalFeature.GetPos() returns value for molecule's default conformer
 (github issue #2530 from greglandrum)
  - Unable to catch RDKit exceptions in linked library when compiling with fvisibility=hidden
 (github issue #2753 from cdvonbargen)
  - Reaction rendering always shows molecules in aromatic form
 (github issue #2976 from greglandrum)
  - Reactions setting unspecified double-bond stereo to STEREOANY
 (github issue #3078 from ricrogz)
  - PDB output flavor&2 documentation change
 (github issue #3089 from adalke)
  - WedgeMolBonds() should prefer degree-1 atoms
 (github issue #3216 from greglandrum)
  - Error in ChemAxon SMILES "parsing"
 (github issue #3320 from IvanChernyshov)
  - Incorrect number of radical electrons calculated for metals
 (github issue #3330 from greglandrum)
  - Problem with lifetime linkage of mols and conformers
 (github issue #3492 from amateurcat)
  - Traceback when pickling ROMol after BCUT descriptors are calculated
 (github issue #3511 from d-b-w)
  - Fix AUTOCORR2D descriptors
 (github pull #3512 from ricrogz)
  - SDMolSupplier requires several attempts to load a SDF file in Python 3.6/3.7
 (github issue #3517 from jaimergp)
  - Remove accidentally included boost header
 (github pull #3518 from ricrogz)
  - legend_height_ should be preserved after drawing the molecule
 (github pull #3520 from greglandrum)
  - remove the include directive for unused <access/tuptoaster.h> header
 (github pull #3525 from rvianello)
  - C++ build fails when configured with RDKIT_USE_BOOST_SERIALIZATION=OFF
 (github issue #3529 from rvianello)
  - Newest RDKIT version allowing chemically invalid smiles
 (github issue #3531 from GintasKam)
  - Behaviour of generate_aligned_coords for erroneous inputs
 (github issue #3539 from dskatov)
  - Drawing artifacts in draw_to_canvas_with_offset
 (github issue #3540 from dskatov)
  - Error adding PNG metadata when kekulize=False
 (github issue #3543 from gayverjr)
  - Add missing methods to remove SubstanceGroup attributes
 (github pull #3547 from greglandrum)
  - Error writing SDF data containing UTF-8 to a StringIO object
 (github issue #3553 from greglandrum)
  - correct handling of amide distances for macrocycles
 (github pull #3559 from hjuinj)
  - rdMolDraw2D, problems during generation of pictures from SMARTS, differences between Cairo and SVG
 (github issue #3572 from wopozka)
  - Fix example of SmilesToMol
 (github pull #3575 from kazuyaujihara)
  - atom/bond notes handle capital letters incorrectly
 (github issue #3577 from greglandrum)
  - Get MolDraw2DQt working again
 (github pull #3592 from greglandrum)
  - Scientific notation in SDF V3000 files
 (github issue #3597 from mark-cresset)
  - Fix: add missing python wrappers for MolDraw2DQt
 (github pull #3613 from greglandrum)
  - V3K mol block parser not saving the chiral flag
 (github issue #3620 from greglandrum)
  - Inconsistent metal disconnectors
 (github issue #3625 from pschwllr)
  - Ring stereochemistry not properly removed from N atoms
 (github issue #3631 from greglandrum)
  - moldraw2djs should not close all polygonal paths
 (github pull #3634 from greglandrum)
  - Unidentifiable C++ Exception with FMCS
 (github issue #3635 from proteneer)
  - Bump catch2 version to allow builds on Apple M1
 (github pull #3641 from naefl)
  - Segmentation fault when parsing InChI
 (github issue #3645 from AnPallo)
  - RDK_BUILD_THREADSAFE_SSS does not work as expected
 (github issue #3646 from pascal-soveaux)
  - Disabling MaeParser and CoordGen Support Breaks the Build
 (github issue #3648 from proteneer)
  - BondStereo info lost in FragmentOnBonds()
 (github pull #3649 from bp-kelley)
  - memory leak when sanitization fails in InChIToMol() 
 (github issue #3655 from greglandrum)
  - Qt GUI libraries being linked into rdmolops.so when Qt support is enabled
 (github issue #3658 from ricrogz)
  - Documentation of Chem.rdmolops.GetMolFrags's frag argument is wrong
 (github issue #3670 from noncomputable)
  - fmcs() + bogus input causes engine crash
 (github issue #3687 from robins)
  - qmol_from_ctab() with NULL crashes engine
 (github issue #3688 from robins)
  - qmol_from_smiles() + bogus input causes engine crash
 (github issue #3689 from robins)
  - Check PIL support for tostring and fromstring
 (github pull #3690 from sirbiscuit)
  - Move MolToSVG() to rdkit.Chem.Draw (Addresses #3694)
 (github pull #3696 from ricrogz)
  - Pandas AttributeError when rendering Molecule in DataFrame
 (github issue #3701 from enricogandini)
  - Memory leak in EnumerateLibrary
 (github issue #3702 from jose-mr)
  - Fix to add ZLIB_INCLUDE_DIRS for Windows build
 (github pull #3714 from kazuyaujihara)
  - Docs/Book: Unexpected unicode character makes pdflatex build fail
 (github issue #3738 from mbanck)
  - Test suite failures if eigen3 is not available
 (github issue #3740 from mbanck)
  - Regression in depiction of double bonds in aromatic rings
 (github issue #3744 from ptosco)
  - RGD with RGroupMatching.GA leaks memory and takes too long
 (github issue #3746 from ptosco)
  - Fix comment to match the code in RemoveHsParameters
 (github pull #3747 from jasondbiggs)
  - Inconsistent canonical tautomer on repeated application
 (github issue #3755 from darintay)
  - bonds no longer highlighted in substruct matches in jupyter
 (github issue #3762 from greglandrum)
  - SubstanceGroup output doesn't correctly quote " symbols
 (github issue #3768 from greglandrum)
  - MolToSmarts inverts direction of dative bond
 (github issue #3774 from IvanChernyshov)
  - Regression in dihedral constraints
 (github issue #3781 from ptosco)
  - Fix pillow error in IPythonConsole.py
 (github pull #3783 from skearnes)
  - lock swig version in MacOS CI builds
 (github pull #3789 from greglandrum)
  - DrawMorganBit errors when useSVG is False
 (github issue #3796 from ncfirth)
  - SubstructLibrary Cached Smiles Holders have bad behavior with bad smiles
 (github issue #3797 from bp-kelley)
  - MolFromSmiles('[He]') produces a diradical helium atom
 (github issue #3805 from jasondbiggs)
  - NaNs from AUTOCORR2D descriptor
 (github issue #3806 from greglandrum)
  - MaeMolSupplier throws an invariant exception on parsing an "undefined" chirality label
 (github issue #3815 from ricrogz)
  - Sanitize molecules when SMILES needs to be produced in PandasTools
 (github pull #3818 from mwojcikowski)
  - Tautomer Query copy constructor is shallow not deep causing segfaults in destructor
 (github issue #3821 from bp-kelley)
  - OptimizeMolecule and OptimizeMoleculeConfs Argument Bug
 (github issue #3824 from xiaoruiDong)
  - rdMolEnumerator.Enumerate() fragile w.r.t. atom ordering
 (github issue #3844 from greglandrum)
  - MinimalLib: Bonds are parallel in SVG but not on an HTML5 Canvas
 (github issue #3852 from dskatov)
  - AddHs creates H atom with nan coordinates on edge case 2D structure
 (github issue #3854 from ricrogz)
  - Build error with static boost libraries (v1.73)
 (github issue #3865 from nielskm)
  - Make sure that added R-groups have non-zero coordinates
 (github pull #3877 from ptosco)
  - Bad H coordinates on fused ring
 (github issue #3879 from greglandrum)
  - SubstructLibrary needs to check bond ring queries as well
 (github issue #3881 from bp-kelley)
  - Fixes Amine.Tertiary.Aromatic definition
 (github pull #3883 from bp-kelley)
  - inconsistency in seedSmarts in FMCS between and GetSubstructureMatches
 (github issue #3886 from proteneer)
  - PandasTools.RGroupDecomposition throws an error when redraw_sidechains is set to True.
 (github pull #3888 from greglandrum)
  - Dev/update glare to py3
 (github pull #3892 from bp-kelley)
  - ConfGen: Macrocycle torsion terms not being used with fused macrocycles
 (github pull #3894 from greglandrum)
  - Broken KNIME link in README
 (github issue #3897 from yurivict)
  - Change class to struct for forward declaration
 (github pull #3906 from bp-kelley)
  - Fixes issues with unlabelled groups on aromatic nitrogens
 (github pull #3908 from ptosco)
  - Fix #3659 regression introduced in #3832
 (github pull #3909 from ricrogz)
  - Error rendering SMARTS queries with atom OR lists
 (github issue #3912 from greglandrum)
  - MoDraw2D: Get tests working without freetype
 (github pull #3923 from greglandrum)
  - RGD default scoring function does not always work as expected
 (github issue #3924 from jones-gareth)
  - MolDraw2D: relative font size changes with bond lengths in molecule
 (github issue #3929 from greglandrum)
  - MolDraw2D: coordinates for reactions not being used
 (github issue #3931 from greglandrum)
  - Follow-on patch for changes in #3899
 (github issue #3932 from greglandrum)
  - Fix MolDraw2DQt exports
 (github pull #3935 from ricrogz)
  - Fix building JavaWrappers on Windows, dynamic linking
 (github pull #3936 from ricrogz)
  - Boost header warnings when compiling
 (github issue #3956 from jasondbiggs)
  - Adds removeAllHydrogenRGroupsAndLabels and fixes kekulization issues
 (github pull #3944 from ptosco)
  - MolToJSONData fails when mol has a property that can't be stringified
 (github issue #3956 from jasondbiggs)
  - RWMol should reset(), not release(), dp_delAtoms and dp_delBonds
 (github pull #3970 from greglandrum)


## New Features and Enhancements:
  - add context managers for writers
 (github issue #2217 from greglandrum)
  - MolToSmiles(kekuleSmiles=True) gives SMILES with aromatic bonds
 (github issue #2788 from adalke)
  - allow specification of color map when drawing similarity maps
 (github issue #2904 from greglandrum)
  - Clean up CMake files
 (github pull #3417 from e-kwsm)
  - Speed up GraphMol/Chirality.cpp/iterateCIPRanks
 (github pull #3482 from jinpan)
  - Removes function which is an exact duplicate of another function
 (github pull #3524 from ptosco)
  - A couple of minor improvements to FindCairo
 (github pull #3535 from ptosco)
  - Give a bit more time to RGD test in debug builds
 (github pull #3536 from ptosco)
  - A couple of fixes to the build system
 (github pull #3538 from ptosco)
  - Modularized WASM module
 (github issue #3561 from dskatov)
  - A couple changes to speedup bulk similarity calculations from Python
 (github pull #3574 from greglandrum)
  - add documentation for the JS wrappers
 (github pull #3583 from greglandrum)
  - add a "comic mode" to MolDraw2D
 (github pull #3584 from greglandrum)
  - Add rendering of SGroup brackets to MolDraw2D
 (github pull #3586 from greglandrum)
  - Update Install.md
 (github pull #3589 from slchan3)
  - Add explicit support for remaining CTAB query bond types
 (github issue #3599 from greglandrum)
  - update Cookbook stereochemistry examples
 (github pull #3604 from vfscalfani)
  - Add support for rendering SGroup data fields to MolDraw2D
 (github pull #3619 from greglandrum)
  - Support rendering the "ABS" flag in MolDraw2D
 (github issue #3623 from greglandrum)
  - Support drawing some query bonds
 (github pull #3624 from greglandrum)
  - Support rendering variable attachment points
 (github pull #3626 from greglandrum)
  - add configuration option to disable atom symbols in the rendering
 (github pull #3630 from greglandrum)
  - Render link nodes in MolDraw2D
 (github issue #3637 from greglandrum)
  - First pass at MolZip (now with bond stereo!)
 (github pull #3644 from bp-kelley)
  - Add molecule annotations/notes to MolDraw2D
 (github pull #3651 from greglandrum)
  - support setting MolDraw2DOptions using JSON from Python
 (github pull #3660 from greglandrum)
  - Make the scope control for Qt more idiomatic
 (github pull #3663 from d-b-w)
  - Expanded MolEnumerator functionality
 (github pull #3664 from greglandrum)
  - add support for generating pattern fps for MolBundles
 (github pull #3665 from greglandrum)
  - Add a callback function to EmbedParameters struct
 (github issue #3667 from jasondbiggs)
  - update SequenceParsers.cpp
 (github pull #3683 from magattaca)
  - MCS: extend completeRingsOnly to cover atoms as well
 (github issue #3693 from driesvr)
  - Add Molbundle search to SWIG
 (github pull #3698 from jones-gareth)
  - Added getMessage method to exceptions
 (github pull #3700 from sroughley)
  - add context manager for MolSuppliers
 (github issue #3703 from greglandrum)
  - Make better use of strictParsing for SGroups
 (github pull #3705 from ptosco)
  - Allow using  POPCNT on big-endian ppc64
 (github pull #3727 from pkubaj)
  - Cleanup: remove fromstring and tostring from functions working with pillow
 (github issue #3730 from greglandrum)
  - Set strictParsing to false in MinimalLib
 (github pull #3737 from ptosco)
  - 3D MCS - Minimal version, no refactoring
 (github pull #3749 from robbason)
  - Include Winsock2.h instead of Windows.h in DebugTrace.h
 (github pull #3756 from dpaoliello)
  - R group match any issue
 (github pull #3767 from jones-gareth)
  - Support new coordgen option to not always make bonds to metals zero-order
 (github pull #3769 from greglandrum)
  - DistanceGeometry: add flag to enforce trans amides
 (github pull #3794 from greglandrum)
  - MolDraw2D: first pass at rendering atom lists
 (github pull #3804 from greglandrum)
  - Issue a warning when embedding a molecule with no Hs
 (github pull #3807 from greglandrum)
  - Add tautomer query to the substructlibrary
 (github pull #3808 from bp-kelley)
  - Enhanced generateDepictionMatching2DStructure functionality
 (github pull #3811 from ptosco)
  - Confgen: add option to use symmetry when doing RMS pruning 
 (github pull #3813 from greglandrum)
  - Remove boost::foreach from public headers
 (github pull #3820 from ricrogz)
  - Adds isotopeLabels and dummyIsotopeLabels MolDrawOptions
 (github pull #3825 from ptosco)
  - Added 2 Cookbook examples
 (github pull #3831 from vfscalfani)
  - Separate MolDraw2DQt into its own library
 (github pull #3832 from d-b-w)
  - Facilities for interactive modification of molecule drawing
 (github pull #3833 from SPKorhonen)
  - cleanup a bunch of compiler warnings
 (github pull #3849 from greglandrum)
  - add a new mol draw option to draw wedge bonds with a single color 
 (github pull #3860 from jasondbiggs)
  - Add Kier Phi descriptor
 (github pull #3864 from greglandrum)
  - Add basic support for hydrogen bonds
 (github pull #3871 from greglandrum)
  - Allow batch editing of molecules: removal only
 (github pull #3875 from greglandrum)
  - Allow retrieving the _ErGAtomTypes property from Python
 (github pull #3878 from ptosco)
  - Exposes InsertMol to python RWMol
 (github pull #3907 from bp-kelley)
  - Use https for Avalon and Inchi downloads
 (github pull #3915 from ptosco)
  - support empty/missing SDT lines for SGroup data
 (github pull #3916 from greglandrum)
  - Cookbook entries should be updated 
 (github issue #3917 from greglandrum)
  - MolDraw2D: support changing annotation colours
 (github pull #3919 from greglandrum)
  - include context managers for the multithreaded suppliers too
 (github pull #3920 from greglandrum)
  - Documentation cleanup and update
 (github pull #3922 from greglandrum)
  - remove an MSVC++ warning caused by #3849
 (github pull #3927 from greglandrum)
  - Adds removeAllHydrogenRGroupsAndLabels and fixes kekulization issues
 (github pull #3944 from ptosco)
  - Remove temporary labels from RGD results
 (github pull #3947 from ptosco)
  - appended a new project depend on RDKit
 (github pull #3955 from kotori-y)
  - Do not add unnecessary R-labels (and an optimization)
 (github pull #3969 from ptosco)
  - Add return codes and make RGroupDecomp less verbose 
 (github pull #3971 from bp-kelley)
  - update to coordgen 2.0.0
 (github pull #3974 from greglandrum)


## Deprecated code (to be removed in a future release):
- The "minimizeOnly" option for coordgen will be removed in the next RDKit release

# Release_2020.09.1
(Changes relative to Release_2020.03.1)


## Backwards incompatible changes
- We've added additional allowed valences for Cl (now 1, 3, 5), Br (now 1, 3,
  5), I (now 1, 3, 5), At (now 1, 3, 5), Xe (now 0, 2, 4, 6), and Po (now 2, 4,
  6). Molecules with atoms in the new valence states will no longer generate
  sanitization errors. Note that this has an impact on the chemistry of
  molecules containing 3-valent I and at least one implict H (present 24 times
  in ChEMBL 27): previously this was incorrectly assigned two implicit Hs, now
  it has no implicit Hs. 
- Aromaticity perception of molecules like `Cc1nnc2n1c1ccccc1n1c(C)nnc12` now
  correctly recognizes the full outer envelope, i.e. the bonds joining the rings
  are now also aromatic.
- FindMCS() may return single atom MCSs, whereas previously it returned an empty
  MCS unless there was at least one commond bond across the input structures.
  So the MCS between molecules `CC` and `CO` is now `[#6]` rather than being null.
- The fontSize()/setFontSize() (FontSize()/SetFontSize()) methods in MolDraw2D
  now work in units of pixels (more or less) instead of the molecule units.
- The Open3DAlign functionality is now in its own separate library - `O3AAlign`
  in cmake. If you are working in C++ and using O3A functionality, you'll need
  to link against this library as well now.
- Due to improvements in the tautomer enumeration code, the method
  `TautomerEnumerator::enumerate` now returns a `TautomerEnumeratorResult`
  object instead of a vector of molecules. Note that if you are iterating over
  the results of a call to `enumerate()` you shouldn't need to change your code.
  If you want to invoke the old (and deprecated, see below) form from C++, call
  `TautomerNumerator::enumerate(mol, nullptr)` or explicitly pass a
  `boost::dynamic_bitset*` to capture the modified atoms.
- The default precision setting for coordgen has been changed. The new default
  was selected to greatly reduce the number of molecules for which it takes a
  very long time to generate coordinates while still producing nice looking
  structures. We may continue to tweak this default value if/when problems
  with it are reported. If you would like to go back to the previous setting, set 
  CoordgenParams.minimizerPrecision to CoordgenParams.sketcherStandardPrecision 
  when you invoke rdCoordGen.AddCoords()
- Uncharger::uncharge() will now neutralize `[Cl,Br,I][O-], [Cl,Br,I](=O)[O-],
  [Cl,Br,I](=O)(=O)[O-], [Cl,Br,I](=O)(=O)(=O)[O-], [O-]N=N[O-], [N,P](=O)[O-],
  [N+](=O)([O-])[O-], P(=O)([O-])[O-], P(=O)([O-])([O-])[O-], S([O-])[O-],
  S(=O)([O-])[O-], S(=O)(=O)([O-])[O-], S(=O)(=O)([O-])OOS(=O)(=O)[O-]`.
  Previously not all of these inorganic acid counterions were consistently
  neutralized.
- The `None` value in the `RGroupCoreAlignment` enum was renamed to `NoAlignment`
  in both C++ and Python, in order to avoid issues when attempting to access it
  from Python.

## Highlights
- There's been another big improvement in the quality of molecule drawings:
  character and font handling is greatly improved thanks to the use of the
  FreeType library
- A new feature has been added to efficiently allow tautomer-insensitive
  substructure search.
- A new, much more accurate, algorithm is available for calculating CIP labels
  on atoms and bonds.
- There's a new rdDeprotect module to allow automatically deprotecting molecules
  before putting them into reactions
- Molecule and reaction metadata can now be added to PNG files generated by
  MolDraw2DCairo

## Acknowledgements
Shrey Aryan, Jinserk Baik, Francois Berenger, Cédric Bouysset, David Cosgrove,
Ivan Chernyshov, Guillaume Godin, Manan Goel, Jan H. Jensen, Gareth Jones, Maria
Kadukova, Eisuke Kawashima, Steven Kearnes, Brian Kelley, Joos Kiener, Kenneth
Lum, Joshua Meyers, Rocco Moretti, Paul R Moses, Dan Nealschneider, Jin Pan,
Joann Prescott-Roy, Matthew Robinson, Jaime Rodríguez-Guerra, Ricardo
Rodriguez-Schmidt, Jeff van Santen, Roger Sayle Vincent F. Scalfani Eric Taw,
Ansgar Schuffenhauer, Paolo Tosco, Ivan Tubert-Brohman, Riccardo Vianello,
Rachel Walker, Maciej Wójcikowski, Christopher Zou, daverona, hjuinj,
intrigus-lgtm, autodataming, paconius, sailfish009

## Bug Fixes:
  - Python tests fail when RDK_BUILD_COMPRESSED_SUPPLIERS is enabled
 (github issue #1888 from greglandrum)
  - ResonanceMolSupplier potentially stuck in infinite loop
 (github issue #2597 from tawe141)
  - ctest pythonTestDirChem failed
 (github issue #2757 from jinserk)
  - Issue with inversion/retention of stereochemistry
 (github issue #2891 from mc-robinson)
  - cannot parse reaction SMILES/SMARTS with dative bonds
 (github issue #2954 from greglandrum)
  - ResonanceMolSupplier can fail with small maxStructs values
 (github issue #3041 from greglandrum)
  - seg fault in ResonanceMolSupplier()
 (github issue #3048 from greglandrum)
  - Bug in image rendering of dative bonds
 (github issue #3056 from IvanChernyshov)
  - Coordinates from coordgen are not centered around the origin
 (github pull #3058 from DavidACosgrove)
  - fix a typo in ScaffoldNetwork/CMakeLists.txt
 (github pull #3060 from greglandrum)
  - Bad double bond placement in polycyclic aromatics
 (github issue #3061 from DavidACosgrove)
  - SGroups with more than one attachment point are now properly parsed
 (github pull #3072 from greglandrum)
  - Reactions should not use strict implicit valence calculations
 (github issue #3097 from mwojcikowski)
  - partial reacting atom detection
 (github issue #3119 from thegodone)
  - DrawMolecules does not center molecules
 (github issue #3126 from JoshuaMeyers)
  - results from coordgen are sometimes not centered
 (github issue #3131 from greglandrum)
  - GCC 10.0.1 compile error
 (github issue #3135 from rvianello)
  - Memory leak when parsing bad SMILES
 (github issue #3139 from intrigus-lgtm)
  - Error breaking StereoBonds in reactions
 (github issue #3147 from mc-robinson)
  - MolOps::removeHs() removes hydrides
 (github issue #3150 from jhjensen2)
  - Kekulization error from CreateScaffoldNetwork
 (github issue #3153 from greglandrum)
  - Fix drawing of N plus
 (github pull #3165 from DavidACosgrove)
  - RWMol::clear() does not explicitly clean up SubstanceGroups or StereoGroups
 (github issue #3167 from greglandrum)
  - Modifying a molecule should not automatically clear SubstanceGroups
 (github issue #3168 from greglandrum)
  - removeHs() should not remove atoms in SubstanceGroups
 (github issue #3169 from greglandrum)
  - fix a memory problem detected in malformed SMILES
 (github pull #3171 from greglandrum)
  - Python wrapper: SetQuery and ExpandQuery for bonds
 (github pull #3172 from i-tub)
  - S-groups: PARENT field should reference index
 (github issue #3175 from greglandrum)
  - rdScaffoldNetwork causes segmenation fault upon None molecule
 (github issue #3177 from AnsgarSchuffenhauer)
  - fix a small inconsistency in the name of the inchi package
 (github pull #3182 from rvianello)
  - Molecule constructed from CXSMILES cannot be translated to SMARTS
 (github issue #3197 from greglandrum)
  - Formatting fix of CalcRMS
 (github pull #3203 from chmnk)
  - fix the CompressedSDMolSupplier python iterator interface
 (github pull #3204 from rvianello)
  - Queries generated from PreprocessReaction cannot be translated to SMARTS
 (github issue #3206 from greglandrum)
  - Attachment point info not being read from V2000 mol blocks
 (github issue #3207 from greglandrum)
  - Memory Sanitizer fails on molFromPickle on empty file
 (github issue #3211 from intrigus-lgtm)
  - Throw exception when reading from stream fails.
 (github pull #3212 from intrigus-lgtm)
  - fix molstogridimage on certain fragments/smarts patterns
 (github pull #3217 from bp-kelley)
  - Lines in wedge bonds being drawn too closely together
 (github issue #3226 from paconius)
  - EnumerateStereochemistry should clear CIP labels
 (github issue #3231 from greglandrum)
  - lock CI cairo version to force an install from the rdkit repo
 (github pull #3240 from greglandrum)
  - XBCORR and XBHEAD in Sgroups no longer cause parse failures
 (github pull #3242 from greglandrum)
  - LINKNODEs are ignored by the CTAB parsers
 (github pull #3247 from greglandrum)
  - add GetStringVectProp() to SubstanceGroup class
 (github pull #3251 from greglandrum)
  - Envelope aromaticity not detected in complex fused system
 (github issue #3256 from greglandrum)
  - Draw.MolsToGridImage repeating atom indices
 (github issue #3258 from greglandrum)
  - Atom indices clash with atom symbols in small pictures.
 (github issue #3262 from DavidACosgrove)
  - MinimalLib Dockerfile is broken at HEAD
 (github issue #3267 from skearnes)
  - Fixes #2757
 (github pull #3268 from greglandrum)
  - RGroupDecomposition restructuring
 (github pull #3270 from bp-kelley)
  - Get PPC builds working
 (github pull #3285 from greglandrum)
  - ScaffoldNetwork not in C# wrappers
 (github pull #3289 from jones-gareth)
  - bonds with "either' stereo cannot be read from JSON
 (github pull #3290 from greglandrum)
  - Small bug fixes and cleanups from fuzz testing
 (github pull #3299 from greglandrum)
  - DrawOptions: bondLineWidth behaving differently since 2020 versions
 (github issue #3305 from kienerj)
  - Not possible to copy SubstanceGroups in Python
 (github issue #3312 from greglandrum)
  - Stereochemistry perception getting confused by a bad drawing.
 (github issue #3314 from greglandrum)
  - SubstanceGroups should not be written with quotes around missing fields
 (github issue #3315 from greglandrum)
  - SetDoubleBondNeighborDirections() not overwriting existing bond directions
 (github issue #3322 from greglandrum)
  - AdjustQueryParameters.adjustSingleBondsBetweenAromaticAtoms does not modify ring bonds
 (github issue #3325 from greglandrum)
  - Fixes for aromatic bond fuzzy queries
 (github pull #3328 from jones-gareth)
  - lock sphinx version in CI due to problem with v3.2.0
 (github pull #3332 from greglandrum)
  - Remove deprecated Sphinx options
 (github pull #3335 from greglandrum)
  - more bug fixes and cleanups from fuzz testing
 (github pull #3339 from greglandrum)
  - unspecified branch bonds in SMARTS don't have aromaticity set
 (github issue #3342 from greglandrum)
  - Incorrect resonance structures in presence of dative bonds
 (github issue #3349 from IvanChernyshov)
  - Converting atoms with high radical counts to InChI generates incorrect results
 (github issue #3365 from greglandrum)
  - Replace fill-opacity= with fill-opacity: in MolDraw2DSVG and tests
 (github pull #3368 from lummyk)
  - Fixes a bug in AddHs() involving sp2 centers with degree 1
 (github pull #3383 from ptosco)
  - Information about charges and isotopes lost when calling AdjustQueryProperties
 (github issue #3388 from greglandrum)
  - prepareMolForDrawing() incorrectly adds chiral Hs if no ring info is present
 (github issue #3392 from greglandrum)
  - CXSMILES parser should not set atom maps for attachment points
 (github issue #3393 from greglandrum)
  - Fixes a couple of query-related bugs
 (github pull #3398 from ptosco)
  - Doing a match of a recursive smarts leaves traces of the previous match
 (github issue #3403 from bp-kelley)
  - Recursive smarts cannot be used in the core for rgroup decomposition
 (github pull #3404 from bp-kelley)
  - Improvements to reaction chirality handling
 (github pull #3412 from greglandrum)
  - V3K mol blocks with no atoms fail to parse
 (github issue #3413 from greglandrum)
  - Problem parsing SGroup data comtaining `""`
 (github issue #3415 from greglandrum)
  - MolEnumerator::enumerate() should call updatePropertyCache()
 (github pull #3420 from greglandrum)
 - Fixed bad draw scale in drawMolecules. Github3391.  Take 2.
 (github pull #3424 from DavidACosgrove)
  - Replace fill-opacity= to fill-opacity: in reaction.out
 (github pull #3426 from daverona)
  - set the ChiralityPossible tag when using the new code with FindMolChiralCenters
 (github pull #3434 from greglandrum)
  - Silence deprecation warning
 (github pull #3439 from ptosco)
  - update minimallib python requirements to python3
 (github pull #3449 from greglandrum)
  - Fix dead links to inchi-trust
 (github pull #3451 from jvansan)
  - ringMatchesRingOnly=True produces a SMARTS query that return no substructure matches
 (github issue #3458 from jaimergp)
  - Normalization rule incorrectly matches sulfones
 (github issue #3460 from greglandrum)
  - BlockLogs was reenabling all logs, not just the ones that were disabled
 (github pull #3466 from bp-kelley)
  - Hydrogen is incorrectly identified as an "early" atom
 (github issue #3470 from greglandrum)
  - Fixes typo that causes the build to fail
 (github pull #3477 from ptosco)
  - Fix a crashing bug with None in rdMolStandardize
 (github pull #3481 from greglandrum)
  - zlib.h not found if not in system directories
 (github issue #3493 from ricrogz)
  - fix paths in ConformerParser tests
 (github pull #3504 from ricrogz)

## New Features and Enhancements:
  - Add GetBestRMS function
 (github issue #1820 from chmnk)
  - Add reorder tautomers function and accompanying tests
 (github pull #3043 from chriswzou)
  - Set RDK_BOOST_VERSION to pass minimum required version to FindBoost
 (github pull #3074 from e-kwsm)
  - bug: the MCS of the molecules CH4 and CH3OH is empty. how to return C? 
 (github issue #3095 from autodataming)
  - start using boost:stacktrace
 (github pull #3124 from greglandrum)
  - Add Fuzzing, fixes #2857
 (github pull #3128 from intrigus-lgtm)
  - Cookbook entry for ETKDG with rings
 (github pull #3129 from hjuinj)
  - Fixes #2795
 (github pull #3134 from manangoel99)
  - Bump Catch2 to v2.12.1
 (github pull #3136 from e-kwsm)
  - Modernize how legacy C headers are included
 (github pull #3137 from e-kwsm)
  - Avoid C preprocessor macros
 (github pull #3138 from e-kwsm)
  - Modernization: use nullptr
 (github pull #3143 from e-kwsm)
  - Update fuzzer dict
 (github pull #3162 from intrigus-lgtm)
  - Add BCUT2D and AUTOCORR2D to desclist
 (github pull #3178 from bp-kelley)
  - Remove usage of the deprecated random_shuffle() function
 (github pull #3187 from greglandrum)
  - clang-tidy modernize-use-default-member-init and modernize-use-emplace
 (github pull #3190 from greglandrum)
  - Tautomer search
 (github pull #3205 from jones-gareth)
  - Add optional timeout to RGroupDecomposition
 (github pull #3223 from greglandrum)
  - Allow symmetrization to be completely disabled in RGD code
 (github issue #3224 from greglandrum)
  - gitignore source and build files from the RingFamilies external lib
 (github pull #3228 from d-b-w)
  - Add new CIP labelling algorithm
 (github pull #3234 from ricrogz)
  - Adds more options to adjustQueryProperties
 (github pull #3235 from greglandrum)
  - Improve SSSR performance for large molecules
 (github pull #3236 from d-b-w)
  - Support using FreeType for text rendering
 (github pull #3237 from DavidACosgrove)
  - Cleanup warnings from clang-10
 (github pull #3238 from greglandrum)
  - DEB packaging: cairo support is needed to generate PNGs
 (github pull #3250 from UnixJunkie)
  - Added call to test legends.
 (github pull #3252 from DavidACosgrove)
  - Improve performance of aromaticity detection for large molecules
 (github pull #3253 from d-b-w)
  - Speed up ring finding by skipping nodes not in rings
 (github pull #3254 from d-b-w)
  - Support enumerating some mol file features into `MolBundles`
 (github pull #3257 from greglandrum)
  - Add cxsmiles query atoms to CTAB parsers and writers
 (github pull #3261 from greglandrum)
  - Update to Coordgen v1.4.1
 (github pull #3265 from ricrogz)
  - ScaffoldNetwork: add feature to count the number of molecules a scaffold originates from
 (github pull #3275 from greglandrum)
  - rgroup speedup
 (github pull #3279 from bp-kelley)
  - Stop trying to assign hybridization to actinides
 (github pull #3281 from greglandrum)
  - Decouple coordgen and maeparser integrations
 (github pull #3286 from greglandrum)
  - Avoid really slow Windows conda builds
 (github pull #3287 from ptosco)
  - Embed default truetype font
 (github pull #3288 from greglandrum)
  - Expanded support for CXSMILES features
 (github pull #3292 from greglandrum)
  - Deprotection Library
 (github pull #3294 from bp-kelley)
  - Use operator() and __call__() consistently across RDKit
 (github pull #3295 from ptosco)
  - Molecule metadata in PNGs
 (github pull #3316 from greglandrum)
  - Cleanup alignment dependencies
 (github pull #3317 from greglandrum)
  - Add the option to minimize structures with coordgen
 (github pull #3319 from greglandrum)
  - Updated code for chirality perception
 (github pull #3324 from greglandrum)
  - Some work on TautomerEnumerator
 (github pull #3327 from ptosco)
  - Add fragmentOnBonds to SWIG wrappers
 (github issue #3329 from greglandrum)
  - Sped up SSSR by not storing every path back to root
 (github pull #3333 from rachelnwalker)
  - Fix Cookbook formatting and added 4 new examples
 (github pull #3345 from vfscalfani)
  - switch to using target_compile_definitions instead of add_definitions
 (github pull #3350 from greglandrum)
  - [GSoC-2020] Generalized and Multithreaded File Reader
 (github pull #3363 from shrey183)
  - support new CIP code and StereoGroups in MolDraw2D_detail::addStereoAnnotation()
 (github issue #3369 from greglandrum)
  - expose additional SubstanceGroup data members to Python
 (github pull #3375 from greglandrum)
  - Add MolDraw2DJS
 (github pull #3376 from greglandrum)
  - Add APK package link for Alpine Linux distribution
 (github pull #3379 from daverona)
  - Add SubstanceGroups to the SWIG Wrappers
 (github pull #3390 from jones-gareth)
  - Add better support for isotopic Hs to removeHs() and addHs()
 (github pull #3396 from ptosco)
  - Add support for abbreviations
 (github pull #3406 from greglandrum)
  - Allow passing explicit removeHs, sanitize and strict flags to the MDL rxn parser
 (github pull #3411 from ptosco)
  - Improvements to reaction chirality handling
 (github pull #3412 from greglandrum)
  - RGD cleanup, optimization and a better fix for #1705
 (github pull #3428 from ptosco)
  - Tautomers with endocyclic double bonds should be preferred over exocyclic ones
 (github issue #3430 from ptosco)
  - RGD: Code modernization and an optimization
 (github pull #3437 from ptosco)
  - expose PNG metadata functions to python
 (github pull #3440 from greglandrum)
  - Replace basestring
 (github pull #3441 from iammosespaulr)
  - Get the Uncharger to deal with a larger set of acids correctly
 (github pull #3448 from ptosco)
  - expose templated coordinate generation to the JS Wrapper
 (github pull #3450 from greglandrum)
  - change default precision for coordgen
 (github pull #3452 from greglandrum)
  - add coordgen support to demo.html
 (github pull #3453 from greglandrum)
  - Two simple MolStandardizer code cleanups
 (github pull #3454 from ptosco)
  - A few improvements to MolStandardize::Normalizer
 (github pull #3455 from ptosco)
  - Add Cookbook entries 30-32
 (github pull #3459 from vfscalfani)
  - A few small tweaks to the drawing code
 (github pull #3464 from greglandrum)
  - Make MetalDisconnector more robust against metallorganics
 (github pull #3465 from greglandrum)
  - Add nocharge algorithm example to cookbook
 (github pull #3467 from vfscalfani)
  - ROMol: add inline impl for common getNumAtoms call
 (github pull #3469 from jinpan)
  - Improve sphinx formatting in rdSubstructLibrary
 (github issue #3471 from cbouy)
  - Cmake config improvements
 (github pull #3478 from rvianello)
  - allow fillColour to be changed from python
 (github pull #3480 from greglandrum)
  - Fix undefined behavior in testCoordGen test
 (github pull #3495 from roccomoretti)
  - Add a version for the pattern fingerprint
 (github pull #3496 from greglandrum)
  - Fixes a number of issues flagged by clang
 (github pull #3498 from ptosco)
  - Update to maeparser v1.2.4
 (github pull #3506 from sailfish009)
  - Fix python invalid escape sequences
 (github pull #3508 from ricrogz)

## Code removed in this release:
- To improve API consistency of the exceptions in RDKit with the default ones in
  the STL, the several `message()` methods and `Invariant::getMessage()` in RDKit's
  exceptions have been removed in favor of `what()`. 
- The old MolHash code has been removed from the C++ code, all wrappers, and the
  PostgreSQL cartridge.

## Deprecated code (to be removed in a future release):
- The function `FileParserUtils::replaceAtomWithQueryAtom()` has been moved to
  the namespace QueryOps. Please use `QueryOps::replaceAtomWithQueryAtom()`
  instead. The version in the `FileParserUtils` namespace will be removed in the
  next release.
- The method `std::vector<ROMOL_SPTR> TautomerEnumerator::enumerate(const ROMol &mol, boost::dynamic_bitset<> *modifiedAtoms, boost::dynamic_bitset<> *modifiedBonds = nullptr)` 
  is deprecated and will be removed in a future release. 
  Please use `TautomerEnumeratorResult TautomerEnumerator::enumerate(const ROMol &mol,bool reassignStereo = true)` 
  instead.
- The `MolDraw2DQt` class is no longer supported since we don't think anyone is
  using it. It will be removed in the 2021.03 release unless we learn otherwise.



# Release_2020.03.1
(Changes relative to Release_2019.09.1)

## Backwards incompatible changes
- Searches for equal molecules (i.e. `mol1 @= mol2`) in the PostgreSQL cartridge
  now use the `do_chiral_sss` option. So if `do_chiral_sss` is false (the
  default), the molecules `CC(F)Cl` and `C[C@H](F)Cl` will be considered to be equal.
  Previously these molecules were always considered to be different.
- Attempting to create a MolSupplier from a filename pointing to an empty file,
  a file that does not exist or sometihing that is not a standard file (i.e.
  something like a directory) now generates an exception.
- The cmake option `RDK_OPTIMIZE_NATIVE` has been renamed to `RDK_OPTIMIZE_POPCNT`

## Highlights:
- The drawings generated by the MolDraw2D objects are now significantly improved
  and can include simple atom and bond annotations (#2931 and #3010)
- An initial implementation of a modified scaffold network algorithm is now
  available (#2911)
- A few new descriptor/fingerprint types are available - BCUTs (#2957), Morse
  atom fingerprints (#1773), Coulomb matrices (#2993), and MHFP and SECFP
  fingerprints (#2643)
- There is a new, and greatly improved, version of the RDKit Cookbook (#2884)
- There is a new version (v3) of the ETKDG conformer generator along with
  optional new terms for handling small rings and macrocyles (http://doi.org/dqnh) (#2999)


## Acknowledgements:
Marcel Baltruschat, Jason Biggs, Eliane Briand, Ben Cornett, David Cosgrove,
Andrew Dalke, Tim Dudgeon, Zhenting Gao, Guillaume Godin, Manan Goel, Gareth
Jones, Zachary Kaplan, Eisuke Kawashima, Steven Kearnes, Brian Kelley, Maxim
Koltsov, Franziska Kruger, Mieszko Manijak, Dan Nealschneider, Daniil
Polykovskiy, Daniel Probst, Sereina Riniker, Matthew Robinson, Steve Roughley,
Kevin Ryan, Vincent F. Scalfani, Ricardo Rodriguez Schmidt, Rim Shayakhmetov,
Aryan Shrey, Nik Stiefl, Matt Swain, Paolo Tosco, Wiep van der Toorn, Riccardo
Vianello, Shuzhe Wang, Piotr Wawrzyniak, Hsiao Yi, 'jasad1', 'luancarvalhomartins'


## Bug Fixes:
  - Mol rendering within DataFrames in a Jupyter Notebook is broken with Pandas 0.25.1
 (github issue #2673 from mrcblt)
  - Removed RDKIT_SIMDIVPICKERS_EXPORT
 (github pull #2740 from ptosco)
  - - enable building RDKitRingDecomposerLib.dll under Windows
 (github pull #2742 from ptosco)
  - Do a windows DLL build as part of the Azure DevOps setup
 (github pull #2743 from greglandrum)
  - Fix data race in sascorer.py
 (github pull #2744 from skearnes)
  - Uncharger not properly setting explicit/implicit H count
 (github issue #2749 from greglandrum)
  - MSVC compile error: MolHash scoped enum cannot be redeclared as unscoped
 (github issue #2752 from mcs07)
  - Molecules whose Y size is very small won't display as SVG
 (github issue #2762 from ptosco)
  - Make the cartridge tests work with PostgreSQL 12
 (github pull #2767 from greglandrum)
  - Salt stripper should consider bond matches as well as atom matches
 (github pull #2768 from greglandrum)
  - Bismuth should count as an early element
 (github issue #2775 from greglandrum)
  - addHs() fails on atoms with "bad" valences
 (github issue #2782 from greglandrum)
  - Element symbol lookup for some transuranics returns incorrect results
 (github issue #2784 from LeanAndMean)
  - [cartridge] molecular equality should use do_chiral_sss setting
 (github issue #2790 from greglandrum)
  - uncharger removes Hs from carbocations instead of adding them
 (github issue #2792 from greglandrum)
  - Fix build without boost serialization library
 (github pull #2796 from maksbotan)
  - Using `SetBoundsMat` significantly slows down conformer generation process.
 (github issue #2800 from hjuinj)
  - rdkit.Ched.rdFMCS.FindMCS generates invalid smarts
 (github issue #2801 from luancarvalhomartins)
  - Remove confId from *FFOptimizeMoleculeConfs Python docs
 (github issue #2805 from ptosco)
  - Hybridization queries on dummy atoms not written properly to SMARTS
 (github issue #2814 from greglandrum)
  - Charge range queries not properly written to SMARTS
 (github issue #2815 from greglandrum)
  - RDKit segfaults in MMFFOptimizeMoleculeConfs()
 (github issue #2820 from ptosco)
  - Trusted Smiles holder doesn't handle ring queries
 (github issue #2830 from bp-kelley)
  - Fix windows substructure crash
 (github pull #2836 from greglandrum)
  - Fix YAeHMOP build
 (github pull #2838 from ptosco)
  - testGithub2245 in testPickers.cpp occasionally fails
 (github issue #2839 from ptosco)
  - add define for RDK_USE_BOOST_SERIALIZATION
 (github pull #2859 from greglandrum)
  - fix start/end atoms when wedging bonds
 (github pull #2861 from greglandrum)
  - Fixes the size of the reduced charge matrix from eHT calculations
 (github pull #2864 from greglandrum)
  - Dev/pvs studio cleanups2
 (github pull #2877 from greglandrum)
  - segfault in MaeMolSupplier
 (github issue #2881 from greglandrum)
  - update maven url in build system
 (github pull #2889 from greglandrum)
  - EnumerateStereoisomers cannot handle STEREOANY bonds
 (github issue #2890 from ricrogz)
  - Update one of the cartridge tests that got missed
 (github pull #2894 from greglandrum)
  - acepentalene aromaticity perception
 (github issue #2895 from adalke)
  - New Similarity Maps drawing code Java Wrappers non-functional
 (github issue #2896 from sroughley)
  - Fix to allow multistructure images in Java/C# and use MCS for c# wrapper
 (github pull #2898 from jones-gareth)
  - Remove bogus URFLib library
 (github pull #2900 from greglandrum)
  - java wrapper build cleanups
 (github pull #2901 from greglandrum)
  - SMARTS parser fails on high-numbered ring closures in branches
 (github issue #2909 from greglandrum)
  - patch to make PandasTools tests pass with pandas v0.22
 (github pull #2913 from greglandrum)
  - fix doctest problem with Pandas v1.0
 (github pull #2918 from greglandrum)
  - Build with -D RDK_BUILD_COORDGEN_SUPPORT=OFF includes a test case that depends on MaeMolSupplier
 (github issue #2929 from rvianello)
  - MinimalLib: get_stereo_tags() should also return unspecified centers
 (github issue #2936 from greglandrum)
  - Fix regression introduced by e245349c
 (github pull #2945 from cornett)
  - Avoid data race warning in SmilesParse.cpp
 (github pull #2946 from skearnes)
  - Empty molecule has non-zero LabuteASA
 (github issue #2948 from jasondbiggs)
  - Fix a problem with aromatic heteroatom tautomer enumeration
 (github pull #2952 from greglandrum)
  - Molecule properties not retained with MolStandardize.rdMolStandardize.Cleanup()
 (github issue #2965 from ZacharyKaplan)
  - Fix build without boost serialization.
 (github pull #2972 from ricrogz)
  - RDKFuncs.chargeParent() core dumps when standardization is skipped
 (bithub issue #2970 from tdudgeon)
  - fix a typo in the scaffold network wrappers and add some tests
 (github pull #2982 from greglandrum)
  - Tautomer enumeration should remove stereo in all tautomers 
 (github issue #2990 from greglandrum)
  - Segmentation fault on EmbedMolecule
 (github issue #3019 from shayakhmetov)
  - Removed dllexport from a function that lives in the anonymous namespace
 (github pull #3027 from ptosco)


## New Features and Enhancements:
  - Morse atom fingerprint
 (github pull #1773 from thegodone)
  - Allow serializing coordinates as doubles
 (github issue #2510 from danpol)
  - Rework MaeMolSupplier, fix #2617
 (github pull #2620 from ricrogz)
  - Implementation of MHFP and SECFP Fingerprints
 (github pull #2643 from daenuprobst)
  - MatchFusedRings does not imply CompleteRingsOnly anymore
 (github pull #2748 from ptosco)
  - Improvements to JS wrappers
 (github pull #2751 from greglandrum)
  - Fix installed header directory structure
 (github pull #2754 from ricrogz)
  - Add doRandom to the header docs
 (github pull #2756 from bp-kelley)
  - Add queryMol data member to MCSResult
 (github pull #2759 from ptosco)
  - Add functions to enable/disable the substructure matching monkey patching in IPythonConsole.py
 (github issue #2786 from greglandrum)
  - Add a function to assign chiral tags from sss atom parity
 (github issue #2823 from ptosco)
  - Support MRV_IMPLICIT_H S groups when reading Mol blocks
 (github issue #2829 from greglandrum)
  - Unset executable flag
 (github pull #2833 from e-kwsm)
  - Remove O(N) behavior of getNumBonds
 (github pull #2847 from bp-kelley)
  - Feature proposal: add remove_stereochemistry=False flag for RemoveHs()
 (github issue #2848 from shayakhmetov)
  - Expose SubstructLibrary serialization stream
 (github pull #2853 from bp-kelley)
  - Fix typo
 (github pull #2862 from e-kwsm)
  - Rename RDK_OPTIMIZE_NATIVE to RDK_OPTIMIZE_POPCNT
 (github pull #2865 from ElianeBriand)
  - Update Draw.MolToImage() and Draw.MolToFile() to use the new drawing code
 (github pull #2866 from greglandrum)
  - Improve PostgreSQL cartridge install documentation
 (github pull #2870 from yellowBirdy)
  - Fixes #2858
 (github pull #2871 from greglandrum)
  - Add a cartridge test to the azure devops config
 (github pull #2873 from greglandrum)
  - Add a new Cookbook v2 to the RDKit docs
 (github pull #2884 from vfscalfani)
  - Add MolVS tautomer canonicalization
 (github pull #2886 from greglandrum)
  - add a convenience function for RGD--Pandas integration
 (github pull #2887 from greglandrum)
  - run clang-tidy with readability-braces-around-statements
 (github pull #2899 from greglandrum)
  - Allow RDProps::clearProp to succeed even if the prop doesn't exist
 (github issue #2910 from greglandrum)
  - Add a scaffold network implementation
 (github pull #2911 from greglandrum)
  - cleanup of the SMILES/SMARTS parsing and writing code
 (github pull #2912 from greglandrum)
  - Add _ctab, _mol2, _pdb to allow direct mol construction from strings
 (github issue #2916 from greglandrum)
  - Parse and handle the stereoCare or STBOX flags in CTABs
 (github pull #2917 from greglandrum)
  - RDKit exceptions do not override the default `what()` method
 (github issue #2920 from ricrogz)
  - Allow custom post-match filters for substructure matching
 (github pull #2927 from greglandrum)
  - Proposed improvements to 2D Drawing Code
 (github issue #2931 from DavidACosgrove)
  - Include running the documentation tests as part of the CI runs
 (github pull #2932 from greglandrum)
  - Add support for phosphine and arsine chirality
 (github issue #2949 from wopozka)
  - A couple additions to the extended Hueckel integration
 (github pull #2955 from greglandrum)
  - Add BCUT 2D descriptors
 (github pull #2957 from bp-kelley)
  - Add multithreaded pattern/fp generator
 (github pull #2973 from bp-kelley)
  - Description for the data files.
 (github pull #2975 from zhentg)
  - Enable larger ring matches in SMARTS expressions
 (github pull #2981 from d-b-w)
  - ScaffoldNetwork rearrangements
 (github pull #2985 from greglandrum)
  - add add_hs() and remove_hs() to JS wrappers
 (github pull #2986 from greglandrum)
  - Add Atom Feature Vectors 
 (github pull #2988 from thegodone)
  - Add CoulombMat calculator
 (github pull #2993 from thegodone)
  - Update azure-pipelines.yml
 (github pull #2997 from greglandrum)
  - Improve Conformational Sampling of Small and Large Ring Molecules
 (github pull #2999 from hjuinj)
  - Fix atom highlighting in notebook PNGs
 (github pull #3000 from greglandrum)
  - adds a one-liner for getting a vector of random smiles for a molecule
 (github pull #3002 from greglandrum)
  - Allow enhanced stereo to be used in substructure search
 (github pull #3003 from d-b-w)
  - Add support for the rest of the v3000 atom properties
 (github pull #3007 from greglandrum)
  - Move jupyter extension logging to the python logger
 (github pull #3008 from bp-kelley)
  - Commit of 2D draw annotation.
 (github pull #3010 from DavidACosgrove)
  - Update Maeparser & Coordgen Dependencies
 (github pull #3011 from ricrogz)
  - Remove unnecessary files
 (github pull #3012 from e-kwsm)
  - allow retrieval of the atoms/bonds modified by the tautomerization
 (github pull #3013 from greglandrum)
  - Add 5 new recipes to Cookbook
 (github pull #3014 from vfscalfani)
  - Turns on cairo support (and testing) in the Azure DevOps CI builds
 (github pull #3022 from greglandrum)
  - Added support for Python FMCS functors
 (github pull #3023 from ptosco)
  - add random seed to docs to get reproducible conformations
 (github pull #3026 from greglandrum)
  - update docs for 2020.03
 (github pull #3028 from greglandrum)
  - update Getting Started in C++ document
 (github pull #3039 from DavidACosgrove)



## Deprecated code (to be removed in a future release):
- To improve API consistency of the exceptions in RDKit with the default ones in
  the STL, the several `message()` methods and `Invariant::getMessage()` in RDKit's
  exceptions are from now on deprecated in favor of `what()`. Both `message()` and
  `Invariant::getMessage()` will be removed in the next release.
- The old MolHash code should be considered deprecated. This release introduces
  a more flexible alternative. Specifically the following pieces will be removed in the next release:
  - The python functionality `rdkit.Chem.rdMolHash.GenerateMoleculeHashString()`
  - The C++ functionality directly present in the header file `GraphMol/MolHash/MolHash.h`

# Release_2019.09.1
(Changes relative to Release_2019.03.1)

## Important
- The atomic van der Waals radii used by the RDKit were corrected/updated in #2154.
  This leads to different results when generating conformations, molecular volumes,
  and molecular shapes. 

## Backwards incompatible changes
- See the note about atomic van der Waals radii above.
- As part of the enhancements to the MolDraw2D class, we changed the type of
  DrawColour from a tuple to be an actual struct. We also added a 4th element to
  capture alpha values. This should have no affect on Python code (the alpha
  value is optional when providing color tuples), but will require changes to C++
  and Java/C# code that is using DrawColour.
- When reading Mol blocks, atoms with the symbol "R" are now converted into
  queries that match any atom when doing a substructure search (analogous to "*"
  in SMARTS). The previous behavior was to only match other dummy atoms
- When loading SDF files using PandasTools.LoadSDF(), we now default to
  producing isomeric smiles in pandas tables.  To reproduce the original
  behavior, use isomericSmiles=False in the call to the function.
- The SMARTS generated by the RDKit no longer contains redundant wildcard
  queries. This means the SMARTS strings generated by this release will generally
  be different from that in previous releases, although the results produced by
  the queries should not change.
- The RGroupDecomposition code now removes Hs from output R groups by default.
  To restore the old behavior create an RGroupDecompositionParameters object and
  set removeHydrogensPostMatch to false.
- The default values for some of the new fingerprint generators have been changed so
  that they more closely resemble the original fingerprinting code. In
  particular most fingerprinters no longer do count simulation by default and
  the RDKit fingerprint now sets two bits per feature by default.
- The SMARTS generated for MCS results using the ringMatchesRingOnly or
  completeRingsOnly options now includes ring-membership queries.

## Highlights:
- The substructure matching code is now about 30% faster. This also improves the
  speed of reaction matching and the FMCS code. (#2500)
- A minimal JavaScript wrapper has been added as part of the core release. (#2444)
- It's now possible to get information about why molecule sanitization failed. (#2587)
- A flexible new molecular hashing scheme has been added. (#2636)

## Acknowledgements:
Patricia Bento, Francois Berenger, Jason Biggs, David Cosgrove, Andrew Dalke,
Thomas Duigou, Eloy Felix, Guillaume Godin, Lester Hedges, Anne Hersey,
Christoph Hillisch, Christopher Ing, Jan Holst Jensen, Gareth Jones, Eisuke
Kawashima, Brian Kelley, Alan Kerstjens, Karl Leswing, Pat Lorton, John
Mayfield, Mike Mazanetz, Dan Nealschneider, Noel O'Boyle, Stephen Roughley,
Roger Sayle, Ricardo Rodriguez Schmidt, Paula Schmiel, Peter St. John, Marvin
Steijaert, Matt Swain, Amol Thakkar Paolo Tosco, Yi-Shu Tu, Ricardo Vianello,
Marc Wittke, '7FeiW', 'c56pony', 'sirbiscuit' 


## Bug Fixes:
  - MCS returning partial rings with completeRingsOnly=True 
 (github issue #945 from greglandrum)
  - Alternating canonical SMILES for fused ring with N
 (github issue #1028 from greglandrum)
  - Atom index out of range error
 (github issue #1868 from A-Thakkar)
  - Incorrect cis/trans stereo symbol for conjugated ring
 (github issue #2023 from baoilleach)
  - Hardcoded max length of SMARTs string cut of input query for FragCatlog
 (github issue #2163 from 7FeiW)
  - VSA_EState {1, ..., 10} calculated by rdkit doesn't seem correct.
 (github issue #2372 from c56pony)
  - MolStandardize: FragmentRemover should not sanitize fragments
 (github issue #2411 from greglandrum)
  - MolStandardize: combinatorial explosion in Normalizer
 (github issue #2414 from greglandrum)
  - MCS code doesn't return envelope MCS if CompleteRingsOnly is true
 (github issue #2420 from greglandrum)
  - RemoveHs() does not remove all hydrogens.
 (github issue #2422 from paulaju)
  - Incorrect assignment of explicit Hs to Al+3 read from mol block
 (github issue #2423 from greglandrum)
  - Cannot set maxProducts > 1000 in RunReactants
 (github issue #2427 from tduigou)
  - Chem.MolStandardize.standardize.Standardizer drops molecular properties
 (github pull #2431 from lilleswing)
  - Canon::rankMolAtoms results in crossed double bonds in rings
 (github issue #2437 from greglandrum)
  - make boost::iostreams optional
 (github pull #2440 from greglandrum)
  - Fix/rgroup sdf isotope
 (github pull #2449 from bp-kelley)
  - Uncharger incorrectly removing charge from boron anions
 (github issue #2452 from greglandrum)
  - Add java builds to azure devops
 (github pull #2460 from greglandrum)
  - Cart fixes
 (github pull #2462 from jones-gareth)
  - Typo in UFF torsion terms for SP2-SP2
 (github issue #2463 from jasad1)
  - Negative atom map values cause depickling to fail
 (github issue #2465 from greglandrum)
  - Deserialization failures crash Java wrapper
 (github issue #2466 from greglandrum)
  - rdkit.six fix and cleanup
 (github pull #2469 from rvianello)
  - dummy atom queries are flagged as complex
 (github issue #2471 from greglandrum)
  - 3D structure display broken in jupyter notebook
 (github issue #2473 from greglandrum)
  - Inconsistent defaults for nonBondedThresh in MMFF optimization
 (github issue #2475 from greglandrum)
  - Fix/rgroup multiple labels
 (github pull #2481 from bp-kelley)
  - 2D Depiction clipped atom highlighting
 (github issue #2486 from DavidACosgrove)
  - BRICSBuild now passes scrambleReagents to children
 (github pull #2488 from greglandrum)
  - Pattern Fingerprint Issues
 (github issue #2501 from jones-gareth)
  - CMake Error: Wrap directories being used when python build is turned off
 (github issue #2516 from jasondbiggs)
  - - fixes ResonanceMolSupplier bug in perceiving conjugated groups
 (github pull #2517 from ptosco)
  - Fix/mmff threadsafety issues
 (github pull #2518 from bp-kelley)
  - update expected SVG output in cartridge tests
 (github pull #2520 from greglandrum)
  - fix to SDWriter docs
 (github pull #2521 from pstjohn)
  - Fix the azure pipelines builds
 (github pull #2522 from greglandrum)
  - Code cleanups from PVS/Studio
 (github pull #2531 from greglandrum)
  - getAtomNeighbors() and getAtomBonds() not in SWIG wrappers.
 (github issue #2532 from greglandrum)
  - Default sanitizerxn doesn't aromatize if possible
 (github issue #2547 from bp-kelley)
  - Add RDKIT_FILEPARSERS_EXPORT to finishMolProcessing
 (github pull #2551 from d-b-w)
  - Chem.rdFMCS.FindMCS hangs for certain ligand pairs
 (github issue #2581 from lohedges)
  - fix the inclusion path for the targets file (#2584)
 (github pull #2589 from rvianello)
  - Fix inocuous typo/bug in Dative bond matching
 (github pull #2593 from ricrogz)
  - E/Z and CIS/TRANS stereo bonds are incorrectly matched
 (github pull #2596 from ricrogz)
  - Uncharger ignores dications
 (github issue #2602 from greglandrum)
  - Possible Garbage Collection Bug in Pharmacophore Generation
 (github issue #2603 from cing)
  - Uncharger incorrectly neutralizes cations when non-neutralizable anions are present.
 (github issue #2605 from greglandrum)
  - Bad valence corrections on Pb, Sn
 (github issue #2606 from greglandrum)
  - Pb and Sn should support valence 2
 (github issue #2607 from greglandrum)
  - Uncharger incorrectly modifying a zwitterion
 (github issue #2610 from greglandrum)
  - CanonicalRankAtomsInFragment breakTies doesn't
 (github issue #2611 from bp-kelley)
  - Pattern fingerprint failing substructure condition in very large molecules
 (github issue #2614 from greglandrum)
  - Memory leak with Chem.Atom() constructor
 (github issue #2639 from AlanKerstjens)
  - Fixes: Atoms in non-standard valences not being properly written to mol blocks
 (github pull #2646 from greglandrum)
  - C++ MCS code returns a null MCS between methylcyclopentane and methylcyclohexane
 (github issue #2662 from ptosco)
  - CXSMILES writer has error if mol comes from v3000 molfile
 (github issue #2666 from d-b-w)
  - MolToCXSmiles generates error for empty molecule
 (github issue #2667 from greglandrum)
  - fix a problem with normalize, ringinfo, and fragments
 (github pull #2685 from greglandrum)
  - Error when a squiggle bond is in an aromatic ring
 (github issue #2695 from greglandrum)
  - Cannot combine multiple range queries on a single atom.
 (github issue #2709 from greglandrum)
  - setBondStereoFromDirections() returning incorrect results.
 (github issue #2712 from greglandrum)
  - update supplier documentation to reflect python 3 iterator syntax
 (github pull #2719 from greglandrum)
  - removeHs messing up double bond stereo in partially sanitized molecules
 (github issue #2721 from greglandrum)
  - seg fault in ReactionRunner
 (github issue #2722 from greglandrum)
  - Intermittent test failures for JavaDistanceGeometryTests
 (github issue #2727 from greglandrum)
  - Fixes a bug in TorsionConstraint
 (github pull #2732 from ptosco)
  - Apply fix for #1592 to _MolsToGridSVG
 (github pull #2737 from yishutu)

## New Features and Enhancements:
  - Added rankAtoms to ROMol wrapper and added Java test case
 (github pull #1540 from sroughley)
  - Use van der Waals radii from blue obelisk
 (github pull #2154 from UnixJunkie)
  - add generateDepictionMatching2DStructure() to SWIG wrappers
 (github issue #2239 from greglandrum)
  - Added OptimizeMoleculeConfs with pre-generated force-field
 (github pull #2401 from ptosco)
  - FreeSASA improvements
 (github pull #2402 from ptosco)
  - Speed up symmetrizeSSSR
 (github issue #2403 from d-b-w)
  - Trim whitespace from mol fragment SMARTS and check SMARTS presence
 (github pull #2406 from ricrogz)
  - Run clang-tidy over the entire codebase
 (github pull #2408 from greglandrum)
  - Enable Azure Pipelines builds for CI
 (github pull #2409 from ricrogz)
  - Add RDProps interface to Conformers
 (github issue #2441 from greglandrum)
  - Add minimal JavaScript wrapper
 (github pull #2444 from greglandrum)
  - Fixes annoying warnings on MSVC
 (github pull #2454 from ptosco)
  - add prepareMolsBeforeDrawing option for drawMols 
 (github pull #2455 from greglandrum)
  - computeGasteigerCharges quality of life improvement for python api
 (github issue #2480 from bp-kelley)
  - Preserve bond direction in fragmentOnBonds
 (github pull #2484 from greglandrum)
  - SanitizeRxn code and docstring cleanup
 (github pull #2491 from greglandrum)
  - Support XYZ format for output
 (github pull #2498 from e-kwsm)
  - vf2 optimisations
 (github pull #2500 from johnmay)
  - Python wrap enhanced stereo setters
 (github pull #2509 from d-b-w)
  - Fix the azure pipelines builds
 (github pull #2522 from greglandrum)
  - add a script for benchmarking fingerprint screenout and substructure performance
 (github pull #2523 from greglandrum)
  - make "R" in CTABs an AtomNull query
 (github pull #2528 from greglandrum)
  - Expose SDF loading options to LoadSDF
 (github pull #2534 from bp-kelley)
  - Remove unused ctest: testCanon
 (github pull #2541 from ricrogz)
  - Update maeparser and coordgen versions
 (github pull #2542 from ricrogz)
  - Improved handling of bond stereo in reactions
 (github pull #2553 from ricrogz)
  - Code simplification for fingerprints to np array
 (github pull #2557 from ChrisHill8)
  - Integrate Unique Ring Families from RingDecomposerLib 
 (github pull #2558 from greglandrum)
  - Allow providing a bounds matrix to EmbedMol
 (github pull #2560 from greglandrum)
  - Enable SimilarityMaps in C++
 (github pull #2562 from greglandrum)
  - Do not run UnitTestMCS.py::TestTimeout in debug builds
 (github pull #2569 from ricrogz)
  - Expose more drawing methods to Python
 (github issue #2571 from greglandrum)
  - Allow Point2D to be constructed from Point3D
 (github pull #2572 from greglandrum)
  - Allows dative bonds to be drawn
 (github pull #2573 from greglandrum)
  - Allow identification of chemistry problems
 (github pull #2587 from greglandrum)
  - Adds MolFragmentToSmarts to generate smarts for a subset of a Molecule
 (github pull #2594 from d-b-w)
  - Removal of redundant wildcards in SMARTS (Null Atom/Bond Query combination)
 (github pull #2595 from ricrogz)
  - Support range-based charge queries from SMARTS
 (github issue #2604 from greglandrum)
  - Keep PDB info from Maestro files if available
 (github pull #2619 from lorton)
  - optimization of the MolStandardize code
 (github pull #2621 from greglandrum)
  - Assign stereo bond labels in molecules read from SMARTS
 (github pull #2623 from ricrogz)
  - Automatically load the ipython extensions running in Jupyter
 (github pull #2626 from bp-kelley)
  - draw zero-order bonds
 (github pull #2630 from greglandrum)
  - Updated cartridge documentation
 (github pull #2635 from msteijaert)
  - Add new mol hashing code
 (github pull #2636 from greglandrum)
  - emolecules link updated
 (github pull #2638 from marcwittke)
  - Update maeparser to 1.2.1 and coorgen to 1.3.1
 (github pull #2652 from ricrogz)
  - Get numpy include dir programmatically
 (github pull #2653 from sirbiscuit)
 - Fix long columns pandas
 (github pull #2655 from sirbiscuit)
  - Added AtomComparator.AtomCompareAnyHeavyAtom and test cases to FMCS code
 (github pull #2656 from sroughley)
  - The C++ MCS code generates ambiguous SMARTS strings
 (github issue #2663 from ptosco)
  - add bond-selector info to SVGs
 (github pull #2664 from greglandrum)
  - support writing CXSMILES from the cartridge
 (github issue #2668 from greglandrum)
  - support the new hashing code in the cartridge
 (github pull #2671 from greglandrum)
  - Adds additional capabilities to the minimal JS wrapper
 (github pull #2676 from greglandrum)
  - Add MolHash to Java Wrappers
 (github issue #2677 from sroughley)
  - A bunch of changes to the new fingerprinter code
 (github pull #2679 from greglandrum)
  - Add viewBox to default SVG output
 (github issue #2680 from bp-kelley)
  - Allow Java to see RGroup labels in the std::map wrapper.
 (github pull #2681 from bp-kelley)
  - Update maeparser to v1.2.2.
 (github pull #2682 from ricrogz)
  - Update coordgen to v1.3.2
 (github pull #2686 from ricrogz)
  - Add a drawOptions object to IPythonConsole
 (github pull #2691 from greglandrum)
  - Make StructureGroups editable from Python
 (github pull #2692 from greglandrum)
  - Update documentation
 (github pull #2697 from greglandrum)
  - Make removeHydrogensPostMatch=true the default for RGD
 (github pull #2701 from greglandrum)
  - Eat our own dogfood, Clone is deprecated so use copy
 (github pull #2711 from bp-kelley)
  - The MCS smartsString may still be ambiguous
 (github issue #2714 from ptosco)
  - Add threaded runner for the filter catalog
 (github pull #2718 from bp-kelley)
  - Add Leader picker implementation
 (github pull #2724 from greglandrum)
  - Add consideration of ring fusion to the MCS algorithm
 (github pull #2731 from ptosco)


## Deprecated code (to be removed in a future release):

- The old MolHash code should be considered deprecated. This release introduces
  a more flexible alternative. Specifically the following pieces will be removed in a future release:
  - The python functionality `rdkit.Chem.rdMolHash.GenerateMoleculeHashString()`
  - The C++ functionality directly present in the header file `GraphMol/MolHash/MolHash.h`


# Release_2019.03.1
(Changes relative to Release_2018.09.1)

## REALLY IMPORTANT ANNOUNCEMENT
- As of this realease (2019.03.1) the RDKit no longer supports Python 2. Please
  read this rdkit-discuss post to learn what your options are if you need to
  keep using Python 2:
  https://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg08354.html

## Backwards incompatible changes
- The fix for github #2245 means that the default behavior of the MaxMinPicker
  is now truly random. If you would like to reproduce the previous behavior,
  provide a seed value of 42.
- The uncharging method in the MolStandardizer now attempts to generate
  canonical results for a given molecule. This may result in different output
  for some molecules.
  
## Highlights:
- There's now a Japanese translation of large parts of the RDKit documentation
- SGroup data can now be read from and written to Mol/SDF files
- The enhanced stereo handling has been improved: the information is now
  accessible from Python, EnumerateStereoisomers takes advantage of it, and it
  can be read from and written to CXSmiles

## Acknowledgements:
Michael Banck, Francois Berenger, Thomas Blaschke, Brian Cole, Andrew Dalke,
Bakary N'tji Diallo, Guillaume Godin, Anne Hersey, Jan Holst Jensen, Sunhwan Jo,
Brian Kelley, Petr Kubat, Karl Leswing, Susan Leung, John Mayfield, Adam Moyer,
Dan Nealschneider, Noel O'Boyle, Stephen Roughley, Takayuki Serizawa, Gianluca
Sforna, Ricardo Rodriguez Schmidt, Gianluca Sforna, Matt Swain, Paolo Tosco, 
Ricardo Vianello, 'John-Videogames', 'magattaca', 'msteijaert', 'paconius', 
'sirbiscuit' 

## Bug Fixes:
  - PgSQL: fix boolean definitions for Postgresql 11
 (github pull #2129 from pkubatrh)
  - update fingerprint tutorial notebook
 (github pull #2130 from greglandrum)
  - Fix typo in RecapHierarchyNode destructor
 (github pull #2137 from iwatobipen)
  - SMARTS roundtrip failure
 (github issue #2142 from mcs07)
  - Error thrown in rdMolStandardize.ChargeParent
 (github issue #2144 from paconius)
  - SMILES parsing inconsistency based on input order
 (github issue #2148 from coleb)
  - MolDraw2D: line width not in python wrapper
 (github issue #2149 from greglandrum)
  - Missing Python API Documentation
 (github issue #2158 from greglandrum)
  - PgSQL: mol_to_svg() changes input molecule.
 (github issue #2174 from janholstjensen)
  - Remove Unicode From AcidBasePair Name
 (github pull #2185 from lilleswing)
  - Inconsistent treatment of `[as]` in SMILES and SMARTS
 (github issue #2197 from greglandrum)
  - RGroupDecomposition fixes, keep userLabels more robust onlyMatchAtRGroups
 (github pull #2202 from bp-kelley)
  - Fix TautomerTransform in operator=
 (github pull #2203 from bp-kelley)
  - testEnumeration hangs/takes where long on 32bit architectures
 (github issue #2209 from mbanck)
  - Silencing some Python 3 warning messages
 (github pull #2223 from coleb)
  - removeHs shouldn't remove atom lists
 (github issue #2224 from rvianello)
  - failure round-tripping mol block with Q atom
 (github issue #2225 from rvianello)
  - problem round-tripping mol files that include bond topology info
 (github issue #2229 from rvianello)
  - aromatic main-group atoms written to SMARTS incorrectly
 (github issue #2237 from greglandrum)
  - findPotentialStereoBonds() stopping too early
 (github issue #2244 from greglandrum)
  - MinMax Diversity picker seeding shows deterministic / non-random behaviour
 (github issue #2245 from sroughley)
  - Fix to serialize binary strings
 (github pull #2264 from bp-kelley)
  - Recognize N in three-membered rings as potentially chiral
 (github issue #2268 from greglandrum)
  - Failure when parsing mol block with M  PXA
 (github issue #2277 from greglandrum)
  - query-query matching failing for atoms constructed from SMARTS
 (github issue #2299 from greglandrum)
  - SMILES parsing fails for dative ring closures
 (github issue #2303 from greglandrum)
  - Missing Dict.h overload: std::string Dict::getVal<std::string>
 (github issue #2308 from greglandrum)
  - fix a problem with the random pickers test
 (github pull #2310 from greglandrum)
  - Some tests still failing on less common platforms.
 (github issue #2311 from giallu)
  - AddHs() using 3D coordinates with 2D conformations
 (github pull #2328 from greglandrum)
  - change to make the SWIG builds work on windows
 (github pull #2340 from greglandrum)
  - uncharger behaves differently on molecules constructed from mol blocks and SMILES
 (github issue #2346 from greglandrum)
  - Memory Error When Writing ToBinary With "AllProps"
 (github issue #2352 from atom-moyer)
  - Seg fault on init if RDBASE is not set
 (github issue #2368 from greglandrum)
  - PandasTools.FrameToGridImage() fails with SVG output
 (github issue #2380 from greglandrum)
  - ClusterMols.GetDistanceMatrix throws a type error in Python 3
 (github issue #2387 from John-Videogames)
  - Uncharging logic reversed: protonate non-acids first 
 (github issue #2392 from Anne Hersey)

## New Features and Enhancements:
  - Allow access to Enhanced Stereochemistry information from Python
 (github issue #2108 from d-b-w)
  - Adopt EnumerateStereoisomers to use extended stereo
 (github issue #2109 from greglandrum)
  - Enable ctest -T memcheck
 (github pull #2113 from ricrogz)
  - Support for parsing/writing SGroups in SD Mol files
 (github pull #2138 from ricrogz)
   - Rename the #define _DEBUG to MMPA_DEBUG in mmpa.cpp
 (github pull #2140 from baoilleach)
  - MolDraw2D: line width should be controlled by MolDrawOptions
 (github issue #2151 from greglandrum)
  - Some refactoring of the distance geometry code
 (github pull #2153 from greglandrum)
  - Less warnings
 (github pull #2155 from UnixJunkie)
  - ShapeTverskyIndex python function
 (github pull #2156 from susanhleung)
  - Skip compound if smiles conversion fails
 (github pull #2168 from msteijaert)
  - Fix #2176: InChI functions should return NULL on un-InChI-able input molecules.
 (github pull #2177 from janholstjensen)
  - Update installation instructions for Linux
 (github pull #2181 from sirbiscuit)
  - Update CMake rules to find external coorgen & maeparser libs
 (github pull #2184 from ricrogz)
  - Update to use the travis Xenial environment
 (github pull #2200 from greglandrum)
  - Do not allow PandasTools to overwrite pandas settings
 (github pull #2206 from sirbiscuit)
  - re-enable (and update) the file parser tests
 (github pull #2208 from greglandrum)
  - Added documentation files written in Japanese into Book directory
 (github pull #2210 from magattaca)
  - Add C++ convenience function for drawing ROMols
 (github issue #2220 from greglandrum)
  - Change boost int types to std types
 (github pull #2233 from bp-kelley)
  - Added exports for SGroup functions
 (github pull #2242 from ricrogz)
  - Use coordMap when starting embedding from random coords
 (github issue #2246 from greglandrum)
  - Improve interactivity of output SVG
 (github pull #2253 from greglandrum)
  - Add options for substructure searching
 (github pull #2254 from greglandrum)
  - keep extra information about bonds from Mol files
 (github pull #2260 from greglandrum)
  - Allow converting mol blocks directly to InChI
 (github pull #2262 from greglandrum)
  - Patch/pains updates
 (github pull #2272 from johnmay)
  - add warning for 2D conformations flagged as 3D
 (github pull #2273 from greglandrum)
  - Store extra CXSMILES data as a property
 (github pull #2281 from ricrogz)
  - Parse enhanced stereo information from CXSMILES
 (github pull #2282 from ricrogz)
  - Robustify parsing of CTABs and SGROUPs
 (github pull #2283 from greglandrum)
  - Write enhanced stereo to cxsmiles
 (github pull #2290 from greglandrum)
  - Allow custom type-handlers in the RDProps interface
 (github pull #2293 from bp-kelley)
  - Add serialization to SubstructLibrary
 (github pull #2295 from bp-kelley)
  - support reading/writing atom props from SD files
 (github pull #2297 from greglandrum)
  - Add test for issue #2285, fix molbundle test
 (github pull #2301 from bp-kelley)
  - Update maeparser & coordgen libraries
 (github pull #2302 from ricrogz)
  - Mem errors clean up
 (github pull #2305 from ricrogz)
  - Add definition of MolFragmentToCXSmiles
 (github pull #2307 from greglandrum)
  - Doc update
 (github pull #2312 from greglandrum)
  - Adds gzstream stream, exposes to swig
 (github pull #2314 from bp-kelley)
  - Remove a bunch of Python2-related warts
 (github pull #2315 from greglandrum)
  - some much-needed optimization work on the new property lists
 (github pull #2317 from greglandrum)
  - Build warnings revisited
 (github pull #2318 from ricrogz)
  - change bogus "3D" to "2D" in a test file
 (github pull #2319 from greglandrum)
  - Allow copying atoms in Python
 (github pull #2322 from d-b-w)
  - fixes an r-group symmetrization problem
 (github pull #2324 from greglandrum)
  - simple docstring fix
 (github pull #2326 from sunhwan)
  - allow using system's catch2 for tests
 (github pull #2327 from giallu)
  - Python wrap DetectAtomStereoChemistry from MolFileStereochem.h
 (github issue #2329 from d-b-w)
  - switch to using cmake to handle the C++ spec
 (github pull #2334 from greglandrum)
  - WIP: optional integration with YAeHMOP
 (github pull #2335 from greglandrum)
  - Exposes substructlibrary to swig
 (github pull #2337 from bp-kelley)
  - Add a skip_all_if_match option to the FragmentRemover
 (github pull #2338 from greglandrum)
  - Dev/general csharp fixes
 (github pull #2341 from bp-kelley)
  - Add a read-only Python wrapper for SGroups
 (github pull #2343 from greglandrum)
  - Expose RGroupDecomposition to SWIG
 (github pull #2345 from greglandrum)
  - update debian build script to python3
 (github pull #2350 from UnixJunkie)
  - add GetStereoIsomerCount() function to EnumerateStereoisomers
 (github pull #2354 from greglandrum)
  - Update coordgenlibs to v1.2.2
 (github pull #2355 from ricrogz)
  - Small fixes to get DLLs to build on Windows
 (github pull #2356 from ptosco)
  - Boost deprecation warning
 (github pull #2357 from d-b-w)
  - Removes an extra debugging cerr statment
 (github pull #2360 from d-b-w)
  - Preserve enhanced stereo in reactions 
 (github pull #2366 from d-b-w)
  - improvements to the Uncharge functionality
 (github pull #2374 from greglandrum)
  - Add ExplicitBitVect prop and query
 (github pull #2384 from bp-kelley)
  - Allow components of the MolStandardize code to be initialized from streams
 (github pull #2385 from greglandrum)


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
 - Dummy atoms now have a default valence of 0 and no maximum
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
 - A segmentation fault that occurred when kekulizing modified
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
 - The replaceSidechains and replaceCore functions have been added
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
