"""
Module containing implementation of SynthonSpace search of
Synthon-based chemical libraries such as Enamine REAL.
  NOTE: This functionality is experimental and the API
and/or results may change in future releases.
"""

from collections.abc import Callable
import os
from typing import overload

import rdkit.Chem.rdEnumerateStereoisomers
import rdkit.Chem.rdFingerprintGenerator
import rdkit.Chem.rdGaussianShape
import rdkit.Chem.rdGeneralizedSubstruct
import rdkit.Chem.rdRascalMCES
import rdkit.Chem.rdchem


class SubstructureResult:
    """Used to return results of SynthonSpace searches."""

    def GetHitMolecules(self) -> list[rdkit.Chem.rdchem.Mol]:
        """A function returning hits from the search"""

    def GetMaxNumResults(self) -> int:
        """
        The upper bound on number of results possible.  There
        may be fewer than this in practice for several reasons
        such as duplicate reagent sets being removed or the
        final product not matching the query even though the
        synthons suggested they would.
        """

    def GetTimedOut(self) -> bool:
        """Returns whether the search timed out or not."""

    def GetCancelled(self) -> bool:
        """Returns whether the search was cancelled or not."""

    def GetBestHit(self) -> object:
        """
        Returns the best hit found in the similarity search, even when none were under the search threshold.  May be empty if not a similarity search or nothing came close to a match.
        """

class SynthonSpaceSearchParams:
    """SynthonSpaceSearch parameters."""

    def __init__(self) -> None: ...

    @property
    def maxHits(self) -> int:
        """
        The maximum number of hits to return.  Default=1000.Use -1 for no maximum.
        """

    @maxHits.setter
    def maxHits(self, arg: int, /) -> None: ...

    @property
    def maxNumFrags(self) -> int:
        """
        The maximum number of fragments the query can be broken into.
          Big molecules will create huge numbers of fragments that may cause
        excessive memory use.  If the number of fragments hits this number,
        fragmentation stops and the search results will likely be incomplete.
          Default=100000.
        """

    @maxNumFrags.setter
    def maxNumFrags(self, arg: int, /) -> None: ...

    @property
    def hitStart(self) -> int:
        """
        The sequence number of the hit to start from.  So that you
        can return the next N hits of a search having already
        obtained N-1.  Default=0
        """

    @hitStart.setter
    def hitStart(self, arg: int, /) -> None: ...

    @property
    def randomSample(self) -> bool:
        """
        If True, returns a random sample of the hits, up to maxHits
        in number.  Default=False.
        """

    @randomSample.setter
    def randomSample(self, arg: bool, /) -> None: ...

    @property
    def randomSeed(self) -> int:
        """
        If using randomSample, this seeds the random number
        generator so as to give reproducible results.  Default=-1
        means use a random seed.
        """

    @randomSeed.setter
    def randomSeed(self, arg: int, /) -> None: ...

    @property
    def buildHits(self) -> bool:
        """
        If false, reports the maximum number of hits that
        the search could produce, but doesn't return them.
        """

    @buildHits.setter
    def buildHits(self, arg: bool, /) -> None: ...

    @property
    def similarityCutoff(self) -> float:
        """
        Similarity cutoff for returning hits by fingerprint similarity.
          At present the fp is hard-coded to be Morgan, bits, radius=2.
          Default=0.5.
        """

    @similarityCutoff.setter
    def similarityCutoff(self, arg: float, /) -> None: ...

    @property
    def fragSimilarityAdjuster(self) -> float:
        """
        Similarities of fragments are generally low due to low bit
        densities.  For the fragment matching, reduce the similarity cutoff
        off by this amount.  Default=0.1.
        """

    @fragSimilarityAdjuster.setter
    def fragSimilarityAdjuster(self, arg: float, /) -> None: ...

    @property
    def approxSimilarityAdjuster(self) -> float:
        """
        The fingerprint search uses an approximate similarity method
        before building a product and doing a final check.  The
        similarityCutoff is reduced by this value for the approximate
        check.  A lower value will give faster run times at the
        risk of missing some hits.  The value you use should have a
        positive correlation with your FOMO.  The default of 0.1 is
        appropriate for Morgan fingerprints.  With RDKit fingerprints,
        0.05 is adequate, and higher than that has been seen to
        produce long run times.
        """

    @approxSimilarityAdjuster.setter
    def approxSimilarityAdjuster(self, arg: float, /) -> None: ...

    @property
    def minHitHeavyAtoms(self) -> int:
        """Minimum number of heavy atoms in a hit.  Default=0."""

    @minHitHeavyAtoms.setter
    def minHitHeavyAtoms(self, arg: int, /) -> None: ...

    @property
    def maxHitHeavyAtoms(self) -> int:
        """Maximum number of heavy atoms in a hit.  Default=-1 means no maximum."""

    @maxHitHeavyAtoms.setter
    def maxHitHeavyAtoms(self, arg: int, /) -> None: ...

    @property
    def minHitMolWt(self) -> float:
        """Minimum molecular weight for a hit.  Default=0.0."""

    @minHitMolWt.setter
    def minHitMolWt(self, arg: float, /) -> None: ...

    @property
    def maxHitMolWt(self) -> float:
        """Maximum molecular weight for a hit.  Default=0.0 mean no maximum."""

    @maxHitMolWt.setter
    def maxHitMolWt(self, arg: float, /) -> None: ...

    @property
    def minHitChiralAtoms(self) -> int:
        """Minimum number of chiral atoms in a hit.  Default=0."""

    @minHitChiralAtoms.setter
    def minHitChiralAtoms(self, arg: int, /) -> None: ...

    @property
    def maxHitChiralAtoms(self) -> int:
        """Maximum number of chiral atoms in a hit.  Default=-1 means no maximum."""

    @maxHitChiralAtoms.setter
    def maxHitChiralAtoms(self, arg: int, /) -> None: ...

    @property
    def numConformers(self) -> int:
        """
        When doing a shape search, the number of conformers to generate for molecules.  Default=100.
        """

    @numConformers.setter
    def numConformers(self, arg: int, /) -> None: ...

    @property
    def confRMSThreshold(self) -> float:
        """
        When doing a shape search, the RMS threshold to use when pruning conformers.  Default=1.0.
        """

    @confRMSThreshold.setter
    def confRMSThreshold(self, arg: float, /) -> None: ...

    @property
    def shapeOverlayOptions(self) -> rdkit.Chem.rdGaussianShape.ShapeOverlayOptions:
        """Options for the shape overlays."""

    @shapeOverlayOptions.setter
    def shapeOverlayOptions(self, arg: rdkit.Chem.rdGaussianShape.ShapeOverlayOptions, /) -> None: ...

    @property
    def bestHit(self) -> bool:
        """
        If True, when doing a shape search it will return the hit conformer with the best shape match to the query conformer.  If False, it just returns the first hit conformer that exceeds the similarity cutoff.  The latter will be faster but the returned hit conformations are likely to be less relevant.
        """

    @bestHit.setter
    def bestHit(self, arg: bool, /) -> None: ...

    @property
    def enumerateUnspecifiedStereo(self) -> bool:
        """
        When doing a shape search, if there is unspecified stereochemistry in either the query or potential hit, enumerate test all possibilities.  Default=False.
        """

    @enumerateUnspecifiedStereo.setter
    def enumerateUnspecifiedStereo(self, arg: bool, /) -> None: ...

    @property
    def stereoEnumOpts(self) -> rdkit.Chem.rdEnumerateStereoisomers.StereoEnumerationOptions:
        """Options for stereoisomer enumeration."""

    @stereoEnumOpts.setter
    def stereoEnumOpts(self, arg: rdkit.Chem.rdEnumerateStereoisomers.StereoEnumerationOptions, /) -> None: ...

    @property
    def timeOut(self) -> int:
        """
        Time limit for search, in seconds.  Default is 600s, 0 means no
        timeout.  Requires an integer
        """

    @timeOut.setter
    def timeOut(self, arg: int, /) -> None: ...

    @property
    def toTryChunkSize(self) -> int:
        """Process possible hits using the given chunk size"""

    @toTryChunkSize.setter
    def toTryChunkSize(self, arg: int, /) -> None: ...

    @property
    def numThreads(self) -> int:
        """
        The number of threads to use for search.  If > 0, will use that
        number.  If <= 0, will use the number of hardware
        threads plus this number.  So if the number of
        hardware threads is 8, and numThreads is -1, it will
        use 7 threads.  Default=1.
        """

    @numThreads.setter
    def numThreads(self, arg: int, /) -> None: ...

    @property
    def useProgressBar(self) -> int:
        """
        Makes a progress bar of given width.  The number given is the number of '*' characters in a full bar.  There will be about another 35 characters or so depending on the size of the job.  Default=0 means no bar.
        """

    @useProgressBar.setter
    def useProgressBar(self, arg: int, /) -> None: ...

    @property
    def excludedVolume(self) -> object:
        """
        Add an excluded volume to use in the shape search.  The volume overlap and mean overlap over clashing atoms will be reported.
        """

    @excludedVolume.setter
    def excludedVolume(self, arg: rdkit.Chem.rdGaussianShape.ShapeInput | None) -> None: ...

    @property
    def maxExcludedVolume(self) -> float:
        """
        Maximum allowed excluded volume for a hit to be accepted.  Default -1.0 means no maximum.
        """

    @maxExcludedVolume.setter
    def maxExcludedVolume(self, arg: float, /) -> None: ...

    @property
    def maxMeanExcludedVolume(self) -> float:
        """
        Maximum mean excluded volume for a hit to be accepted.  The mean is the total excluded volume divided by the number of clashing atoms (within 2 CARBON_RAD of an excluded volume atom).  To try and distinguish between a mild clash over the whole hit and a few atoms having a really bad clash.
        """

    @maxMeanExcludedVolume.setter
    def maxMeanExcludedVolume(self, arg: float, /) -> None: ...

    @property
    def possibleHitsFile(self) -> str:
        """
        Name of a file to save the possible hits to. These are the combinations of synthons that might match the query but need building and final checking.  Each line has a space-separated list of the synthons and the hit's name.  The file will be emptied and re-filled if it already exists.
        """

    @possibleHitsFile.setter
    def possibleHitsFile(self, arg: str, /) -> None: ...

    @property
    def maxPossibleHitsToWrite(self) -> int:
        """
        Maximum number of lines to write to possibleHitsFile.  When dealing with huge synthon spaces it's very easy to fill a disk.  Default=10M.
        """

    @maxPossibleHitsToWrite.setter
    def maxPossibleHitsToWrite(self, arg: int, /) -> None: ...

    @property
    def writePossibleHitsAndStop(self) -> bool:
        """
        If True, creates the possibleHitsFile and stops without doing the final building and checking.  Default is False.
        """

    @writePossibleHitsAndStop.setter
    def writePossibleHitsAndStop(self, arg: bool, /) -> None: ...

    def setUserConformerGenerator(self, func: Callable[[str, int], rdkit.Chem.rdchem.Mol | None]) -> None:
        """
        Allows you to provide a function that will be called instead of the default
         conformer generator to generate conformers for the synthons.  The function should
         take a SMILES string and the maximum number of conformers to generated and
         return a molecule object.
        """

    def __setattr__(self, name: str, value: object | None) -> None: ...

class ShapeBuildParams:
    """Parameters for building shape objects for SynthonSpaceSearch."""

    def __init__(self) -> None: ...

    @property
    def numConfs(self) -> int:
        """Maximum number of conformers per synthon or query.  Default=10"""

    @numConfs.setter
    def numConfs(self, arg: int, /) -> None: ...

    @property
    def rmsThreshold(self) -> float:
        """RMS threshold to use when pruning conformations.  Default=1.0."""

    @rmsThreshold.setter
    def rmsThreshold(self, arg: float, /) -> None: ...

    @property
    def shapeSimThreshold(self) -> float:
        """
        When generating shapes, similarity threshold for pruning.  No 2 shapes for each synthon or query will be more similar than this threshold.  Default=1.9.
        """

    @shapeSimThreshold.setter
    def shapeSimThreshold(self, arg: float, /) -> None: ...

    @property
    def numThreads(self) -> int:
        """
        The number of threads to use for shape building.  If > 0, will use that number.  If <= 0, will use the number of hardware threads plus this number.Default=1.
        """

    @numThreads.setter
    def numThreads(self, arg: int, /) -> None: ...

    @property
    def randomSeed(self) -> int:
        """
        Seed for random number generator.  Default=-1 means use system random seed.
        """

    @randomSeed.setter
    def randomSeed(self, arg: int, /) -> None: ...

    @property
    def stereoEnumOpts(self) -> rdkit.Chem.rdEnumerateStereoisomers.StereoEnumerationOptions:
        """Options for stereoisomer enumeration."""

    @stereoEnumOpts.setter
    def stereoEnumOpts(self, arg: rdkit.Chem.rdEnumerateStereoisomers.StereoEnumerationOptions, /) -> None: ...

    @property
    def useProgressBar(self) -> int:
        """
        Makes a progress bar of given width.  The number given is the number of '*' characters in a full bar.  There will be about another 35 characters or so depending on the size of the job.  Default=0 means no bar.
        """

    @useProgressBar.setter
    def useProgressBar(self, arg: int, /) -> None: ...

    @property
    def maxSynthonAtoms(self) -> int:
        """
        If >0, sets a maximum number of heavy atoms, excluding dummies, for synthon to have a shape made.  Default=0.
        """

    @maxSynthonAtoms.setter
    def maxSynthonAtoms(self, arg: int, /) -> None: ...

    @property
    def maxEmbedAttempts(self) -> int:
        """
        Maximum number of attempts for embedding a single synthon.  Default=10.
        """

    @maxEmbedAttempts.setter
    def maxEmbedAttempts(self, arg: int, /) -> None: ...

    @property
    def timeOut(self) -> int:
        """
        Maximum time in seconds to spend on each synthon when generating conformers.  Default=600 means no timeout.
        """

    @timeOut.setter
    def timeOut(self, arg: int, /) -> None: ...

    @property
    def interimFile(self) -> str:
        """
        Interim file to write the SynthonSpace to during shape generation.  In the event of a failure, a restart from this file will be possible.
        """

    @interimFile.setter
    def interimFile(self, arg: str, /) -> None: ...

    @property
    def interimWrites(self) -> int:
        """
        If an interim file has been given, every this many shapes write a new version of the file.  Default=1000.
        """

    @interimWrites.setter
    def interimWrites(self, arg: int, /) -> None: ...

    def setUserConformerGenerator(self, func: Callable[[str, int], rdkit.Chem.rdchem.Mol | None]) -> None:
        """
        Allows you to provide a function that will be called instead of the default
         conformer generator to generate conformers for the synthons.  The function should
         take a SMILES string and the maximum number of conformers to generated and
         return a molecule object.
        """

    def __setattr__(self, name: str, value: object | None) -> None: ...

class SynthonSpace:
    """SynthonSpaceSearch object."""

    def __init__(self) -> None: ...

    def ReadTextFile(self, inFile: str | os.PathLike) -> None:
        """Reads text file of the sort used by ChemSpace/Enamine."""

    def ReadDBFile(self, inFile: str | os.PathLike, numThreads: int = 1) -> None:
        """
        Reads binary database file.  Takes optional number of threads,default=1.
        """

    def WriteDBFile(self, outFile: str | os.PathLike) -> None:
        """Writes binary database file."""

    def WriteEnumeratedFile(self, outFile: str | os.PathLike) -> None:
        """Writes enumerated library to file."""

    def GetNumReactions(self) -> int:
        """Returns number of reactions in the SynthonSpace."""

    def GetNumProducts(self) -> int:
        """
        Returns number of products in the SynthonSpace, with multiple
        counting of any duplicates.
        """

    def GetNumSynthons(self) -> int:
        """Returns number of synthons in the SynthonSpace."""

    def GetNumSynthonsWithShapes(self) -> int:
        """Returns the number of synthons in the SynthonSpace that have a shape."""

    def Summarise(self) -> None:
        """Writes a summary of the SynthonSpace to stdout."""

    def ReportSynthonUsage(self) -> None:
        """Writes a summary of the synthon usage in the SynthonSpace to stdout."""

    def GetSynthonFingerprintType(self) -> str:
        """
        Returns the information string for the fingerprint generator
        used to create this space.
        """

    @overload
    def SubstructureSearch(self, query: rdkit.Chem.rdchem.Mol, substructMatchParams: rdkit.Chem.rdchem.SubstructMatchParameters | None = None, params: SynthonSpaceSearchParams | None = None) -> SubstructureResult:
        """Does a substructure search in the SynthonSpace."""

    @overload
    def SubstructureSearch(self, query: rdkit.Chem.rdGeneralizedSubstruct.ExtendedQueryMol, substructMatchParams: rdkit.Chem.rdchem.SubstructMatchParameters | None = None, params: SynthonSpaceSearchParams | None = None) -> SubstructureResult:
        """
        Does a substructure search in the SynthonSpace using an
        extended query.
        """

    @overload
    def SubstructureSearch(self, query: rdkit.Chem.rdchem.Mol, substructMatchParams: rdkit.Chem.rdchem.SubstructMatchParameters | None, params: SynthonSpaceSearchParams | None, startLine: int, finishLine: int) -> SubstructureResult: ...

    @overload
    def SubstructureSearch(self, query: rdkit.Chem.rdchem.Mol, substructMatchParams: rdkit.Chem.rdchem.SubstructMatchParameters | None = None, params: SynthonSpaceSearchParams | None = None, *, startLine: int, finishLine: int) -> SubstructureResult: ...

    @overload
    def SubstructureSearch(self, query: rdkit.Chem.rdGeneralizedSubstruct.ExtendedQueryMol, substructMatchParams: rdkit.Chem.rdchem.SubstructMatchParameters | None, params: SynthonSpaceSearchParams | None, startLine: int, finishLine: int) -> SubstructureResult:
        """
        Take the contents of params.possibleHitsFile, which is assumed to have
        been written by an earlier search, and extract those that are indeed
        hits.  It makes sense that params is the same as the one used to
        generate the possible hits, but this is not essential.  You could search
        at a higher similarity threshold than used to create the possible hits,
        for example.
        """

    @overload
    def SubstructureSearch(self, query: rdkit.Chem.rdGeneralizedSubstruct.ExtendedQueryMol, substructMatchParams: rdkit.Chem.rdchem.SubstructMatchParameters | None = None, params: SynthonSpaceSearchParams | None = None, *, startLine: int, finishLine: int) -> SubstructureResult: ...

    def SubstructureSearchIncremental(self, query: rdkit.Chem.rdchem.Mol, callback: Callable[[list[rdkit.Chem.rdchem.Mol]], bool | None], substructMatchParams: rdkit.Chem.rdchem.SubstructMatchParameters | None = None, params: SynthonSpaceSearchParams | None = None) -> None:
        """
        Does a substructure search in the SynthonSpace returning results in the callback.
        """

    @overload
    def FingerprintSearch(self, query: rdkit.Chem.rdchem.Mol, fingerprintGenerator: rdkit.Chem.rdFingerprintGenerator.FingerprintGenerator64, params: SynthonSpaceSearchParams | None = None) -> SubstructureResult:
        """
        Does a fingerprint search in the SynthonSpace using the
        FingerprintGenerator passed in.
        """

    @overload
    def FingerprintSearch(self, query: rdkit.Chem.rdchem.Mol, fingerprintGenerator: rdkit.Chem.rdFingerprintGenerator.FingerprintGenerator64, params: SynthonSpaceSearchParams | None, startLine: int, finishLine: int) -> SubstructureResult:
        """
        Take the contents of params.possibleHitsFile, which is assumed to have
        been written by an earlier search, and extract those that are indeed
        hits.  It makes sense that params is the same as the one used to
        generate the possible hits, but this is not essential.  You could search
        at a higher similarity threshold than used to create the possible hits,
        for example.
        Duplicate SMILES strings produced by different reactions will
        be returned.
        """

    def FingerprintSearchIncremental(self, query: rdkit.Chem.rdchem.Mol, fingerprintGenerator: rdkit.Chem.rdFingerprintGenerator.FingerprintGenerator64, callback: Callable[[list[rdkit.Chem.rdchem.Mol]], bool | None], params: SynthonSpaceSearchParams | None = None) -> None:
        """
        Does a fingerprint search in the SynthonSpace using the
        FingerprintGenerator passed in, returning results the callback.
        """

    @overload
    def RascalSearch(self, query: rdkit.Chem.rdchem.Mol, rascalOptions: rdkit.Chem.rdRascalMCES.RascalOptions, params: SynthonSpaceSearchParams | None = None) -> SubstructureResult:
        """
        Does a search using the Rascal similarity score.  The similarity
        threshold used is provided by rascalOptions, and the one in
        params is ignored.
        """

    @overload
    def RascalSearch(self, query: rdkit.Chem.rdchem.Mol, rascalOptions: rdkit.Chem.rdRascalMCES.RascalOptions, params: SynthonSpaceSearchParams | None, startLine: int, finishLine: int) -> SubstructureResult:
        """
        Take the contents of params.possibleHitsFile, which is assumed to have
        been written by an earlier search, and extract those that are indeed
        hits.  It makes sense that params is the same as the one used to
        generate the possible hits, but this is not essential.  You could search
        at a higher similarity threshold than used to create the possible hits,
        for example.
        Duplicate SMILES strings produced by different reactions will
        be returned.
        """

    def RascalSearchIncremental(self, query: rdkit.Chem.rdchem.Mol, rascalOptions: rdkit.Chem.rdRascalMCES.RascalOptions | None, callback: Callable[[list[rdkit.Chem.rdchem.Mol]], bool | None], params: SynthonSpaceSearchParams | None = None) -> None:
        """
        Does a search using the Rascal similarity score.  The similarity
        threshold used is provided by rascalOptions, and the one in
        params is ignored.  Returns results iteratively in the callback.
        """

    @overload
    def ShapeSearch(self, query: rdkit.Chem.rdchem.Mol, params: SynthonSpaceSearchParams | None = None) -> SubstructureResult:
        """
        Perform a shape similarity search with the given query molecule
        across the synthonspace library.  Duplicate SMILES strings produced by
        different reactions will be returned.  Requires a query with at least
        1 3D conformer.  Only the first conformer will be used in the search.
        """

    @overload
    def ShapeSearch(self, query: rdkit.Chem.rdchem.Mol, params: SynthonSpaceSearchParams, startLine: int, finishLine: int) -> SubstructureResult:
        """
        Take the contents of params.possibleHitsFile, which is assumed to have
        been written by an earlier search, and extract those that are indeed
        hits.  It makes sense that params is the same as the one used to
        generate the possible hits, but this is not essential.  You could search
        at a higher similarity threshold than used to create the possible hits,
        for example.
        Duplicate SMILES strings produced by different reactions will
        be returned.  Requires a query with at least 1 3D conformer.  Only
        the first conformer will be used in the search.
        """

    def BuildSynthonFingerprints(self, fingerprintGenerator: rdkit.Chem.rdFingerprintGenerator.FingerprintGenerator64, progressBarWidth: int = 0) -> None:
        """
        Build the synthon fingerprints ready for similarity searching.  This
        is done automatically when the first similarity search is done, but if
        converting a text file to binary format it might need to be done
        explicitly.  If progressBarWidth is > 0, a progress bar of that width
        plus about 35 characters is displayed.
        """

    def BuildSynthonShapes(self, py_params: ShapeBuildParams | None = None) -> None:
        """
        Build shapes for the synthons.  The conformations are generated, pruned with the given threshold, which is passed directly to EmbedMultipleConfs.
        """

def ConvertTextToDBFile(inFilename: str | os.PathLike, outFilename: str | os.PathLike, fpGen: rdkit.Chem.rdFingerprintGenerator.FingerprintGenerator64 | None = None, py_shapeParams: ShapeBuildParams | None = None) -> None:
    """
    Convert the text file into the binary DB file in our format.
      Assumes that all synthons from a reaction are contiguous in the input file.
      This uses a lot less memory than using ReadTextFile() followed by
      WriteDBFile().
    - inFilename the name of the text file
    - outFilename the name of the binary file
    - optional fingerprint generator
    """

def FormattedIntegerString(value: int) -> str:
    """Format an integer with spaces every 3 digits for ease of reading"""
