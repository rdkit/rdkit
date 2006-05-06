# $Id$
#
# Copyright (C) 2005-2006 Rational Discovery LLC
#  All Rights Reserved
#

# name of the application (used as the window caption)
applicationName = 'RDPharm3D'
#applicationName = 'CCQuickAlign'

# value (in Angstroms) for the 3D distance slop
distance3DSlop=0.3
# multiplier (as a fraction) for the topological distance
# slop:
distance2DSlopFactor=0.5

# neighborhood radius in PyMol (in Angstroms):
neighborhoodRadius=5.0

# radius of features in 2D sketches:
featHighlightRadius=0.75

# Embeddings with a relative energy larger than this value
# will not be active (i.e. the user won't be able to select
# or display them):
embeddingEnergyTol=50.0

# name of the protein object in PyMol
pymolProteinName='Protein'
# name of the template object in PyMol
pymolTemplateName='Template'
# name used for the reference pharmacophore in PyMol:
refPcophoreName="Ref_Pharmacophore"
# name of the h-bond object in PyMol
pymolHBondName='hbnd'
# name of the collision object in PyMol
pymolCollisionName='close_contacts'
# distance cutoff for things to be considered collisions:
pymolCollisionCutoff=3.0
# text used to select the molecule when displaying h bonds and collisions:
pymolMolSelText='(%(molName)s)'
# text used to select the protein when displaying h bonds and collisions:
pymolProteinSelText='(%(proteinName)s and not het)'



# message displayed in results tree when a molecule
# doesn't match every feature:
resultsMissingFeatsMessage='Some features did not match'
# message displayed in results tree when a molecule
# doesn't have a valid mapping:
resultsNoMappingMessage='No compatible mappings found'
# message displayed in results tree when a mapping
# doesn't generate embeddings:
resultsNoEmbeddingMessage='No valid conformations were found'


# used to construct the names of un-named molecules
# in the results tree:
resultsMoleculeNameBase='Mol'

# if this is true, we'll actually do the pharmacophore search
# on the results from doing a grab from CDX (instead of just adding
# the molecules:
searchOnCDXGrab=False

# if this is false, we'll disable lots of functionality
allowExpertMode=True

# maximum number of points the non-expert user can add to the pre-defined
# pharmacophore:
maxAddablePoints=2

# name of the file used to collect log messages (for debugging):
logFilename='%s_log.txt'%applicationName

# if this is set, you'll always be able to enable/disable display
# of h bonds between the ligand and the protein (assuming that
# there is a protein)
hBondsAlwaysEnabled=True

# This is the default number of aligments the program tries to
# generate when a mapping item is expanded. It's also the number of
# alignments generated when "Generate More Alignments" is selected.
numAlignmentsDesired=5

# This is the maximum number of random distance matrices to generate
# while trying to generate the alignments:
maxAlignmentsToTry=40

# this is the maximum distance allowed between atoms for the
# automatic mol-mol correspondence algorithm to consider them
# as possibilities
maxAtomAtomCorrespondenceDistance=1.5

# should we forbid the automatic correspondence algorithm to
# suggest duplicate correspondences? (i.e. one reference atom
# mapping to multiple probe atoms)
noDupeCorrespondences=True

# force constant for the springs used to connect atoms when
# refining Alignments
refineForceConstant=25

# name used for refined alignments in pymol:
refinedAlignmentName='Refined_Alignment'

pymolCorrespondenceName='correspondence'
