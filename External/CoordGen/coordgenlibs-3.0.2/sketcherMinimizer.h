#pragma once
/*
 *  sketcherMinimizer.h
 *
 *  Created by Nicola Zonta on 13/04/2010.
 *   Copyright Schrodinger, LLC. All rights reserved.
 *
 */
#include <map>
#include <vector>
#include <iostream>

#include "CoordgenConfig.hpp"

#include "CoordgenFragmentBuilder.h"
#include "CoordgenMinimizer.h"
#include "sketcherMinimizerFragment.h"
#include "sketcherMinimizerMarchingSquares.h"
#include "sketcherMinimizerMolecule.h"
#include "sketcherMinimizerResidue.h"
#include "sketcherMinimizerResidueInteraction.h"
#include "sketcherMinimizerRing.h"

static const float SKETCHER_STANDARD_PRECISION = 1.f;
static const float SKETCHER_QUICK_PRECISION = 0.2f;
static const float SKETCHER_BEST_PRECISION = 3.f;

class sketcherMinimimizerInteraction;

namespace schrodinger
{
namespace mae
{
class Block;
}
} // namespace schrodinger

typedef struct {
    std::vector<sketcherMinimizerPointF> additionVectors;
    std::vector<sketcherMinimizerPointF> centers;
    std::vector<int> counters;
} proximityData;

/* class to handle templates of common difficult ring structures */
class CoordgenTemplates
{
  public:
    CoordgenTemplates() = default;
    ~CoordgenTemplates()
    {
        for (auto molecule : m_templates) {
            for (auto atom : molecule->_atoms) {
                delete atom;
            }
            for (auto bond : molecule->_bonds) {
                delete bond;
            }
            delete molecule;
        }
        m_templates.clear();
    }
    std::vector<sketcherMinimizerMolecule*>& getTemplates()
    {
        return m_templates;
    }
    void setTemplateDir(std::string&& dir)
    {
        m_templateDir = std::move(dir);
        if (m_templateDir.back() != '/') {
            m_templateDir += "/";
        }
    }
    std::string getTemplateDir() { return m_templateDir; }

  private:
    std::vector<sketcherMinimizerMolecule*> m_templates;
    std::string m_templateDir = "";
};

/* main class. Creates 2d coordinates for molecular inputs */
class EXPORT_COORDGEN sketcherMinimizer
{

  public:
    sketcherMinimizer(float precision = SKETCHER_STANDARD_PRECISION);
    ~sketcherMinimizer();

    /* fields used in rdkit */
    std::vector<sketcherMinimizerFragment*> _fragments;
    CoordgenMinimizer m_minimizer;

    /* run coordinates generation and return true if the pose is considered
     * optimal */
    bool runGenerateCoordinates();

    /* initialize data from given molecule */
    void initialize(sketcherMinimizerMolecule* minMol);

    /* run a force-field based minimization on the given molecule */
    void minimizeMolecule(sketcherMinimizerMolecule* molecule);

    // void initializeFromMolecule(ChmMol& mol);
    void writeMinimizationData();

    /* run a force-field based minimization */
    void forceFieldMinimize();

    /* if mol contains separate molecules, split them into a vector */
    void splitIntoMolecules(sketcherMinimizerMolecule* mol,
                            std::vector<sketcherMinimizerMolecule*>& mols);

    /* split molecules into rigid fragments */
    void findFragments();

    /* constrain coordinates on all atoms */
    void constrainAllAtoms();

    /* constrain coordinates on atoms corresponding to true */
    void constrainAtoms(const std::vector<bool>& constrained);

    /* fix cooordinates (i.e. guarantee they will not change) on atoms marked as
     * true */
    void fixAtoms(const std::vector<bool>& fixed);

    /* set a flag to enable/disable the scoring of interactions with residues */
    void setScoreResidueInteractions(bool b);

    /* add bond to m_extrabonds */
    void addExtraBond(sketcherMinimizerBond* bond);

    /* getters */
    const std::vector<sketcherMinimizerAtom*>& getAtoms() const { return m_atoms; }
    const std::vector<sketcherMinimizerMolecule*>& getMolecules() const { return m_molecules; }
    const std::vector<sketcherMinimizerAtom*>& getReferenceAtoms() const { return m_referenceAtoms; }
    const std::vector<sketcherMinimizerBond*>& getReferenceBonds() const { return m_referenceBonds; }
    const std::vector<sketcherMinimizerFragment*>& getIndependentFragments() const { return m_independentFragments; }
    bool getTreatNonterminalBondsToMetalAsZOBs() const {return m_treatNonterminalBondsToMetalAsZOBs;}
    const std::vector<sketcherMinimizerBendInteraction*>& getBendInteractions() const;
    const std::vector<sketcherMinimizerStretchInteraction*>& getStretchInteractions() const;
    const std::vector<sketcherMinimizerInteraction*>& getInteractions() const;

    /* setters */
    void setFragments(std::vector<sketcherMinimizerFragment*> fragments) { _fragments = fragments; }
    void setTreatNonterminalBondsToMetalAsZOBs(bool b) {m_treatNonterminalBondsToMetalAsZOBs = b;}
    void setEvenAngles(bool b);
    void setSkipMinimization(bool b);
    void setForceOpenMacrocycles(bool b);

    /// Interactions with CoordgenMinimizer and CoordgenFragmentBuilder
    /* find a list of carbons from the backbone C=O of a peptide */
    std::set<sketcherMinimizerAtom*>
    getChetoCs(const std::vector<sketcherMinimizerAtom*>& allAtoms) const;

    /* find a list of nitrogens from the backbon NH of a peptide */
    std::set<sketcherMinimizerAtom*>
    getAminoNs(const std::vector<sketcherMinimizerAtom*>& allAtoms) const;

    /* find a list of alpha carbons of a peptide */
    std::set<sketcherMinimizerAtom*>
    getAlphaCs(const std::vector<sketcherMinimizerAtom*>& allAtoms,
               const std::set<sketcherMinimizerAtom*>& chetoCs,
               const std::set<sketcherMinimizerAtom*>& aminoNs) const;

    /* clear all the interactions loaded in the minimizer and free memory */
    void clearInteractions();

    /* setup all constraints of given molecule */
    void addInteractionsOfMolecule(sketcherMinimizerMolecule* molecule,
                                   bool intrafragmentClashes = false);

    void addExtraInteraction(sketcherMinimizerMolecule* molecule,
                             sketcherMinimizerInteraction* interaction);

    void buildFromFragments(bool b);
    bool avoidClashesOfMolecule(
        sketcherMinimizerMolecule* molecule,
        const std::vector<sketcherMinimizerInteraction*>& extraInteractions =
            std::vector<sketcherMinimizerInteraction*>());

    /// static members
    /* put atoms in a canonical order to reduce dependency from order in the
     * input vector */
    static void canonicalOrdering(sketcherMinimizerMolecule* minMol);

    /* pick one atom out of the vector. Arbitrary criteria such as atomic number
     * and connectivity are used */
    static sketcherMinimizerAtom*
    pickBestAtom(std::vector<sketcherMinimizerAtom*>& atoms);

    /* return a score of the alignment between direction and templat.first,
     * weight on the angle between the two and templat.second */
    static float
    testAlignment(const sketcherMinimizerPointF& direction,
                  const std::pair<sketcherMinimizerPointF, float>& templat);

    /* find the best alignment of a fragment to its parent and set invert in
     * case the fragment needs to be flipped */
    static sketcherMinimizerPointF scoreDirections(
        sketcherMinimizerFragment* fragment, float angle,
        const std::vector<std::pair<sketcherMinimizerPointF, float>>&
            directions,
        bool& invert);

    /* align the fragment to its parent */
    static void
    alignWithParentDirection(sketcherMinimizerFragment* f,
                             const sketcherMinimizerPointF& position,
                             float angle);

    /* align the fragment to its parent in the case of constrained coordinates
     */
    static bool
    alignWithParentDirectionConstrained(sketcherMinimizerFragment* fragment,
                                        const sketcherMinimizerPointF& position,
                                        float angle);

    /* align the fragment to its parent in the case of unconstrained coordinates
     */
    static bool
    alignWithParentDirectionUnconstrained(sketcherMinimizerFragment* fragment,
                                          float angle);

    /* get all bonds to a terminal atom */
    static std::vector<sketcherMinimizerBond*>
    getAllTerminalBonds(sketcherMinimizerFragment* fragment);

    /* return a list of vectors the given fragment can be aligned with and a
     * score of the importance of each */
    static std::vector<std::pair<sketcherMinimizerPointF, float>>
    findDirectionsToAlignWith(sketcherMinimizerFragment* fragment);

    /* if the three atoms share a ring, return it */
    static sketcherMinimizerRing* sameRing(const sketcherMinimizerAtom* at1,
                                           const sketcherMinimizerAtom* at2,
                                           const sketcherMinimizerAtom* at3);

    /* if the two atoms share a ring, return it */
    static sketcherMinimizerRing* sameRing(const sketcherMinimizerAtom* at1,
                                           const sketcherMinimizerAtom* at2);

    /* if the two atoms share a bond, return it */
    static sketcherMinimizerBond* getBond(const sketcherMinimizerAtom* a1,
                                          const sketcherMinimizerAtom* a2);

    /* calculate root mean square deviation between templates and points */
    static float RMSD(const std::vector<sketcherMinimizerPointF>& templates,
                      const std::vector<sketcherMinimizerPointF>& points);

    /* singular value decomposition for 2x2 matrices.
     used for 2D alignment. */
    static void svd(float* a, float* U, float* Sig, float* V);

    /* set m to a rotation matrix to align ref to points */
    static void
    alignmentMatrix(const std::vector<sketcherMinimizerPointF>& ref,
                    const std::vector<sketcherMinimizerPointF>& points,
                    float* m);

    static void
    checkIdentity(std::vector<unsigned int> solution, int newSol,
                  std::vector<bool>& matrix,
                  std::vector<sketcherMinimizerPointF>& templateCoordinates,
                  std::vector<std::vector<size_t>>& molBonds,
                  std::vector<std::vector<size_t>>& templateBonds,
                  std::vector<std::vector<size_t>>& molCisTransChains,
                  std::vector<bool>& molIsCis, size_t size, bool& found,
                  std::vector<unsigned int>& mapping);

    /* compare atoms and bonds to template and map which atom is which in case
     * of a positive match */
    static bool compare(const std::vector<sketcherMinimizerAtom*>& atoms,
                        const std::vector<sketcherMinimizerBond*>& bonds,
                        sketcherMinimizerMolecule* templ,
                        std::vector<unsigned int>& mapping);

    /* calculate morgan scores for the given input */
    static int morganScores(const std::vector<sketcherMinimizerAtom*>& atoms,
                            const std::vector<sketcherMinimizerBond*>& bonds,
                            std::vector<int>& scores);

    /* load the templates from the template file */
    static void setTemplateFileDir(std::string dir);
    static void loadTemplates();
    static CoordgenTemplates m_templates;

private:
    /*all non-terminal bonds to a metal atom are treated as if they were zero order bonds (this usually results
     in a longer bond*/
    bool m_treatNonterminalBondsToMetalAsZOBs = true;

    CoordgenFragmentBuilder m_fragmentBuilder;

    std::vector<sketcherMinimizerAtom*> m_atoms;
    std::vector<sketcherMinimizerAtom*> m_referenceAtoms;
    std::vector<sketcherMinimizerResidue*> m_residues;
    std::vector<sketcherMinimizerResidueInteraction*> m_residueInteractions;

    std::vector<sketcherMinimizerFragment*> m_independentFragments;

    std::vector<sketcherMinimizerBond*> m_bonds;
    std::vector<sketcherMinimizerBond*> m_referenceBonds;
    std::vector<sketcherMinimizerBond*> m_proximityRelations;
    std::vector<sketcherMinimizerBond*> m_extraBonds;
    std::vector<sketcherMinimizerMolecule*> m_molecules;

    /* return the position of res, which is part of SSE, given that the first
     * residue of SSE is placed at startF and consecutive residues are placed
     * increment away from each other. All distances are expressed in floats,
     * where 0.f is an arbitrary starting point, 0.5 is the opposite side of the
     * curve and 1.0 is again the starting point */
    float getResidueDistance(float startF, float increment,
                             sketcherMinimizerResidue* res,
                             const std::vector<sketcherMinimizerResidue*>& SSE) const;

    /* return the vector index corresponding to floatPosition */
    int getShapeIndex(const std::vector<sketcherMinimizerPointF>& shape,
                      float floatPosition) const;

    /* return true if the molecules structure is reasonable (e.g. reasonable
     * amount of fused rings) */
    bool structurePassSanityCheck() const;

    /* clear data and free memory */
    void clear();

    /* flag atoms that will be drawn with 90° angles (e.g. phosphate P) */
    void flagCrossAtoms();

    /* assign coordinates to all molecules and residues */
    void minimizeAll();

    /* find the best angle to rotate each molecule */
    void bestRotation();

    /* add info to choose the best angle so that, if present, peptide chains are
     * horizontal */
    void addBestRotationInfoForPeptides(
        std::vector<std::pair<float, float>>& angles,
        const std::vector<sketcherMinimizerAtom*>& atoms);

    /* if a peptide chain is present make sure that the N term is on the left */
    void maybeFlipPeptides(const std::vector<sketcherMinimizerAtom*>& atoms,
                           float& scoreX);

    /* consider flipping the molecule horizontally and/or vertically */
    void maybeFlip();

    /* mark atoms with wedges as above or below the plane to correctly draw
     * crossing bonds */
    void assignPseudoZ();

    /* write wedges and dashed bonds to mark stereochemistry */
    void writeStereoChemistry();

    /* arrange multiple molecules next to each other */
    void arrangeMultipleMolecules();

    /* arrange molecules that have parts that interact with each other so that
     * they are close */
    void placeMoleculesWithProximityRelations(
        std::vector<sketcherMinimizerMolecule*> proximityMols);

    /* place molecules so that one is in the middle and other are around */
    void placeMolResidueLigandStyle(sketcherMinimizerMolecule* mol,
                                    sketcherMinimizerMolecule* parent);

    /* if the molecule has more than one interaction and they cross mirror its
     coordinates so they don't cross anymore */
    void flipIfCrossingInteractions(sketcherMinimizerMolecule* mol);

    /* build data vectors to place molecules with proximity relations */
    std::vector<proximityData> buildProximityDataVector(
        std::vector<sketcherMinimizerMolecule*>& proximityMols,
        std::map<sketcherMinimizerMolecule*, sketcherMinimizerAtom*>& molMap);

    /* translate molecules with proximity relations */
    void translateMoleculesWithProximityRelations(
        std::vector<sketcherMinimizerMolecule*>& proximityMols,
        std::map<sketcherMinimizerMolecule*, sketcherMinimizerAtom*>& molMap,
        std::map<sketcherMinimizerMolecule*, sketcherMinimizerPointF>&
            templateCenters,
        std::vector<proximityData>& proximityDataVecto);

    /* rotate molecules with proximity relations */
    void rotateMoleculesWithProximityRelations(
        std::vector<sketcherMinimizerMolecule*>& proximityMols,
        std::map<sketcherMinimizerMolecule*, sketcherMinimizerAtom*>& molMap,
        std::vector<proximityData>& proximityDataVector);

    /* place residues in concentric contours around the ligand */
    void placeResiduesInCrowns();

    /* place residues from the given strands of consecutive residues to fill
     the given path */
    bool fillShape(std::vector<std::vector<sketcherMinimizerResidue*>>& SSEs,
                   const std::vector<sketcherMinimizerPointF>& shape,
                   int shapeN);

    /* place a single strand of consecutive residues */
    void placeSSE(const std::vector<sketcherMinimizerResidue*>& SSE,
                  const std::vector<sketcherMinimizerPointF>& shape, int shapeN,
                  std::vector<bool>& penalties,
                  std::set<sketcherMinimizerResidue*>& outliers,
                  bool placeOnlyInteracting = false);

    /* score the position of the given strands */
    float scoreSSEPosition(const std::vector<sketcherMinimizerResidue*>& SSE,
                           const std::vector<sketcherMinimizerPointF>& shape,
                           int shapeN, std::vector<bool>& penalties, float startingPosition,
                           float increment);

    /* score the distance between the two given points of connected residues */
    float scoreSSEBondStretch(const sketcherMinimizerPointF& coordinates1,
                              const sketcherMinimizerPointF& coordinates2);

    /* solution represent the placement chosen for residues in SSE. Mark the
     * corresponding sections of the crown to prevent other residues to be
     * placed
     * there */
    void markSolution(const std::pair<float, float>& solution,
                      const std::vector<sketcherMinimizerResidue*>& SSE,
                      const std::vector<sketcherMinimizerPointF>& shape,
                      std::vector<bool>& penalties,
                      std::set<sketcherMinimizerResidue*>& outliers);

    /* return a concentric shape around the ligand. CrownN controls how far away
     from the ligand the shape is */
    std::vector<sketcherMinimizerPointF> shapeAroundLigand(int crownN);

    /* group residues in strands of consecutive residues */
    std::vector<std::vector<sketcherMinimizerResidue*>>
    groupResiduesInSSEs(const std::vector<sketcherMinimizerResidue*>& residues);

    /* score the position of given residues */
    float
    scoreResiduePosition(int index,
                         const std::vector<sketcherMinimizerPointF>& shape,
                         int shapeN, std::vector<bool>& penalties,
                         sketcherMinimizerResidue* residue);

    /* assign coordinates to residues */
    void placeResidues(const std::vector<sketcherMinimizerAtom*>& atoms =
                           std::vector<sketcherMinimizerAtom*>(0));

    /* assign coordinates to residues in the context of a protein-protein
     interaction diagram */
    void placeResiduesProteinOnlyMode();

    /* assign coordinates to residues in a protein-protein interaction
     diagram shaped as a circle */
    void placeResiduesProteinOnlyModeCircleStyle(
        const std::map<std::string, std::vector<sketcherMinimizerResidue*>>&
            chains);

    /* assign coordinates to residues in a protein-protein interaction
     diagram shaped as a LID */
    void placeResiduesProteinOnlyModeLIDStyle(
        const std::map<std::string, std::vector<sketcherMinimizerResidue*>>&
            chains);

    /* order residues for drawing, so that residues interacting together are
     * drawn one after the other and residues with more interactions are drawn
     * first */
    std::vector<sketcherMinimizerResidue*> orderResiduesOfChains(
        const std::map<std::string, std::vector<sketcherMinimizerResidue*>>&
            chains);

    /* find center position for each chain of residues using a meta-molecule
     * approach, building a molecule where each atom represents a chain and each
     * bond connects two interacting chains */
    std::map<std::string, sketcherMinimizerPointF>
    computeChainsStartingPositionsMetaMol(
        const std::map<std::string, std::vector<sketcherMinimizerResidue*>>&
            chains);

    /* place interacting residues closer to each other, so they end up at the
     * periphery of the chain */
    void shortenInteractions(
        const std::map<std::string, std::vector<sketcherMinimizerResidue*>>&
            chains);

    /* explore positions in a grid around the current one to ease clashes */
    sketcherMinimizerPointF exploreGridAround(
        const sketcherMinimizerPointF& centerOfGrid, unsigned int levels,
        float gridD, float dx = 0.f, float dy = 0.f,
        float distanceFromAtoms = -1.f, bool watermap = false,
        sketcherMinimizerResidue* residueForInteractions = nullptr,
        const sketcherMinimizerPointF& direction = sketcherMinimizerPointF(0,
                                                                           1));

    sketcherMinimizerPointF exploreMolPosition(sketcherMinimizerMolecule* mol,
                                               unsigned int levels, float gridD,
                                               float distanceFromAtoms = -1.f);

    /* add given angle to a vector of angles, clustering them if a close angle
     * is already in the vector */
    void addToVector(float weight, float angle,
                     std::vector<std::pair<float, float>>& angles);

    void assignLongestChainFromHere(sketcherMinimizerFragment* f);
    void assignNumberOfChildrenAtomsFromHere(sketcherMinimizerFragment* f);

    /* for each residue, find the closest atom among atoms, or all atoms
     if none are given */
    void
    findClosestAtomToResidues(const std::vector<sketcherMinimizerAtom*>& atoms =
                                  std::vector<sketcherMinimizerAtom*>(0));

    /* initialize data and coordinates for each fragment */
    void initializeFragments();

    /* for tracking coordinates */
    float sin_flip = 0.f;
    float cos_flip = 0.f;
    float centerX = 0.f;
    float centerY = 0.f;
    int flipX = 1;
    int flipY = 1;
};
