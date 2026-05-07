/*
 Contributors: Nicola Zonta
 Copyright Schrodinger, LLC. All rights reserved
 */

#ifndef COORDGEN_MINIMIZER_H
#define COORDGEN_MINIMIZER_H

#include "CoordgenConfig.hpp"
#include <iostream>
#include <map>
#include <set>
#include <vector>

class sketcherMinimizerInteraction;
class sketcherMinimizerStretchInteraction;
class sketcherMinimizerBendInteraction;
class sketcherMinimizerClashInteraction;

class sketcherMinimizerMolecule;
class sketcherMinimizerResidue;
class sketcherMinimizerResidueInteraction;
class sketcherMinimizerAtom;
class sketcherMinimizerBond;
class sketcherMinimizerRing;
class sketcherMinimizerFragment;
class sketcherMinimizerPointF;

class CoordgenFragmentDOF;
class CoordgenMinimizer;

/* class to manage the solutions to a minimization problem. Each solution is a
 vector that contains values for each of the degrees of freedom of the problem
 */
class CoordgenDOFSolutions
{
  public:
    CoordgenDOFSolutions(const CoordgenMinimizer* minimizer,
                         sketcherMinimizerMolecule* molecule,
                         std::vector<CoordgenFragmentDOF*> allDofs)
        : m_minimizer(minimizer), m_molecule(molecule),
          m_allDofs(std::move(allDofs))
    {
    }
    /*
     calculate the value of the scoring function on the currently loaded
     solution
     */
    float scoreCurrentSolution();

    /*
     get the solution that is currently loaded
     */
    std::vector<short unsigned int> getCurrentSolution();

    /*
     load the given solution (i.e. set each degree of freedom in the molecule to
     the  given value)
     */
    void loadSolution(const std::vector<short unsigned int>& solution);

    /*
     return the best scoring solution that has been found so far
     */
    std::pair<std::vector<short unsigned int>, float> findBestSolution() const;

    /*
     check if the given solution has already been scored
     */
    bool hasSolution(const std::vector<short unsigned int>& solution);

    std::vector<CoordgenFragmentDOF*>& getAllDofs() { return m_allDofs; }

  private:
    const CoordgenMinimizer* m_minimizer;
    sketcherMinimizerMolecule* m_molecule;
    std::map<std::vector<short unsigned int>, float> m_solutions;
    std::vector<CoordgenFragmentDOF*> m_allDofs;
};

/*
 minimizer class that resolves clashes in a molecule. It can explore degrees of
 freedom conformations (e.g. flip around single bonds) and perform a force-field
 based free atom minimization
 */
class EXPORT_COORDGEN CoordgenMinimizer
{
  public:
    CoordgenMinimizer();

    ~CoordgenMinimizer();

    /* clear all the interactions loaded in the minimizer and free memory */
    void clearInteractions();

    /* run a force-field based minimization */
    void run();

    /* Apply forces and take a step in the minimization. Returns false if
     * converged, true if not. */
    bool applyForces(float maxd = 3);

    /* run a force-field based minimization on residues only */
    void minimizeResidues();

    /* run a force-field based minimization on the given molecule */
    void minimizeMolecule(sketcherMinimizerMolecule* molecule);

    /* solve clashes of residues in a protein-protein interaction LID */
    void minimizeProteinOnlyLID(
        const std::map<std::string, std::vector<sketcherMinimizerResidue*>>&
            chains);

    /* setup constraints and run a force-field based minimization */
    void minimizeAll();

    /* setup constraints on residues only */
    void setupInteractionsOnlyResidues();

    /* setup constraints on residues in a protein-protein interaction scenario
     */
    void setupInteractionsProteinOnly(
        const std::map<std::string, std::vector<sketcherMinimizerResidue*>>&
            chains);

    /* setup all constraints */
    void setupInteractions(bool intrafragmentClashes = false);

    /* setup all constraints of given molecule */
    void addInteractionsOfMolecule(sketcherMinimizerMolecule* molecule,
                                   bool intrafragmentClashes = false);

    /* score the forcefield value of the current conformation */
    float scoreInteractions();

    /* add a list of intermolecular clash constraints between two given
     * molecules
     */
    bool findIntermolecularClashes(sketcherMinimizerMolecule* mol1,
                                   sketcherMinimizerMolecule* mol2,
                                   float threshold);

    /* add a list of intermolecular clash constraints between given molecules */
    bool findIntermolecularClashes(
        const std::vector<sketcherMinimizerMolecule*>& mols, float threshold);


    /* precision of the minimization. Higher values result in higher times and
     * better results */
    float getPrecision() const;
    void setPrecision(float f);

    /* score all clashes of the given molecule */
    float
    scoreClashes(sketcherMinimizerMolecule* molecule,
                 bool residueInteractions = false,
                 bool scoreProximityRelationsOnOppositeSides = true) const;

    /* score the penalty for intersecting bonds */
    float scoreCrossBonds(sketcherMinimizerMolecule* molecule,
                          bool residueInteractions = false) const;

    /* score the penalty for degrees of freedom set to non-ideal values */
    float scoreDofs(sketcherMinimizerMolecule* molecule) const;

    /* score the penalty of atoms placed inside rings */
    float scoreAtomsInsideRings() const;

    /* assign a penalty if a molecule (A) is bound to another (B) with more than
     * a proximity relation and they are on different sides of A and involving
     * different fragments. This forces the algorithm to look for poses where
     * all the atoms of A that have proximity relations with B are on the same
     * side of A. */
    float scoreProximityRelationsOnOppositeSides() const;

    /*
     run the minimization and return true if the result is considered ideal
     */
    bool avoidClashes();

    /*
     run the minimization on the given molecule. Returns true if the result is
     considered ideal
     */
    bool avoidClashesOfMolecule(
        sketcherMinimizerMolecule* molecule,
        const std::vector<sketcherMinimizerInteraction*>& extraInteractions =
            std::vector<sketcherMinimizerInteraction*>());

    /*
     set up the DoF list and run a minimization on the given molecule
     */
    bool flipFragments(sketcherMinimizerMolecule* molecule, float& clashE);

    /*
     run a search on the degrees of freedom, exploring combinations of degrees
     of freedom
     */
    bool runLocalSearch(sketcherMinimizerMolecule* molecule,
                        const std::vector<CoordgenFragmentDOF*>& dofs,
                        int levels, float& clashE,
                        CoordgenDOFSolutions& solutions);

    /*
     iteratively grow the pool of solutions by mutating the best scoring one by
     one degree of freedom
     */
    bool growSolutions(
        std::set<std::vector<short unsigned int>>& allScoredSolutions,
        int& currentTier,
        std::map<std::vector<short unsigned int>, float>& growingSolutions,
        CoordgenDOFSolutions& solutions, float& bestScore);

    /*
     run the search to find good scoring solutions to the problem. Each degree
     of freedom has a tier value that ensures that simpler and more aesthetic
     pleasing ones (e.g. flipping fragments) are searched before more complex
     ones (i.e. putting substituents inside macrocycles). This is alternative
     and preferred to runExhaustiveSearch and runLocalSearch
     */
    bool runSearch(int tier, CoordgenDOFSolutions& solutions);

    /* build a list of tuples of the given order representing all combinations
     * of dofs */
    std::vector<std::vector<CoordgenFragmentDOF*>>
    buildTuplesOfDofs(const std::vector<CoordgenFragmentDOF*>& dofs,
                      unsigned int order) const;

    /*
     run en exhaustive search on all combinations of DoFs.
     */
    bool runExhaustiveSearch(sketcherMinimizerMolecule* molecule,
                             std::vector<CoordgenFragmentDOF*> dofs,
                             float& clashE, CoordgenDOFSolutions& solutions);

    void runExhaustiveSearchLevel(
        sketcherMinimizerMolecule* molecule,
        const std::vector<CoordgenFragmentDOF*>::iterator& iterator,
        std::vector<CoordgenFragmentDOF*>& dofs, float& bestResult, bool& abort,
        CoordgenDOFSolutions& solutions);

    /* return true if the given bonds clash */
    bool bondsClash(sketcherMinimizerBond* bond,
                    sketcherMinimizerBond* bond2) const;

    /*
     quick function to avoid clashes of terminal atoms without running a
     minimization
     */
    void avoidTerminalClashes(sketcherMinimizerMolecule* molecule,
                              float& clashE);

    /*
     check if the ring system cannot be drown with regular polygons and needs a
     FF based minimization
     */
    static void
    maybeMinimizeRings(const std::vector<sketcherMinimizerRing*>& rings);

    /*
     avoid clashes of terminal atoms of the same fragment without running a
     minimization
     */
    static void avoidInternalClashes(sketcherMinimizerFragment* fragment);

    /*
     assign coordinates to each atom from the current values of DoFs of the
     fragments
     */
    void buildFromFragments(bool firstTime = false) const;

    /*
     assign coordinates to each atom from the current values of DoFs of the
     fragments
     */
    void buildMoleculeFromFragments(sketcherMinimizerMolecule* molecule,
                                    bool firstTime = false) const;

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

    /* check the atom for clashes with other atoms */
    static void checkForClashes(sketcherMinimizerAtom* a);

    /*
     return true if atoms have NaN coordinates
     */
    static bool
    hasNaNCoordinates(const std::vector<sketcherMinimizerAtom*>& atoms);
    bool hasNaNCoordinates();

    /*
     return true if the atom has valid 3d coordinates
     */
    static bool
    hasValid3DCoordinates(const std::vector<sketcherMinimizerAtom*>& atoms);

    /* use 3d coordinates in 2d (e.g. when a reasonable 2d structure cannot be
     * found) */
    static void
    fallbackOn3DCoordinates(const std::vector<sketcherMinimizerAtom*>& atoms);

    /*
     add the given constraint to the minimizer
     */
    void addExtraInteraction(sketcherMinimizerMolecule* molecule,
                             sketcherMinimizerInteraction* interaction);

    std::vector<sketcherMinimizerAtom*> m_atoms;
    std::vector<sketcherMinimizerBond*> m_bonds;
    bool m_evenAngles;
    std::vector<sketcherMinimizerResidue*> m_residues;
    std::vector<sketcherMinimizerResidueInteraction*> m_residueInteractions;
    std::vector<sketcherMinimizerFragment*> m_fragments;
    std::vector<sketcherMinimizerMolecule*> m_molecules;
    std::vector<float> energy_list;
    std::vector<std::vector<sketcherMinimizerPointF>> all_coordinates;

    bool skipMinimization, skipAvoidClashes, skipFlipFragments,
        m_scoreResidueInteractions;

    const std::vector<sketcherMinimizerBendInteraction*>& getBendInteractions() const {return _bendInteractions;};
    const std::vector<sketcherMinimizerStretchInteraction*>& getStretchInteractions() const {return _stretchInteractions;};
    const std::vector<sketcherMinimizerInteraction*>& getInteractions() const {return _interactions;};


    /*
     add clash constraints of the given molecule
     */
    void addClashInteractionsOfMolecule(sketcherMinimizerMolecule* molecule,
                                        bool intrafragmentClashes);

    /*
     add stretch constraints of the given molecule
     */
    void addStretchInteractionsOfMolecule(sketcherMinimizerMolecule* molecule);

    /*
     add angle bend constraints of the given molecule
     */
    void addBendInteractionsOfMolecule(sketcherMinimizerMolecule* molecule);

    /*
     add constraints to avoid deviating from constrained coordinates (e.g. for
     alignment)
    */
    void
    addConstrainedInteractionsOfMolecule(sketcherMinimizerMolecule* molecule);

    /*
     add constraints to avoid chiral inversion of the given molecule
     */
    void addChiralInversionConstraintsOfMolecule(
        sketcherMinimizerMolecule* molecule);

    /*
     add constraints to keep peptide chains linear
     */
    void addPeptideBondInversionConstraintsOfMolecule(
        sketcherMinimizerMolecule* molecule);

private:
    /*
     get lists of four atoms that form a chain and are each present in one of
     the four sets respectively
     */
    void getFourConsecutiveAtomsThatMatchSequence(
        std::vector<std::vector<sketcherMinimizerAtom*>>&
            consecutiveAtomsGroups,
        const std::set<sketcherMinimizerAtom*>& firstSet,
        const std::set<sketcherMinimizerAtom*>& secondSet,
        const std::set<sketcherMinimizerAtom*>& thirdSet,
        const std::set<sketcherMinimizerAtom*>& fourthSet) const;

    std::vector<sketcherMinimizerInteraction*> _interactions;
    std::vector<sketcherMinimizerStretchInteraction*> _stretchInteractions;
    std::vector<sketcherMinimizerBendInteraction*> _bendInteractions;

    std::vector<sketcherMinimizerInteraction*> _intramolecularClashInteractions;
    std::vector<sketcherMinimizerInteraction*> _extraInteractions;
    std::map<sketcherMinimizerMolecule*,
             std::vector<sketcherMinimizerInteraction*>>
        _extraInteractionsOfMolecule;

    unsigned int m_maxIterations;
    float m_precision;
};

#endif /* defined(COORDGEN_MINIMIZER_H) */
