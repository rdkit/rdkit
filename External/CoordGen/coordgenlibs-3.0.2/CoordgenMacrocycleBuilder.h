/*
 Contributors: Nicola Zonta
 Copyright Schrodinger, LLC. All rights reserved
 */

#ifndef COORDGEN_MACROCYCLE_BUILDER_H
#define COORDGEN_MACROCYCLE_BUILDER_H

#include "CoordgenConfig.hpp"
#include <cstdlib>
#include <iostream>
#include <vector>

class sketcherMinimizerAtom;
class sketcherMinimizerRing;
class sketcherMinimizerBond;
class sketcherMinimizerPointF;

struct doubleBondConstraint {
    bool trans;
    int previousAtom, atom1, atom2, followingAtom;
};

struct ringConstraint {
    ringConstraint(int a, sketcherMinimizerRing* r, bool fo)
    {
        ring = r;
        atom = a;
        forceOutside = fo;
    }
    bool forceOutside;
    int atom;
    sketcherMinimizerRing* ring;
};

struct vertexCoords {
    bool operator!=(const vertexCoords& rhs) const
    {
        if (x != rhs.x) {
            return true;
        }
        if (y != rhs.y) {
            return true;
        }
        if (z != rhs.z) {
            return true;
        }
        return false;
    }
    bool operator<(const vertexCoords& rhs) const
    {
        if (x < rhs.x) {
            return true;
        }
        if (y < rhs.y) {
            return true;
        }
        if (z < rhs.z) {
            return true;
        }
        return false;
    }
    friend const vertexCoords operator+(const vertexCoords& v1,
                                        const vertexCoords& v2)
    {
        return vertexCoords(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
    }
    friend const vertexCoords operator-(const vertexCoords& v1,
                                        const vertexCoords& v2)
    {
        return vertexCoords(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
    }

    bool operator==(const vertexCoords& rhs) const
    {
        return x == rhs.x && y == rhs.y && z == rhs.z;
    }

    vertexCoords(int ix, int iy, int iz)
    {
        x = ix;
        y = iy;
        z = iz;
    }
    int x, y, z;

  private:
    friend std::ostream& operator<<(std::ostream& os, const vertexCoords& v);
};

struct pathRestraints {
    std::vector<int> heteroAtoms;
    std::vector<std::pair<int, int>> substitutedAtoms;
};

struct pathConstraints {
    std::vector<doubleBondConstraint> doubleBonds;
    std::vector<ringConstraint> ringConstraints;
    std::vector<int> forceOutside;
};

struct hexCoords {
    hexCoords(int ix, int iy)
    {
        x = ix;
        y = iy;
    }
    bool operator==(const hexCoords& rhs) const
    {
        return x == rhs.x && y == rhs.y;
    }
    int x, y;
    int distanceFrom(const hexCoords& origin) const
    {
        int dx = abs(x - origin.x);
        int dy = abs(y - origin.y);
        int dz = abs((-x - y) - (-origin.x - origin.y));
        int max = (dx > dy) ? dx : dy;
        return (dz > max) ? dz : max;
    }
    hexCoords rotate30Degrees()
    {
        int z = -x - y;
        return {-z, -x};
    }
    vertexCoords toVertexCoords() const { return {x, y, -x - y}; }

  private:
    friend std::ostream& operator<<(std::ostream& os, const hexCoords& h);
};

/*
 unit hexagon used to build polyominoes for macrocycle shapes
 */
struct Hex {
    Hex(hexCoords coords) : m_coords(coords) {}
    void setCoords(hexCoords coords) { m_coords = coords; }
    int x() const { return m_coords.x; }
    int y() const { return m_coords.y; }
    int z() const { return -x() - y(); }
    hexCoords coords() const { return m_coords; }
    hexCoords m_coords;
    std::vector<hexCoords> neighbors() const;
    static std::vector<hexCoords> neighboringPositions(hexCoords h);
    vertexCoords followingVertex(vertexCoords v) const;
};

/*
hex polyomino (geometrical figure built on a hexagon lattice). All functions
assume that the polyomino has no holes
 */
class EXPORT_COORDGEN Polyomino
{
  public:
    Polyomino();
    Polyomino(const Polyomino& rhs);
    ~Polyomino();
    Polyomino& operator=(const Polyomino& rhs);

    /*
     explore the topology of the polyominoes and returns true if they have the
     same. Takes into account translations, rotations and mirroring
     */
    bool isTheSameAs(Polyomino& p) const;

    /* returns the number of hexagons in the polyomino */
    size_t size() const;

    /* empties the polyomino */
    void clear();

    /* marks one hexagon to be a pentagon */
    void
    markOneVertexAsPentagon(); // to get a path with an odd number of vertices

    /* returns all hexagons that share the given vertex */
    std::vector<Hex*> vertexNeighbors(vertexCoords v) const;

    /*
     returns neighboring positions that are not shared with other hexagons
     */
    std::vector<hexCoords> freeVertexNeighborPositions(vertexCoords v) const;

    /*
     return the contour path of the polyomino
     */
    std::vector<vertexCoords> getPath() const;

    /* return the hexagon at the given position */
    Hex* getHex(hexCoords coords) const;

    /* find an outer vertex. Used to start path around polyomino. */
    vertexCoords findOuterVertex() const;

    /* number of Hexs present at the given vertex */
    size_t hexagonsAtVertex(vertexCoords v) const;

    /* build a round-ish polyomino with the given number of vertices */
    void buildWithVerticesN(int totVertices);

    /* build a box-like polyomino of the given size */
    void buildSkewedBoxShape(int x, int y,
                             bool pentagon = false); // squared shape

    /*
     build a zig-zag box polyomino of the given size
     */
    void buildRaggedBoxShape(int x, int y, bool pentagon = false);

    /*
     build a zig-zag box polyomino of the given size
     */
    void buildRaggedSmallerBoxShape(
        int x, int y,
        bool pentagon = false); // box alternating rows of length x and x-1

    /*
     build a zig-zag box polyomino of the given size
     */
    void buildRaggedBiggerBoxShape(
        int x, int y,
        bool pentagon = false); // box alternating rows of length x and x+1

    /*
     return vector of coordinates of unoccupied hexagons bordering the polyomino
     */
    std::vector<hexCoords> allFreeNeighbors() const;

    /*
     return number of neighbors of hexagon at given coordinates
     */
    int countNeighbors(hexCoords) const;

    /* add an hexagon at given coordinates */
    void addHex(hexCoords coords);

    /* remove the hexagon at given coordinates */
    void removeHex(hexCoords coords);

    /* does removing this hexagon yield another polyomino with the same number
     of vertices? true if the hexagon has 3 neighbors all next to each other
     */
    bool isEquivalentWithout(hexCoords c) const;

    /* give the coordinates of an hypotetical substituent bound to an atom in
     * position pos */
    vertexCoords coordinatesOfSubstituent(vertexCoords pos) const;

    // holds pointers to all hexagons
    std::vector<Hex*> m_list;

    // holds coordinates to vertices marked as pentagonal
    std::vector<vertexCoords> pentagonVertices;

  private:
    /* mark vertex as pentagon (i.e. fuse it with neighbor vertex) */
    void setPentagon(vertexCoords c);

    /* resize the grid to new size i */
    void resizeGrid(int i) const;

    /* reassign hexagons to the grid */
    void reassignHexs() const;

    /* get index of hexagon at the given coordinates */
    int getIndexInList(hexCoords coords) const;

    /* hold pointers to all positions in the
    grid, NULL pointers for empty positions */
    std::vector<Hex*> mutable m_grid;

    /* size of the grid */
    mutable int m_gridSize;
};

/*
 class that builds coordinates for macrocycles
 */
class EXPORT_COORDGEN CoordgenMacrocycleBuilder
{
  public:
    CoordgenMacrocycleBuilder() = default;
    ~CoordgenMacrocycleBuilder() = default;

    /* assign coordinates to macrocycle */
    std::vector<sketcherMinimizerPointF>
    newMacrocycle(sketcherMinimizerRing* r,
                  std::vector<sketcherMinimizerAtom*> atoms) const;

    /* assign coordinates by removing one bond, generating the coordinates of
     the resulting molecule and re-adding the missing bond */
    bool openCycleAndGenerateCoords(sketcherMinimizerRing* ring) const;

    /* find most suitable bond to be broken */
    sketcherMinimizerBond* findBondToOpen(sketcherMinimizerRing* ring) const;

    /* Add constraints to stereoactive double bonds to avoid inversions
     public to be accessed by tests */
    std::vector<doubleBondConstraint>
    getDoubleBondConstraints(std::vector<sketcherMinimizerAtom*>& atoms) const;

    /* Skip the polyomino approach and fall back to opening the macrocycle when
     * generating coordinates */
    bool m_forceOpenMacrocycles = false;

    float getPrecision() const;

    void setPrecision(float f);

  private:
    /*
     precision of calculation. Higher values result in better results but longer
     calculation times
     */
    float m_precision;

    /* build polyominos with the given number of vertices */
    std::vector<Polyomino> buildSquaredShapes(int totVertices) const;

    /* remove duplicate polyominos from the list */
    std::vector<Polyomino> removeDuplicates(std::vector<Polyomino>& pols) const;

    std::vector<ringConstraint>
    getRingConstraints(std::vector<sketcherMinimizerAtom*>& atoms) const;

    /* get total number of atoms bound to a that are not on parent size */
    int getNumberOfChildren(sketcherMinimizerAtom* a,
                            sketcherMinimizerAtom* parent) const;

    /* get restraints to try to avoid things like having heteroatoms or
     * substituents pointing inside of the macrocycle */
    pathRestraints
    getPathRestraints(std::vector<sketcherMinimizerAtom*>& atoms) const;

    /* a shape is not allowed to break a constraint (i.e. invert
     stereochemistry on a double bond or put a fused ring inside
     the macrocycle) */
    pathConstraints
    getPathConstraints(std::vector<sketcherMinimizerAtom*>& atoms) const;

    /* build a list of polyominoes with the same number of vertices by removing
     * hexagons with 3 neighbors */
    std::vector<Polyomino> listOfEquivalent(const Polyomino& p) const;

    std::vector<Polyomino>
    listOfEquivalents(const std::vector<Polyomino>& l) const;

    /* check that no ring constraints are violated */
    bool checkRingConstraints(std::vector<ringConstraint>& ringConstraints,
                              Polyomino& p, std::vector<vertexCoords>& path,
                              std::vector<int>& neighborNs, int& startI) const;

    /* check that no double bond constraints are violated */
    bool checkDoubleBoundConstraints(
        std::vector<doubleBondConstraint>& dbConstraints,
        std::vector<vertexCoords>& vertices, int& startI) const;

    /* score the path restraints */
    int scorePathRestraints(pathRestraints& pr, Polyomino& p,
                            std::vector<vertexCoords>& path,
                            std::vector<int>& neighborNs, int& startI) const;

    /* check that no path constraints are violated */
    bool scorePathConstraints(pathConstraints& pc, Polyomino& p,
                              std::vector<vertexCoords>& path,
                              std::vector<int>& neighborNs, int& startI) const;

    /* score restraints and constraints of the given path */
    int scorePath(Polyomino& p, std::vector<vertexCoords>& path,
                  std::vector<int>& neighborNs, int& startI,
                  pathConstraints& pc, pathRestraints& pr) const;

    /* acceptable score to consider a shape appropriate */
    int acceptableShapeScore(int numberOfAtoms) const;
    std::vector<int> getVertexNeighborNs(Polyomino& p,
                                         std::vector<vertexCoords>& path) const;

    /* return lowest period after which there is a rotation symmetry */
    int getLowestPeriod(std::vector<int>& neighbors) const;

    /* explore the given polyominoes and score them */
    bool matchPolyominoes(std::vector<Polyomino>& pols, pathConstraints& pc,
                          pathRestraints& pr, int& bestP, int& bestScore,
                          int& bestStart, int& checkedMacrocycles) const;

    /* explore the best starting point to place atoms around the given polyomino
     */
    bool matchPolyomino(Polyomino& p, pathConstraints& pc, pathRestraints& pr,
                        int& bestStart, int& bestScore) const;

    /* write the coordinates of the given path */
    void
    writePolyominoCoordinates(std::vector<vertexCoords>& path,
                              const std::vector<sketcherMinimizerAtom*>& atoms,
                              int startI) const;

    /* return the coordinates of the given vertex */
    sketcherMinimizerPointF coordsOfVertex(vertexCoords& v) const;
};

#endif /* defined(COORDGEN_MACROCYCLE_BUILDER_H) */
