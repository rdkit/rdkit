/*
 Contributors: Nicola Zonta
 Copyright Schrodinger, LLC. All rights reserved
 */
#include <algorithm>
#include <queue>

#include "CoordgenMacrocycleBuilder.h"
#include "CoordgenMinimizer.h"
#include "sketcherMinimizer.h"
#include "sketcherMinimizerMaths.h"
#include "sketcherMinimizerRing.h"
#include "sketcherMinimizerStretchInteraction.h"


#define MAX_MACROCYCLES 40
#define PATH_FAILED -1000
#define HETEROATOM_RESTRAINT 1
#define SUBSTITUTED_ATOM_RESTRAINT 10
#define SUBSTITUENTS_TO_SAME_VERTEX_RESTRAINT 200
#define SUBSTITUENT_CLASH_WITH_PATH_RESTRAINT 400
#define SQRT3HALF 0.8660254037844386

using namespace std;

std::ostream& operator<<(std::ostream& os, const vertexCoords& v)
{
    os << "(" << v.x << "," << v.y << "," << v.z << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, const hexCoords& h)
{
    os << "(" << h.x << "," << h.y << ")";
    return os;
}

Polyomino::Polyomino()
{
    resizeGrid(1);
}

Polyomino::Polyomino(const Polyomino& rhs)
{

    clear();
    pentagonVertices = rhs.pentagonVertices;
    resizeGrid(1);
    for (auto i : rhs.m_list) {
        addHex(i->coords());
    }
    reassignHexs();
}

Polyomino& Polyomino::operator=(const Polyomino& rhs)
{
    clear();
    resizeGrid(1);
    pentagonVertices = rhs.pentagonVertices;
    for (auto i : rhs.m_list) {
        addHex(i->coords());
    }
    reassignHexs();
    return *this;
}

Polyomino::~Polyomino()
{
    clear();
}

void Polyomino::clear()
{
    for (auto& i : m_list) {
        delete i;
    }
    m_list.clear();
}

size_t Polyomino::size() const
{
    return m_list.size();
}

void Polyomino::resizeGrid(int i) const
{
    m_grid.resize((2 * i + 1) *
                  (2 * i + 1)); // grid can hold coordinates from -i to i
    m_gridSize = i;
    reassignHexs();
}

void Polyomino::reassignHexs() const
{
    for (auto& j : m_grid) {
        j = nullptr;
    }
    for (auto hex : m_list) {
        m_grid[getIndexInList(hex->coords())] = hex;
    }
}

bool Polyomino::isTheSameAs(Polyomino& p) const
{
    // mirrored polyominoes are considered to be different
    if (size() != p.size()) {
        return false;
    }
    vector<hexCoords> targetCoords;
    for (Hex* hex : p.m_list) {
        targetCoords.push_back(hex->coords());
    }
    if (targetCoords.empty()) {
        return true; // both polyominoes are empty
    }
    int lowestx = m_list[0]->coords().x;
    int lowesty = m_list[0]->coords().y;
    for (Hex* hex : m_list) {
        int x = hex->coords().x;
        int y = hex->coords().y;
        if (x < lowestx) {
            lowestx = x;
        }
        if (y < lowesty) {
            lowesty = y;
        }
    }
    for (unsigned int i = 0; i < 6;
         i++) { // explore 6 possible positions rotating the coordinates 30
                // degrees
        int lowestTargetX = 0, lowestTargetY = 0;
        for (unsigned int j = 0; j < targetCoords.size(); j++) {
            if (j == 0 || targetCoords[j].x < lowestTargetX) {
                lowestTargetX = targetCoords[j].x;
            }
            if (j == 0 || targetCoords[j].y < lowestTargetY) {
                lowestTargetY = targetCoords[j].y;
            }
        }
        // translate
        for (auto& targetCoord : targetCoords) {
            targetCoord = hexCoords(targetCoord.x + lowestx - lowestTargetX,
                                    targetCoord.y + lowesty - lowestTargetY);
        }
        // check

        bool same = true;
        for (auto targetCoord : targetCoords) {
            if (!getHex(targetCoord)) {
                same = false;
                break;
            }
        }
        if (same) {
            return true;
        }
        // rotate
        for (auto& targetCoord : targetCoords) {
            targetCoord = targetCoord.rotate30Degrees();
        }
    }
    return false;
}

Hex* Polyomino::getHex(hexCoords coords) const
{
    return m_grid[getIndexInList(coords)];
}

void Polyomino::buildSkewedBoxShape(int x, int y, bool pentagon)
{
    clear();
    for (int yy = 0; yy < y; yy++) {
        for (int xx = 0; xx < x; xx++) {
            addHex(hexCoords(xx, yy));
        }
    }
    if (pentagon) {
        markOneVertexAsPentagon();
    }
}

void Polyomino::buildRaggedSmallerBoxShape(
    int x, int y,
    bool pentagon) // box alternating rows of length x and x-1
{
    clear();
    int startx = 0;
    for (int yy = 0; yy < y; yy++) {
        // bigger row
        for (int xx = 0; xx < x; xx++) {
            addHex(hexCoords(startx + xx, yy));
        }
        yy++;
        if (yy >= y) {
            break;
        }
        for (int xx = 0; xx < x - 1; xx++) {
            addHex(hexCoords(startx + xx, yy));
        }
        startx--;
    }
    if (pentagon) {
        markOneVertexAsPentagon();
    }
}

void Polyomino::buildRaggedBiggerBoxShape(
    int x, int y,
    bool pentagon) // box alternating rows of length x and x+1
{
    clear();
    int startx = 0;
    for (int yy = 0; yy < y; yy++) {
        // smaller row
        for (int xx = 0; xx < x; xx++) {
            addHex(hexCoords(startx + xx, yy));
        }
        yy++;
        if (yy >= y) {
            break;
        }
        for (int xx = 0; xx < x + 1; xx++) {
            addHex(hexCoords(startx + xx - 1, yy));
        }
        startx--;
    }
    if (pentagon) {
        markOneVertexAsPentagon();
    }
}

void Polyomino::buildRaggedBoxShape(int x, int y, bool pentagon)
{
    clear();
    int startx = 0;
    for (int yy = 0; yy < y; yy++) {
        for (int xx = 0; xx < x; xx++) {
            addHex(hexCoords(startx + xx, yy));
        }
        yy++;
        if (yy >= y) {
            break;
        }
        for (int xx = 0; xx < x; xx++) {
            addHex(hexCoords(startx + xx, yy));
        }
        startx--;
    }
    if (pentagon) {
        markOneVertexAsPentagon();
    }
}

void Polyomino::buildWithVerticesN(int totVertices)
{
    clear();
    addHex(hexCoords(0, 0));
    addHex(hexCoords(1, 0));

    int vertices = 10;
    while (vertices < totVertices) { // find a hexagon with 2 neighbors to add.
        std::vector<hexCoords> nextCoords = allFreeNeighbors();
        int lowestDistance = -1;
        unsigned int bestI = 0;
        for (unsigned int i = 0; i < nextCoords.size(); i++) {
            hexCoords c = nextCoords[i];
            if (countNeighbors(c) != 2) {
                continue;
            }
            int distance = c.distanceFrom(hexCoords(0, 0));

            if (lowestDistance == -1 || distance < lowestDistance) {
                lowestDistance = distance;
                bestI = i;
            }
        }
        assert(lowestDistance != -1);

        addHex(nextCoords[bestI]);
        for (unsigned int i = 0; i < nextCoords.size(); i++) {
            if (i == bestI) {
                continue;
            }
            hexCoords c = nextCoords[i];
            if (countNeighbors(c) == 3) {
                addHex(c);
            }
        }
        vertices += 2;
    }
    if ((vertices - totVertices) == 1) {
        markOneVertexAsPentagon();
    }
}

void Polyomino::markOneVertexAsPentagon() // TO DO: should check if one vertex
                                          // is already marked to avoid marking
                                          // it again. Not a problem while we
                                          // only mark one vertex per polyomino
{
    vector<vertexCoords> pathV = getPath();
    size_t previousMultiplicity = hexagonsAtVertex(pathV[pathV.size() - 1]);
    size_t thisMultiplicity = hexagonsAtVertex(pathV[0]);
    size_t nextMultiplicity = 0;
    for (size_t i = 0; i < pathV.size(); i++) {
        size_t nextI = (i >= pathV.size() - 1) ? 0 : i + 1;
        nextMultiplicity = hexagonsAtVertex(pathV[nextI]);
        if (previousMultiplicity == 2 && thisMultiplicity == 1 &&
            nextMultiplicity == 2) {
            setPentagon(pathV[i]);
            return;
        }
        previousMultiplicity = thisMultiplicity;
        thisMultiplicity = nextMultiplicity;
    }
    previousMultiplicity = hexagonsAtVertex(pathV[pathV.size() - 1]);
    thisMultiplicity = hexagonsAtVertex(pathV[0]);
    nextMultiplicity = 0;
    for (size_t i = 0; i < pathV.size(); i++) {
        size_t nextI = (i >= pathV.size() - 1) ? 0 : i + 1;
        nextMultiplicity = hexagonsAtVertex(pathV[nextI]);
        if (previousMultiplicity == 1 && thisMultiplicity == 2 &&
            nextMultiplicity == 1) {
            setPentagon(pathV[i]);
            return;
        }
        previousMultiplicity = thisMultiplicity;
        thisMultiplicity = nextMultiplicity;
    }
}

void Polyomino::setPentagon(vertexCoords c) // the marked vertex and the
                                            // following one will count as a
                                            // single vertex to turn the hexagon
                                            // into a pentagon
{
    pentagonVertices.push_back(c);
}

size_t Polyomino::hexagonsAtVertex(vertexCoords v) const
{
    return vertexNeighbors(v).size();
}

vector<hexCoords> Polyomino::freeVertexNeighborPositions(vertexCoords v) const
{
    vector<hexCoords> out;
    int direction = v.x + v.y + v.z;
    if (direction != 1 && direction != -1) {
        cerr << "wrong input to free vertex neighbor positions " << v << endl;
        return out;
    }
    Hex* h = getHex(hexCoords(v.x - direction, v.y));
    if (!h) {
        out.emplace_back(v.x - direction, v.y);
    }
    h = getHex(hexCoords(v.x, v.y - direction));
    if (!h) {
        out.emplace_back(v.x, v.y - direction);
    }
    h = getHex(hexCoords(v.x, v.y)); // z - direction
    if (!h) {
        out.emplace_back(v.x, v.y);
    }
    return out;
}

vector<Hex*> Polyomino::vertexNeighbors(vertexCoords v) const
{
    vector<Hex*> out;
    int direction = v.x + v.y + v.z;
    if (direction != 1 && direction != -1) {
        cerr << "wrong input to vertex Neighbors " << v << endl;
        return out;
    }
    Hex* h = getHex(hexCoords(v.x - direction, v.y));
    if (h) {
        out.push_back(h);
    }
    h = getHex(hexCoords(v.x, v.y - direction));
    if (h) {
        out.push_back(h);
    }
    h = getHex(hexCoords(v.x, v.y)); // z - direction
    if (h) {
        out.push_back(h);
    }
    return out;
}

vertexCoords Polyomino::findOuterVertex()
    const // find an outer vertex. Used to start path around polyomino.
{
    // find a hexagon with a free vertex in direction (1, 0, 0). Such a hexagon
    // is guarateed to exist in a polyomino.
    for (auto h : m_list) {
        vertexCoords vert(h->x() + 1, h->y(), h->z());
        if (hexagonsAtVertex(vert) == 1) {
            return vert;
        }
    }
    cerr << "something went wrong in finding the outer vertex" << endl;
    return {0, 0, 0};
}

int Polyomino::countNeighbors(hexCoords h) const
{
    int out = 0;
    vector<hexCoords> neighs = Hex::neighboringPositions(h);
    for (auto neigh : neighs) {
        if (getHex(neigh) != nullptr) {
            out++;
        }
    }
    return out;
}

std::vector<hexCoords> Polyomino::allFreeNeighbors() const
{
    for (auto i : m_list) { // make sure that the grid is big enough to store
                            // every neighbor
        getIndexInList(hexCoords(i->x() + 1, i->y() + 1));
        getIndexInList(hexCoords(i->x() - 1, i->y() - 1));
    }
    std::vector<hexCoords> out;
    std::vector<bool> visited(
        m_grid.size(),
        false); // keep track if a neighbors has already been checked or not
    for (auto i : m_list) {
        std::vector<hexCoords> neighborsCoords = i->neighbors();
        for (auto neighborsCoord : neighborsCoords) {
            bool isPresent = (getHex(neighborsCoord) != nullptr);
            if (isPresent) {
                continue;
            }
            int index = getIndexInList(neighborsCoord);
            if (visited[index]) {
                continue;
            }
            visited[index] = true;
            out.push_back(neighborsCoord);
        }
    }
    return out;
}

int Polyomino::getIndexInList(hexCoords coords) const
{
    int x = coords.x;
    int y = coords.y;
    if (abs(x) > m_gridSize) {
        resizeGrid(abs(x));
    }
    if (abs(y) > m_gridSize) {
        resizeGrid(abs(y));
    }
    int a = m_gridSize + x;
    int b = m_gridSize + y;
    int i = (2 * m_gridSize + 1) * a + b;
    return i;
}

void Polyomino::addHex(hexCoords coords)
{
    int index = getIndexInList(coords);
    assert(m_grid[index] == nullptr);
    Hex* h = new Hex(coords);
    m_list.push_back(h);
    m_grid[index] = h;
}

void Polyomino::removeHex(hexCoords coords)
{
    int index = getIndexInList(coords);
    Hex* hex = m_grid[getIndexInList(coords)];
    assert(hex != nullptr);
    for (unsigned int i = 0; i < m_list.size(); i++) {
        if (m_list[i] == hex) {
            m_list.erase(m_list.begin() + i);
            break;
        }
    }
    delete hex;
    m_grid[index] = nullptr;
}

bool Polyomino::isEquivalentWithout(hexCoords c) const
{
    /* does removing this hexagon yield another polyomino with the same number
     of vertices? true if hte hexagon has 3 neighbors all next to each other */
    if (countNeighbors(c) != 3) {
        return false;
    }
    vector<hexCoords> neighs = Hex::neighboringPositions(c);
    for (unsigned int i = 0; i < neighs.size(); i++) {
        int i2 = (i - 1 + 6) % 6;
        int i3 = (i - 2 + 6) % 6;
        if (getHex(neighs[i]) != nullptr && getHex(neighs[i2]) != nullptr &&
            getHex(neighs[i3]) != nullptr) {
            return true;
        }
    }
    return false;
}

vertexCoords Polyomino::coordinatesOfSubstituent(const vertexCoords pos) const
{
    vector<Hex*> neighbors = vertexNeighbors(pos);
    assert(!neighbors.empty());
    assert(neighbors.size() < 3);
    vertexCoords out = pos;

    if (neighbors.size() == 1) {
        /*
         the vertex coordinates differ from its parent's by a single +1 or -1
         (d). The substituent will be along the same direction, so the other two
         coordinates will be incremented by -d. e.g. (1, 0, 0)-> (1, -1, -1)
         or (0, -1, 0)-> (1, -1, 1)
         */
        vertexCoords parentCoords = neighbors[0]->coords().toVertexCoords();
        vertexCoords v = pos - parentCoords;
        int d = -1;
        if (v.x + v.y + v.z > 0) {
            d = 1;
        }
        if (v.x == 0) {
            v.x -= d;
        }
        if (v.y == 0) {
            v.y -= d;
        }
        if (v.z == 0) {
            v.z -= d;
        }
        out = parentCoords + v;
    } else if (neighbors.size() == 2) {
        /*
         the two parents share two vertices. One is pos, the other one is the
         neighbor we are looking for
         */
        vertexCoords parent1Coords = neighbors[0]->coords().toVertexCoords();
        vertexCoords v = pos - parent1Coords;
        vertexCoords parent2Coords = neighbors[1]->coords().toVertexCoords();
        out = parent2Coords - v;
    }
    return out;
}

vector<vertexCoords> Polyomino::getPath() const
{
    vector<vertexCoords> out;
    vertexCoords firstVertex = findOuterVertex();
    vertexCoords currentVertex = firstVertex;
    vector<Hex*> neighbors = vertexNeighbors(firstVertex);
    assert(neighbors.size() == 1);
    Hex* currentHex = neighbors[0];
    vertexCoords nextVertex = currentHex->followingVertex(currentVertex);
    do {
        bool skip = false;
        if (!pentagonVertices.empty()) {
            for (auto pentagonVertice : pentagonVertices) {
                if (pentagonVertice == currentVertex) {
                    skip = true;
                    break;
                }
            }
        }
        if (!skip) {
            out.push_back(currentVertex);
        }

        currentVertex = nextVertex;

        neighbors = vertexNeighbors(currentVertex);
        assert(neighbors.size() <= 2);
        if (neighbors.size() == 2) {
            currentHex =
                (neighbors[0] != currentHex ? neighbors[0] : neighbors[1]);
        }

        nextVertex = currentHex->followingVertex(currentVertex);

    } while (currentVertex != firstVertex);
    return out;
}

vector<hexCoords> Hex::neighboringPositions(hexCoords h)
{
    int xx = h.x;
    int yy = h.y;
    vector<hexCoords> out;
    out.emplace_back(xx + 1, yy); // z-1
    out.emplace_back(xx + 1, yy - 1);
    out.emplace_back(xx, yy - 1); // z+1
    out.emplace_back(xx - 1, yy); // z+1
    out.emplace_back(xx - 1, yy + 1);
    out.emplace_back(xx, yy + 1); // z-1
    return out;
}

vector<hexCoords> Hex::neighbors() const
{
    return neighboringPositions(m_coords);
}
vertexCoords Hex::followingVertex(vertexCoords v) const
{
    int dx = v.x - x();
    int dy = v.y - y();
    int dz = v.z - z();
    if (dx + dy + dz != 1 && dx + dy + dz != -1) {
        cerr << "wrong input to transform to following vertex" << endl;
    }
    if (dx == 0 && dy == 0) {
        dx = -dz;
        dy = 0;
        dz = 0;
    } else if (dx == 0 && dz == 0) {
        dz = -dy;
        dx = 0;
        dy = 0;
    } else if (dy == 0 && dz == 0) {
        dy = -dx;
        dx = 0;
        dz = 0;
    } else {
        cerr << "wrong input to transform to following vertex" << endl;
    }
    return {x() + dx, y() + dy, z() + dz};
}

float CoordgenMacrocycleBuilder::getPrecision() const
{
    return m_precision;
}

void CoordgenMacrocycleBuilder::setPrecision(float f)
{
    m_precision = f;
}

vector<sketcherMinimizerPointF> CoordgenMacrocycleBuilder::newMacrocycle(
    sketcherMinimizerRing* ring, vector<sketcherMinimizerAtom*> atoms) const
{
    // TODO the coordinates should be built on the returned vector, never saved
    // to the atoms in atoms.
    int natoms = static_cast<int>(atoms.size());
    Polyomino p;
    p.buildWithVerticesN(natoms);
    vector<Polyomino> pols;
    pols.push_back(p);
    vector<Polyomino> squarePols = buildSquaredShapes(natoms);
    for (unsigned int i = 0; i < squarePols.size(); i++) {
        assert(squarePols[i].getPath().size() == atoms.size());
    }
    pols.reserve(pols.size() + squarePols.size());
    pols.insert(pols.end(), squarePols.begin(), squarePols.end());
    pols = removeDuplicates(pols);
    bool found = false;
    pathRestraints pr = getPathRestraints(atoms);
    pathConstraints pc = getPathConstraints(atoms);
    int scoreOfChosen = PATH_FAILED;
    int startOfChosen = 0;
    int bestStart = 0;
    int bestScore = PATH_FAILED;
    Polyomino chosenP = pols[0];
    int acceptableScore = acceptableShapeScore(natoms);
    int checkedMacrocycles = 0;
    if (!m_forceOpenMacrocycles) {
        do {
            int bestP = 0;
            found = matchPolyominoes(pols, pc, pr, bestP, bestScore, bestStart,
                                     checkedMacrocycles);
            if (bestScore > scoreOfChosen) {
                startOfChosen = bestStart;
                scoreOfChosen = bestScore;
                chosenP = pols[bestP];
                if (bestScore > acceptableScore) {
                    break;
                }
            }
            if (checkedMacrocycles > MAX_MACROCYCLES) {
                break;
            }

            pols = listOfEquivalents(pols);
            pols = removeDuplicates(pols);

        } while (!pols.empty());
    }
    if (found) {
        vector<vertexCoords> path = chosenP.getPath();

        writePolyominoCoordinates(path, atoms, startOfChosen);
        if (!chosenP.pentagonVertices.empty()) {
            atoms.at(0)->molecule->requireMinimization();
        }
    } else { // could not find a shape. fallback methods

        if (!openCycleAndGenerateCoords(ring)) {
            vector<sketcherMinimizerPointF> coords =
                CoordgenFragmentBuilder::listOfCoordinatesFromListofRingAtoms(
                    atoms);
            int i = 0;
            for (sketcherMinimizerAtom* atom : atoms) {
                atom->setCoordinates(coords[i]);
                ++i;
            }
        }
        atoms.at(0)->molecule->requireMinimization();
    }
    vector<sketcherMinimizerPointF> coordinates;
    coordinates.reserve(atoms.size());
    for (sketcherMinimizerAtom* atom : atoms) {
        coordinates.push_back(atom->getCoordinates());
    }
    return coordinates;
}

sketcherMinimizerBond*
CoordgenMacrocycleBuilder::findBondToOpen(sketcherMinimizerRing* ring) const
{
    sketcherMinimizerBond* bestBond = nullptr;
    float bestScore = 0.f;
    for (sketcherMinimizerBond* bond : ring->_bonds) {
        float score = 0.f;
        if (ring->isMacrocycle()) {
            if (bond->getBondOrder() != 1) {
                continue;
            }
            bool nextToMultipleBond = false;
            for (auto otherBond : bond->getStartAtom()->bonds) {
                if (otherBond->isStereo()) {
                    nextToMultipleBond = true;
                    break;
                }
            }
            for (auto otherBond : bond->getEndAtom()->bonds) {
                if (otherBond->isStereo()) {
                    nextToMultipleBond = true;
                    break;
                }
            }
            if (nextToMultipleBond) {
                continue;
            }
        }
        score += bond->rings.size() * 10;
        score += bond->getStartAtom()->neighbors.size();
        score += bond->getEndAtom()->neighbors.size();
        score /= bond->crossingBondPenaltyMultiplier;
        if (bestBond == nullptr || score < bestScore) {
            bestScore = score;
            bestBond = bond;
        }
    }
    return bestBond;
}

bool CoordgenMacrocycleBuilder::openCycleAndGenerateCoords(
    sketcherMinimizerRing* ring) const // generate coordinates of whole fragment
                                       // by opening one bond of the macrocycle
{
    map<sketcherMinimizerAtom*, sketcherMinimizerAtom*> atomMap;
    sketcherMinimizer min(getPrecision());
    min.setSkipMinimization(true);
    min.setForceOpenMacrocycles(true);
    sketcherMinimizerBond* bondToBreak = findBondToOpen(ring);
    if (!bondToBreak) {
        return false;
    }
    auto* minMol = new sketcherMinimizerMolecule;
    sketcherMinimizerAtom* a = ring->_atoms[0];
    // vector<sketcherMinimizerAtom*> atoms = a->getFragment()->getAtoms();
    vector<sketcherMinimizerAtom*> atoms = a->molecule->getAtoms();
    for (sketcherMinimizerAtom* at : atoms) {
        if (at->isResidue()) {
            continue;
        }
        auto* new_at = new sketcherMinimizerAtom;
        atomMap[at] = new_at;
        new_at->templateCoordinates = at->coordinates;
        new_at->coordinates = at->coordinates;
        new_at->atomicNumber = at->atomicNumber;
        new_at->_generalUseN = at->_generalUseN;
        new_at->molecule = minMol;
        new_at->_implicitHs = at->_implicitHs;
        minMol->_atoms.push_back(new_at);
    }

    //  vector<sketcherMinimizerBond*> bonds = a->getFragment()->getBonds();
    vector<sketcherMinimizerBond*> bonds = a->molecule->getBonds();

    for (sketcherMinimizerBond* bo : bonds) {
        if (bo == bondToBreak || bo->isResidueInteraction()) {
            continue;
        }
        auto* new_bo = new sketcherMinimizerBond;
        new_bo->bondOrder = bo->bondOrder;
        new_bo->startAtom = atomMap[bo->startAtom];
        new_bo->endAtom = atomMap[bo->endAtom];
        new_bo->isZ = bo->isZ;
        new_bo->m_ignoreZE = bo->m_ignoreZE;
        minMol->_bonds.push_back(new_bo);
    }
    min.initialize(minMol);
    min.findFragments();
    min.buildFromFragments(true);
    auto* brokenBond = new sketcherMinimizerBond;
    brokenBond->bondOrder = bondToBreak->bondOrder;
    brokenBond->startAtom = atomMap[bondToBreak->startAtom];
    brokenBond->endAtom = atomMap[bondToBreak->endAtom];
    brokenBond->isZ = bondToBreak->isZ;
    minMol->_bonds.push_back(brokenBond);
    minMol->forceUpdateStruct(minMol->_atoms, minMol->_bonds, minMol->_rings);
    std::vector<sketcherMinimizerInteraction*> extraInteractions;
    extraInteractions.push_back(new sketcherMinimizerStretchInteraction(
        brokenBond->startAtom, brokenBond->endAtom));
    min.avoidClashesOfMolecule(minMol, extraInteractions);

    for (sketcherMinimizerAtom* atom : atoms) {
        sketcherMinimizerAtom* otherAtom = atomMap[atom];
        if (otherAtom->rigid) {
            atom->rigid = true;
        }
        atom->setCoordinates(otherAtom->getCoordinates());
    }
    for (auto r : ring->getAtoms()[0]->fragment->getRings()) {
        r->coordinatesGenerated = true;
    }
    for (auto bondRing : bondToBreak->rings) {
        if (!bondRing->isMacrocycle() && bondRing != ring) {
            bondRing->coordinatesGenerated = false;
        }
    }
    delete brokenBond;
    return true;
}

vector<Polyomino>
CoordgenMacrocycleBuilder::listOfEquivalents(const vector<Polyomino>& l) const
{
    vector<Polyomino> out;
    for (const auto& i : l) {
        vector<Polyomino> newV = listOfEquivalent(i);
        out.reserve(out.size() + newV.size());
        out.insert(out.end(), newV.begin(), newV.end());
    }
    return out;
}

vector<Polyomino>
CoordgenMacrocycleBuilder::listOfEquivalent(const Polyomino& p)
    const // build a list of polyominoes with the same number of
          // vertices by removing hexagons with 3 neighbors
{
    vector<Polyomino> out;
    vector<Hex*> l = p.m_list;
    size_t pentagonVs = p.pentagonVertices.size();
    for (auto& i : l) {
        hexCoords c = i->coords();
        if (p.isEquivalentWithout(c)) {
            Polyomino newP = p;
            newP.pentagonVertices.clear();
            newP.removeHex(c);
            for (size_t i = 0; i < pentagonVs; i++) {
                newP.markOneVertexAsPentagon();
            }
            out.push_back(newP);
        }
    }
    return out;
}

pathConstraints CoordgenMacrocycleBuilder::getPathConstraints(
    vector<sketcherMinimizerAtom*>& atoms) const
{
    pathConstraints pc;
    pc.doubleBonds = getDoubleBondConstraints(atoms);
    pc.ringConstraints = getRingConstraints(atoms);
    return pc;
}

vector<ringConstraint> CoordgenMacrocycleBuilder::getRingConstraints(
    vector<sketcherMinimizerAtom*>& atoms) const
{
    vector<ringConstraint> out;
    for (int i = 0; i < static_cast<int>(atoms.size()); i++) {
        sketcherMinimizerAtom* a = atoms[i];
        if (a->rings.size() > 1) {
            for (unsigned int rr = 0; rr < a->rings.size(); rr++) {
                sketcherMinimizerRing* r = a->rings[rr];
                if (r->_atoms.size() < MACROCYCLE) { // this excludes the
                                                     // current cycle and all
                                                     // fused macrocycles
                    bool forceOutside = false;
                    for (auto n : a->neighbors) {
                        if (find(atoms.begin(), atoms.end(), n) ==
                            atoms.end()) {
                            if (r->containsAtom(n)) {
                                forceOutside = true;
                            }
                            break;
                        }
                    }
                    out.emplace_back(i, r, forceOutside);
                }
            }
        }
    }
    return out;
}

vector<doubleBondConstraint>
CoordgenMacrocycleBuilder::getDoubleBondConstraints(
    vector<sketcherMinimizerAtom*>& atoms) const
{
    vector<doubleBondConstraint> out;
    if (atoms.size() >= MACROCYCLE) {
        for (unsigned int i = 0; i < atoms.size(); i++) {
            unsigned int index = i;
            unsigned int index2 = (i + 1) % atoms.size();
            sketcherMinimizerBond* b =
                sketcherMinimizer::getBond(atoms[i], atoms[index2]);
            if (!b) {
                cerr << "bad input to get double bond constraints" << endl;
                break;
            }
            if (b->bondOrder != 2) {
                continue;
            }
            bool smallRingBond = false;
            if (b->rings.size() > 1) {
                for (auto& ring : b->rings) {
                    if (ring->_atoms.size() < MACROCYCLE) {
                        smallRingBond = true;
                        break;
                    }
                }
            }
            if (smallRingBond) {
                continue;
            }
            int previousI =
                static_cast<int>((i + atoms.size() - 1) % atoms.size());
            int followingI = (i + 2) % atoms.size();

            bool isTrans = !b->isZ; // is the atom trans in the ring? (isZ
                                    // stores the absolute chirality)

            if (b->startAtom !=
                atoms[i]) { // atoms[i] could be the end atom of b
                int swap = previousI;
                previousI = followingI;
                followingI = swap;
                index = index2;
                index2 = i;
            }
            if (b->startAtomCIPFirstNeighbor() != atoms[previousI]) {
                isTrans = !isTrans;
            }
            if (b->endAtomCIPFirstNeighbor() != atoms[followingI]) {
                isTrans = !isTrans;
            }
            doubleBondConstraint constraint;
            constraint.trans = isTrans;
            constraint.atom1 = index;
            constraint.atom2 = index2;
            constraint.previousAtom = previousI;
            constraint.followingAtom = followingI;
            out.push_back(constraint);
        }
    }
    return out;
}

int CoordgenMacrocycleBuilder::getNumberOfChildren(
    sketcherMinimizerAtom* a, sketcherMinimizerAtom* parent) const
{
    int n = 0;
    map<sketcherMinimizerAtom*, bool> visited;
    queue<sketcherMinimizerAtom*> q;
    q.push(a);
    visited[parent] = true;
    while (!q.empty()) {
        sketcherMinimizerAtom* thisA = q.front();
        q.pop();
        visited[thisA] = true;
        n++;
        for (auto n : thisA->neighbors) {
            if (visited[n]) {
                continue;
            }
            q.push(n);
        }
    }
    return n;
}

pathRestraints CoordgenMacrocycleBuilder::getPathRestraints(
    vector<sketcherMinimizerAtom*>& atoms) const
{
    pathRestraints pr;
    vector<int> heteroAtoms;
    vector<pair<int, int>> substitutedAtoms;
    for (unsigned int i = 0; i < atoms.size(); i++) {
        if (atoms[i]->atomicNumber != 6) {
            heteroAtoms.push_back(i);
        }
        if (atoms[i]->neighbors.size() != 2) {
            int totN = 0;
            int prevI = static_cast<int>((i + atoms.size() - 1) % atoms.size());
            int postI = (i + 1) % atoms.size();
            for (unsigned int j = 0; j < atoms[i]->neighbors.size(); j++) {
                sketcherMinimizerAtom* n = atoms[i]->neighbors[j];
                if (n == atoms[prevI]) {
                    continue;
                }
                if (n == atoms[postI]) {
                    continue;
                }
                totN += getNumberOfChildren(n, atoms[i]);
            }
            substitutedAtoms.emplace_back(i, totN);
        }
    }
    pr.heteroAtoms = heteroAtoms;
    pr.substitutedAtoms = substitutedAtoms;
    return pr;
}

bool CoordgenMacrocycleBuilder::checkDoubleBoundConstraints(
    vector<doubleBondConstraint>& dbConstraints, vector<vertexCoords>& vertices,
    int& startI) const
{
    for (auto& dbConstraint : dbConstraints) {
        size_t counter = (startI + dbConstraint.previousAtom) % vertices.size();
        sketcherMinimizerPointF p1 = coordsOfVertex(vertices[counter]);
        counter = (startI + dbConstraint.atom1) % vertices.size();
        sketcherMinimizerPointF p2 = coordsOfVertex(vertices[counter]);
        counter = startI + dbConstraint.atom2;
        if (counter >= vertices.size()) {
            counter -= vertices.size();
        }
        sketcherMinimizerPointF p3 = coordsOfVertex(vertices[counter]);
        counter = startI + dbConstraint.followingAtom;
        if (counter >= vertices.size()) {
            counter -= vertices.size();
        }
        sketcherMinimizerPointF p4 = coordsOfVertex(vertices[counter]);
        if (sketcherMinimizerMaths::sameSide(p1, p4, p2, p3) ==
            dbConstraint.trans) {
            return false;
        }
    }
    return true;
}

int CoordgenMacrocycleBuilder::scorePath(Polyomino& p,
                                         vector<vertexCoords>& path,
                                         vector<int>& neighborNs, int& startI,
                                         pathConstraints& pc,
                                         pathRestraints& pr) const
{
    if (!scorePathConstraints(pc, p, path, neighborNs, startI)) {
        return PATH_FAILED;
    }
    return scorePathRestraints(pr, p, path, neighborNs, startI);
}

int CoordgenMacrocycleBuilder::scorePathRestraints(pathRestraints& pr,
                                                   Polyomino& p,
                                                   vector<vertexCoords>& path,
                                                   vector<int>& neighborNs,
                                                   int& startI) const
{
    int score = 0;
    set<vertexCoords> usedSubstituentCoordinates;
    for (int heteroAtom : pr.heteroAtoms) {
        int counter = (heteroAtom + startI) % neighborNs.size();
        if (neighborNs[counter] == 1) {
            score -= HETEROATOM_RESTRAINT; // heteroatom placed towards outside
        }
    }
    for (unsigned int i = 0; i < pr.substitutedAtoms.size(); i++) {
        int counter =
            (pr.substitutedAtoms[i].first + startI) % neighborNs.size();
        if (neighborNs[counter] == 2) {
            score -= SUBSTITUTED_ATOM_RESTRAINT *
                     pr.substitutedAtoms[i]
                         .second; // substituted atoms placed towards inside
            vertexCoords substituentsCoordinates =
                p.coordinatesOfSubstituent(path[i]);
            if (usedSubstituentCoordinates.find(substituentsCoordinates) !=
                usedSubstituentCoordinates.end()) {
                score -= SUBSTITUENTS_TO_SAME_VERTEX_RESTRAINT;
            }
            if (find(path.begin(), path.end(), substituentsCoordinates) !=
                path.end()) {
                // substituent would clash with path
                score -= SUBSTITUENT_CLASH_WITH_PATH_RESTRAINT;
            }
            usedSubstituentCoordinates.insert(substituentsCoordinates);
        }
    }
    return score;
}

bool CoordgenMacrocycleBuilder::scorePathConstraints(pathConstraints& pc,
                                                     Polyomino& p,
                                                     vector<vertexCoords>& path,
                                                     vector<int>& neighborNs,
                                                     int& startI) const
{

    if (!checkRingConstraints(pc.ringConstraints, p, path, neighborNs,
                              startI)) {
        return false;
    }
    if (!checkDoubleBoundConstraints(pc.doubleBonds, path, startI)) {
        return false;
    }
    return true;
}

bool CoordgenMacrocycleBuilder::checkRingConstraints(
    vector<ringConstraint>& ringConstraints, Polyomino& p,
    vector<vertexCoords>& path, vector<int>& neighborNs,
    int& startI) const // make sure that all shared atoms with a fused small
                       // ring map to the same hexagon
{
    std::map<sketcherMinimizerRing*, vector<hexCoords>> allowedHexs;
    for (auto& ringConstraint : ringConstraints) {
        unsigned int counter = (ringConstraint.atom + startI) % path.size();
        if (ringConstraint.forceOutside) {
            if (neighborNs[counter] != 1) {
                return false;
            }
        }
        sketcherMinimizerRing* r = ringConstraint.ring;
        vector<hexCoords> newPos = p.freeVertexNeighborPositions(path[counter]);
        vector<hexCoords> oldPos = allowedHexs[r];
        vector<hexCoords> nextPos;
        if (oldPos.empty()) {
            nextPos = newPos;
        } else {
            for (auto toCheck : newPos) {
                if (find(oldPos.begin(), oldPos.end(), toCheck) !=
                    oldPos.end()) {
                    // hexCoords was already marked as allowed and is now
                    // confirmed
                    nextPos.push_back(toCheck);
                }
            }
        }
        if (nextPos.empty()) {
            return false;
        }
        allowedHexs[r] = nextPos;
    }
    return true;
}
vector<int>
CoordgenMacrocycleBuilder::getVertexNeighborNs(Polyomino& p,
                                               vector<vertexCoords>& path) const
{
    vector<int> out;
    out.reserve(path.size());
    for (auto i : path) {
        out.push_back(static_cast<int>(p.hexagonsAtVertex(i)));
    }
    return out;
}

int CoordgenMacrocycleBuilder::getLowestPeriod(
    std::vector<int>& neighbors) const
{
    for (unsigned int period = 1; period < neighbors.size(); period++) {
        bool hasDifference = false;
        for (unsigned int i = 0; i < neighbors.size(); i++) {
            unsigned int i2 = (i + period) % neighbors.size();
            if (neighbors[i] != neighbors[i2]) {
                hasDifference = true;
                break;
            }
        }
        if (!hasDifference) {
            return period;
        }
    }
    return static_cast<int>(neighbors.size());
}

bool CoordgenMacrocycleBuilder::matchPolyominoes(vector<Polyomino>& pols,
                                                 pathConstraints& pc,
                                                 pathRestraints& pr, int& bestP,
                                                 int& bestScore, int& bestStart,
                                                 int& checkedMacrocycles) const
{
    bestStart = 0;
    bestP = 0;
    bestScore = PATH_FAILED;
    int bestScoreOfThis = PATH_FAILED;
    int bestStartOfThis = 0;
    bool matched = false;
    for (unsigned int i = 0; i < pols.size(); i++) {
        if (matchPolyomino(pols[i], pc, pr, bestStartOfThis, bestScoreOfThis)) {
            matched = true;
            if (bestScoreOfThis > bestScore) {
                bestScore = bestScoreOfThis;
                bestStart = bestStartOfThis;
                bestP = i;
                if (bestScore == 0) {
                    return true;
                }
            }
        }
        if (checkedMacrocycles++ > MAX_MACROCYCLES) {
            break;
        }
    }
    return matched;
}

bool CoordgenMacrocycleBuilder::matchPolyomino(Polyomino& p,
                                               pathConstraints& pc,
                                               pathRestraints& pr,
                                               int& bestStart,
                                               int& bestScore) const
{
    vector<vertexCoords> path = p.getPath();
    vector<int> neighborNs = getVertexNeighborNs(p, path);
    bestStart = 0;
    bestScore = PATH_FAILED;
    for (int startI = 0; startI < getLowestPeriod(neighborNs);
         startI++) { // test starting from every vertex.
        int score = scorePath(p, path, neighborNs, startI, pc, pr);
        if (score > bestScore) {
            bestScore = score;
            bestStart = startI;
            if (bestScore == 0) {
                break;
            }
        }
    }
    return bestScore > PATH_FAILED;
}

sketcherMinimizerPointF
CoordgenMacrocycleBuilder::coordsOfVertex(vertexCoords& v) const
{
    return sketcherMinimizerPointF(
        static_cast<float>(BONDLENGTH * SQRT3HALF * v.x -
                           BONDLENGTH * SQRT3HALF * v.z),
        static_cast<float>(-BONDLENGTH * 0.5 * v.x + BONDLENGTH * v.y +
                           -BONDLENGTH * 0.5 * v.z));
}

void CoordgenMacrocycleBuilder::writePolyominoCoordinates(
    vector<vertexCoords>& path, const vector<sketcherMinimizerAtom*>& atoms,
    int startI) const
{

    for (unsigned int n = 0; n < atoms.size(); n++) {
        unsigned int counter = (n + startI) % path.size();
        vertexCoords hCoords = path[counter];
        if (!atoms[n]->coordinatesSet) {
            atoms[n]->setCoordinates(coordsOfVertex(hCoords));
        }
    }
}

vector<Polyomino>
CoordgenMacrocycleBuilder::removeDuplicates(vector<Polyomino>& pols) const
{
    // This algorithm has potential to be n^2 Could be a potential bottleneck
    // for performance

    vector<Polyomino> out;
    for (auto& pol : pols) {
        bool duplicate = false;
        for (auto& j : out) {
            if (pol.isTheSameAs(j)) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) {
            out.push_back(pol);
        }
    }
    return out;
}

vector<Polyomino>
CoordgenMacrocycleBuilder::buildSquaredShapes(int totVertices) const
{
    // a square like this has 4x + 4 y - 2 vertices. Only do if number of
    // vertices (or number of vertices+1) is divisible by 2 but not by 4.  x + y
    // = (N + 2)/ 4

    vector<Polyomino> out;
    bool pentagon = false;
    if ((totVertices % 2) > 0) {
        pentagon = true;
        totVertices++;
    }
    if (totVertices % 4 == 0) {
        if (totVertices >= 12) {
            // we can build one with ragged borders and even number of y. We'll
            // only build the ones alternating rows with length x, x+1, because
            // the x, x-1 one is equivant
            int xandy = totVertices / 4;
            int mid = xandy / 2; // truncate
            for (int i = 1; i < mid; ++i) {
                int i2 = xandy - i;
                if ((i2 % 2) == 0 && i > 1) {
                    Polyomino p;
                    p.buildRaggedBiggerBoxShape(i, i2, pentagon);
                    out.insert(out.begin(), p);
                }
                if ((i % 2) == 0 && i2 > 1) {
                    Polyomino p;
                    p.buildRaggedBiggerBoxShape(i2, i, pentagon);
                    out.insert(out.begin(), p);
                }
            }
        }

    } else {
        // we'll build the skewed box, the ragged box, the boxes with odd y
        // alternating x and x+1, and the boxes with odd y alternating x and x-1
        int xandy = (totVertices + 2) / 4;
        int mid = xandy / 2; // truncate
        for (int i = 1; i < mid + 1;
             i++) { // build all possible combinations of shapes
            Polyomino p;
            p.buildSkewedBoxShape(i, xandy - i, pentagon);
            out.insert(out.begin(),
                       p); // add at the beginning because we want rounder
                           // shapes to be evaluated first (they have more space
                           // inside and less risk of substituents clashes
            int i2 = xandy - i;
            if (i < 2 || i2 < 2) {
                continue;
            }
            {
                Polyomino p;
                p.buildRaggedBoxShape(i, i2, pentagon);
                out.insert(out.begin(), p);
            }
            {
                Polyomino p;
                p.buildRaggedBoxShape(i2, i, pentagon);
                out.insert(out.begin(), p);
            }
            if (i2 % 2 > 0) {
                Polyomino p;
                p.buildRaggedBiggerBoxShape(i, i2, pentagon);
                out.insert(out.begin(), p);
            }
            if (i % 2 > 0) {
                Polyomino p;
                p.buildRaggedBiggerBoxShape(i2, i, pentagon);
                out.insert(out.begin(), p);
            }
            if (i > 2 && i2 % 2 > 0) {
                Polyomino p;
                p.buildRaggedSmallerBoxShape(i, i2, pentagon);
                out.insert(out.begin(), p);
            }
            if (i2 > 2 && i % 2 > 0) {
                Polyomino p;
                p.buildRaggedSmallerBoxShape(i2, i, pentagon);
                out.insert(out.begin(), p);
            }
        }
    }
    return out;
}

int CoordgenMacrocycleBuilder::acceptableShapeScore(int numberOfAtoms) const
{
    if (numberOfAtoms < 10) {
        return 0;
    } else {
        return numberOfAtoms * SUBSTITUTED_ATOM_RESTRAINT / 2;
    }
}
