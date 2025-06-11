//
//  Copyright (C) 2003-2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file ROMol.h

  \brief Defines the primary molecule class \c ROMol as well as associated
  typedefs

*/

#include <RDGeneral/export.h>
#ifndef RD_ROMOL_H
#define RD_ROMOL_H

/// Std stuff
#include <cstddef>
#include <iterator>
#include <utility>
#include <map>

// boost stuff
#include <RDGeneral/BoostStartInclude.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/dynamic_bitset.hpp>

#ifdef RDK_USE_BOOST_SERIALIZATION
#include <boost/serialization/split_member.hpp>
#endif
#include <RDGeneral/BoostEndInclude.h>

// our stuff
#include <RDGeneral/types.h>
#include <RDGeneral/RDProps.h>
#include "Atom.h"
#include "Bond.h"
#include "Conformer.h"
#include "RDMol.h"
#include "SubstanceGroup.h"
#include "StereoGroup.h"
#include "RingInfo.h"
#include <GraphMol/rdmol_throw.h>

namespace RDKit {
class SubstanceGroup;
class Atom;
class Bond;
class MolPickler;
class RWMol;
class QueryAtom;
class QueryBond;
class RingInfo;

template <class T1, class T2>
class AtomIterator_;
class BondIterator_;
class ConstBondIterator_;

template <class T1, class T2>
class AromaticAtomIterator_;
template <class T1, class T2>
class HeteroatomIterator_;
template <class T1, class T2>
class QueryAtomIterator_;
template <class T1, class T2>
class MatchingAtomIterator_;

RDKIT_GRAPHMOL_EXPORT extern const int ci_RIGHTMOST_ATOM;
RDKIT_GRAPHMOL_EXPORT extern const int ci_LEADING_BOND;
RDKIT_GRAPHMOL_EXPORT extern const int ci_ATOM_HOLDER;

//! ROMol is a molecule class that is intended to have a fixed topology
/*!
  This is the primary class for most molecule operations.

  If you need to be manipulating the molecule (e.g. adding or deleting
  atoms or bonds, use an RWMol instead.

  <b>Notes:</b>
    - each ROMol maintains a Dict of \c properties:
        - Each \c property is keyed by name and can store an
          arbitrary type.
        - \c Properties can be marked as \c calculated, in which case
          they will be cleared when the \c clearComputedProps() method
          is called.
        - Because they have no impact upon chemistry, all \c property
          operations are \c const, this allows extra flexibility for
          clients who need to store extra data on ROMol objects.

    - each ROMol has collections of \c bookmarks for Atoms and Bonds:
        - the Atom bookmarks and Bond bookmarks are stored separately
          from each other
        - each \c bookmark, an integer, can map to more than one
          Atom or Bond
        - these are currently used in molecule construction, but
          could also be useful for reaction mapping and the like

    - information about rings (SSSR and the like) is stored in the
      molecule's RingInfo pointer.

 */

//! \name C++11 Iterators

template <class Graph, class Vertex,
          class Iterator = typename Graph::vertex_iterator>
struct CXXAtomIterator {
  Graph *graph;
  Iterator vstart, vend;

  struct CXXAtomIter {
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = Vertex;
    using pointer = Vertex *;
    using reference = Vertex &;

    Graph *graph;
    Iterator pos;
    value_type current;

    CXXAtomIter(Graph *graph, Iterator pos)
        : graph(graph), pos(pos), current(nullptr) {}

    const value_type &operator*() {
      // The const_cast is to maintain the same behaviour as the previous
      // interface.
      current = const_cast<value_type>((*graph)[*pos]);
      return current;
    }
    CXXAtomIter &operator++() {
      ++pos;
      return *this;
    }
    bool operator==(const CXXAtomIter &it) const { return pos == it.pos; }
    bool operator!=(const CXXAtomIter &it) const { return pos != it.pos; }
  };

  CXXAtomIterator(Graph *graph) : graph(graph) {
    auto vs = graph->getVertices();
    vstart = vs.first;
    vend = vs.second;
  }
  CXXAtomIterator(Graph *graph, Iterator start, Iterator end)
      : graph(graph), vstart(start), vend(end) {};
  CXXAtomIter begin() { return {graph, vstart}; }
  CXXAtomIter end() { return {graph, vend}; }
};

template <class Graph, class Edge,
          class Iterator = typename Graph::edge_iterator>
struct CXXBondIterator {
  Graph *graph;
  Iterator vstart, vend;

  struct CXXBondIter {
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = Edge;
    using pointer = Edge *;
    using reference = Edge &;

    Graph *graph;
    Iterator pos;
    value_type current;

    CXXBondIter(Graph *graph, Iterator pos)
        : graph(graph), pos(pos), current(nullptr) {}

    const value_type &operator*() {
      // The const_cast is to maintain the same behaviour as the previous
      // interface.
      current = const_cast<value_type>((*graph)[*pos]);
      return current;
    }
    CXXBondIter &operator++() {
      ++pos;
      return *this;
    }
    bool operator==(const CXXBondIter &it) const { return pos == it.pos; }
    bool operator!=(const CXXBondIter &it) const { return pos != it.pos; }
  };

  CXXBondIterator(Graph *graph) : graph(graph) {
    auto vs = graph->getEdges();
    vstart = vs.first;
    vend = vs.second;
  }
  CXXBondIterator(Graph *graph, Iterator start, Iterator end)
      : graph(graph), vstart(start), vend(end) {};
  CXXBondIter begin() { return {graph, vstart}; }
  CXXBondIter end() { return {graph, vend}; }
};

class RDKIT_GRAPHMOL_EXPORT ROMol {
  RDMol *dp_mol;
  bool d_ownedBySelf;

  friend class MolPickler;
  friend class RDMol;
  friend class RWMol;

 public:
  // MolGraph compatibility descriptor and iterator definitions

  // For compatibility, since vertex_descriptor is often passed by reference,
  // this is kept as an exact match of size_t.
  using vertex_descriptor = size_t;
  // edge_descriptor needs to be incompatible with vertex_descriptor, without
  // an implicit conversion chain, to avoid ambiguity compile errors.
  class edge_descriptor {
    uint32_t index;

    friend class ROMol;

   public:
    edge_descriptor() : index(uint32_t(-1)) {}
    explicit edge_descriptor(uint32_t i) : index(i) {}
    edge_descriptor(const edge_descriptor &that) = default;
    edge_descriptor &operator=(const edge_descriptor &that) = default;
    bool operator==(const edge_descriptor &that) const {
      return index == that.index;
    }
    bool operator!=(const edge_descriptor &that) const {
      return index != that.index;
    }
    uint32_t operator*() const { return index; }
    edge_descriptor &operator++() {
      ++index;
      return *this;
    }
    edge_descriptor operator++(int) {
      edge_descriptor prev = *this;
      ++index;
      return prev;
    }
    edge_descriptor &operator--() {
      --index;
      return *this;
    }
    edge_descriptor operator--(int) {
      edge_descriptor prev = *this;
      --index;
      return prev;
    }
    uint32_t operator-(const edge_descriptor &that) const {
      return index - that.index;
    }

    edge_descriptor operator+(size_t offset) const {
      return edge_descriptor(index + offset);
    }
    edge_descriptor operator-(size_t offset) const {
      return edge_descriptor(index - offset);
    }
    edge_descriptor &operator+=(size_t offset) {
      index += offset;
      return *this;
    }
    edge_descriptor &operator-=(size_t offset) {
      index -= offset;
      return *this;
    }
  };

  // This iterator is for iterating over edges of the graph in order.
  class edge_iterator {
    edge_descriptor i;

   public:
    edge_iterator() = default;
    explicit edge_iterator(edge_descriptor index) : i(index) {}
    edge_iterator(const edge_iterator &) = default;
    edge_iterator &operator=(const edge_iterator &) = default;
    bool operator==(const edge_iterator &that) const { return i == that.i; }
    bool operator!=(const edge_iterator &that) const { return i != that.i; }
    edge_descriptor operator*() const { return i; }
    edge_iterator &operator++() {
      ++i;
      return *this;
    }
    edge_iterator operator++(int) {
      edge_iterator prev = *this;
      ++i;
      return prev;
    }
    edge_iterator &operator--() {
      --i;
      return *this;
    }
    edge_iterator operator--(int) {
      edge_iterator prev = *this;
      --i;
      return prev;
    }
    size_t operator-(const edge_iterator &that) const {
      return i - that.i;
    }

    edge_iterator operator+(size_t offset) const {
      return edge_iterator(i + offset);
    }
    edge_iterator operator-(size_t offset) const {
      return edge_iterator(i - offset);
    }
    edge_iterator &operator+=(size_t offset) {
      i += offset;
      return *this;
    }
    edge_iterator &operator-=(size_t offset) {
      i -= offset;
      return *this;
    }

    using iterator_category = std::forward_iterator_tag;
    using difference_type = size_t;
    using value_type = edge_descriptor;
    using pointer = edge_descriptor *;
    using reference = edge_descriptor &;
  };

  // This iterator is for iterating over the edges of a vertex.
  // This is similar to adjacency_iterator, but yielding the bond,
  // instead of the opposite vertex.
  class out_edge_iterator {
    const uint32_t *p;

   public:
    out_edge_iterator() : p(nullptr) {}
    explicit out_edge_iterator(const uint32_t *edges) : p(edges) {}
    out_edge_iterator(const out_edge_iterator &) = default;
    out_edge_iterator &operator=(const out_edge_iterator &) = default;
    bool operator==(const out_edge_iterator &that) const { return p == that.p; }
    bool operator!=(const out_edge_iterator &that) const { return p != that.p; }
    edge_descriptor operator*() const { return edge_descriptor{*p}; }
    out_edge_iterator &operator++() {
      ++p;
      return *this;
    }
    out_edge_iterator operator++(int) {
      out_edge_iterator prev = *this;
      ++p;
      return prev;
    }
    out_edge_iterator &operator--() {
      --p;
      return *this;
    }
    out_edge_iterator operator--(int) {
      out_edge_iterator prev = *this;
      --p;
      return prev;
    }
    std::ptrdiff_t operator-(const out_edge_iterator &that) const {
      return p - that.p;
    }

    out_edge_iterator operator+(size_t offset) const {
      return out_edge_iterator(p + offset);
    }
    out_edge_iterator operator-(size_t offset) const {
      return out_edge_iterator(p - offset);
    }
    out_edge_iterator &operator+=(size_t offset) {
      p += offset;
      return *this;
    }
    out_edge_iterator &operator-=(size_t offset) {
      p -= offset;
      return *this;
    }

    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = edge_descriptor;
    using pointer = edge_descriptor *;
    using reference = edge_descriptor &;
  };

  // This iterator is for iterating over vertices of the graph in order.
  class vertex_iterator {
    vertex_descriptor i;

   public:
    vertex_iterator() : i(size_t(-1)) {}
    explicit vertex_iterator(vertex_descriptor index) : i(index) {}
    vertex_iterator(const vertex_iterator &) = default;
    vertex_iterator &operator=(const vertex_iterator &) = default;
    bool operator==(const vertex_iterator &that) const { return i == that.i; }
    bool operator!=(const vertex_iterator &that) const { return i != that.i; }
    vertex_descriptor operator*() const { return i; }
    vertex_iterator &operator++() {
      ++i;
      return *this;
    }
    vertex_iterator operator++(int) {
      vertex_iterator prev = *this;
      ++i;
      return prev;
    }
    vertex_iterator &operator--() {
      --i;
      return *this;
    }
    vertex_iterator operator--(int) {
      vertex_iterator prev = *this;
      --i;
      return prev;
    }
    size_t operator-(const vertex_iterator &that) const {
      return i - that.i;
    }

    vertex_iterator operator+(size_t offset) const {
      return vertex_iterator(i + offset);
    }
    vertex_iterator operator-(size_t offset) const {
      return vertex_iterator(i - offset);
    }
    vertex_iterator &operator+=(size_t offset) {
      i += offset;
      return *this;
    }
    vertex_iterator &operator-=(size_t offset) {
      i -= offset;
      return *this;
    }

    using iterator_category = std::forward_iterator_tag;
    using difference_type = size_t;
    using value_type = vertex_descriptor;
    using pointer = vertex_descriptor *;
    using reference = vertex_descriptor &;
  };

  // This iterator is for iterating over the vertices adjacent to a vertex.
  // This is similar to out_edge_iterator, but yielding the opposite vertex,
  // instead of the bond.
  class adjacency_iterator {
    const uint32_t *p;

   public:
    adjacency_iterator() : p(nullptr) {}
    explicit adjacency_iterator(const uint32_t *neighbors) : p(neighbors) {}
    adjacency_iterator(const adjacency_iterator &) = default;
    adjacency_iterator &operator=(const adjacency_iterator &) = default;
    bool operator==(const adjacency_iterator &that) const {
      return p == that.p;
    }
    bool operator!=(const adjacency_iterator &that) const {
      return p != that.p;
    }
    vertex_descriptor operator*() const { return vertex_descriptor{*p}; }
    adjacency_iterator &operator++() {
      ++p;
      return *this;
    }
    adjacency_iterator operator++(int) {
      adjacency_iterator prev = *this;
      ++p;
      return prev;
    }
    adjacency_iterator &operator--() {
      --p;
      return *this;
    }
    adjacency_iterator operator--(int) {
      adjacency_iterator prev = *this;
      --p;
      return prev;
    }
    std::ptrdiff_t operator-(const adjacency_iterator& that) const {
      return p - that.p;
    }

    adjacency_iterator operator+(size_t offset) const {
      return adjacency_iterator(p + offset);
    }
    adjacency_iterator operator-(size_t offset) const {
      return adjacency_iterator(p - offset);
    }
    adjacency_iterator &operator+=(size_t offset) {
      p += offset;
      return *this;
    }
    adjacency_iterator &operator-=(size_t offset) {
      p -= offset;
      return *this;
    }

    bool operator<(const adjacency_iterator &that) const {
      return p < that.p;
    }

    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = vertex_descriptor;
    using pointer = vertex_descriptor *;
    using reference = vertex_descriptor &;
  };

  //! \cond TYPEDEFS

  //! \name typedefs
  //! @{
  typedef edge_iterator EDGE_ITER;
  typedef out_edge_iterator OEDGE_ITER;
  typedef vertex_iterator VERTEX_ITER;
  typedef adjacency_iterator ADJ_ITER;
  typedef std::pair<EDGE_ITER, EDGE_ITER> BOND_ITER_PAIR;
  typedef std::pair<OEDGE_ITER, OEDGE_ITER> OBOND_ITER_PAIR;
  typedef std::pair<VERTEX_ITER, VERTEX_ITER> ATOM_ITER_PAIR;
  typedef std::pair<ADJ_ITER, ADJ_ITER> ADJ_ITER_PAIR;

  typedef std::vector<Atom *> ATOM_PTR_VECT;
  typedef ATOM_PTR_VECT::iterator ATOM_PTR_VECT_I;
  typedef ATOM_PTR_VECT::const_iterator ATOM_PTR_VECT_CI;
  typedef std::vector<Bond *> BOND_PTR_VECT;
  typedef BOND_PTR_VECT::iterator BOND_PTR_VECT_I;
  typedef BOND_PTR_VECT::const_iterator BOND_PTR_VECT_CI;

  typedef std::list<Atom *> ATOM_PTR_LIST;
  typedef ATOM_PTR_LIST::iterator ATOM_PTR_LIST_I;
  typedef ATOM_PTR_LIST::const_iterator ATOM_PTR_LIST_CI;
  typedef std::list<Bond *> BOND_PTR_LIST;
  typedef BOND_PTR_LIST::iterator BOND_PTR_LIST_I;
  typedef BOND_PTR_LIST::const_iterator BOND_PTR_LIST_CI;

  // list of conformations
  typedef std::list<CONFORMER_SPTR> CONF_SPTR_LIST;
  typedef CONF_SPTR_LIST::iterator CONF_SPTR_LIST_I;
  typedef CONF_SPTR_LIST::const_iterator CONF_SPTR_LIST_CI;
  typedef std::pair<CONF_SPTR_LIST_I, CONF_SPTR_LIST_I> CONFS_I_PAIR;

  // ROFIX: these will need to be readonly somehow?
  typedef std::map<int, ATOM_PTR_LIST> ATOM_BOOKMARK_MAP;
  typedef std::map<int, BOND_PTR_LIST> BOND_BOOKMARK_MAP;

  typedef class AtomIterator_<Atom, ROMol> AtomIterator;
  typedef class AtomIterator_<const Atom, const ROMol> ConstAtomIterator;
  typedef class BondIterator_ BondIterator;
  typedef class ConstBondIterator_ ConstBondIterator;
  typedef class AromaticAtomIterator_<Atom, ROMol> AromaticAtomIterator;
  typedef class AromaticAtomIterator_<const Atom, const ROMol>
      ConstAromaticAtomIterator;
  typedef class HeteroatomIterator_<Atom, ROMol> HeteroatomIterator;
  typedef class HeteroatomIterator_<const Atom, const ROMol>
      ConstHeteroatomIterator;
  typedef class QueryAtomIterator_<Atom, ROMol> QueryAtomIterator;
  typedef class QueryAtomIterator_<const Atom, const ROMol>
      ConstQueryAtomIterator;
  typedef class MatchingAtomIterator_<Atom, ROMol> MatchingAtomIterator;
  typedef class MatchingAtomIterator_<const Atom, const ROMol>
      ConstMatchingAtomIterator;

  typedef CONF_SPTR_LIST_I ConformerIterator;
  typedef CONF_SPTR_LIST_CI ConstConformerIterator;

  //! @}
  //! \endcond

  //! C++11 Range iterator
  /*!
    <b>Usage</b>
    \code
      for(auto atom : mol.atoms()) {
         atom->getIdx();
      };
    \endcode
   */

  CXXAtomIterator<ROMol, Atom *> atoms() { return {this}; }

  // NOTE: This yields pointers to non-const Atoms, because the previous
  // interface did.
  CXXAtomIterator<const ROMol, Atom *> atoms() const { return {this}; }

  // NOTE: This yields pointers to non-const Atoms, because the previous
  // interface did.
  CXXAtomIterator<const ROMol, Atom *, adjacency_iterator>
  atomNeighbors(Atom const *at) const {
    auto pr = getAtomNeighbors(at);
    return {this, pr.first, pr.second};
  }

  CXXAtomIterator<ROMol, Atom *, adjacency_iterator> atomNeighbors(
      Atom const *at) {
    auto pr = getAtomNeighbors(at);
    return {this, pr.first, pr.second};
  }

  // NOTE: This yields pointers to non-const Bonds, because the previous
  // interface did.
  CXXBondIterator<const ROMol, Bond *, out_edge_iterator>
  atomBonds(Atom const *at) const {
    auto pr = getAtomBonds(at);
    return {this, pr.first, pr.second};
  }

  CXXBondIterator<ROMol, Bond *, out_edge_iterator> atomBonds(
      Atom const *at) {
    auto pr = getAtomBonds(at);
    return {this, pr.first, pr.second};
  }

  /*!
  <b>Usage</b>
  \code
    for(auto bond : mol.bonds()) {
       bond->getIdx();
    };
  \endcode
 */

  CXXBondIterator<ROMol, Bond *> bonds() { return {this}; }

  // NOTE: This yields pointers to non-const Bonds, because the previous
  // interface did.
  CXXBondIterator<const ROMol, Bond *> bonds() const { return {this}; }

 protected:
  ROMol(RDMol *mol) : dp_mol(mol), d_ownedBySelf(false) {}

 public:
  ROMol();

  //! copy constructor with a twist
  /*!
    \param other     the molecule to be copied
    \param quickCopy (optional) if this is true, the resulting ROMol will not
         copy any of the properties or bookmarks and conformers from \c other.
    This can
         make the copy substantially faster (thus the name).
    \param confId (optional) if this is >=0, the resulting ROMol will contain
    only
         the specified conformer from \c other.
  */
  ROMol(const ROMol &other, bool quickCopy = false, int confId = -1);
  //! construct a molecule from a pickle string
  ROMol(const std::string &binStr);
  //! construct a molecule from a pickle string
  ROMol(const std::string &binStr, unsigned int propertyFlags);

  ROMol(ROMol &&o) noexcept;
  ROMol &operator=(ROMol &&o) noexcept;
  ROMol &operator=(const ROMol &) =
      delete;  // disable assignment, RWMol's support assignment

  RDMol &asRDMol() { return *dp_mol; }
  const RDMol &asRDMol() const { return *dp_mol; }

  virtual ~ROMol();

  //! @}
  //! \name Atoms
  //! @{

  //! returns our number of atoms
  unsigned int getNumAtoms() const;
  unsigned int getNumAtoms(bool onlyExplicit) const;
  //! returns our number of heavy atoms (atomic number > 1)
  unsigned int getNumHeavyAtoms() const;
  //! returns a pointer to a particular Atom
  Atom *getAtomWithIdx(unsigned int idx);
  //! \overload
  const Atom *getAtomWithIdx(unsigned int idx) const;
  //! \overload
  template <class U>
  Atom *getAtomWithIdx(const U idx) {
    return getAtomWithIdx(rdcast<unsigned int>(idx));
  }
  //! \overload
  template <class U>
  const Atom *getAtomWithIdx(const U idx) const {
    return getAtomWithIdx(rdcast<unsigned int>(idx));
  }
  //! returns the degree (number of neighbors) of an Atom in the graph
  unsigned int getAtomDegree(const Atom *at) const;
  //! @}

  //! \name Bonds
  //! @{

  //! returns our number of Bonds
  unsigned int getNumBonds(bool onlyHeavy = true) const;
  //! returns a pointer to a particular Bond
  Bond *getBondWithIdx(unsigned int idx);
  //! \overload
  const Bond *getBondWithIdx(unsigned int idx) const;
  //! \overload
  template <class U>
  Bond *getBondWithIdx(const U idx) {
    return getBondWithIdx(rdcast<unsigned int>(idx));
  }
  //! \overload
  template <class U>
  const Bond *getBondWithIdx(const U idx) const {
    return getBondWithIdx(rdcast<unsigned int>(idx));
  }
  //! returns a pointer to the bond between two atoms, Null on failure
  Bond *getBondBetweenAtoms(unsigned int idx1, unsigned int idx2);
  //! \overload
  const Bond *getBondBetweenAtoms(unsigned int idx1, unsigned int idx2) const;
  //! \overload
  template <class U, class V>
  Bond *getBondBetweenAtoms(const U idx1, const V idx2) {
    return getBondBetweenAtoms(rdcast<unsigned int>(idx1),
                               rdcast<unsigned int>(idx2));
  }
  //! \overload
  template <class U, class V>
  const Bond *getBondBetweenAtoms(const U idx1, const V idx2) const {
    return getBondBetweenAtoms(rdcast<unsigned int>(idx1),
                               rdcast<unsigned int>(idx2));
  }

  //! @}

  //! \name Topology
  //! @{

  //! returns a pointer to our RingInfo structure
  //! <b>Note:</b> the client should not delete this.
  RingInfo *getRingInfo() const;

  //! provides access to all neighbors around an Atom
  /*!
    \param at the atom whose neighbors we are looking for

    <b>Usage</b>
    \code
      ... mol is a const ROMol & ...
      ... atomPtr is a const Atom * ...
      ... requires #include <boost/range/iterator_range.hpp>
      for (const auto &nbri :
           boost::make_iterator_range(m.getAtomNeighbors(atomPtr))) {
        const auto &nbr = (*m)[nbri];
        // nbr is an atom pointer
      }

    \endcode
  */
  ADJ_ITER_PAIR getAtomNeighbors(Atom const *at) const;

  //! provides access to all Bond objects connected to an Atom
  /*!
    \param at the atom whose neighbors we are looking for

    <b>Usage</b>
    \code
      ... mol is a const ROMol & ...
      ... atomPtr is a const Atom * ...
      ... requires #include <boost/range/iterator_range.hpp>
      for (const auto &nbri :
           boost::make_iterator_range(m.getAtomBonds(atomPtr))) {
        const auto &nbr = (*m)[nbri];
        // nbr is a bond pointer
      }
    \endcode
    or, if you need a non-const Bond *:
    \code
      ... mol is a const ROMol & ...
      ... atomPtr is a const Atom * ...
      ... requires #include <boost/range/iterator_range.hpp>
      for (const auto &nbri :
           boost::make_iterator_range(m.getAtomBonds(atomPtr))) {
        auto nbr = (*m)[nbri];
        // nbr is a bond pointer
      }
    \endcode
  */
  OBOND_ITER_PAIR getAtomBonds(Atom const *at) const;

  //! returns an iterator pair for looping over all Atoms
  /*!

    <b>Usage</b>
    \code

      ROMol::VERTEX_ITER atBegin,atEnd;
      boost::tie(atBegin,atEnd) = mol.getVertices();
      while(atBegin!=atEnd){
        ATOM_SPTR at2=mol[*atBegin];
        ... do something with the Atom ...
        ++atBegin;
      }
    \endcode
  */
  ATOM_ITER_PAIR getVertices() const;
  //! returns an iterator pair for looping over all Bonds
  /*!

    <b>Usage</b>
    \code

      ROMol::EDGE_ITER firstB,lastB;
      boost::tie(firstB,lastB) = mol.getEdges();
      while(firstB!=lastB){
        BOND_SPTR bond = mol[*firstB];
        ... do something with the Bond ...
        ++firstB;
      }
    \endcode
  */
  BOND_ITER_PAIR getEdges() const;

  //! brief returns a pointer to our underlying BGL object
  /*!
      This can be useful if you need to call other BGL algorithms:

      Here's an example:
      \code
         ... mol is a const ROMol ...
         ... mapping is an INT_VECT ...
         mapping.resize(mol.getNumAtoms());
         const MolGraph &G_p = mol.getTopology();
         int res = boost::connected_components(G_p,&mapping[0]);
      \endcode
   */
  ROMol const &getTopology() const { return *this; }
  //! @}

  //! \name Iterators
  //! @{

  //! get an AtomIterator pointing at our first Atom
  AtomIterator beginAtoms();
  //! \overload
  ConstAtomIterator beginAtoms() const;
  //! get an AtomIterator pointing at the end of our Atoms
  AtomIterator endAtoms();
  //! \overload
  ConstAtomIterator endAtoms() const;
  //! get a BondIterator pointing at our first Bond
  BondIterator beginBonds();
  //! \overload
  ConstBondIterator beginBonds() const;
  //! get a BondIterator pointing at the end of our Bonds
  BondIterator endBonds();
  //! \overload
  ConstBondIterator endBonds() const;

  //! get an AtomIterator pointing at our first aromatic Atom
  AromaticAtomIterator beginAromaticAtoms();
  //! \overload
  ConstAromaticAtomIterator beginAromaticAtoms() const;
  //! get an AtomIterator pointing at the end of our Atoms
  AromaticAtomIterator endAromaticAtoms();
  //! \overload
  ConstAromaticAtomIterator endAromaticAtoms() const;

  //! get an AtomIterator pointing at our first hetero Atom
  HeteroatomIterator beginHeteros();
  //! \overload
  ConstHeteroatomIterator beginHeteros() const;
  //! get an AtomIterator pointing at the end of our Atoms
  HeteroatomIterator endHeteros();
  //! \overload
  ConstHeteroatomIterator endHeteros() const;

  //! get an AtomIterator pointing at our first Atom that matches \c query
  QueryAtomIterator beginQueryAtoms(QueryAtom const *query);
  //! \overload
  ConstQueryAtomIterator beginQueryAtoms(QueryAtom const *) const;
  //! get an AtomIterator pointing at the end of our Atoms
  QueryAtomIterator endQueryAtoms();
  //! \overload
  ConstQueryAtomIterator endQueryAtoms() const;

  bool hasQuery() const;

  //! get an AtomIterator pointing at our first Atom that matches \c query
  MatchingAtomIterator beginMatchingAtoms(bool (*query)(Atom *));
  //! \overload
  ConstMatchingAtomIterator beginMatchingAtoms(
      bool (*query)(const Atom *)) const;
  //! get an AtomIterator pointing at the end of our Atoms
  MatchingAtomIterator endMatchingAtoms();
  //! \overload
  ConstMatchingAtomIterator endMatchingAtoms() const;

  inline ConformerIterator beginConformers() {
    return dp_mol->getConformersCompat().begin();
  }
  inline ConformerIterator endConformers() {
    return dp_mol->getConformersCompat().end();
  }
  inline ConstConformerIterator beginConformers() const {
    return dp_mol->getConformersCompat().begin();
  }
  inline ConstConformerIterator endConformers() const {
    return dp_mol->getConformersCompat().end();
  }

  //! @}

  //! \name Properties
  //! @{

  //! clears all of our \c computed \c properties
  /*!
    <b>Notes:</b>
    - This must be allowed on a const ROMol for backward compatibility,
    but DO NOT call this while this molecule might be being modified in another
    thread, NOR while any properties on the same molecule or its atoms or bonds
    might be being read or written.
  */
  void clearComputedProps(bool includeRings = true) const;
  //! calculates any of our lazy \c properties
  /*!
    <b>Notes:</b>
       - this calls \c updatePropertyCache() on each of our Atoms and Bonds
  */
  void updatePropertyCache(bool strict = true);

  bool needsUpdatePropertyCache() const;
  void clearPropertyCache();

  //! @}

  //! \name Misc
  //! @{
  //! sends some debugging info to a stream
  void debugMol(std::ostream &str) const;
  //! @}

  Atom *operator[](const vertex_descriptor &v) { return getAtomWithIdx(v); }
  const Atom *operator[](const vertex_descriptor &v) const {
    return getAtomWithIdx(v);
  }

  Bond *operator[](const edge_descriptor &e) { return getBondWithIdx(*e); }
  const Bond *operator[](const edge_descriptor &e) const {
    return getBondWithIdx(*e);
  }

  //! Gets a reference to the groups of atoms with relative stereochemistry
  /*!
    Stereo groups are also called enhanced stereochemistry in the SDF/Mol3000
    file format.
  */
  const std::vector<StereoGroup> &getStereoGroups() const;

  //! Sets groups of atoms with relative stereochemistry
  /*!
    \param stereo_groups the new set of stereo groups. All will be replaced.

    Stereo groups are also called enhanced stereochemistry in the SDF/Mol3000
    file format. stereo_groups should be std::move()ed into this function.
  */
  void setStereoGroups(std::vector<StereoGroup> stereo_groups);


  //! returns a list with the names of our \c properties
  STR_VECT getPropList(bool includePrivate = true,
                       bool includeComputed = true) const;

  //! sets a \c property value
  /*!
    \param key the name under which the \c property should be stored.
    If a \c property is already stored under this name, it will be
    replaced.
    \param val the value to be stored
    \param computed (optional) allows the \c property to be flagged
    \c computed.

    <b>Notes:</b>
    - This must be allowed on a const ROMol for backward compatibility,
    but DO NOT call this while this molecule might be being modified in another
    thread, NOR while any properties on the same molecule or its atoms or bonds
    might be being read or written.
  */
  //! \overload
  template <typename T>
  void setProp(const std::string &key, T val, bool computed = false) const {
    dp_mol->setMolProp(PropToken(key), val, computed);
  }
  //! allows retrieval of a particular property value
  /*!
    \param key the name under which the \c property should be stored.
    If a \c property is already stored under this name, it will be
    replaced.
    \param res a reference to the storage location for the value.

    <b>Notes:</b>
    - if no \c property with name \c key exists, a KeyErrorException will be
    thrown.
    - the \c boost::lexical_cast machinery is used to attempt type
    conversions.
    If this fails, a \c boost::bad_lexical_cast exception will be thrown.
  */
  //! \overload
  template <typename T>
  void getProp(const std::string &key, T &res) const {
    PropToken token(key);
    if constexpr (std::is_same_v<T, STR_VECT>) {
      if (token == detail::computedPropNameToken) {
        dp_mol->getComputedPropList(res);
        return;
      }
    }
    bool found = dp_mol->getMolPropIfPresent(token, res);
    if (!found) {
      throw KeyErrorException(key);
    }
  }

  //! \overload
  template <typename T>
  T getProp(const std::string &key) const {
    PropToken token(key);
    if constexpr (std::is_same_v<T, STR_VECT>) {
      if (token == detail::computedPropNameToken) {
        STR_VECT res;
        dp_mol->getComputedPropList(res);
        return res;
      }
    }
    return dp_mol->getMolProp<T>(token);
  }

  //! returns whether or not we have a \c property with name \c key
  //!  and assigns the value if we do
  //! \overload
  template <typename T>
  bool getPropIfPresent(const std::string &key, T &res) const {
    PropToken token(key);
    if constexpr (std::is_same_v<T, STR_VECT>) {
      if (token == detail::computedPropNameToken) {
        dp_mol->getComputedPropList(res);
        return true;
      }
    }
    return dp_mol->getMolPropIfPresent(token, res);
  }

  bool hasProp(const std::string &key) const;

  //! clears the value of a \c property
  /*!
    <b>Notes:</b>
    - if no \c property with name \c key exists, a KeyErrorException
    will be thrown.
    - if the \c property is marked as \c computed, it will also be removed
    from our list of \c computedProperties

    <b>Notes:</b>
    - This must be allowed on a const ROMol for backward compatibility,
    but DO NOT call this while this molecule might be being modified in another
    thread, NOR while any properties on the same molecule or its atoms or bonds
    might be being read or written.
  */
  //! \overload
  void clearProp(const std::string &key) const;

  //! update the properties from another
  /*
    \param source    Source to update the properties from
    \param preserve  Existing If true keep existing data, else override from
    the source
  */

  void updateProps(const ROMol& source, bool preserveExisting = false);
  void updateProps(const RDProps& source, bool preserveExisting = false);

  [[noreturn]] Dict& getDict() {
    raiseNonImplementedFunction("GetDict");
  }
  [[noreturn]] const Dict& getDict() const {
    raiseNonImplementedFunction("GetDict");
  }

  //! Clear all properties.
  void clear();

  //! \name Bookmarks
  //! @{

  //! associates an Atom pointer with a bookmark
  void setAtomBookmark(Atom *at, int mark);
  //! associates an Atom pointer with a bookmark
  void replaceAtomBookmark(Atom *at, int mark);
  //! returns the first Atom associated with the \c bookmark provided
  Atom *getAtomWithBookmark(int mark);
  //! returns the Atom associated with the \c bookmark provided
  //! a check is made to ensure it is the only atom with that bookmark
  Atom *getUniqueAtomWithBookmark(int mark);
  //! returns all Atoms associated with the \c bookmark provided
  ATOM_PTR_LIST &getAllAtomsWithBookmark(int mark);
  //! removes a \c bookmark from our collection
  void clearAtomBookmark(int mark);
  //! removes a particular Atom from the list associated with the \c bookmark
  void clearAtomBookmark(int mark, const Atom *atom);
  //! blows out all atomic \c bookmarks
  void clearAllAtomBookmarks();
  //! queries whether or not any atoms are associated with a \c bookmark
  bool hasAtomBookmark(int mark) const;
  //! returns a pointer to all of our atom \c bookmarks
  ATOM_BOOKMARK_MAP *getAtomBookmarks();

  //! associates a Bond pointer with a bookmark
  void setBondBookmark(Bond *bond, int mark);
  //! returns the first Bond associated with the \c bookmark provided
  Bond *getBondWithBookmark(int mark);
  //! returns the Bond associated with the \c bookmark provided
  //! a check is made to ensure it is the only bond with that bookmark
  Bond *getUniqueBondWithBookmark(int mark);
  //! returns all bonds associated with the \c bookmark provided
  BOND_PTR_LIST &getAllBondsWithBookmark(int mark);
  //! removes a \c bookmark from our collection
  void clearBondBookmark(int mark);
  //! removes a particular Bond from the list associated with the \c bookmark
  void clearBondBookmark(int mark, const Bond *bond);
  //! blows out all bond \c bookmarks
  void clearAllBondBookmarks();
  //! queries whether or not any bonds are associated with a \c bookmark
  bool hasBondBookmark(int mark) const;
  //! returns a pointer to all of our bond \c bookmarks
  BOND_BOOKMARK_MAP *getBondBookmarks();

  //! @}

  //! \name Conformers
  //! @{

  //! return the conformer with a specified ID
  //! if the ID is negative the first conformation will be returned
  const Conformer &getConformer(int id = -1) const;

  //! return the conformer with a specified ID
  //! if the ID is negative the first conformation will be returned
  Conformer &getConformer(int id = -1);

  //! Delete the conformation with the specified ID
  void removeConformer(unsigned int id);

  //! Clear all the conformations on the molecule
  void clearConformers();

  //! Add a new conformation to the molecule
  /*!
    \param conf - conformation to be added to the molecule, this molecule takes
    ownership
                  of the conformer
    \param assignId - a unique ID will be assigned to the conformation if
    true
                      otherwise it is assumed that the conformation already has
    an (unique) ID set
  */
  unsigned int addConformer(Conformer *conf, bool assignId = false);

  unsigned int getNumConformers() const;

  friend RDKIT_GRAPHMOL_EXPORT std::vector<SubstanceGroup> &getSubstanceGroups(
      ROMol &);
  friend RDKIT_GRAPHMOL_EXPORT const std::vector<SubstanceGroup> &
  getSubstanceGroups(const ROMol &);
  void clearSubstanceGroups();


  #ifdef RDK_USE_BOOST_SERIALIZATION
    //! \name boost::serialization support
    //! @{
    template <class Archive>
    void save(Archive &ar, const unsigned int version) const;
    template <class Archive>
    void load(Archive &ar, const unsigned int version);
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    //! @}
  #endif

 protected:
#ifndef WIN32
 private:
#endif

  //! adds an Atom to our collection
  /*!
    \param atom          pointer to the Atom to add
    \param updateLabel   (optional) if this is true, the new Atom will be
                         our \c activeAtom
    \param takeOwnership (optional) if this is true, we take ownership of \c
    atom
                         instead of copying it.

    \return the index of the new atom
  */
  unsigned int addAtom(Atom *atom, bool updateLabel = true,
                       bool takeOwnership = false);
  //! adds a Bond to our collection
  /*!
    \param bond          pointer to the Bond to add
    \param takeOwnership (optional) if this is true, we take ownership of \c
    bond
                         instead of copying it.

    \return the new number of bonds
  */
  unsigned int addBond(Bond *bond, bool takeOwnership = false);

  //! adds a Bond to our collection
  /*!
    \param bond          pointer to the Bond to add

    \return the new number of bonds

    <b>Note:</b> since this is using a smart pointer, we don't need to worry
    about
    issues of ownership.
  */
  void initFromOther(const ROMol &other, bool quickCopy, int confId);
};

using MolGraph = ROMol;

typedef std::vector<ROMol> MOL_VECT;
typedef boost::shared_ptr<ROMol> ROMOL_SPTR;
typedef std::vector<ROMol *> MOL_PTR_VECT;
typedef std::vector<ROMOL_SPTR> MOL_SPTR_VECT;

typedef MOL_PTR_VECT::const_iterator MOL_PTR_VECT_CI;
typedef MOL_PTR_VECT::iterator MOL_PTR_VECT_I;
//! Number of vertices for boost VF2
RDKIT_GRAPHMOL_EXPORT uint32_t num_vertices(const ROMol &mol);
//! Number of edges for boost VF2
RDKIT_GRAPHMOL_EXPORT uint32_t num_edges(const ROMol &mol);

};  // namespace RDKit
#endif
