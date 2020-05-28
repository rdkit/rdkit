//  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <RDGeneral/export.h>
#ifndef __RD_FILTER_MATCHER_H__
#define __RD_FILTER_MATCHER_H__
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include "FilterMatcherBase.h"
#include <GraphMol/MolPickler.h>

namespace RDKit {

namespace {
std::string getArgName(const boost::shared_ptr<FilterMatcherBase> &arg) {
  if (arg.get()) return arg->getName();
  return "<nullmatcher>";
}
}  // namespace

namespace FilterMatchOps {
class RDKIT_FILTERCATALOG_EXPORT And : public FilterMatcherBase {
  boost::shared_ptr<FilterMatcherBase> arg1;
  boost::shared_ptr<FilterMatcherBase> arg2;

 public:
  // !Default Constructor for serialization
  And() : FilterMatcherBase("And"), arg1(), arg2() {}

  //! Constructs an Ander
  //! True if arg1 and arg2 FilterMatchers are true

  And(const FilterMatcherBase &arg1, const FilterMatcherBase &arg2)
      : FilterMatcherBase("And"), arg1(arg1.copy()), arg2(arg2.copy()) {}

  And(const boost::shared_ptr<FilterMatcherBase> &arg1,
      const boost::shared_ptr<FilterMatcherBase> &arg2)
      : FilterMatcherBase("And"), arg1(arg1), arg2(arg2) {}

  And(const And &rhs)
      : FilterMatcherBase(rhs), arg1(rhs.arg1), arg2(rhs.arg2) {}

  virtual std::string getName() const {
    return "(" + getArgName(arg1) + " " + FilterMatcherBase::getName() + " " +
           getArgName(arg2) + ")";
  }

  bool isValid() const {
    return arg1.get() && arg2.get() && arg1->isValid() && arg2->isValid();
  }

  bool hasMatch(const ROMol &mol) const {
    PRECONDITION(isValid(),
                 "FilterMatchOps::And is not valid, null arg1 or arg2");
    return arg1->hasMatch(mol) && arg2->hasMatch(mol);
  }

  bool getMatches(const ROMol &mol, std::vector<FilterMatch> &matchVect) const {
    PRECONDITION(isValid(),
                 "FilterMatchOps::And is not valid, null arg1 or arg2");
    std::vector<FilterMatch> matches;
    if (arg1->getMatches(mol, matches) && arg2->getMatches(mol, matches)) {
      matchVect = matches;
      return true;
    }
    return false;
  }

  boost::shared_ptr<FilterMatcherBase> copy() const {
    return boost::shared_ptr<FilterMatcherBase>(new And(*this));
  }

 private:
#ifdef RDK_USE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    RDUNUSED_PARAM(version);
    ar &boost::serialization::base_object<FilterMatcherBase>(*this);

    ar &arg1;
    ar &arg2;
  }
#endif
};

class RDKIT_FILTERCATALOG_EXPORT Or : public FilterMatcherBase {
  boost::shared_ptr<FilterMatcherBase> arg1;
  boost::shared_ptr<FilterMatcherBase> arg2;

 public:
  // !Default Constructor for serialization
  Or() : FilterMatcherBase("Or"), arg1(), arg2() {}

  //! Constructs or Ander
  //! true if arg1 or arg2 are true
  Or(const FilterMatcherBase &arg1, const FilterMatcherBase &arg2)
      : FilterMatcherBase("Or"), arg1(arg1.copy()), arg2(arg2.copy()) {}

  Or(const boost::shared_ptr<FilterMatcherBase> &arg1,
     const boost::shared_ptr<FilterMatcherBase> &arg2)
      : FilterMatcherBase("Or"), arg1(arg1), arg2(arg2) {}

  Or(const Or &rhs) : FilterMatcherBase(rhs), arg1(rhs.arg1), arg2(rhs.arg2) {}

  virtual std::string getName() const {
    return "(" + getArgName(arg1) + " " + FilterMatcherBase::getName() + " " +
           getArgName(arg2) + ")";
  }

  bool isValid() const {
    return arg1.get() && arg2.get() && arg1->isValid() && arg2->isValid();
  }

  bool hasMatch(const ROMol &mol) const {
    PRECONDITION(isValid(), "Or is not valid, null arg1 or arg2");
    return arg1->hasMatch(mol) || arg2->hasMatch(mol);
  }

  bool getMatches(const ROMol &mol, std::vector<FilterMatch> &matchVect) const {
    PRECONDITION(isValid(),
                 "FilterMatchOps::Or is not valid, null arg1 or arg2");
    // we want both matches to run in order to accumulate all matches
    //  into matchVect, otherwise the or can be arbitrary...
    bool res1 = arg1->getMatches(mol, matchVect);
    bool res2 = arg2->getMatches(mol, matchVect);
    return res1 || res2;
  }

  boost::shared_ptr<FilterMatcherBase> copy() const {
    return boost::shared_ptr<FilterMatcherBase>(new Or(*this));
  }

#ifdef RDK_USE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    RDUNUSED_PARAM(version);
    ar &boost::serialization::base_object<FilterMatcherBase>(*this);
    ar &arg1;
    ar &arg2;
  }
#endif
};

class RDKIT_FILTERCATALOG_EXPORT Not : public FilterMatcherBase {
  boost::shared_ptr<FilterMatcherBase> arg1;

 public:
  // !Default Constructor for serialization
  Not() : FilterMatcherBase("Not"), arg1() {}

  //! Constructs a Noter
  //! true if arg1 is false (note, never returns matches
  //  from getMatches since a false internal match matches
  //  nothing!
  Not(const FilterMatcherBase &arg1)
      : FilterMatcherBase("Not"), arg1(arg1.copy()) {}

  Not(const boost::shared_ptr<FilterMatcherBase> &arg1)
      : FilterMatcherBase("Not"), arg1(arg1) {}

  Not(const Not &rhs) : FilterMatcherBase(rhs), arg1(rhs.arg1) {}

  virtual std::string getName() const {
    return "(" + FilterMatcherBase::getName() + " " + getArgName(arg1) + ")";
  }

  bool isValid() const { return arg1.get() && arg1->isValid(); }

  bool hasMatch(const ROMol &mol) const {
    PRECONDITION(isValid(), "FilterMatchOps::Not: arg1 is null");
    return !arg1->hasMatch(mol);
  }

  bool getMatches(const ROMol &mol, std::vector<FilterMatch> &) const {
    PRECONDITION(isValid(), "FilterMatchOps::Not: arg1 is null");
    // If we are a not, we really can't hold the match for
    //  this query since by definition it won't exist!
    std::vector<FilterMatch> matchVect;
    return !arg1->getMatches(mol, matchVect);
  }

  boost::shared_ptr<FilterMatcherBase> copy() const {
    return boost::shared_ptr<FilterMatcherBase>(new Not(*this));
  }

 private:
#ifdef RDK_USE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    RDUNUSED_PARAM(version);
    ar &boost::serialization::base_object<FilterMatcherBase>(*this);
    ar &arg1;
  }
#endif
};
}  // namespace FilterMatchOps

RDKIT_FILTERCATALOG_EXPORT extern const char *SMARTS_MATCH_NAME_DEFAULT;
class RDKIT_FILTERCATALOG_EXPORT SmartsMatcher : public FilterMatcherBase {
  ROMOL_SPTR d_pattern;
  unsigned int d_min_count{0};
  unsigned int d_max_count;

 public:
  //! Construct a SmartsMatcher
  SmartsMatcher(const std::string &name = SMARTS_MATCH_NAME_DEFAULT)
      : FilterMatcherBase(name),
        d_pattern(),
        
        d_max_count(UINT_MAX) {}

  //! Construct a SmartsMatcher from a query molecule
  /*
    \param pattern  query molecule used as the substructure search
    \param unsigned int minCount  minimum number of times the pattern needs to
    appear
    \param maxCount the maximum number of times the pattern should appear
    a value of UINT_MAX indicates the pattern can exist any number of times.
    [default UINT_MAX]

  */
  SmartsMatcher(const ROMol &pattern, unsigned int minCount = 1,
                unsigned int maxCount = UINT_MAX);

  //! Construct a SmartsMatcher
  /*
    \param name     name for the smarts pattern
    \param pattern  query molecule used as the substructure search
    \param unsigned int minCount  minimum number of times the pattern needs to
    appear
    \param maxCount the maximum number of times the pattern should appear
    a value of UINT_MAX indicates the pattern can exist any number of times.
    [default UINT_MAX]

  */

  SmartsMatcher(const std::string &name, const ROMol &pattern,
                unsigned int minCount = 1, unsigned int maxCount = UINT_MAX);

  //! Construct a SmartsMatcher from a smarts pattern
  /*
    \param name     name for the smarts pattern
    \param smarts   smarts pattern to use for the filter
    \param unsigned int minCount  minimum number of times the pattern needs to
    appear
    \param maxCount the maximum number of times the pattern should appear
    a value of UINT_MAX indicates the pattern can exist any number of times.
    [default UINT_MAX]
  */

  SmartsMatcher(const std::string &name, const std::string &smarts,
                unsigned int minCount = 1, unsigned int maxCount = UINT_MAX);

  //! Construct a SmartsMatcher from a shared_ptr
  /*
    \param name     name for the smarts pattern
    \param pattern  shared_ptr query molecule used as the substructure search
    \param unsigned int minCount  minimum number of times the pattern needs to
    appear
    \param maxCount the maximum number of times the pattern should appear
    a value of UINT_MAX indicates the pattern can exist any number of times.
    [default UINT_MAX]
  */

  SmartsMatcher(const std::string &name, ROMOL_SPTR onPattern,
                unsigned int minCount = 1, unsigned int maxCount = UINT_MAX);

  SmartsMatcher(const SmartsMatcher &rhs);

  //! Returns True if the Smarts pattern is valid
  bool isValid() const { return d_pattern.get(); }

  //! Return the shared_ptr to the underlying query molecule
  const ROMOL_SPTR &getPattern() const { return d_pattern; }
  //! Set the smarts pattern for the matcher
  void setPattern(const std::string &smarts);
  //! Set the query molecule for the matcher
  void setPattern(const ROMol &mol);
  //! Set the shared query molecule for the matcher
  void setPattern(const ROMOL_SPTR &pat) { d_pattern = pat; }

  //! Get the minimum match count for the pattern to be true
  unsigned int getMinCount() const { return d_min_count; }
  //! Set the minimum match count for the pattern to be true
  void setMinCount(unsigned int val) { d_min_count = val; }
  //! Get the maximum match count for the pattern to be true
  unsigned int getMaxCount() const { return d_max_count; }
  //! Set the maximum match count for the pattern to be true
  void setMaxCount(unsigned int val) { d_max_count = val; }

  virtual bool getMatches(const ROMol &mol,
                          std::vector<FilterMatch> &matchVect) const;
  virtual bool hasMatch(const ROMol &mol) const;
  virtual boost::shared_ptr<FilterMatcherBase> copy() const {
    return boost::shared_ptr<FilterMatcherBase>(new SmartsMatcher(*this));
  }

 private:
#ifdef RDK_USE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template <class Archive>
  void save(Archive &ar, const unsigned int version) const {
    RDUNUSED_PARAM(version);
    ar &boost::serialization::base_object<FilterMatcherBase>(*this);
    std::string res;
    MolPickler::pickleMol(*d_pattern.get(), res);
    ar &res;
    ar &d_min_count;
    ar &d_max_count;
  }
  template <class Archive>
  void load(Archive &ar, const unsigned int version) {
    ar &boost::serialization::base_object<FilterMatcherBase>(*this);
    {
      RDUNUSED_PARAM(version);
      std::string res;
      ar &res;
      d_pattern = boost::shared_ptr<ROMol>(new ROMol(res));
    }
    ar &d_min_count;
    ar &d_max_count;
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER();
#endif
};

// ------------------------------------------------------------------
// Syntactic sugar for the following style patterns
// Add exclusion patterns
//   using FilterMatchOps;
//   And(new SmartsMatcher(pat1),
//                  new Not(SmartsMatcher(pat2)))
// The exclusion match never adds any FilterMatches when getMatches
//  is called, the main intent is for it to be used with an
//  And construct, such as:
//    And(SmartsMatcher(..), ExclusionList(...))
//
//  which will return the SmartsMatcher FilterMatch only if no patterns
//    in the exclusion list are found.
class RDKIT_FILTERCATALOG_EXPORT ExclusionList : public FilterMatcherBase {
  std::vector<boost::shared_ptr<FilterMatcherBase>> d_offPatterns;

 public:
  ExclusionList() : FilterMatcherBase("Not any of"), d_offPatterns() {}

  //! Constructs an ExclusionList
  //! true if non of the FilterMatcherBases are true
  //! Syntactic sugar for
  //!  using FilterMatchOps;
  //!  And(Not(SmartsMatcher(pat1),
  //!                 And(Not(SmartsMatcher(pat2)),
  //!                                And(Not(Single...

  ExclusionList(
      const std::vector<boost::shared_ptr<FilterMatcherBase>> &offPatterns)
      : FilterMatcherBase("Not any of"), d_offPatterns(offPatterns) {}

  virtual std::string getName() const {
    std::string res;
    res = "(" + FilterMatcherBase::getName();
    for (size_t i = 0; i < d_offPatterns.size(); ++i) {
      res += " " + d_offPatterns[i]->getName();
    }
    res += ")";
    return res;
  }

  bool isValid() const {
    for (size_t i = 0; i < d_offPatterns.size(); ++i)
      if (!d_offPatterns[i]->isValid()) return false;
    return true;
  }

  void addPattern(const FilterMatcherBase &base) {
    PRECONDITION(base.isValid(), "Invalid FilterMatcherBase");
    d_offPatterns.push_back(base.copy());
  }

  void setExclusionPatterns(
      const std::vector<boost::shared_ptr<FilterMatcherBase>> &offPatterns) {
    d_offPatterns = offPatterns;
  }

  virtual bool getMatches(const ROMol &mol, std::vector<FilterMatch> &) const {
    PRECONDITION(isValid(),
                 "ExclusionList: one of the exclusion pattens is invalid");
    bool result = true;
    for (size_t i = 0; i < d_offPatterns.size() && result; ++i) {
      result &= !d_offPatterns[i]->hasMatch(mol);
    }

    return result;
  }

  virtual bool hasMatch(const ROMol &mol) const {
    PRECONDITION(isValid(),
                 "ExclusionList: one of the exclusion pattens is invalid");
    bool result = true;
    for (size_t i = 0; i < d_offPatterns.size() && result; ++i) {
      result &= !d_offPatterns[i]->hasMatch(mol);
    }

    return result;
  }

  virtual boost::shared_ptr<FilterMatcherBase> copy() const {
    return boost::shared_ptr<FilterMatcherBase>(new ExclusionList(*this));
  }

 private:
#ifdef RDK_USE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    RDUNUSED_PARAM(version);
    ar &boost::serialization::base_object<FilterMatcherBase>(*this);
    ar &d_offPatterns;
  }
#endif
};

class RDKIT_FILTERCATALOG_EXPORT FilterHierarchyMatcher
    : public FilterMatcherBase {
  std::vector<boost::shared_ptr<FilterHierarchyMatcher>> d_children;
  boost::shared_ptr<FilterMatcherBase> d_matcher;

 public:
  // !Default Constructor for serialization
  FilterHierarchyMatcher() : FilterMatcherBase(), d_matcher() {}
  //! Constructs a FilterHierarchyMatcher from a FilterMatchBase
  //!  A FilterHierarchyMatcher is a tree hierarchy where to
  //!  match a child node, one needs to match the parent first.
  //!  For each branch, the lowest nodes are returned when
  //!   getting the filter matches.
  /*
      \param matcher FilterMatcherBase to match this node against
  */
  FilterHierarchyMatcher(const FilterMatcherBase &matcher)
      : FilterMatcherBase(), d_matcher(matcher.copy()) {}

  //! Return the name for this node (from the underlying FilterMatcherBase)
  virtual std::string getName() const {
    if (d_matcher.get()) {
      return d_matcher->getName();
    }
    return "FilterMatcherHierarchy root";
  }

  //! returns true if this node has a valid matcher
  bool isValid() const { return d_matcher->isValid(); }

  //! Set a new FilterMatcherBase for this node
  /*
    \param matcher The new FilterMatcherBase
  */
  void setPattern(const FilterMatcherBase &matcher) {
    PRECONDITION(matcher.isValid(), "Adding invalid patterns is not allowed.");
    d_matcher = matcher.copy();
    PRECONDITION(getName() == d_matcher->getName(), "Opps");
  }

  //! add a FilterHierarchy as a child.
  //!  returns the FilterHierarchy pointer used in the tree (this is a
  //!   shallow copy of the original)
  /*
    \param hierarchy The new FilterHierarchyMatcher child for this node
  */
  boost::shared_ptr<FilterHierarchyMatcher> addChild(
      const FilterHierarchyMatcher &hierarchy) {
    PRECONDITION(hierarchy.d_matcher.get() && hierarchy.d_matcher->isValid(),
                 "Only one root node is allowed in a FilterHierarchyMatcher");

    d_children.push_back(boost::shared_ptr<FilterHierarchyMatcher>(
        new FilterHierarchyMatcher(hierarchy)));
    return d_children.back();
  }

  //! returns the FilterMatches against the given molecule
  /*
    \param mol The molecule to match against
    \param matches The vector of FilterMatch objects that match
  */
  virtual bool getMatches(const ROMol &mol,
                          std::vector<FilterMatch> &matches) const;

  //! Does this node match the molecule
  /*
    \param mol The molecule to match against
  */
  virtual bool hasMatch(const ROMol &mol) const {
    std::vector<FilterMatch> temp;
    return getMatches(mol, temp);
  }

  //! copys the FilterHierarchyMatcher into a FilterMatcherBase
  virtual boost::shared_ptr<FilterMatcherBase> copy() const {
    return boost::shared_ptr<FilterMatcherBase>(
        new FilterHierarchyMatcher(*this));
  }

 private:
#ifdef RDK_USE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    RDUNUSED_PARAM(version);
    ar &boost::serialization::base_object<FilterMatcherBase>(*this);
    ar &d_children;
    ar &d_matcher;
  }
#endif
};

#ifdef RDK_USE_BOOST_SERIALIZATION
// Register all known filter matcher types for serialization
template <class Archive>
void registerFilterMatcherTypes(Archive &ar) {
  ar.register_type(static_cast<FilterMatchOps::And *>(nullptr));
  ar.register_type(static_cast<FilterMatchOps::Or *>(nullptr));
  ar.register_type(static_cast<FilterMatchOps::Not *>(nullptr));
  ar.register_type(static_cast<SmartsMatcher *>(nullptr));
  ar.register_type(static_cast<ExclusionList *>(nullptr));
  ar.register_type(static_cast<FilterHierarchyMatcher *>(nullptr));
}
#endif
}  // namespace RDKit

#ifdef RDK_USE_BOOST_SERIALIZATION
BOOST_CLASS_VERSION(RDKit::SmartsMatcher, 1)
BOOST_CLASS_VERSION(RDKit::ExclusionList, 1)
BOOST_CLASS_VERSION(RDKit::FilterHierarchyMatcher, 1)
#endif

#endif
