//
//  Copyright (C) 2014 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include <GraphMol/new_canon.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/hanoiSort.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>

#include <iostream>
#include <vector>
#include <random>
#include <cstdlib>

using namespace RDKit;

int pcmp(const void *a, const void *b) {
  if ((*(int *)a) < (*(int *)b)) {
    return -1;
  } else if ((*(int *)a) > (*(int *)b)) {
    return 1;
  }
  return 0;
}
int icmp(int a, int b) {
  if (a < b) {
    return -1;
  } else if (a > b) {
    return 1;
  }
  return 0;
}

class int_compare_ftor {
  const int *dp_ints{nullptr};

 public:
  int_compare_ftor(){};
  int_compare_ftor(const int *ints) : dp_ints(ints){};
  int operator()(int i, int j) const {
    PRECONDITION(dp_ints, "no ints");
    unsigned int ivi = dp_ints[i];
    unsigned int ivj = dp_ints[j];
    if (ivi < ivj) {
      return -1;
    } else if (ivi > ivj) {
      return 1;
    } else {
      return 0;
    }
  }
};

void qs1(const std::vector<std::vector<int>> &vects) {
  BOOST_LOG(rdInfoLog) << "sorting (qsort) vectors" << std::endl;
  for (auto tv : vects) {
    int *data = &tv.front();
    qsort(data, tv.size(), sizeof(int), pcmp);
    for (unsigned int j = 1; j < tv.size(); ++j) {
      TEST_ASSERT(tv[j] >= tv[j - 1]);
    }
  }
  BOOST_LOG(rdInfoLog) << "done: " << vects.size() << std::endl;
}

void hs1(const std::vector<std::vector<int>> &vects) {
  BOOST_LOG(rdInfoLog) << "sorting (hanoi sort) vectors" << std::endl;
  for (const auto &vect : vects) {
    const int *data = &vect.front();
    int_compare_ftor icmp(data);
    int *indices = (int *)malloc(vect.size() * sizeof(int));
    for (unsigned int j = 0; j < vect.size(); ++j) {
      indices[j] = j;
    }
    int *count = (int *)malloc(vect.size() * sizeof(int));
    int *changed = (int *)malloc(vect.size() * sizeof(int));
    memset(changed, 1, vect.size() * sizeof(int));
    RDKit::hanoisort(indices, vect.size(), count, changed, icmp);
    for (unsigned int j = 1; j < vect.size(); ++j) {
      TEST_ASSERT(data[indices[j]] >= data[indices[j - 1]]);
    }
    free(count);
    free(indices);
    free(changed);
  }
  BOOST_LOG(rdInfoLog) << "done: " << vects.size() << std::endl;
}

void test1() {
  BOOST_LOG(rdInfoLog) << "Testing the hanoi sort" << std::endl;

  typedef boost::random::mersenne_twister<std::uint32_t, 32, 4, 2, 31,
                                          0x9908b0df, 11, 7, 0x9d2c5680, 15,
                                          0xefc60000, 18, 3346425566U>
      rng_type;
  typedef boost::uniform_int<> distrib_type;
  typedef boost::variate_generator<rng_type &, distrib_type> source_type;
  rng_type generator(42u);

  const unsigned int nVects = 500000;
  const unsigned int vectSize = 50;
  const unsigned int nClasses = 15;

  distrib_type dist(0, nClasses);
  source_type randomSource(generator, dist);

  BOOST_LOG(rdInfoLog) << "populating vectors" << std::endl;
  std::vector<std::vector<int>> vects(nVects);
  for (unsigned int i = 0; i < nVects; ++i) {
    vects[i] = std::vector<int>(vectSize);
    for (unsigned int j = 0; j < vectSize; ++j) {
      vects[i][j] = randomSource();
    }
  }

  // qs1(vects);
  hs1(vects);
  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
};

class atomcomparefunctor {
  Canon::canon_atom *d_atoms{nullptr};

 public:
  atomcomparefunctor(){};
  atomcomparefunctor(Canon::canon_atom *atoms) : d_atoms(atoms){};
  int operator()(int i, int j) const {
    PRECONDITION(d_atoms, "no atoms");
    unsigned int ivi, ivj;

    // always start with the current class:
    ivi = d_atoms[i].index;
    ivj = d_atoms[j].index;
    if (ivi < ivj) {
      return -1;
    } else if (ivi > ivj) {
      return 1;
    }

    ivi = d_atoms[i].atom->getAtomicNum();
    ivj = d_atoms[j].atom->getAtomicNum();
    if (ivi < ivj) {
      return -1;
    } else if (ivi > ivj) {
      return 1;
    }

    return 0;
  }
};
class atomcomparefunctor2 {
  Canon::canon_atom *d_atoms{nullptr};

 public:
  atomcomparefunctor2(){};
  atomcomparefunctor2(Canon::canon_atom *atoms) : d_atoms(atoms){};
  int operator()(int i, int j) const {
    PRECONDITION(d_atoms, "no atoms");
    unsigned int ivi, ivj;

    // always start with the current class:
    ivi = d_atoms[i].index;
    ivj = d_atoms[j].index;
    if (ivi < ivj) {
      return -1;
    } else if (ivi > ivj) {
      return 1;
    }

    // start by comparing degree
    ivi = d_atoms[i].atom->getDegree();
    ivj = d_atoms[j].atom->getDegree();
    if (ivi < ivj) {
      return -1;
    } else if (ivi > ivj) {
      return 1;
    }

    // move onto atomic number
    ivi = d_atoms[i].atom->getAtomicNum();
    ivj = d_atoms[j].atom->getAtomicNum();
    if (ivi < ivj) {
      return -1;
    } else if (ivi > ivj) {
      return 1;
    }

    return 0;
  }
};

void test2() {
  BOOST_LOG(rdInfoLog) << "Testing hanoi with a functor." << std::endl;
  // make sure that hanoi works with a functor and "molecule data"
  {
    std::string smi = "FC1C(Cl)C1C";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    std::vector<Canon::canon_atom> atoms(m->getNumAtoms());
    std::vector<int> indices(m->getNumAtoms());
    for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
      atoms[i].atom = m->getAtomWithIdx(i);
      atoms[i].index = 0;
      indices[i] = i;
    }
    atomcomparefunctor ftor(&atoms.front());

    int *data = &indices.front();
    int *count = (int *)malloc(atoms.size() * sizeof(int));
    int *changed = (int *)malloc(atoms.size() * sizeof(int));
    memset(changed, 1, atoms.size() * sizeof(int));
    RDKit::hanoisort(data, atoms.size(), count, changed, ftor);

    for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
      // std::cerr<<indices[i]<<" "<<" index: "<<atoms[indices[i]].index<<"
      // count: "<<count[indices[i]]<<std::endl;
      if (i > 0) {
        TEST_ASSERT(atoms[indices[i]].atom->getAtomicNum() >=
                    atoms[indices[i - 1]].atom->getAtomicNum());
        if (atoms[indices[i]].atom->getAtomicNum() !=
            atoms[indices[i - 1]].atom->getAtomicNum()) {
          TEST_ASSERT(count[indices[i]] != 0);
        } else {
          TEST_ASSERT(count[indices[i]] == 0);
        }
      } else {
        TEST_ASSERT(count[indices[i]] != 0);
      }
    }
    delete m;
    free(count);
    free(changed);
  }
  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
};

void test3() {
  BOOST_LOG(rdInfoLog) << "Testing basic partition refinement." << std::endl;
  // basic partition refinement
  {
    std::string smi = "FC1C(Cl)CCC1C";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    std::vector<Canon::canon_atom> atoms(m->getNumAtoms());
    initCanonAtoms(*m, atoms, true);
    atomcomparefunctor ftor(&atoms.front());

    RDKit::Canon::canon_atom *data = &atoms.front();
    int *count = (int *)malloc(atoms.size() * sizeof(int));
    int *order = (int *)malloc(atoms.size() * sizeof(int));
    int activeset;
    int *next = (int *)malloc(atoms.size() * sizeof(int));
    int *changed = (int *)malloc(atoms.size() * sizeof(int));
    memset(changed, 1, atoms.size() * sizeof(int));
    char *touched = (char *)malloc(atoms.size() * sizeof(char));
    memset(touched, 0, atoms.size() * sizeof(char));

    RDKit::Canon::CreateSinglePartition(atoms.size(), order, count, data);
    RDKit::Canon::ActivatePartitions(atoms.size(), order, count, activeset,
                                     next, changed);

    // std::cerr<<"----------------------------------"<<std::endl;
    // for(unsigned int i=0;i<m->getNumAtoms();++i){
    //   std::cerr<<i<<" "<<atoms[i].index<<" "<<count[i]<<" "<<next[i]<<"
    //   "<<order[i]<<std::endl;
    // }

    RDKit::Canon::RefinePartitions(*m, data, ftor, false, order, count,
                                   activeset, next, changed, touched);

    // std::cerr<<"----------------------------------"<<std::endl;
    // for(unsigned int i=0;i<m->getNumAtoms();++i){
    //   std::cerr<<i<<" "<<atoms[i].index<<" "<<count[i]<<" "<<next[i]<<"
    //   "<<order[i]<<std::endl;
    // }

    // std::cerr<<"----------------------------------"<<std::endl;
    // for(unsigned int i=0;i<m->getNumAtoms();++i){
    //   std::cerr<<order[i]<<" "<<atoms[order[i]].invar<<" index:
    //   "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
    // }

    TEST_ASSERT(order[0] == 1);
    TEST_ASSERT(order[6] == 0);
    TEST_ASSERT(order[7] == 3);
    TEST_ASSERT(count[order[0]] == 6);
    TEST_ASSERT(count[order[1]] == 0);
    TEST_ASSERT(count[order[6]] == 1);
    TEST_ASSERT(count[order[7]] == 1);

    delete m;
    free(count);
    free(order);
    free(next);
    free(changed);
    free(touched);
  }
  {
    // this time with smarter invariants
    std::string smi = "FC1C(Cl)CCC1C";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    std::vector<Canon::canon_atom> atoms(m->getNumAtoms());
    initCanonAtoms(*m, atoms, true);
    atomcomparefunctor2 ftor(&atoms.front());

    RDKit::Canon::canon_atom *data = &atoms.front();
    int *count = (int *)malloc(atoms.size() * sizeof(int));
    int *order = (int *)malloc(atoms.size() * sizeof(int));
    int activeset;
    int *next = (int *)malloc(atoms.size() * sizeof(int));
    int *changed = (int *)malloc(atoms.size() * sizeof(int));
    memset(changed, 1, atoms.size() * sizeof(int));
    char *touched = (char *)malloc(atoms.size() * sizeof(char));
    memset(touched, 0, atoms.size() * sizeof(char));

    RDKit::Canon::CreateSinglePartition(atoms.size(), order, count, data);
    RDKit::Canon::ActivatePartitions(atoms.size(), order, count, activeset,
                                     next, changed);

    RDKit::Canon::RefinePartitions(*m, data, ftor, false, order, count,
                                   activeset, next, changed, touched);

    // std::cerr<<"----------------------------------"<<std::endl;
    // for(unsigned int i=0;i<m->getNumAtoms();++i){
    //   std::cerr<<order[i]<<" "<<" index: "<<atoms[order[i]].index<<" count:
    //   "<<count[order[i]]<<std::endl;
    // }

    TEST_ASSERT(order[0] == 7);
    TEST_ASSERT(order[1] == 0);
    TEST_ASSERT(order[2] == 3);
    TEST_ASSERT(order[3] == 4);
    TEST_ASSERT(order[5] == 1);
    TEST_ASSERT(count[order[0]] == 1);
    TEST_ASSERT(count[order[1]] == 1);
    TEST_ASSERT(count[order[2]] == 1);
    TEST_ASSERT(count[order[3]] == 2);
    TEST_ASSERT(count[order[4]] == 0);
    TEST_ASSERT(count[order[5]] == 3);
    TEST_ASSERT(count[order[6]] == 0);
    delete m;
    free(count);
    free(order);
    free(next);
    free(changed);
    free(touched);
  }
  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
};

class atomcomparefunctor3 {
  Canon::canon_atom *dp_atoms{nullptr};
  const ROMol *dp_mol{nullptr};
  unsigned int getAtomNeighborhood(unsigned int i) const {
    unsigned int res = 0;
    const Atom *at = dp_mol->getAtomWithIdx(i);
    std::vector<unsigned int> nbrs(at->getDegree());
    unsigned int nbridx = 0;
    ROMol::OEDGE_ITER beg, end;
    boost::tie(beg, end) = dp_mol->getAtomBonds(at);
    while (beg != end) {
      const Bond *bond = (*dp_mol)[*beg];
      nbrs[nbridx] =
          static_cast<unsigned int>(100 * bond->getBondTypeAsDouble()) +
          dp_atoms[bond->getOtherAtomIdx(i)].index;
      ++beg;
      ++nbridx;
    }
    std::sort(nbrs.begin(), nbrs.end());
    for (nbridx = 0; nbridx < at->getDegree(); ++nbridx) {
      res += (nbridx + 1) * 1000 + nbrs[nbridx];
    }
    return res;
  }
  int basecomp(int i, int j) const {
    PRECONDITION(dp_atoms, "no atoms");
    unsigned int ivi, ivj;

    // always start with the current class:
    ivi = dp_atoms[i].index;
    ivj = dp_atoms[j].index;
    if (ivi < ivj) {
      return -1;
    } else if (ivi > ivj) {
      return 1;
    }

    // start by comparing degree
    ivi = dp_atoms[i].atom->getDegree();
    ivj = dp_atoms[j].atom->getDegree();
    if (ivi < ivj) {
      return -1;
    } else if (ivi > ivj) {
      return 1;
    }

    // move onto atomic number
    ivi = dp_atoms[i].atom->getAtomicNum();
    ivj = dp_atoms[j].atom->getAtomicNum();
    if (ivi < ivj) {
      return -1;
    } else if (ivi > ivj) {
      return 1;
    }

    return 0;
  }

 public:
  bool df_useNbrs{false};
  atomcomparefunctor3(){};
  atomcomparefunctor3(Canon::canon_atom *atoms, const ROMol &m)
      : dp_atoms(atoms), dp_mol(&m), df_useNbrs(false){};
  int operator()(int i, int j) const {
    PRECONDITION(dp_atoms, "no atoms");
    PRECONDITION(dp_mol, "no molecule");
    int v = basecomp(i, j);
    if (v) {
      return v;
    }
    unsigned int ivi, ivj;
    if (df_useNbrs) {
      ivi = dp_atoms[i].index + 1 + getAtomNeighborhood(i);
      ivj = dp_atoms[j].index + 1 + getAtomNeighborhood(j);
      // std::cerr<<"               "<<i<<"-"<<j<<": "<<ivi<<"
      // "<<ivj<<std::endl;
      if (ivi < ivj) {
        return -1;
      } else if (ivi > ivj) {
        return 1;
      }
    }
    return 0;
  }
};

void test4() {
  BOOST_LOG(rdInfoLog) << "Testing partition refinement with neighbors."
                       << std::endl;
  // partition refinement with neighbors
  {
    std::string smi = "FC1C(Cl)CCC1C";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    std::vector<Canon::canon_atom> atoms(m->getNumAtoms());
    initCanonAtoms(*m, atoms, true);
    atomcomparefunctor3 ftor(&atoms.front(), *m);
    RDKit::Canon::canon_atom *data = &atoms.front();
    int *count = (int *)malloc(atoms.size() * sizeof(int));
    int *order = (int *)malloc(atoms.size() * sizeof(int));
    int activeset;
    int *next = (int *)malloc(atoms.size() * sizeof(int));
    int *changed = (int *)malloc(atoms.size() * sizeof(int));
    memset(changed, 1, atoms.size() * sizeof(int));
    char *touched = (char *)malloc(atoms.size() * sizeof(char));
    memset(touched, 0, atoms.size() * sizeof(char));

    RDKit::Canon::CreateSinglePartition(atoms.size(), order, count, data);
    RDKit::Canon::ActivatePartitions(atoms.size(), order, count, activeset,
                                     next, changed);
    // std::cerr<<"1----------------------------------"<<std::endl;
    //  for(unsigned int i=0;i<m->getNumAtoms();++i){
    //    std::cerr<<order[i]<<" "<<" index: "<<atoms[order[i]].index<<" count:
    //    "<<count[order[i]]<<" next: "<<next[order[i]]<<" changed:
    //    "<<changed[order[i]]<<std::endl;
    // }
    RDKit::Canon::RefinePartitions(*m, data, ftor, false, order, count,
                                   activeset, next, changed, touched);

    // std::cerr<<"2----------------------------------"<<std::endl;
    //  for(unsigned int i=0;i<m->getNumAtoms();++i){
    //    std::cerr<<order[i]<<" "<<" index: "<<atoms[order[i]].index<<" count:
    //    "<<count[order[i]]<<" next: "<<next[order[i]]<<" changed:
    //    "<<changed[order[i]]<<std::endl;
    // }

    // std::cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
    ftor.df_useNbrs = true;
    RDKit::Canon::ActivatePartitions(atoms.size(), order, count, activeset,
                                     next, changed);
    // std::cerr<<"3----------------------------------"<<std::endl;
    // for(unsigned int i=0;i<m->getNumAtoms();++i){
    //    std::cerr<<order[i]<<" "<<" index: "<<atoms[order[i]].index<<" count:
    //    "<<count[order[i]]<<" next: "<<next[order[i]]<<" changed:
    //    "<<changed[order[i]]<<std::endl;
    // }
    RDKit::Canon::RefinePartitions(*m, data, ftor, true, order, count,
                                   activeset, next, changed, touched);

    // std::cerr<<"----------------------------------"<<std::endl;
    for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
      // std::cerr<<order[i]<<" "<<" index: "<<atoms[order[i]].index<<" count:
      // "<<count[order[i]]<<std::endl;
      TEST_ASSERT(count[order[i]] == 1);
      if (i > 0) {
        TEST_ASSERT(ftor(order[i], order[i - 1]) >= 0);
      }
    }
    delete m;
    free(count);
    free(order);
    free(next);
    free(changed);
    free(touched);
  }

  {
    std::string smi = "FC1C(CO)CCC1CC";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    std::vector<Canon::canon_atom> atoms(m->getNumAtoms());
    initCanonAtoms(*m, atoms, true);
    atomcomparefunctor3 ftor(&atoms.front(), *m);

    RDKit::Canon::canon_atom *data = &atoms.front();
    int *count = (int *)malloc(atoms.size() * sizeof(int));
    int *order = (int *)malloc(atoms.size() * sizeof(int));
    int activeset;
    int *next = (int *)malloc(atoms.size() * sizeof(int));
    int *changed = (int *)malloc(atoms.size() * sizeof(int));
    memset(changed, 1, atoms.size() * sizeof(int));
    char *touched = (char *)malloc(atoms.size() * sizeof(char));
    memset(touched, 0, atoms.size() * sizeof(char));

    RDKit::Canon::CreateSinglePartition(atoms.size(), order, count, data);
    RDKit::Canon::ActivatePartitions(atoms.size(), order, count, activeset,
                                     next, changed);

    RDKit::Canon::RefinePartitions(*m, data, ftor, false, order, count,
                                   activeset, next, changed, touched);
    // std::cerr<<"----------------------------------"<<std::endl;
    // for(unsigned int i=0;i<m->getNumAtoms();++i){
    //   std::cerr<<order[i]<<" "<<" index: "<<atoms[order[i]].index<<" count:
    //   "<<count[order[i]]<<std::endl;
    // }

    ftor.df_useNbrs = true;
    RDKit::Canon::ActivatePartitions(atoms.size(), order, count, activeset,
                                     next, changed);
    RDKit::Canon::RefinePartitions(*m, data, ftor, true, order, count,
                                   activeset, next, changed, touched);

    // std::cerr<<"----------------------------------"<<std::endl;
    for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
      // std::cerr<<order[i]<<" "<<" index: "<<atoms[order[i]].index<<" count:
      // "<<count[order[i]]<<std::endl;
      TEST_ASSERT(count[order[i]] == 1);
      if (i > 0) {
        // std::cerr<<"  ftor: "<<ftor(order[i],order[i-1])<<std::endl;
        TEST_ASSERT(ftor(order[i], order[i - 1]) >= 0);
      }
    }
    delete m;
    free(count);
    free(order);
    free(next);
    free(changed);
    free(touched);
  }

  {
    std::string smi = "FC1C(CC)CCC1CC";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    std::vector<Canon::canon_atom> atoms(m->getNumAtoms());
    initCanonAtoms(*m, atoms, true);
    atomcomparefunctor3 ftor(&atoms.front(), *m);

    RDKit::Canon::canon_atom *data = &atoms.front();
    int *count = (int *)malloc(atoms.size() * sizeof(int));
    int *order = (int *)malloc(atoms.size() * sizeof(int));
    int activeset;
    int *next = (int *)malloc(atoms.size() * sizeof(int));
    int *changed = (int *)malloc(atoms.size() * sizeof(int));
    memset(changed, 1, atoms.size() * sizeof(int));
    char *touched = (char *)malloc(atoms.size() * sizeof(char));
    memset(touched, 0, atoms.size() * sizeof(char));

    RDKit::Canon::CreateSinglePartition(atoms.size(), order, count, data);
    RDKit::Canon::ActivatePartitions(atoms.size(), order, count, activeset,
                                     next, changed);

    // std::cerr<<"----------------------------------"<<std::endl;
    // for(unsigned int i=0;i<m->getNumAtoms();++i){
    //   std::cerr<<order[i]<<" "<<atoms[order[i]].invar<<" index:
    //   "<<atoms[order[i]].index<<std::endl;
    // }

    RDKit::Canon::RefinePartitions(*m, data, ftor, false, order, count,
                                   activeset, next, changed, touched);

    // std::cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
    ftor.df_useNbrs = true;

    RDKit::Canon::ActivatePartitions(atoms.size(), order, count, activeset,
                                     next, changed);
    RDKit::Canon::RefinePartitions(*m, data, ftor, true, order, count,
                                   activeset, next, changed, touched);
    // std::cerr<<"----------------------------------"<<std::endl;

    for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
      // std::cerr<<order[i]<<" "<<" index: "<<atoms[order[i]].index<<" count:
      // "<<count[order[i]]<<std::endl;
      if (i > 0) {
        // std::cerr<<"  ftor: "<<ftor(order[i],order[i-1])<<std::endl;
        TEST_ASSERT(ftor(order[i], order[i - 1]) >= 0);
      }
    }

    // here we can't manage to get everything  unique
    TEST_ASSERT(order[0] == 4 && count[4] == 2);
    TEST_ASSERT(order[1] == 9 && count[9] == 0);
    TEST_ASSERT(order[2] == 0 && count[0] == 1);
    TEST_ASSERT(order[3] == 3 && count[3] == 2);
    TEST_ASSERT(order[4] == 8 && count[8] == 0);
    TEST_ASSERT(order[5] == 5 && count[5] == 2);
    TEST_ASSERT(order[6] == 6 && count[6] == 0);
    TEST_ASSERT(order[7] == 2 && count[2] == 2);
    TEST_ASSERT(order[8] == 7 && count[7] == 0);
    TEST_ASSERT(order[9] == 1 && count[1] == 1);

    delete m;
    free(count);
    free(order);
    free(next);
    free(changed);
    free(touched);
  }

  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
};

void test5() {
  BOOST_LOG(rdInfoLog) << "testing canonicalization via tie breaking."
                       << std::endl;
  // canonicalization via tie breaking
  {
    std::string smi = "FC1C(CC)CCC1CC";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    std::vector<Canon::canon_atom> atoms(m->getNumAtoms());
    initCanonAtoms(*m, atoms, true);
    atomcomparefunctor3 ftor(&atoms.front(), *m);

    RDKit::Canon::canon_atom *data = &atoms.front();
    int *count = (int *)malloc(atoms.size() * sizeof(int));
    int *order = (int *)malloc(atoms.size() * sizeof(int));
    int activeset;
    int *next = (int *)malloc(atoms.size() * sizeof(int));
    int *changed = (int *)malloc(atoms.size() * sizeof(int));
    memset(changed, 1, atoms.size() * sizeof(int));
    char *touched = (char *)malloc(atoms.size() * sizeof(char));
    memset(touched, 0, atoms.size() * sizeof(char));

    RDKit::Canon::CreateSinglePartition(atoms.size(), order, count, data);
    RDKit::Canon::ActivatePartitions(atoms.size(), order, count, activeset,
                                     next, changed);

    // std::cerr<<"----------------------------------"<<std::endl;
    // for(unsigned int i=0;i<m->getNumAtoms();++i){
    //   std::cerr<<order[i]<<" "<<atoms[order[i]].invar<<" index:
    //   "<<atoms[order[i]].index<<std::endl;
    // }

    RDKit::Canon::RefinePartitions(*m, data, ftor, false, order, count,
                                   activeset, next, changed, touched);

    // std::cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
    ftor.df_useNbrs = true;
    RDKit::Canon::ActivatePartitions(atoms.size(), order, count, activeset,
                                     next, changed);
    RDKit::Canon::RefinePartitions(*m, data, ftor, true, order, count,
                                   activeset, next, changed, touched);

    // std::cerr<<"----------------------------------"<<std::endl;
    // for(unsigned int i=0;i<m->getNumAtoms();++i){
    //   std::cerr<<order[i]<<" "<<" index: "<<atoms[order[i]].index<<" count:
    //   "<<count[order[i]]<<std::endl;
    // }

    // here we can't manage to get everything  unique
    TEST_ASSERT(order[0] == 4 && count[4] == 2);
    TEST_ASSERT(order[1] == 9 && count[9] == 0);
    TEST_ASSERT(order[2] == 0 && count[0] == 1);
    TEST_ASSERT(order[3] == 3 && count[3] == 2);
    TEST_ASSERT(order[4] == 8 && count[8] == 0);
    TEST_ASSERT(order[5] == 5 && count[5] == 2);
    TEST_ASSERT(order[6] == 6 && count[6] == 0);
    TEST_ASSERT(order[7] == 2 && count[2] == 2);
    TEST_ASSERT(order[8] == 7 && count[7] == 0);
    TEST_ASSERT(order[9] == 1 && count[1] == 1);

    RDKit::Canon::BreakTies(*m, data, ftor, true, order, count, activeset, next,
                            changed, touched);
    for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
      // std::cerr<<order[i]<<" "<<" index: "<<atoms[order[i]].index<<" count:
      // "<<count[order[i]]<<std::endl;
      TEST_ASSERT(count[order[i]] == 1);
    }
    delete m;
    free(count);
    free(order);
    free(next);
    free(changed);
    free(touched);
  }
  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
};

void test6() {
  BOOST_LOG(rdInfoLog) << "testing canonicalization using the wrapper."
                       << std::endl;
// canonicalization using the wrapper
#if 1
  {
    std::string smi = "FC1C(CC)CCC1CC";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);

    std::vector<unsigned int> atomRanks;
    RDKit::Canon::rankMolAtoms(*m, atomRanks);
    boost::dynamic_bitset<> seen(m->getNumAtoms());
    for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
      TEST_ASSERT(!seen[atomRanks[i]]);
      seen.set(atomRanks[i], 1);
    }
    // std::copy(atomRanks.begin(),atomRanks.end(),std::ostream_iterator<unsigned
    // int>(std::cerr," "));
    // std::cerr<<std::endl;
    TEST_ASSERT(atomRanks[0] == 2);
    TEST_ASSERT(atomRanks[1] == 9);
    TEST_ASSERT(atomRanks[2] == 7);
    TEST_ASSERT(atomRanks[3] == 3);
    TEST_ASSERT(atomRanks[4] == 0);
    TEST_ASSERT(atomRanks[5] == 5);
    TEST_ASSERT(atomRanks[6] == 6);
    TEST_ASSERT(atomRanks[7] == 8);
    TEST_ASSERT(atomRanks[8] == 4);
    TEST_ASSERT(atomRanks[9] == 1);
    delete m;
  }

  {
    std::string smi = "CC[C@@H]1CCC[C@@H](O1)C(=O)O";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);

    std::vector<unsigned int> atomRanks;
    RDKit::Canon::rankMolAtoms(*m, atomRanks);
    boost::dynamic_bitset<> seen(m->getNumAtoms());
    for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
      // std::cerr<<i<<" "<<atomRanks[i]<<std::endl;
      TEST_ASSERT(!seen[atomRanks[i]]);
      seen.set(atomRanks[i], 1);
    }

    // for(unsigned int ii=0;ii<atomRanks.size();++ii){
    //    std::cerr<<ii<<":"<<atomRanks[ii]<<std::endl;
    // }
    TEST_ASSERT(atomRanks[0] == 0);
    TEST_ASSERT(atomRanks[1] == 4);
    TEST_ASSERT(atomRanks[2] == 9);
    TEST_ASSERT(atomRanks[3] == 5);
    TEST_ASSERT(atomRanks[4] == 3);
    TEST_ASSERT(atomRanks[5] == 6);
    TEST_ASSERT(atomRanks[6] == 10);
    TEST_ASSERT(atomRanks[7] == 7);
    TEST_ASSERT(atomRanks[8] == 8);
    TEST_ASSERT(atomRanks[9] == 1);
    TEST_ASSERT(atomRanks[10] == 2);

    delete m;
  }

  {
    std::string smi =
        "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](CCCN=C(N)N)C(=O)N[C@@H](Cc3c["
        "nH]c4ccccc34)C(=O)OCc5ccccc5";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);

    std::vector<unsigned int> atomRanks;
    RDKit::Canon::rankMolAtoms(*m, atomRanks);
    boost::dynamic_bitset<> seen(m->getNumAtoms());
    for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
      // std::cerr<<i<<" "<<atomRanks[i]<<std::endl;
      TEST_ASSERT(!seen[atomRanks[i]]);
      seen.set(atomRanks[i], 1);
    }

    // for(unsigned int ii=0;ii<atomRanks.size();++ii){
    //   std::cerr<<ii<<":"<<atomRanks[ii]<<std::endl;
    // }
    delete m;
  }
#endif
  {
    std::string smi = "BrC=C1CCC(C(=O)O1)c2cccc3ccccc23";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);

    std::vector<unsigned int> atomRanks;
    RDKit::Canon::rankMolAtoms(*m, atomRanks);
    boost::dynamic_bitset<> seen(m->getNumAtoms());
    for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
      // std::cerr<<i<<" "<<atomRanks[i]<<std::endl;
      TEST_ASSERT(!seen[atomRanks[i]]);
      seen.set(atomRanks[i], 1);
    }

    // for(unsigned int ii=0;ii<atomRanks.size();++ii){
    //   std::cerr<<ii<<":"<<atomRanks[ii]<<std::endl;
    // }
    delete m;
  }

  {
    std::string smi = "CC12CCCC1CCCC2";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);

    // start w/o tie breaking here; we shouldn't need it.
    std::vector<unsigned int> atomRanks;
    RDKit::Canon::rankMolAtoms(*m, atomRanks, false);
    boost::dynamic_bitset<> seen(m->getNumAtoms());
    for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
      //      std::cerr<<"      "<<i<<" "<<atomRanks[i]<<std::endl;
      TEST_ASSERT(!seen[atomRanks[i]]);
      seen.set(atomRanks[i], 1);
    }
    delete m;
  }

  {
    std::string smi = "CC12CCCC1C1CCC3CC(O)CCC3(C)C1CC2";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);

    // start w/o tie breaking here; we shouldn't need it.
    std::vector<unsigned int> atomRanks;
    RDKit::Canon::rankMolAtoms(*m, atomRanks, false);
    boost::dynamic_bitset<> seen(m->getNumAtoms());
    for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
      TEST_ASSERT(!seen[atomRanks[i]]);
      seen.set(atomRanks[i], 1);
    }
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
};

namespace {

ROMol *_renumber(const ROMol *m, std::vector<unsigned int> &nVect,
                 const std::string & /*inSmiles*/) {
  ROMol *nm = MolOps::renumberAtoms(*m, nVect);
  TEST_ASSERT(nm);
  TEST_ASSERT(nm->getNumAtoms() == m->getNumAtoms());
  TEST_ASSERT(nm->getNumBonds() == m->getNumBonds());
  // MolOps::assignStereochemistry(*nm, true, true);
  // for (unsigned int ii = 0; ii < nm->getNumAtoms(); ++ii) {
  //   if (nm->getAtomWithIdx(ii)->hasProp("_CIPCode")) {
  //     TEST_ASSERT(m->getAtomWithIdx(nVect[ii])->hasProp("_CIPCode"));
  //     std::string ocip =
  //         m->getAtomWithIdx(nVect[ii])->getProp<std::string>("_CIPCode");
  //     std::string ncip =
  //         nm->getAtomWithIdx(ii)->getProp<std::string>("_CIPCode");
  //     if (ocip != ncip) {
  //       std::cerr << "  cip mismatch: " << inSmiles << std::endl;
  //       std::cerr << "      " << nVect[ii] << ": " << ocip << " -> " << ii
  //                 << ": " << ncip << std::endl;
  //       std::cerr << "      " << MolToSmiles(*nm, true) << std::endl;
  //     }
  //     TEST_ASSERT(ocip == ncip);
  //   }
  // }
  return nm;
}

void _renumberTest(const ROMol *m, std::string inSmiles,
                   unsigned int numRenumbers) {
  PRECONDITION(m, "no molecule");
  //    std::cerr<<">>>>>>>>>>>>>>>>>>>>>>>>>>>"<<std::endl;
  std::string osmi = MolToSmiles(*m, true);
  std::vector<unsigned int> idxV(m->getNumAtoms());
  for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
    idxV[i] = i;
  }

  for (unsigned int i = 0; i < numRenumbers; ++i) {
    //      std::cerr<<"---------------------------------------------------"<<std::endl;
    std::vector<unsigned int> nVect(idxV);
    std::shuffle(nVect.begin(), nVect.end(), std::mt19937(0xf00d));
    //      for(unsigned int j=0;j<m->getNumAtoms();++j){
    //        std::cerr<<"Renumber: "<<nVect[j]<<"->"<<j<<std::endl;
    //      }
    ROMol *nm = _renumber(m, nVect, inSmiles);
    nm->setProp(common_properties::_StereochemDone, 1);

    std::string smi = MolToSmiles(*nm, true);
    if (smi != osmi) {
      std::cerr << "  input: " << inSmiles << ", Renumbering round: " << i
                << std::endl;
      std::cerr << osmi << std::endl;
      std::cerr << smi << std::endl;
      m->setProp("_Name", "orig");
      std::cerr << MolToMolBlock(*m) << std::endl;
      nm->setProp("_Name", "renumber");
      std::cerr << MolToMolBlock(*nm) << std::endl;
      for (unsigned int j = 0; j < m->getNumAtoms(); ++j) {
        std::cerr << "Renumber: " << nVect[j] << "->" << j << std::endl;
      }
    }
    delete nm;
    TEST_ASSERT(smi == osmi);
  }
}

void _renumberTest2(const ROMol *m, std::string inSmiles,
                    unsigned int numRenumbers) {
  PRECONDITION(m, "no molecule");

  unsigned int nAtoms = m->getNumAtoms();
  std::vector<unsigned int> idxV(m->getNumAtoms());
  for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
    idxV[i] = i;
  }

  for (unsigned int i = 0; i < numRenumbers; ++i) {
    std::vector<unsigned int> nVect(idxV);
    std::shuffle(nVect.begin(), nVect.end(), std::mt19937(0xf00d));

    ROMol *nm = _renumber(m, nVect, inSmiles);

    UINT_VECT ranks(nAtoms);
    Canon::rankMolAtoms(*nm, ranks, true);
    char *ranksSet = (char *)malloc(nAtoms * sizeof(char));
    memset(ranksSet, 0, nAtoms * sizeof(char));
    for (unsigned int rank : ranks) {
      ranksSet[rank] = 1;
    }
    for (unsigned int i = 0; i < nAtoms; i++) {
      if (ranksSet[i] != 1) {
        std::cerr << "Molecule has non unique ranks: " << MolToSmiles(*nm, true)
                  << ", Renumbering round: " << i << std::endl;
        for (unsigned int i = 0; i < nAtoms; i++) {
          std::cerr << "AtomIdx: " << i << " Rank: " << ranks[i] << std::endl;
        }
      }
      TEST_ASSERT(ranksSet[i] == 1);
    }
    delete nm;
    free(ranksSet);
  }
}
}  // namespace

void test7a() {
  BOOST_LOG(rdInfoLog) << "testing some specific ordering problems"
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  std::string smi1, smi2;
  {
    std::string fName = rdbase + "/Code/GraphMol/test_data/canon_reorder1.mol";
    RWMol *m = MolFileToMol(fName, false, false);
    TEST_ASSERT(m);
    MolOps::sanitizeMol(*m);
    std::vector<unsigned int> atomRanks;
    //    std::cerr <<"\n\n\n\n\n\n\n\n\n\n\n\n>--------------" << std::endl;
    RDKit::Canon::rankMolAtoms(*m, atomRanks, false);
    //    std::cerr <<"---------------" << std::endl;
    //    for(unsigned int i=0;i<m->getNumAtoms();++i){
    //      std::cerr<<" "<<i+1<<" "<<atomRanks[i]<<std::endl;
    //    }
    //    std::cerr <<"---------------" << std::endl;
    smi1 = MolToSmiles(*m, true);
    delete m;
  }
  {
    std::string fName = rdbase + "/Code/GraphMol/test_data/canon_reorder2.mol";
    RWMol *m = MolFileToMol(fName, false, false);
    TEST_ASSERT(m);
    MolOps::sanitizeMol(*m);
    std::vector<unsigned int> atomRanks;
    //    std::cerr <<">--------------" << std::endl;
    RDKit::Canon::rankMolAtoms(*m, atomRanks, false);
    //    std::cerr <<"---------------" << std::endl;
    //    for(unsigned int i=0;i<m->getNumAtoms();++i){
    //      std::cerr<<" "<<i+1<<" "<<atomRanks[i]<<std::endl;
    //    }
    //    std::cerr <<"---------------" << std::endl;
    smi2 = MolToSmiles(*m, true);
    delete m;
  }
  if (smi1 != smi2) {
    std::cerr << smi1 << "\n" << smi2 << std::endl;
  }
  TEST_ASSERT(smi1 == smi2);
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

std::string smis[] = {
    "C[C@@H]1CCC[C@H](C)[C@H]1C", "N[C@@]1(C[C@H]([18F])C1)C(=O)O",
    "CC12CCCC1C1CCC3CC(O)CCC3(C)C1CC2",
    "CC(C)CCCC[C@@H]1C[C@H](/C=C/[C@]2(C)CC[C@H](O)CC2)[C@@H](O)[C@H]1O",
    "C[C@@]12CCC[C@H]1[C@@H]1CC[C@H]3C[C@@H](O)CC[C@]3(C)[C@H]1CC2",
    "CCCN[C@H]1CC[C@H](NC)CC1",
    "O=S(=O)(NC[C@H]1CC[C@H](CNCc2ccc3ccccc3c2)CC1)c1ccc2ccccc2c1",
    "CC(C)[C@H]1CC[C@H](C(=O)N[C@H](Cc2ccccc2)C(=O)O)CC1",
    "O=[N+]([O-])c1ccccc1S(=O)(=O)NC[C@H]1CC[C@H](CNCC2Cc3ccccc3CC2)CC1",
    "Oc1ccc2c(Cc3ccc(OCCN4CCCCC4)cc3)c([C@H]3CC[C@H](O)CC3)sc2c1",
    "O=C(c1ccc(OCCN2CCCCC2)cc1)c1c2ccc(O)cc2sc1[C@H]1CC[C@H](O)CC1",
    "N#Cc1ccc2c(c1)CCN(CC[C@@H]1CC[C@@H](NC(=O)c3ccnc4ccccc34)CC1)C2",
    "COCCOC[C@H](CC1(C(=O)N[C@H]2CC[C@@H](C(=O)O)CC2)CCCC1)C(=O)O",
    "c1ccc(CN[C@H]2CC[C@H](Nc3ccc4[nH]ncc4c3)CC2)cc1",
    "CCC1=C(C)CN(C(=O)NCCc2ccc(S(=O)(=O)NC(=O)N[C@H]3CC[C@H](C)CC3)cc2)C1=O",
    "C[C@H]1C[C@H](C1)N1CCC1", "C[C@H]1C[C@H](C1)N1CCN(C)CC1",
    "CN1CCN(CC1)[C@H]1C[C@H](C1)c1ncc2c(N)nccn12",
    "CN1CCN(CC1)[C@H]1C[C@H](C1)c1nc(-c2ccc3ccc(nc3c2)-c2ccccc2)c2c(N)nccn12",
    "C12C3C1C3C4C5C4C52", "N[C@H]1C2CC3CC1C[C@](O)(C3)C2",
    "O=C(CN1CCN(c2ccc(C(F)(F)F)cn2)CC1)N[C@H]1C2CC3CC1C[C@](O)(C3)C2",
    "COc1cc([C@H]2[C@H](C)[C@H](C)[C@H]2c2ccc(O)c(OC)c2)ccc1O",
    "N[C@@H]1[C@H]2CNC[C@@H]12",
    "N[C@@H]1[C@H]2CN(c3nc4c(cc3F)c(=O)c(C(=O)O)cn4C3CC3)C[C@@H]12",
    // some examples that came up while doing a torture test in ZINC
    "CN1CCCNCCN(C)CCC[NH2+]CC1", "CN1CCC[NH2+]CCN(C)CCC[NH+](C)CC1",
    "O=P([O-])([O-])C[NH+]1CCCN(CP(=O)([O-])O)CC[NH+](CP(=O)([O-])[O-])CCC[NH+]"
    "(CP(=O)([O-])[O-])CC1",
    "O=C(CCNC(=O)Cn1cnc2ccccc2c1=O)NCCc1c[nH]c2ccccc12",
    "C1CNCC[NH2+]CC(C2CNCC[NH2+]CCC[NH2+]CCNC2)CNCC[NH2+]C1",
    // C60
    "C12=C3C4=C5C6=C1C7=C8C9=C1C%10=C%11C(=C29)C3=C2C3=C4C4=C5C5=C9C6=C7C6="
    "C7C8=C1C1=C8C%10=C%10C%11=C2C2=C3C3=C4C4=C5C5=C%11C%12=C(C6=C95)C7=C1C1=C%"
    "12C5=C%11C4=C3C3=C5C(=C81)C%10=C23",
    // C70
    "C12=C3C4=C5C6=C7C8=C9C%10=C%11C%12=C%13C%10=C%10C8=C5C1=C%10C1=C%13C5="
    "C8C1=C2C1=C3C2=C3C%10=C%13C%14=C3C1=C8C1=C3C5=C%12C5=C8C%11=C%11C9=C7C7="
    "C9C6=C4C2=C2C%10=C4C(=C29)C2=C6C(=C8C8=C9C6=C4C%13=C9C(=C%141)C3=C85)C%11="
    "C27",
    // Bernd's example1
    "C12C3C4C3C3C1C2C43",
    // Bernd's example2
    "C12C3C1C1C4C5C4C4C6C(C6C3C3C1C43)C25",
    // doubled house
    "C12C3C45C67C8C9C66C14C21C35C78C961", "C12C3C45C6C7C11C27C41C356",
    // Problematic round-tripping
    "COC(=O)/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(\\C)/C=C/C(=O)[O-]",
    "COC(=O)/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(\\C)/C=C/C=C(\\C)/C=C/C(=O)[O-]",
    "c1cc2ccc(ccc3ccc1cc3)cc2", "C13C6C1C2C4C2C3C5C4C56",
    "C45C1C6C3C6C5C4C2C3C12", "C45C2C6C3C6C5C4C1C3C12",
    "Cl[C@H]1[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H]1Cl",
    "N[C@]1(C(=O)O)C[C@H](n2oc(=O)[nH]c2=O)C1", "CC1CC(C)CC(C)C1",
    "C[C@H]1C[C@@H](C)C[C@H](C)C1", "C[C@H]1C[C@@H](C)C[C@@H](C)C1",
    // Stereochemistry in large rings
    "C1[C@@H](C(C)=O)[C@H]2[C@H](C(C)=C1)[C@@H]1O[C@H]2[C@@H](O)C/C=C\\CC1",
    // Chembl 20 examples
    "COC(=O)CC[C@]12C[C@]13CC[C@]1(C)[C@@H]4[C@@H](C[C@@]1(C)[C@@H]3CC[C@H]2C("
    "C)(C)O)O[C@@]1(C[C@@H](C)C(=O)O1)C[C@H]4C",
    "C/C=C/"
    "[C@H](O)C(C)(C)[C@@H]1CC=CC=CC=C[C@H](OC)Cc2nc(co2)C(=O)O[C@H](C(C)(C)[C@@"
    "H](O)/C=C/C)CC=C[C@H]2O[C@H]2C=CC=Cc2nc(co2)C(=O)O1",
    "N[C@]1(C(=O)O)C[C@@H](n2oc(=O)[nH]c2=O)C1",
    "COc1ccc([C@H]2C3(CO)C4[N@](C)C5C2(CO)C2[N@](C)C3C4(CO)[C@H](c3ccc(OC)cc3)"
    "C52CO)cc1",
    "NCCNC(=O)[C@]1(O)C[C@@H](O)[C@H](O)[C@H](O)C1",
    // CIP renumber problems
    //  "O=P(O)(O)O[C@@H]1[C@H](OP(=O)(O)O)[C@@H](OP(=O)(O)O)[C@H](OP(=O)(O)O)[C@@H](OP(=O)(O)O[32P](=O)(O)O)[C@@H]1OP(=O)(O)O",
    //  "O=C(NCCO[C@H]1[C@@H](O)[C@H](OP(=O)(O)O)[C@@H](OP(=O)(O)O)[C@H](O)[C@H]1OP(=O)(O)O)NCCO[C@H]1[C@@H](O)[C@H](OP(=O)(O)O)[C@@H](O[PH](O)(O)O)[C@H](O)[C@H]1OP(=O)(O)O",
    //  "O=C(NCCO[C@@H]1[C@H](OP(=O)(O)O)[C@@H](OP(=O)(O)O)[C@H](O[PH](O)(O)O)[C@@H](OP(=O)(O)O)[C@@H]1OP(=O)(O)O)NCCO[C@H]1[C@@H](OP(=O)(O)O)[C@H](O)[C@@H](OP(=O)(O)O)[C@H](OP(=O)(O)O)[C@H]1O",
    //  "O=C(O)[C@H]1[C@@H](C(=O)O)[C@@H](C(=O)O)[C@@H]1C(=O)O",
    "C[C@H]1[C@H](C)[C@@H](C)[C@H]1C",
    "COc1cc([C@H]2[C@@](NC(=O)c3ccc(NC(=O)C4CCCC4)cc3)(C(=O)O)[C@@H](c3ccc(OC(="
    "O)c4cccs4)c(OC)c3)[C@]2(NC(=O)c2ccc(NC(=O)C3CCCC3)cc2)C(=O)O)ccc1OC(=O)"
    "c1cccs1",
    "COc1cc([C@H]2[C@@](NC(=O)c3ccc(NC(=O)OC(C)(C)C)cc3)(C(=O)O)[C@@H](c3ccc("
    "OC(=O)c4cccs4)c(OC)c3)[C@]2(NC(=O)c2ccc(NC(=O)OC(C)(C)C)cc2)C(=O)O)ccc1OC("
    "=O)c1cccs1",
    "COc1cc([C@H]2[C@](NC(=O)c3ccc(NC(=O)OC(C)(C)C)cc3)(C(=O)O)[C@H](c3ccc(OC(="
    "O)c4cccs4)c(OC)c3)[C@]2(NC(=O)c2ccc(NC(=O)OC(C)(C)C)cc2)C(=O)O)ccc1OC(=O)"
    "c1cccs1",
    "CCC[C@H]1CC[C@H]([C@H]2CC[C@H](OC(=O)[C@H]3[C@@H](c4ccc(O)cc4)[C@H](C(=O)"
    "O[C@H]4CC[C@H]([C@H]5CC[C@H](CCC)CC5)CC4)[C@@H]3c3ccc(O)cc3)CC2)CC1",
    // test molecules with atom-mapping numbers
    "[O:1]=[C:2]([CH2:3][C:4]1=[CH:5][CH:6]=[CH:7][CH:8]=[CH:9]1)[NH2:10]",
    // chembl molecules with multiple fragments
    // CHEMBL439119
    "CCCCCCC1C23C4=c5c6c7c8c9c%10c%11c%12c%13c%14c%15c%16c%17c%18c%19c%20c%21c%"
    "22c%23c(c5c5c6c6c8c%11c8c%11c%12c%15c%12c(c%20%16)c%21c%15c%23c5c(c68)c%"
    "15c%12%11)C2(C[N+]1(C)C)C%22C%19c1c-%18c2c5c(c13)C4C7C9=C5C1(C2C%17%14)C("
    "CCCCCC)[N+](C)(C)CC%10%131.[I-].[I-]",
    // CHEMBL1203199
    "C[C@H](NC(=O)[C@H]1Cc2c(sc3ccccc23)CN1)c1ccccc1.Cl",
    // CHEMBL501667
    "CCn1c2ccc3cc2c2cc(ccc21)C(=O)c1ccc(cc1)Cn1cc[n+](c1)Cc1ccc(cc1)-c1cccc(-"
    "c2ccc(cc2)C[n+]2ccn(c2)Cc2ccc(cc2)C3=O)c1C(=O)O.[Br-].[Br-]",
    // CHEMBL12438
    "CCCCCCCCCCCCCCCCCCNC(=O)OC[C@H]1C[C@H]([C@@H2]OC(=O)N(Cc2cccc[n+]2CC)C(C)="
    "O)C1.[I-]",
    // CHEMBL1172371
    "CC.CCCCCCCCCC(C(=O)NCCc1ccc(OP(=S)(Oc2ccc(CCNC(=O)C(CCCCCCCCC)P(=O)([O-])"
    "O)cc2)N(C)/N=C/c2ccc(OP3(Oc4ccc(/C=N/"
    "N(C)P(=S)(Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)([O-])O)cc5)Oc5ccc(CCNC(=O)C("
    "CCCCCCCCC)P(=O)([O-])O)cc5)cc4)=NP(Oc4ccc(/C=N/"
    "N(C)P(=S)(Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)([O-])O)cc5)Oc5ccc(CCNC(=O)C("
    "CCCCCCCCC)P(=O)([O-])O)cc5)cc4)(Oc4ccc(/C=N/"
    "N(C)P(=S)(Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)([O-])O)cc5)Oc5ccc(CCNC(=O)C("
    "CCCCCCCCC)P(=O)([O-])O)cc5)cc4)=NP(Oc4ccc(/C=N/"
    "N(C)P(=S)(Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)([O-])O)cc5)Oc5ccc(CCNC(=O)C("
    "CCCCCCCCC)P(=O)([O-])O)cc5)cc4)(Oc4ccc(/C=N/"
    "N(C)P(=S)(Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)([O-])O)cc5)Oc5ccc(CCNC(=O)C("
    "CCCCCCCCC)P(=O)([O-])O)cc5)cc4)=N3)cc2)cc1)P(=O)([O-])O.CCCCCCCCCCCCCCCC["
    "NH2+]OC(CO)C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)CO.CCCCCCCCCCCCCCCC[NH2+]OC(CO)"
    "C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)CO.CCCCCCCCCCCCCCCC[NH2+]OC(CO)C(O)C(OC1OC("
    "CO)C(O)C(O)C1O)C(O)CO.CCCCCCCCCCCCCCCC[NH2+]OC(CO)C(O)C(OC1OC(CO)C(O)C(O)"
    "C1O)C(O)CO.CCCCCCCCCCCCCCCC[NH2+]OC(CO)C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)CO."
    "CCCCCCCCCCCCCCCC[NH2+]OC(CO)C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)CO."
    "CCCCCCCCCCCCCCCC[NH2+]OC(CO)C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)CO."
    "CCCCCCCCCCCCCCCC[NH2+]OC(CO)C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)CO."
    "CCCCCCCCCCCCCCCC[NH2+]OC(CO)C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)CO."
    "CCCCCCCCCCCCCCCC[NH2+]OC(CO)C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)CO."
    "CCCCCCCCCCCCCCCC[NH2+]OC(CO)C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)CO."
    "CCCCCCCCCCCCCCCC[NH2+]OC(CO)C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)CO",
    // Examples first reviewer
    "C12C3C4C1C3C1C3C4C1C23",  // does not initially work
    "C12C3C1C1C4C5C(C23)C5C14", "C12C3C4C5C1C1C4C5C3C21",
    "C12C3C4C1C2C1C2C4C2C31", "C12C3C4C5C2C2C6C1C(C5C36)C42",
    "C12C3C4C5C1C1C3C3C5C1C2C43",
    "C12C3C4C5C6C1C1C7C3C3C5C1C1C6C3C2C7C41",  // does not initially work
    "C12C3C4C5C6C1C1C7C3C3C6C6C4C7C2C3C5C16",
    "C12C3C4C5C6C7C8C9C1C6C1C(C37)C9C5C2C8C41",
    "C12C3C4C5C6C7C1C1C8C4C4C9C2C2C5C5C1C1C3C(C7C45)C2C8C6C91",  // does not
                                                                 // initially
                                                                 // work
    "C12C3C4C5C6C7C1C1C8C4C4C7C7C3C3C8C8C2C2C5C3C4C(C1C72)C68",
    "C12C3C4C5C6C7C8C1C1C9C5C5C%10C2C2C%11C%12C%13C3C3C7C%10C7C4C%11C1C3C(C5C8%"
    "12)C(C62)C7C9%13",  // does not initially work
    // drawn examples first reviewer
    "C12C3C4C1CC5C46C7C5C1C57C6C53C1C2", "C1C2C3C4CC5C6C1C17C8C61C5C48C3C27",
    // part of github #3490
    "C[C@H](C1)C[C@]12CCN2", "EOS"};

void test7() {
  BOOST_LOG(rdInfoLog) << "testing stability w.r.t. renumbering." << std::endl;
  unsigned int i = 0;
  while (smis[i] != "EOS") {
    std::string smiles = smis[i++];
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    MolOps::assignStereochemistry(*m, true);
    //_renumberTest(m, smiles, 1000);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

std::string molbl1 =
    "CHEMBL1950780                                                        \n"
    "     RDKit          2D                                               \n"
    "                                                                     \n"
    " 12 12  0  0  0  0  0  0  0  0999 V2000                              \n"
    "   16.2083   -6.1750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   16.2083   -7.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   16.9204   -7.4083    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   17.6324   -7.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   17.6324   -6.1750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   16.9204   -5.7583    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   16.9204   -4.9333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   18.3480   -5.7646    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   18.3462   -7.4135    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   16.9204   -8.2333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   15.4945   -7.4135    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   15.4927   -5.7646    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  5  6  1  0                                                         \n"
    "  6  7  1  6                                                         \n"
    "  1  2  1  0                                                         \n"
    "  5  8  1  6                                                         \n"
    "  1  6  1  0                                                         \n"
    "  4  9  1  6                                                         \n"
    "  2  3  1  0                                                         \n"
    "  3 10  1  1                                                         \n"
    "  3  4  1  0                                                         \n"
    "  2 11  1  6                                                         \n"
    "  4  5  1  0                                                         \n"
    "  1 12  1  6                                                         \n"
    "M  END";

std::string molbl2 =
    "CHEMBL1874247                                                        \n"
    "     RDKit          2D                                               \n"
    "                                                                     \n"
    " 12 12  0  0  0  0  0  0  0  0999 V2000                              \n"
    "    0.0000    1.6500    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.4289    0.8250    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -1.4289    0.8250    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.4289   -0.8250    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -1.4289   -0.8250    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.0000   -1.6500    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.0000    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.7145    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.7145    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.7145   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.7145   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.0000   -0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  7  1  1  1                                                         \n"
    "  8  2  1  1                                                         \n"
    "  9  3  1  1                                                         \n"
    " 10  4  1  6                                                         \n"
    " 11  5  1  6                                                         \n"
    " 12  6  1  1                                                         \n"
    "  7  8  1  0                                                         \n"
    "  7  9  1  0                                                         \n"
    "  8 10  1  0                                                         \n"
    "  9 11  1  0                                                         \n"
    " 10 12  1  0                                                         \n"
    " 11 12  1  0                                                         \n"
    "M  END";

void test8() {
  BOOST_LOG(rdInfoLog) << "testing smiles round-tripping." << std::endl;
  std::string rdbase = getenv("RDBASE");
  {
    std::string fName = rdbase + "/Code/GraphMol/test_data/iChi1b.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    std::string smi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2 = MolToSmiles(*m, true);
    if (smi1 != smi2) {
      std::cerr << smi1 << "\n" << smi2 << std::endl;
    }
    TEST_ASSERT(smi1 == smi2);
    delete m;
  }

  {
    unsigned int i = 0;
    while (smis[i] != "EOS") {
      std::string smiles = smis[i++];
      ROMol *m = SmilesToMol(smiles);
      TEST_ASSERT(m);
      //      std::cerr<<"MolToSMILES 1"<<std::endl;
      std::string smi1 = MolToSmiles(*m, true);
      delete m;

      m = SmilesToMol(smi1);
      TEST_ASSERT(m);
      //      std::cerr<<"MolToSMILES 2"<<std::endl;
      std::string smi2 = MolToSmiles(*m, true);
      delete m;
      if (smi1 != smi2) {
        std::cerr << "Input smiles: " << smiles << "\n1. Iter: " << smi1
                  << "\n2. Iter: " << smi2 << std::endl;
      }
      TEST_ASSERT(smi1 == smi2);
    }
    {
      ROMol *m = MolBlockToMol(molbl1);
      TEST_ASSERT(m);
      std::string smiles = MolToSmiles(*m, true);
      delete m;

      m = SmilesToMol(smiles);
      TEST_ASSERT(m);
      //      std::cerr<<"MolToSMILES 1"<<std::endl;
      std::string smi1 = MolToSmiles(*m, true);
      delete m;

      if (smiles != smi1) {
        std::cerr << smiles << "\n" << smi1 << std::endl;
      }
      TEST_ASSERT(smiles == smi1);
    }
    {
      ROMol *m = MolBlockToMol(molbl2);
      TEST_ASSERT(m);
      std::string smiles = MolToSmiles(*m, true);
      delete m;

      m = SmilesToMol(smiles);
      TEST_ASSERT(m);
      //      std::cerr<<"MolToSMILES 1"<<std::endl;
      std::string smi1 = MolToSmiles(*m, true);
      delete m;

      if (smiles != smi1) {
        std::cerr << smiles << "\n" << smi1 << std::endl;
      }
      TEST_ASSERT(smiles == smi1);
    }
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void test9() {
  BOOST_LOG(rdInfoLog) << "testing chiral invariants." << std::endl;
  std::string rdbase = getenv("RDBASE");

  {
    std::string smi = "C[C@](F)(Cl)I";
    RWMol *m = SmilesToMol(smi, 0, 0);
    TEST_ASSERT(m);
    MolOps::sanitizeMol(*m);
    std::vector<unsigned int> atomRanks;
    // std::cerr<<smi<<std::endl;
    RDKit::Canon::chiralRankMolAtoms(*m, atomRanks);
    // std::copy(atomRanks.begin(),atomRanks.end(),std::ostream_iterator<unsigned
    // int>(std::cerr," "));
    // std::cerr<<std::endl;
    TEST_ASSERT(atomRanks[0] < atomRanks[2]);
    TEST_ASSERT(atomRanks[0] < atomRanks[3]);
    TEST_ASSERT(atomRanks[0] < atomRanks[4]);
    TEST_ASSERT(atomRanks[2] < atomRanks[3]);
    TEST_ASSERT(atomRanks[2] < atomRanks[4]);
    TEST_ASSERT(atomRanks[3] < atomRanks[4]);
    delete m;
  }

  {
    std::string smi = "CC[C@](F)(Cl)C=C";
    RWMol *m = SmilesToMol(smi, 0, 0);
    TEST_ASSERT(m);
    MolOps::sanitizeMol(*m);
    std::vector<unsigned int> atomRanks;
    // std::cerr<<smi<<std::endl;
    RDKit::Canon::chiralRankMolAtoms(*m, atomRanks);
    // std::copy(atomRanks.begin(),atomRanks.end(),std::ostream_iterator<unsigned
    // int>(std::cerr," "));
    // std::cerr<<std::endl;
    TEST_ASSERT(atomRanks[1] < atomRanks[3]);
    TEST_ASSERT(atomRanks[1] < atomRanks[4]);
    TEST_ASSERT(atomRanks[1] < atomRanks[5]);
    TEST_ASSERT(atomRanks[3] < atomRanks[4]);
    TEST_ASSERT(atomRanks[4] > atomRanks[5]);
    TEST_ASSERT(atomRanks[4] > atomRanks[5]);
    delete m;
  }

  {
    // make sure we aren't breaking ties
    std::string smi = "C[C@](C)(Cl)I";
    RWMol *m = SmilesToMol(smi, 0, 0);
    TEST_ASSERT(m);
    MolOps::sanitizeMol(*m);
    std::vector<unsigned int> atomRanks;
    // std::cerr<<smi<<std::endl;
    RDKit::Canon::chiralRankMolAtoms(*m, atomRanks);
    // std::copy(atomRanks.begin(),atomRanks.end(),std::ostream_iterator<unsigned
    // int>(std::cerr," "));
    // std::cerr<<std::endl;
    TEST_ASSERT(atomRanks[0] == atomRanks[2]);
    TEST_ASSERT(atomRanks[0] < atomRanks[3]);
    TEST_ASSERT(atomRanks[0] < atomRanks[4]);
    TEST_ASSERT(atomRanks[2] < atomRanks[3]);
    TEST_ASSERT(atomRanks[2] < atomRanks[4]);
    TEST_ASSERT(atomRanks[3] < atomRanks[4]);
    delete m;
  }

  {
    std::string smi = "N[C@H]1C2CC3CC1C[C@](O)(C3)C2";
    RWMol *m = SmilesToMol(smi, 0, 0);
    TEST_ASSERT(m);
    MolOps::sanitizeMol(*m);
    std::vector<unsigned int> atomRanks;
    // std::cerr<<smi<<std::endl;
    RDKit::Canon::chiralRankMolAtoms(*m, atomRanks);
    // std::copy(atomRanks.begin(),atomRanks.end(),std::ostream_iterator<unsigned
    // int>(std::cerr," "));
    // std::cerr<<std::endl;
    TEST_ASSERT(atomRanks[0] > atomRanks[1]);
    TEST_ASSERT(atomRanks[0] < atomRanks[9]);
    TEST_ASSERT(atomRanks[2] == atomRanks[6]);
    TEST_ASSERT(atomRanks[7] == atomRanks[11]);
    TEST_ASSERT(atomRanks[3] == atomRanks[5]);
    TEST_ASSERT(atomRanks[2] > atomRanks[3]);
    TEST_ASSERT(atomRanks[2] > atomRanks[11]);
    TEST_ASSERT(atomRanks[3] < atomRanks[11]);
    delete m;
  }

  {
    // this one was a chiral ranking problem
    std::string smi = "COC(C)CC(C)(C)O";
    RWMol *m = SmilesToMol(smi, 0, 0);
    TEST_ASSERT(m);
    MolOps::sanitizeMol(*m);
    std::vector<unsigned int> atomRanks;
    // std::cerr<<smi<<std::endl;
    RDKit::Canon::chiralRankMolAtoms(*m, atomRanks);
    // std::copy(atomRanks.begin(),atomRanks.end(),std::ostream_iterator<unsigned
    // int>(std::cerr," "));
    // std::cerr<<std::endl;
    TEST_ASSERT(atomRanks[1] > atomRanks[8]);
    TEST_ASSERT(atomRanks[5] > atomRanks[2]);
    delete m;
  }

  {
    // are double bonds being handled correctly?
    std::string smi = "OC[C@H](F)C=O";
    RWMol *m = SmilesToMol(smi, 0, 0);
    TEST_ASSERT(m);
    MolOps::sanitizeMol(*m);
    std::vector<unsigned int> atomRanks;
    // std::cerr<<smi<<std::endl;
    RDKit::Canon::chiralRankMolAtoms(*m, atomRanks);
    // std::copy(atomRanks.begin(),atomRanks.end(),std::ostream_iterator<unsigned
    // int>(std::cerr," "));
    // std::cerr<<std::endl;
    TEST_ASSERT(atomRanks[0] < atomRanks[5]);
    TEST_ASSERT(atomRanks[1] < atomRanks[4]);
    delete m;
  }

  {
    // are double bonds being handled correctly?
    std::string smi = "O=C[C@H](F)CO";
    RWMol *m = SmilesToMol(smi, 0, 0);
    TEST_ASSERT(m);
    MolOps::sanitizeMol(*m);
    std::vector<unsigned int> atomRanks;
    // std::cerr<<smi<<std::endl;
    RDKit::Canon::chiralRankMolAtoms(*m, atomRanks);
    // std::copy(atomRanks.begin(),atomRanks.end(),std::ostream_iterator<unsigned
    // int>(std::cerr," "));
    // std::cerr<<std::endl;
    TEST_ASSERT(atomRanks[0] > atomRanks[5]);
    TEST_ASSERT(atomRanks[1] > atomRanks[4]);
    delete m;
  }

  {
    // are double bonds being handled correctly?
    std::string smi = "CC[C@](C)(CF)C=O";
    RWMol *m = SmilesToMol(smi, 0, 0);
    TEST_ASSERT(m);
    MolOps::sanitizeMol(*m);
    std::vector<unsigned int> atomRanks;
    // std::cerr<<smi<<std::endl;
    RDKit::Canon::chiralRankMolAtoms(*m, atomRanks);
    // std::copy(atomRanks.begin(),atomRanks.end(),std::ostream_iterator<unsigned
    // int>(std::cerr," "));
    // std::cerr<<std::endl;
    TEST_ASSERT(atomRanks[4] > atomRanks[6]);
    TEST_ASSERT(atomRanks[1] < atomRanks[4]);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void test10() {
  BOOST_LOG(rdInfoLog) << "testing unique ranks in w.r.t. renumbering."
                       << std::endl;
  unsigned int i = 0;
  while (smis[i] != "EOS") {
    std::string smiles = smis[i++];
    //    std::cerr<< ">>>Molecule: " << smiles << std::endl;
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    MolOps::assignStereochemistry(*m, true);
    _renumberTest2(m, smiles, 1);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void test11() {
  BOOST_LOG(rdInfoLog) << "testing mol fragments." << std::endl;
  {
    std::string smi =
        "C[C@H]([C@H](c1ccccc1)O)N2CCCCC2.C[C@@H]([C@H](c1ccccc1)O)N2CCCCC2";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    std::vector<std::string> vfragsmi;
    std::vector<std::vector<int>> frags;
    unsigned int numFrag = MolOps::getMolFrags(*m, frags);
    for (unsigned i = 0; i < numFrag; ++i) {
      std::string smii =
          MolFragmentToSmiles(*m, frags[i], nullptr, nullptr, nullptr, true);
      // std::cout << "Test "<< smii << std::endl;
      vfragsmi.push_back(smii);
    }
    std::string smi1 = MolToSmiles(*m, true);
    delete m;

    smi = "C[C@@H]([C@H](c1ccccc1)O)N2CCCCC2.C[C@H]([C@H](c1ccccc1)O)N2CCCCC2";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    std::string smi2 = MolToSmiles(*m, true);
    delete m;
    // std::cout << smi1 << "\n" << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void test12() {
  BOOST_LOG(rdInfoLog) << "testing protein round-tripping." << std::endl;
  std::string rdbase = getenv("RDBASE");
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/2FVD.pdb";
    ROMol *m = PDBFileToMol(fName);
    TEST_ASSERT(m);
    std::string smi1 = MolToSmiles(*m, true);
    delete m;

    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2 = MolToSmiles(*m, true);
    delete m;
    // std::cout << smi1 << "\n" << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testGithub1567() {
  BOOST_LOG(rdInfoLog)
      << "testing github #1567: Non-canonical result from MolFragmentToSmiles()"
      << std::endl;
  {
    ROMol *m1 = SmilesToMol("CC1CN(Cc2cccc(C)c2)C1");
    TEST_ASSERT(m1);
    int m1Ats_a[6] = {1, 12, 3, 4, 5, 11};
    std::vector<int> m1Ats(m1Ats_a, m1Ats_a + 6);
    int m1Bnds_a[5] = {12, 11, 3, 4, 13};
    std::vector<int> m1Bnds(m1Bnds_a, m1Bnds_a + 5);
    std::string smi1 = MolFragmentToSmiles(*m1, m1Ats, &m1Bnds);

    ROMol *m2 = SmilesToMol("CN(CCC)Cc1cccc(C)c1");
    TEST_ASSERT(m2);
    int m2Ats_a[6] = {3, 2, 1, 5, 6, 12};
    std::vector<int> m2Ats(m2Ats_a, m2Ats_a + 6);
    int m2Bnds_a[5] = {2, 1, 4, 5, 12};
    std::vector<int> m2Bnds(m2Bnds_a, m2Bnds_a + 5);
    std::string smi2 = MolFragmentToSmiles(*m2, m2Ats, &m2Bnds);

    TEST_ASSERT(smi1 == smi2);

    delete m1;
    delete m2;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testCanonicalDiastereomers() {
  // FIX: this is another one that we dno't currently handle properly
#if 0
  BOOST_LOG(rdInfoLog) << "testing diastereomer problem." << std::endl;

  auto m1 = "F[C@@H](Cl)[C@H](F)Cl"_smiles;
  auto m2 = "F[C@H](Cl)[C@@H](F)Cl"_smiles;
  auto smi1 = MolToSmiles(*m1);
  auto smi2 = MolToSmiles(*m2);
  TEST_ASSERT(smi1 != smi2);
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
#endif
}

void testRingsAndDoubleBonds() {
// FIX: we don't currently handle this case properly
#if 0
  BOOST_LOG(rdInfoLog)
      << "testing some particular ugly para-stereochemistry examples."
      << std::endl;
  std::vector<std::string> smis = {"C/C=C/C=C/C=C/C=C/C", "C/C=C1/C[C@H](O)C1",
                                   "C/C=C1/CC[C@H](O)CC1"};
  for (const auto smi : smis) {
    SmilesParserParams ps;
    ps.sanitize = false;
    ps.removeHs = false;
    std::unique_ptr<ROMol> mol(SmilesToMol(smi, ps));
    TEST_ASSERT(mol);
    mol->setProp(common_properties::_StereochemDone, 1);
    mol->updatePropertyCache();
    MolOps::setBondStereoFromDirections(*mol);
    std::cerr << "   " << MolToSmiles(*mol) << std::endl;
    _renumberTest(mol.get(), smi, 500);
    std::cerr << "   " << MolToSmiles(*mol) << std::endl;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
#endif
}

int main() {
  RDLog::InitLogs();
#if 1
  test1();
  test2();
  test3();
  test4();
  test5();
  test6();
  test7a();
  test9();
  test10();
  test11();
  test12();
  test7();
  test8();
  testGithub1567();
#endif
  testRingsAndDoubleBonds();
  testCanonicalDiastereomers();
  return 0;
}
