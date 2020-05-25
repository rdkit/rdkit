/*
 * This is an extensive modification by Greg Landrum of
 * pieces from several files in the vflib-2.0 distribution
 *
 * The initial version of the modifications was completed
 *   in April 2009.
 *
 * the original author of the vflib files is:
 *    Author: P. Foggia
 *  http://amalfi.dis.unina.it/graph/db/vflib-2.0/doc/vflib.html
 *
 */
#include <boost/graph/adjacency_list.hpp>
#include <vector>
#include <algorithm>
#include <cstring>

#ifndef __BGL_VF2_SUB_STATE_H__
#define __BGL_VF2_SUB_STATE_H__
//#define RDK_VF2_PRUNING
#define RDK_ADJ_ITER typename Graph::adjacency_iterator

namespace boost {
namespace detail {
typedef std::uint32_t node_id;
const node_id NULL_NODE = 0xFFFF;
struct NodeInfo {
  node_id id;
  node_id in;
  node_id out;
};

template <class Graph>
struct Pair {
  node_id n1, n2;
  bool hasiter{false};
  RDK_ADJ_ITER nbrbeg, nbrend;

  Pair() : n1(NULL_NODE), n2(NULL_NODE) {}
};

/**
 * The ordering by in/out degree
 */
static bool nodeInfoComp1(const NodeInfo &a, const NodeInfo &b) {
  if (a.out < b.out) return true;
  if (a.out > b.out) return false;
  if (a.in < b.in) return true;
  if (a.in > b.in) return false;
  return false;
}

/**
 * The ordering by frequency/valence.
 * The frequency is in the out field, the valence in `in'.
 */
static int nodeInfoComp2(const NodeInfo &a, const NodeInfo &b) {
  if (!a.in && b.in) return 1;
  if (a.in && !b.in) return -1;
  if (a.out < b.out) return -1;
  if (a.out > b.out) return 1;
  if (a.in < b.in) return -1;
  if (a.in > b.in) return 1;
  return 0;
}

template <class Graph, class VertexDescr, class EdgeDescr>
VertexDescr getOtherIdx(const Graph &g, const EdgeDescr &edge,
                        const VertexDescr &vertex) {
  VertexDescr tmp = boost::source(edge, g);
  if (tmp == vertex) {
    tmp = boost::target(edge, g);
  }
  return tmp;
}

/*----------------------------------------------------
 * Sorts the nodes of a graphs, returning a
 * heap-allocated vector (using new) with the node ids
 * in the proper orders.
 * The sorting criterion takes into account:
 *    1 - The number of nodes with the same in/out
 *        degree.
 *    2 - The valence of the nodes.
 * The nodes at the beginning of the vector are
 * the most singular, from which the matching should
 * start.
 *--------------------------------------------------*/
template <class Graph>
node_id *SortNodesByFrequency(const Graph *g) {
  std::vector<NodeInfo> vect;
  vect.reserve(boost::num_vertices(*g));
  typename Graph::vertex_iterator bNode, eNode;
  boost::tie(bNode, eNode) = boost::vertices(*g);
  while (bNode != eNode) {
    NodeInfo t;
    t.id = vect.size();
    t.in = boost::out_degree(*bNode, *g);  // <- assuming undirected graph
    t.out = boost::out_degree(*bNode, *g);
    vect.push_back(t);
    ++bNode;
  }
  std::sort(vect.begin(), vect.end(), nodeInfoComp1);

  unsigned int run = 1;
  for (unsigned int i = 0; i < vect.size(); i += run) {
    for (run = 1; i + run < vect.size() && vect[i + run].in == vect[i].in &&
                  vect[i + run].out == vect[i].out;
         ++run)
      ;
    for (unsigned int j = 0; j < run; ++j) {
      vect[i + j].in += vect[i + j].out;
      vect[i + j].out = run;
    }
  }
  std::sort(vect.begin(), vect.end(), nodeInfoComp2);

  node_id *nodes = new node_id[vect.size()];
  for (unsigned int i = 0; i < vect.size(); ++i) {
    nodes[i] = vect[i].id;
  }

  return nodes;
}

/*----------------------------------------------------------
 * class VF2SubState
 * A representation of the SSS current state
 ---------------------------------------------------------*/
template <class Graph, class VertexCompatible, class EdgeCompatible,
          class MatchChecking>
class VF2SubState {
 private:
  Graph *g1, *g2;
  VertexCompatible &vc;
  EdgeCompatible &ec;
  MatchChecking &mc;
  unsigned int n1, n2;

  unsigned int core_len;
  unsigned int t1_len;
  unsigned int t2_len;  // Core nodes are also counted by these...
  node_id *core_1;
  node_id *core_2;
  node_id *term_1;
  node_id *term_2;

  node_id *order;

  long *share_count;
  int *vs_compared;

 public:
  VF2SubState(Graph *ag1, Graph *ag2, VertexCompatible &avc,
              EdgeCompatible &aec, MatchChecking &amc, bool sortNodes = false)
      : g1(ag1),
        g2(ag2),
        vc(avc),
        ec(aec),
        mc(amc),
        n1(num_vertices(*ag1)),
        n2(num_vertices(*ag2)) {
    if (sortNodes) {
      order = SortNodesByFrequency(ag1);
    } else {
      order = nullptr;
    }

    core_len = 0;
    t1_len = 0;
    t2_len = 0;

    core_1 = new node_id[n1];
    core_2 = new node_id[n2];
    term_1 = new node_id[n1];
    term_2 = new node_id[n2];
    share_count = new long;

    for (unsigned int i = 0; i < n1; i++) {
      core_1[i] = NULL_NODE;
      term_1[i] = 0;
    }
    for (unsigned int i = 0; i < n2; i++) {
      core_2[i] = NULL_NODE;
      term_2[i] = 0;
    }
    vs_compared = nullptr;
    // vs_compared = new int[n1*n2];
    // memset((void *)vs_compared,0,n1*n2*sizeof(int));

    // es_compared = new std::map<unsigned int,bool>();
    *share_count = 1;
  };

  VF2SubState(const VF2SubState &state)
      : g1(state.g1),
        g2(state.g2),
        vc(state.vc),
        ec(state.ec),
        mc(state.mc),
        n1(state.n1),
        n2(state.n2),
        order(state.order),
        vs_compared(state.vs_compared)
  // es_compared(state.es_compared)
  {
    core_len = state.core_len;
    t1_len = state.t1_len;
    t2_len = state.t2_len;

    core_1 = state.core_1;
    core_2 = state.core_2;
    term_1 = state.term_1;
    term_2 = state.term_2;
    share_count = state.share_count;

    ++(*share_count);
  };

  ~VF2SubState() {
    if (--*share_count == 0) {
      delete[] core_1;
      delete[] core_2;
      delete[] term_1;
      delete[] term_2;
      delete share_count;
      delete[] order;
      // delete [] vs_compared;
      // delete es_compared;
    }
  };

  bool IsGoal() { return core_len == n1; };
  bool MatchChecks(const node_id c1[], const node_id c2[]) {
    return mc(c1, c2);
  };
  bool IsDead() { return n1 > n2 || t1_len > t2_len; };
  unsigned int CoreLen() { return core_len; }
  Graph *GetGraph1() { return g1; }
  Graph *GetGraph2() { return g2; }

  bool NextPair(Pair<Graph> &pair) {
    if (pair.n1 == NULL_NODE) pair.n1 = 0;
    if (pair.n2 == NULL_NODE)
      pair.n2 = 0;
    else
      pair.n2++;

#if 0
    std::cerr<<" **** np: "<< prev_n1<<","<<prev_n2<<std::endl;
    std::cerr<<"in_1 ";
    for(unsigned int i=0;i<n1;++i){
      std::cerr<<"("<<in_1[i]<<","<<out_1[i]<<"), ";
    } 
    std::cerr<<std::endl;
    std::cerr<<"in_2 ";
    for(unsigned int i=0;i<n2;++i){
      std::cerr<<"("<<in_2[i]<<","<<out_2[i]<<"), ";
    } 
    std::cerr<<std::endl;
#endif
    if (t1_len > core_len && t2_len > core_len) {
      while (pair.n1 < n1 &&
             (core_1[pair.n1] != NULL_NODE || term_1[pair.n1] == 0)) {
        pair.n1++;
        pair.n2 = 0;
      }

      /* Initialize VF2 Plus neighbor iterator.
       * The next query node (pair.n1) has been selected from the terminal
       * set and is therefore adjacent to an already mapped atom (in
       * core_1). Rather than select pair.n2 from all atoms (0...n2) we can
       * select it from the neighbors of this mapped atom (0...deg(nbor))
       * since it must also be adajcent to this mapped atom!
       */
      if (!pair.hasiter) {
        RDK_ADJ_ITER n1iter_beg, n1iter_end;
        boost::tie(n1iter_beg, n1iter_end) =
            boost::adjacent_vertices(pair.n1, *g1);

        while (n1iter_beg != n1iter_end && core_1[*n1iter_beg] == NULL_NODE)
          ++n1iter_beg;

        assert(n1iter_beg != n1iter_end);

        boost::tie(pair.nbrbeg, pair.nbrend) =
            boost::adjacent_vertices(core_1[*n1iter_beg], *g2);
        pair.hasiter = true;
      }
    } else if (pair.n1 == 0 && order != nullptr) {
      // Optimisation: if the order vector is laid out in a DFS/BFS then this
      // loop can be replaced with:
      //   pair.n1=order[core_len];
      // :)
      unsigned int i = 0;
      while (i < n1 && core_1[pair.n1 = order[i]] != NULL_NODE) i++;
      if (i == n1) pair.n1 = n1;
    } else {
      while (pair.n1 < n1 && core_1[pair.n1] != NULL_NODE) {
        pair.n1++;
        pair.n2 = 0;
      }
    }

    /* VF2 Plus iterator available? */
    if (pair.hasiter) {
      while (pair.nbrbeg < pair.nbrend && core_2[*pair.nbrbeg] != NULL_NODE) {
        ++pair.nbrbeg;
      }

      if (pair.nbrbeg < pair.nbrend) {
        pair.n2 = *pair.nbrbeg;
        ++pair.nbrbeg;
      } else {
        pair.n2 = n2;
      }
    } else if (t1_len > core_len && t2_len > core_len) {
      while (pair.n2 < n2 &&
             (core_2[pair.n2] != NULL_NODE || term_2[pair.n2] == 0)) {
        pair.n2++;
      }
    } else {
      while (pair.n2 < n2 && core_2[pair.n2] != NULL_NODE) {
        pair.n2++;
      }
    }
    return pair.n1 < n1 && pair.n2 < n2;
  };
  bool IsFeasiblePair(node_id node1, node_id node2) {
    assert(node1 < n1);
    assert(node2 < n2);
    assert(core_1[node1] == NULL_NODE);
    assert(core_2[node2] == NULL_NODE);

    // std::cerr<<"  ifp:"<<node1<<"-"<<node2<<"
    // "<<vs_compared->size()<<std::endl;
    // int &isCompat=vs_compared[node1*n2+node2];
    // if(isCompat==0){
    //   isCompat=vc(node1,node2)?1:-1;
    // }
    // if( isCompat<0 ){
    //   //std::cerr<<"  short1"<<std::endl;
    //   return false;
    // }

    // O(1) check for adjacency list
    if (boost::out_degree(node1, *g1) > boost::out_degree(node2, *g2))
      return false;
    if (!vc(node1, node2)) return false;

    unsigned int other1, other2;
#ifdef RDK_VF2_PRUNING
    unsigned int term1 = 0, term2 = 0;
    unsigned int new1 = 0, new2 = 0;
#endif

    // Check the out edges of node1
    typename Graph::out_edge_iterator bNbrs, eNbrs;
    boost::tie(bNbrs, eNbrs) = boost::out_edges(node1, *g1);
    while (bNbrs != eNbrs) {
      other1 = getOtherIdx(*g1, *bNbrs, node1);
      if (core_1[other1] != NULL_NODE) {
        other2 = core_1[other1];
        typename Graph::edge_descriptor oEdge;
        bool found;
        boost::tie(oEdge, found) = boost::edge(node2, other2, *g2);
        if (!found || !ec(*bNbrs, oEdge)) {
          // std::cerr<<"  short2"<<std::endl;
          return false;
        }
      }
#ifdef RDK_VF2_PRUNING
      else {
        if (term_1[other1]) ++term1;
        if (!term_1[other1]) ++new1;
      }
#endif
      ++bNbrs;
    }

#ifdef RDK_VF2_PRUNING
    // Check the out edges of node2
    boost::tie(bNbrs, eNbrs) = boost::out_edges(node2, *g2);
    while (bNbrs != eNbrs) {
      other2 = getOtherIdx(*g2, *bNbrs, node2);
      if (core_2[other2] != NULL_NODE) {
        // do nothing
      } else {
        if (term_2[other2]) ++term2;
        if (!term_2[other2]) ++new2;
      }
      ++bNbrs;
    }
    // std::cerr<<(termin1 <= termin2 && termout1 <= termout2 &&
    // (termin1+termout1+new1)<=(termin2+termout2+new2))<<std::endl;

    // n.b. term1+new1 == boost::out_degree(node1) and
    //      term2+new2 == boost::out_degree(node2)
    return term1 <= term2 && (term1 + new1) <= (term2 + new2);
#else
    return true;
#endif
  };
  void AddPair(node_id node1, node_id node2) {
    assert(node1 < n1);
    assert(node2 < n2);
    assert(core_len < n1);
    assert(core_len < n2);

    ++core_len;
    if (!term_1[node1]) {
      term_1[node1] = core_len;
      ++t1_len;
    }

    if (!term_2[node2]) {
      term_2[node2] = core_len;
      ++t2_len;
    }

    core_1[node1] = node2;
    core_2[node2] = node1;

    typename Graph::out_edge_iterator bNbrs, eNbrs;
    // FIX: this is explicitly ignoring directionality
    boost::tie(bNbrs, eNbrs) = boost::out_edges(node1, *g1);
    while (bNbrs != eNbrs) {
      unsigned int other = getOtherIdx(*g1, *bNbrs, node1);
      if (!term_1[other]) {
        term_1[other] = core_len;
        ++t1_len;
      }
      ++bNbrs;
    }

    // FIX: this is explicitly ignoring directionality
    boost::tie(bNbrs, eNbrs) = boost::out_edges(node2, *g2);
    while (bNbrs != eNbrs) {
      unsigned int other = getOtherIdx(*g2, *bNbrs, node2);
      if (!term_2[other]) {
        term_2[other] = core_len;
        ++t2_len;
      }
      ++bNbrs;
    }
  };
  void GetCoreSet(node_id c1[], node_id c2[]) {
    unsigned int i, j;
    for (i = 0, j = 0; i < n1; ++i) {
      if (core_1[i] != NULL_NODE) {
        c1[j] = i;
        c2[j] = core_1[i];
        ++j;
      }
    }
  };
  VF2SubState *Clone() { return new VF2SubState(*this); };
  void BackTrack(node_id node1, node_id node2) {
    if (term_1[node1] == core_len) {
      term_1[node1] = 0;
      --t1_len;
    }

    typename Graph::out_edge_iterator bNbrs, eNbrs;
    boost::tie(bNbrs, eNbrs) = boost::out_edges(node1, *g1);
    while (bNbrs != eNbrs) {
      unsigned int other = getOtherIdx(*g1, *bNbrs, node1);
      if (term_1[other] == core_len) {
        term_1[other] = 0;
        --t1_len;
      }
      ++bNbrs;
    }

    if (term_2[node2] == core_len) {
      term_2[node2] = 0;
      --t2_len;
    }

    boost::tie(bNbrs, eNbrs) = boost::out_edges(node2, *g2);
    while (bNbrs != eNbrs) {
      unsigned int other = getOtherIdx(*g2, *bNbrs, node2);
      if (term_2[other] == core_len) {
        term_2[other] = 0;
        --t2_len;
      }
      ++bNbrs;
    }

    core_1[node1] = NULL_NODE;
    core_2[node2] = NULL_NODE;
    --core_len;
  };
  bool Match(node_id c1[], node_id c2[]) {
    if (IsGoal()) {
      GetCoreSet(c1, c2);
      if (MatchChecks(c1, c2)) return true;
    }

    if (IsDead()) return false;

    Pair<Graph> pair;
    while (NextPair(pair)) {
      if (IsFeasiblePair(pair.n1, pair.n2)) {
        AddPair(pair.n1, pair.n2);
        if (Match(c1, c2))  // recurse
          return true;
        BackTrack(pair.n1, pair.n2);
      }
    }
    return false;
  }

  template <class DoubleBackInsertionSequence>
  bool MatchAll(node_id c1[], node_id c2[], DoubleBackInsertionSequence &res,
                unsigned int lim = 0) {
    if (IsGoal()) {
      GetCoreSet(c1, c2);
      if (MatchChecks(c1, c2)) {
        typename DoubleBackInsertionSequence::value_type newSeq;
        for (unsigned int i = 0; i < core_len; ++i) {
          newSeq.push_back(std::pair<int, int>(c1[i], c2[i]));
        }
        res.push_back(newSeq);
        return lim && res.size() >= lim;
      }
    }

    if (IsDead()) return false;

    Pair<Graph> pair;
    while (NextPair(pair)) {
      if (IsFeasiblePair(pair.n1, pair.n2)) {
        AddPair(pair.n1, pair.n2);
        if (MatchAll(c1, c2, res, lim))  // recurse
          return true;
        BackTrack(pair.n1, pair.n2);
      }
    }
    return false;
  }
};

/*-------------------------------------------------------------
 * static bool match(pn, c1, c2, s)
 * Finds a matching between two graphs, if it exists, starting
 * from state s.
 * Returns true a match has been found.
 * *pn is assigned the numbero of matched nodes, and
 * c1 and c2 will contain the ids of the corresponding nodes
 * in the two graphs.
 ------------------------------------------------------------*/
template <class SubState>
bool match(int *pn, node_id c1[], node_id c2[], SubState &s) {
  if (s.Match(c1, c2)) {
    // not needed, pn = num query atoms (n1)...
    *pn = s.CoreLen();
    return true;
  }
  return false;
}

/*-------------------------------------------------------------
 * static bool match(c1, c2, vis, usr_data, pcount)
 * Visits all the matchings between two graphs,  starting
 * from state s.
 * Returns true if the caller must stop the visit.
 * Stops when max_results is reached, set max_results to 0 to
 * keep going until there are no more matches
 *
 ------------------------------------------------------------*/
template <class SubState, class DoubleBackInsertionSequence>
bool match(node_id c1[], node_id c2[], SubState &s,
           DoubleBackInsertionSequence &res, unsigned int max_results) {
  s.MatchAll(c1, c2, res, max_results);
  return !res.empty();
}
};  // end of namespace detail

template <
    class Graph, class VertexLabeling  // binary predicate
    ,
    class EdgeLabeling  // binary predicate
    ,
    class MatchChecking  // binary predicate
    ,
    class
    BackInsertionSequence  // contains
                           // std::pair<vertex_descriptor,vertex_descriptor>
    >
bool vf2(const Graph &g1, const Graph &g2, VertexLabeling &vertex_labeling,
         EdgeLabeling &edge_labeling, MatchChecking &match_checking,
         BackInsertionSequence &F) {
  detail::VF2SubState<const Graph, VertexLabeling, EdgeLabeling, MatchChecking>
      s0(&g1, &g2, vertex_labeling, edge_labeling, match_checking, false);
  detail::node_id *ni1 = new detail::node_id[num_vertices(g1)];
  detail::node_id *ni2 = new detail::node_id[num_vertices(g2)];
  int n = 0;

  F.clear();
  F.resize(0);
  if (match(&n, ni1, ni2, s0)) {
    for (unsigned int i = 0; i < num_vertices(g1); i++) {
      F.push_back(std::pair<int, int>(ni1[i], ni2[i]));
    }
  }
  delete[] ni1;
  delete[] ni2;

  return !F.empty();
};
template <class Graph, class VertexLabeling  // binary predicate
          ,
          class EdgeLabeling  // binary predicate
          ,
          class MatchChecking  // binary predicate
          ,
          class DoubleBackInsertionSequence  // contains a back insertion
                                             // sequence
          >
bool vf2_all(const Graph &g1, const Graph &g2, VertexLabeling &vertex_labeling,
             EdgeLabeling &edge_labeling, MatchChecking &match_checking,
             DoubleBackInsertionSequence &F, unsigned int max_results = 1000) {
  detail::VF2SubState<const Graph, VertexLabeling, EdgeLabeling, MatchChecking>
      s0(&g1, &g2, vertex_labeling, edge_labeling, match_checking, false);
  detail::node_id *ni1 = new detail::node_id[num_vertices(g1)];
  detail::node_id *ni2 = new detail::node_id[num_vertices(g2)];

  F.clear();
  F.resize(0);

  match(ni1, ni2, s0, F, max_results);

  delete[] ni1;
  delete[] ni2;

  return !F.empty();
};
}  // end of namespace boost
#endif

#undef RDK_VF2_PRUNING
#undef RDK_ADJ_ITER
