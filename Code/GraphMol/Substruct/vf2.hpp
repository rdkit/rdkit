// $Id$
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

namespace boost{
  namespace detail {
    typedef unsigned short node_id;
    const node_id NULL_NODE=0xFFFF;
    struct NodeInfo {
      node_id id;
      node_id in;
      node_id out;
    };

    /**
     * The ordering by in/out degree
     */
    static bool nodeInfoComp1(const NodeInfo &a, const NodeInfo &b) {
      if(a.out < b.out) return true;
      if(a.out > b.out) return false;
      if(a.in < b.in) return true;
      if(a.in > b.in) return false;
      return false;
    }

    /**
     * The ordering by frequency/valence.
     * The frequency is in the out field, the valence in `in'.
     */
    static int nodeInfoComp2(const NodeInfo &a, const NodeInfo &b) {
      if (!a.in && b.in ) return 1;
      if (a.in && !b.in) return -1;
      if (a.out < b.out) return -1;
      if (a.out > b.out) return 1;
      if( a.in < b.in ) return -1;
      if (a.in > b.in) return 1;
      return 0;
    }

    template <class Graph,class VertexDescr,class EdgeDescr> 
    VertexDescr getOtherIdx(const Graph &g,const EdgeDescr &edge,const VertexDescr &vertex) {
      VertexDescr tmp=boost::source(edge,g);
      if(tmp==vertex){
        tmp=boost::target(edge,g);
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
    node_id* SortNodesByFrequency(const Graph *g) {
      std::vector<NodeInfo> vect;
      vect.reserve(boost::num_vertices(*g));
      typename Graph::vertex_iterator bNode,eNode;
      boost::tie(bNode,eNode) = boost::vertices(*g);
      while(bNode!=eNode){
        NodeInfo t;
        t.id=vect.size();
        t.in=boost::out_degree(*bNode,*g);// <- assuming undirected graph
        t.out=boost::out_degree(*bNode,*g); 
        vect.push_back(t);
        ++bNode;
      }
      std::sort(vect.begin(),vect.end(),nodeInfoComp1);
    
      unsigned int run=1;
      for(unsigned int i=0; i<vect.size(); i+=run){
        for(run=1; i+run<vect.size() && 
              vect[i+run].in==vect[i].in && 
              vect[i+run].out==vect[i].out;
            ++run) 
          ;
        for(unsigned int j=0; j<run; ++j) {
          vect[i+j].in += vect[i+j].out;
          vect[i+j].out=run;
        }
      }
      std::sort(vect.begin(),vect.end(),nodeInfoComp2);
    
      node_id *nodes=new node_id[vect.size()];
      for(unsigned int i=0; i<vect.size(); ++i){
        nodes[i]=vect[i].id;
      }
    
      return nodes;
    }

    /*----------------------------------------------------------
     * class VF2SubState
     * A representation of the SSS current state
     ---------------------------------------------------------*/
    template <class Graph,class VertexCompatible,class EdgeCompatible,class MatchChecking >
    class VF2SubState
    { 
    private:
      Graph *g1, *g2;
      VertexCompatible &vc;
      EdgeCompatible &ec;
      MatchChecking &mc;
      unsigned int n1, n2;

      unsigned int core_len, orig_core_len;
      unsigned int added_node1;
      unsigned int t1both_len, t2both_len;
      unsigned int t1in_len, t1out_len; 
      unsigned int t2in_len, t2out_len; // Core nodes are also counted by these...
      node_id *core_1;
      node_id *core_2;
      node_id *in_1;
      node_id *in_2;
      node_id *out_1;
      node_id *out_2;

      node_id *order;

      long *share_count;
      int *vs_compared;
    
    public:
      VF2SubState(Graph *ag1, Graph *ag2,
                  VertexCompatible &avc,
                  EdgeCompatible &aec,
                  MatchChecking &amc,
                  bool sortNodes=false) : g1(ag1), g2(ag2), vc(avc), ec(aec), mc(amc),
                                          n1(num_vertices(*ag1)),n2(num_vertices(*ag2)) {
        if (sortNodes){
          order = SortNodesByFrequency(ag1);
        } else {
          order = NULL;
        }

        core_len=orig_core_len=0;
        t1both_len=t1in_len=t1out_len=0;
        t2both_len=t2in_len=t2out_len=0;

        added_node1=NULL_NODE;

        core_1=new node_id[n1];
        core_2=new node_id[n2];
        in_1=new node_id[n1];
        in_2=new node_id[n2];
        out_1=new node_id[n1];
        out_2=new node_id[n2];
        share_count = new long;

        for(unsigned int i=0; i<n1; i++){
          core_1[i]=NULL_NODE;
          in_1[i]=0;
          out_1[i]=0;
        }
        for(unsigned int i=0; i<n2; i++){
          core_2[i]=NULL_NODE;
          in_2[i]=0;
          out_2[i]=0;
        }
        vs_compared=0;
        //vs_compared = new int[n1*n2];
        //memset((void *)vs_compared,0,n1*n2*sizeof(int));
        
        //es_compared = new std::map<unsigned int,bool>();
        *share_count = 1;
      };

      VF2SubState(const VF2SubState &state) :
        g1(state.g1), g2(state.g2), vc(state.vc), ec(state.ec), mc(state.mc),
        n1(state.n1),n2(state.n2), order(state.order),vs_compared(state.vs_compared)
        //es_compared(state.es_compared)
      {

        core_len=orig_core_len=state.core_len;
        t1in_len=state.t1in_len;
        t1out_len=state.t1out_len;
        t1both_len=state.t1both_len;
        t2in_len=state.t2in_len;
        t2out_len=state.t2out_len;
        t2both_len=state.t2both_len;

        added_node1=NULL_NODE;

        core_1=state.core_1;
        core_2=state.core_2;
        in_1=state.in_1;
        in_2=state.in_2;
        out_1=state.out_1;
        out_2=state.out_2;
        share_count=state.share_count;

        ++(*share_count);
      };

      ~VF2SubState(){
        if (-- *share_count == 0) {
          delete [] core_1;
          delete [] core_2;
          delete [] in_1;
          delete [] out_1;
          delete [] in_2;
          delete [] out_2;
          delete share_count;
          delete [] order;
          //delete [] vs_compared;
          //delete es_compared;
        }
      }; 

      bool IsGoal() { return core_len==n1 ; };
      bool MatchChecks(const node_id c1[],const node_id c2[]){
        return mc(c1,c2);
      };
      bool IsDead() { return n1>n2  || 
          t1both_len>t2both_len ||
          t1out_len>t2out_len ||
          t1in_len>t2in_len;
      };
      unsigned int CoreLen() { return core_len; }
      Graph *GetGraph1() { return g1; }
      Graph *GetGraph2() { return g2; }

      bool NextPair(node_id *pn1, node_id *pn2,
                    node_id prev_n1=NULL_NODE, node_id prev_n2=NULL_NODE){
        if (prev_n1==NULL_NODE)
          prev_n1=0;
        if (prev_n2==NULL_NODE)
          prev_n2=0;
        else
          prev_n2++;

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
        if (t1both_len>core_len && t2both_len>core_len) {
          while (prev_n1<n1 &&
                 (core_1[prev_n1]!=NULL_NODE || out_1[prev_n1]==0
                  || in_1[prev_n1]==0) ) {
            prev_n1++;    
            prev_n2=0;
          }
        }
        else if (t1out_len>core_len && t2out_len>core_len) {
          while (prev_n1<n1 &&
                 (core_1[prev_n1]!=NULL_NODE || out_1[prev_n1]==0) ){
            prev_n1++;    
            prev_n2=0;
          }
        }
        else if (t1in_len>core_len && t2in_len>core_len) {
          while (prev_n1<n1 &&
                 (core_1[prev_n1]!=NULL_NODE || in_1[prev_n1]==0) ) {
            prev_n1++;    
            prev_n2=0;
          }
        }
        else if (prev_n1==0 && order!=NULL) {
          unsigned int i=0;
          while (i<n1 && core_1[prev_n1=order[i]] != NULL_NODE)
            i++;
          if (i==n1)
            prev_n1=n1;
        }
        else {
          while (prev_n1<n1 && core_1[prev_n1]!=NULL_NODE ){
            prev_n1++;    
            prev_n2=0;
          }
        }

        if (t1both_len>core_len && t2both_len>core_len) {
          while (prev_n2<n2 &&
                 (core_2[prev_n2]!=NULL_NODE || out_2[prev_n2]==0
                  || in_2[prev_n2]==0) ) {
            prev_n2++;    
          }
        }
        else if (t1out_len>core_len && t2out_len>core_len) {
          while (prev_n2<n2 &&
                 (core_2[prev_n2]!=NULL_NODE || out_2[prev_n2]==0) ) {
            prev_n2++;    
          }
        }
        else if (t1in_len>core_len && t2in_len>core_len) {
          while (prev_n2<n2 &&
                 (core_2[prev_n2]!=NULL_NODE || in_2[prev_n2]==0) ) {
            prev_n2++;    
          }
        }
        else {
          while (prev_n2<n2 && core_2[prev_n2]!=NULL_NODE ){
            prev_n2++;    
          }
        }
        //std::cerr<<" "<< prev_n1<<"<"<<n1<<" "<<prev_n2<<"<"<<n2;
        if (prev_n1<n1 && prev_n2<n2) {
          *pn1=prev_n1;
          *pn2=prev_n2;
          //std::cerr<<"  Found"<<std::endl;
          return true;
        }
        //std::cerr<<"  nope"<< std::endl;
        return false;
      };
      bool IsFeasiblePair(node_id node1, node_id node2){
        assert(node1 < n1);
        assert(node2 < n2);
        assert(core_1[node1] == NULL_NODE);
        assert(core_2[node2] == NULL_NODE);

        //std::cerr<<"  ifp:"<<node1<<"-"<<node2<<" "<<vs_compared->size()<<std::endl;
        // int &isCompat=vs_compared[node1*n2+node2];
        // if(isCompat==0){
        //   isCompat=vc(node1,node2)?1:-1;
        // }
        // if( isCompat<0 ){
        //   //std::cerr<<"  short1"<<std::endl;
        //   return false;
        // }
        if(!vc(node1,node2)) return false;

        unsigned int other1, other2;
        unsigned int termout1 = 0, termout2 = 0, termin1 = 0, termin2 = 0;
        unsigned int new1 = 0, new2 = 0;

        // Check the out edges of node1
        typename Graph::out_edge_iterator bNbrs,eNbrs;
        boost::tie(bNbrs,eNbrs) = boost::out_edges(node1,*g1);
        while(bNbrs!=eNbrs){
          other1=getOtherIdx(*g1,*bNbrs,node1);
          if (core_1[other1] != NULL_NODE) {
            other2 = core_1[other1];
            typename Graph::edge_descriptor oEdge;
            bool found;
            boost::tie(oEdge,found) = boost::edge(node2,other2,*g2);
            if(!found || !ec(*bNbrs,oEdge) ){
              //std::cerr<<"  short2"<<std::endl;
              return false;
            }
          } else {
            if (in_1[other1]) ++termin1;
            if (out_1[other1]) ++termout1;
            if (!in_1[other1] && !out_1[other1]) ++new1;
          }
          ++bNbrs;
        }

        // Check the out edges of node2
        boost::tie(bNbrs,eNbrs) = boost::out_edges(node2,*g2);
        while(bNbrs!=eNbrs){
          other2=getOtherIdx(*g2,*bNbrs,node2);
          if (core_2[other2] != NULL_NODE) {
            // do nothing
          } else {
            if (in_2[other2]) ++termin2;
            if (out_2[other2]) ++termout2;
            if (!in_2[other2] && !out_2[other2]) ++new2;
          }
          ++bNbrs;
        }
        //std::cerr<<(termin1 <= termin2 && termout1 <= termout2 && (termin1+termout1+new1)<=(termin2+termout2+new2))<<std::endl;

        return termin1 <= termin2 && termout1 <= termout2 && (termin1+termout1+new1)<=(termin2+termout2+new2);
      };
      void AddPair(node_id node1, node_id node2){
        assert(node1 < n1);
        assert(node2 < n2);
        assert(core_len < n1);
        assert(core_len < n2);

        ++core_len;
        added_node1 = node1;

        if (!in_1[node1]) {
          in_1[node1] = core_len;
          ++t1in_len;
          if (out_1[node1]) ++t1both_len;
        }
        if (!out_1[node1]) {
          out_1[node1] = core_len;
          ++t1out_len;
          if (in_1[node1]) ++t1both_len;
        }

        if (!in_2[node2]) {
          in_2[node2] = core_len;
          ++t2in_len;
          if (out_2[node2]) ++t2both_len;
        }
        if (!out_2[node2]) {
          out_2[node2] = core_len;
          ++t2out_len;
          if (in_2[node2]) ++t2both_len;
        }

        core_1[node1] = node2;
        core_2[node2] = node1;

        typename Graph::out_edge_iterator bNbrs,eNbrs;
        // FIX: this is explicitly ignoring directionality
        boost::tie(bNbrs,eNbrs) = boost::out_edges(node1,*g1);
        while(bNbrs!=eNbrs){
          unsigned int other = getOtherIdx(*g1,*bNbrs,node1);
          if (!in_1[other]) {
            in_1[other] = core_len;
            ++t1in_len;
            if (out_1[other])  ++t1both_len;
          }
          if (!out_1[other]) {
            out_1[other] = core_len;
            ++t1out_len;
            if (in_1[other])  ++t1both_len;
          }
          ++bNbrs;
        }

        // FIX: this is explicitly ignoring directionality
        boost::tie(bNbrs,eNbrs) = boost::out_edges(node2,*g2);
        while(bNbrs!=eNbrs){
          unsigned int other = getOtherIdx(*g2,*bNbrs,node2);
          if (!in_2[other]) {
            in_2[other] = core_len;
            ++t2in_len;
            if (out_2[other]) ++t2both_len;
          }
          if (!out_2[other]) {
            out_2[other] = core_len;
            ++t2out_len;
            if (in_2[other]) ++t2both_len;
          }
          ++bNbrs;
        }
      };
      void GetCoreSet(node_id c1[], node_id c2[]){
        unsigned int i, j;
        for (i = 0, j = 0; i < n1; ++i){
          if (core_1[i] != NULL_NODE) {
            c1[j] = i;
            c2[j] = core_1[i];
            ++j;
          }
        }
    
      };
      VF2SubState *Clone(){
        return new VF2SubState(*this);
      };
      void BackTrack(){
        assert(core_len - orig_core_len <= 1);
        assert(added_node1 != NULL_NODE);

        if (orig_core_len < core_len) {

          if (in_1[added_node1] == core_len) in_1[added_node1] = 0;
          if (out_1[added_node1] == core_len)  out_1[added_node1] = 0;

          typename Graph::out_edge_iterator bNbrs,eNbrs;
          boost::tie(bNbrs,eNbrs) = boost::out_edges(added_node1,*g1);
          while(bNbrs!=eNbrs){
            unsigned int other = getOtherIdx(*g1,*bNbrs,added_node1);
            if (out_1[other] == core_len) out_1[other] = 0;
            if (in_1[other] == core_len) in_1[other] = 0;
            ++bNbrs;
          }

          unsigned int node2 = core_1[added_node1];
          if (in_2[node2] == core_len) in_2[node2] = 0;
          if (out_2[node2] == core_len) out_2[node2] = 0;

          boost::tie(bNbrs,eNbrs) = boost::out_edges(node2,*g2);
          while(bNbrs!=eNbrs){
            unsigned int other = getOtherIdx(*g2,*bNbrs,node2);
            if (out_2[other] == core_len) out_2[other] = 0;
            if (in_2[other] == core_len) in_2[other] = 0;
            ++bNbrs;
          }

          core_1[added_node1] = NULL_NODE;
          core_2[node2] = NULL_NODE;

          core_len = orig_core_len;
          added_node1 = NULL_NODE;
        }
      };
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
    bool match(int *pn, node_id c1[], node_id c2[], SubState &s)
    {
      if (s.IsGoal() ) { 
        s.GetCoreSet(c1, c2);
        if(s.MatchChecks(c1,c2)) {
          *pn=s.CoreLen();
          return true;
        }
      }

      if (s.IsDead())
        return false;
      //std::cerr<<"  > match: "<<*pn<<" "<<&s<<std::endl;
      node_id n1=NULL_NODE, n2=NULL_NODE;
      bool found=false;
      while (!found && s.NextPair(&n1, &n2, n1, n2)) {
        //std::cerr<<"           "<<n1<<","<<n2<<std::endl;
        if (s.IsFeasiblePair(n1, n2)){
          SubState *s1=s.Clone();
          s1->AddPair(n1, n2);
          found=match(pn, c1, c2, *s1);
          s1->BackTrack();
          delete s1;
        }
      }
      //std::cerr<<"  < returning: "<<found<<" "<<*pn<<" "<<&s<<std::endl;
      return found;
    }

    /*-------------------------------------------------------------
     * static bool match(c1, c2, vis, usr_data, pcount)
     * Visits all the matchings between two graphs,  starting
     * from state s.
     * Returns true if the caller must stop the visit.
     * Stops when there are no more matches
     *
     ------------------------------------------------------------*/
    template <class SubState,class DoubleBackInsertionSequence>
    bool match(node_id c1[], node_id c2[], SubState &s, DoubleBackInsertionSequence &res) {
      if (s.IsGoal()){
        s.GetCoreSet(c1, c2);
        if(s.MatchChecks(c1,c2)) {
          typename DoubleBackInsertionSequence::value_type newSeq;
          for(unsigned int i=0;i<s.CoreLen();++i){
            newSeq.push_back(std::pair<int,int>(c1[i],c2[i]));
          }
          res.push_back(newSeq);
        }
        return false;
      }

      if (s.IsDead())
        return false;

      node_id n1=NULL_NODE, n2=NULL_NODE;
      while (s.NextPair(&n1, &n2, n1, n2)) {
        if (s.IsFeasiblePair(n1, n2)){
          SubState *s1=s.Clone();
          s1->AddPair(n1, n2);
          if (match(c1, c2, *s1,res)){
            s1->BackTrack(); 
            delete s1;
            return true;
          }
          else {
            s1->BackTrack(); 
            delete s1;
          }
        }
      }
      return false;
    }
  }; //end of namespace detail

  template <  class Graph
              , class VertexLabeling    // binary predicate
              , class EdgeLabeling      // binary predicate
              , class MatchChecking      // binary predicate
              , class BackInsertionSequence   // contains std::pair<vertex_descriptor,vertex_descriptor>
              >
  bool vf2(const Graph &g1,const Graph &g2,
           VertexLabeling& vertex_labeling,
           EdgeLabeling& edge_labeling,
           MatchChecking& match_checking,
           BackInsertionSequence& F){
    detail::VF2SubState<const Graph,VertexLabeling,EdgeLabeling,MatchChecking> s0(&g1,&g2,vertex_labeling,
                                                                                  edge_labeling,match_checking,false);
    detail::node_id *ni1 = new detail::node_id[num_vertices(g1)];
    detail::node_id *ni2 = new detail::node_id[num_vertices(g2)];
    int n=0;
    
    F.clear();
    F.resize(0);
    if(match(&n,ni1,ni2,s0)){
      for(unsigned int i=0;i<num_vertices(g1);i++){
        F.push_back(std::pair<int,int>(ni1[i],ni2[i]));
      }
    }
    delete [] ni1;
    delete [] ni2;
    
    return !F.empty();
  };
  template <  class Graph
              , class VertexLabeling    // binary predicate
              , class EdgeLabeling      // binary predicate
              , class MatchChecking      // binary predicate
              , class DoubleBackInsertionSequence   // contains a back insertion sequence
              >
  bool vf2_all(const Graph& g1, const Graph& g2,
               VertexLabeling& vertex_labeling,
               EdgeLabeling& edge_labeling,
               MatchChecking& match_checking,
               DoubleBackInsertionSequence& F) {
    detail::VF2SubState<const Graph,VertexLabeling,EdgeLabeling,MatchChecking> s0(&g1,&g2,vertex_labeling,
                                                                                  edge_labeling,match_checking,false);
    detail::node_id *ni1 = new detail::node_id[num_vertices(g1)];
    detail::node_id *ni2 = new detail::node_id[num_vertices(g2)];
    
    F.clear();
    F.resize(0);

    match(ni1,ni2,s0,F);

    delete [] ni1;
    delete [] ni2;
    
    return !F.empty();
  };
} // end of namespace boost
#endif

