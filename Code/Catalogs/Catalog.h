//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef __RD_CATALOG_H__
#define __RD_CATALOG_H__

// Boost graph stuff
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION >= 104000
#include <boost/property_map/property_map.hpp>
#else
#include <boost/property_map.hpp>
#endif


// for some typedefs
#include <RDGeneral/types.h>
#include <RDGeneral/StreamOps.h>

namespace RDCatalog {
  const int versionMajor=1;
  const int versionMinor=0;
  const int versionPatch=0;
  const int endianId=0xDEADBEEF;
  
  //-----------------------------------------------------------------------------
  //! abstract base class for a catalog object
  template <class entryType, class paramType>
  class Catalog {
  public:
    //------------------------------------
    Catalog() : d_fpLength(0), dp_cParams(0) {};

    //------------------------------------
    virtual ~Catalog(){
      delete dp_cParams;
    }
    
    //------------------------------------
    //! return a serialized form of the Catalog as an std::string
    virtual std::string Serialize() const = 0;
    
    //------------------------------------
    //! adds an entry to the catalog
    /*!

      \param entry          the entry to be added
      \param updateFPLength (optional) if this is true, our internal
      fingerprint length will also be updated.
      
    */
    virtual unsigned int addEntry(entryType *entry, bool updateFPLength = true) = 0;
    
    //------------------------------------
    //! returns a particular entry in the Catalog
    virtual const entryType* getEntryWithIdx(unsigned int idx) const = 0;
    
    //------------------------------------
    //! returns the number of entries
    virtual unsigned int getNumEntries() const = 0;
    
    //------------------------------------
    //! returns the length of our fingerprint
    unsigned int getFPLength() const {return d_fpLength;}
    
    //------------------------------------
    //! sets our fingerprint length
    void setFPLength(unsigned int val) {d_fpLength = val;}
    
    //------------------------------------
    //! sets our parameters by copying the \c params argument
    void setCatalogParams(paramType *params) {
      PRECONDITION(params,"bad parameter object");
      //if we already have a paramter object throw an exception
      PRECONDITION(!dp_cParams,"A parameter object already exists on the catalog" );
      /*
        if (dp_cParams) {
        // we already have parameter object on the catalog
        // can't overwrite it
        PRECONDITION(0, "A parameter object already exist on the catalog");
        }*/
      dp_cParams = new paramType(*params);
    }
    
    //------------------------------------
    //! returns a pointer to our parameters
    const paramType *getCatalogParams() const { return dp_cParams;}

  protected:
    // this is the ID that will be assigned to the next entry 
    // added to the catalog - need not be same as the number of entries 
    // in the catalog and does not correspond with the
    // id of the entry in the catalog.
    // this is more along the lines of bitId
    unsigned int d_fpLength; //!< the length of our fingerprint
    paramType *dp_cParams;   //!< our params object

  };



  //-----------------------------------------------------------------------------
  //! A Catalog with a hierarchical structure
  /*!

    The entries of a HierarchCatalog are arranged in a directed graph

    <b>The difference between <i>Indices</i> and <i>Bit Ids</i></b>
    
    A HierarchCatalog may contain more entries than the user is actually
    interested in.  For example a HierarchCatalog constructed to contain
    orders 5 through 8 may well contain information about orders 1-5,
    in order to facilitate some search optimizations.

    - <i>Bit Ids</i> refer to the "interesting" bits.  
    So, in the above example, Bit Id \c 0 will be the first entry
    with order 5.
    - <i>Indices</i> refer to the underlying structure of the catalog.
    So, in the above example, the entry with index \c 0 will be
    the first entry with order 1.

  */
  template <class entryType, class paramType, class orderType> 
  class HierarchCatalog : public Catalog <entryType, paramType> {
    // the entries in the catalog can be traversed using the edges
    // in a desired order
  public:
    //! used by the BGL to set up the node properties in our graph
    struct vertex_entry_t {
      enum { num=1003 };
      typedef boost::vertex_property_tag kind;
    };
    typedef boost::property<vertex_entry_t, entryType *> EntryProperty;
    
    //! the type of the graph itself:
    typedef boost::adjacency_list<boost::vecS,
					 boost::vecS, // FIX: should be using setS for edges so that parallel edges are never added (page 225 BGL book)
      // but that seems result in compile errors
      boost::bidirectionalS,
	     EntryProperty> CatalogGraph;
    
    typedef boost::graph_traits<CatalogGraph> CAT_GRAPH_TRAITS;
    typedef typename CAT_GRAPH_TRAITS::vertex_iterator VER_ITER;
    typedef std::pair<VER_ITER, VER_ITER> ENT_ITER_PAIR;
    typedef typename CAT_GRAPH_TRAITS::adjacency_iterator DOWN_ENT_ITER;
    typedef std::pair<DOWN_ENT_ITER, DOWN_ENT_ITER> DOWN_ENT_ITER_PAIR;
    
    //------------------------------------
    HierarchCatalog<entryType, paramType, orderType>() {};
    
    //------------------------------------
    //! Construct by making a copy of the input \c params object
    HierarchCatalog<entryType, paramType, orderType>(paramType *params) : Catalog<entryType,paramType>() {
      this->setCatalogParams(params);
    }

    //------------------------------------
    //! Construct from a \c pickle (a serialized form of the HierarchCatalog)
    HierarchCatalog<entryType, paramType, orderType>(const std::string &pickle) {
      this->initFromString(pickle);
    }
    
    //------------------------------------
    ~HierarchCatalog() {
      destroy();
    }

    //------------------------------------
    //! serializes this object to a stream
    void toStream(std::ostream &ss) const {
      PRECONDITION(this->getCatalogParams(),"NULL parameter object");
      
      // the i/o header:
      RDKit::streamWrite(ss,endianId);
      RDKit::streamWrite(ss,versionMajor);
      RDKit::streamWrite(ss,versionMinor);
      RDKit::streamWrite(ss,versionPatch);

      // information about the catalog itself:
      int tmpUInt;
      tmpUInt = this->getFPLength();
      RDKit::streamWrite(ss,tmpUInt);
      tmpUInt = this->getNumEntries();
      RDKit::streamWrite(ss,tmpUInt);

      //std::cout << ">>>>-------------------------------" << std::endl;
      //std::cout << "\tlength: " << getFPLength() << " " << getNumEntries() << std::endl;

      // add the params object:
      this->getCatalogParams()->toStream(ss);
      //std::cout << "\tparams: " << getCatalogParams()->getLowerFragLength();
      //std::cout << " " << getCatalogParams()->getUpperFragLength();
      //std::cout << " " << getCatalogParams()->getNumFuncGroups();
      //std::cout << std::endl;
      
      // write the entries in order:
      for(unsigned int i=0;i<getNumEntries();i++){
        this->getEntryWithIdx(i)->toStream(ss);
      }

      // finally the adjacency list:
      for(unsigned int i=0;i<getNumEntries();i++){
        RDKit::INT_VECT children=this->getDownEntryList(i);
        tmpUInt = children.size();
        RDKit::streamWrite(ss,tmpUInt);
        for(RDKit::INT_VECT::const_iterator ivci=children.begin();
            ivci!=children.end();
            ivci++){
	  RDKit::streamWrite(ss,*ivci);
	}
      }
    }

    //------------------------------------
    //! serializes this object and returns the resulting \c pickle
    std::string Serialize() const {
      std::stringstream ss(std::ios_base::binary|std::ios_base::out|std::ios_base::in);
      this->toStream(ss);
      return ss.str();
    }

    //------------------------------------
    //! fills the contents of this object from a stream containing a \c pickle
    void initFromStream(std::istream &ss) {
      int tmpInt;
      // FIX: at the moment we ignore the header info:
      RDKit::streamRead(ss,tmpInt);
      RDKit::streamRead(ss,tmpInt);
      RDKit::streamRead(ss,tmpInt);
      RDKit::streamRead(ss,tmpInt);
      
      unsigned int tmpUInt;
      RDKit::streamRead(ss,tmpUInt);// fp length
      this->setFPLength(tmpUInt);
      
      unsigned int numEntries;
      RDKit::streamRead(ss,numEntries);
      //std::cout << "<<<-------------------------------" << std::endl;
      //std::cout << "\tlength: " << getFPLength() << " " << numEntries << std::endl;


      // grab the params:
      paramType *params = new paramType();
      params->initFromStream(ss);
      this->setCatalogParams(params);

      //std::cout << "\tparams: " << getCatalogParams()->getLowerFragLength();
      //std::cout << " " << getCatalogParams()->getUpperFragLength();
      //std::cout << " " << getCatalogParams()->getNumFuncGroups();
      //std::cout << std::endl;

      // now all of the entries:
      for(unsigned int i=0;i<numEntries;i++){
        entryType *entry = new entryType();
        entry->initFromStream(ss);
        this->addEntry(entry,false);
      }

      // and, finally, the adjacency list:
      for(unsigned int i=0;i<numEntries;i++){
        unsigned int nNeighbors;
        RDKit::streamRead(ss,nNeighbors);
        for(unsigned int j=0;j<nNeighbors;j++){
          RDKit::streamRead(ss,tmpInt);
          this->addEdge(i,tmpInt);
	}
      }
    }
    
    //------------------------------------
    unsigned int getNumEntries() const {
      return boost::num_vertices(d_graph);
    }

    //------------------------------------
    //! fills the contents of this object from a string containing a \c pickle
    void initFromString(const std::string &text){
      std::stringstream ss(std::ios_base::binary|std::ios_base::out|std::ios_base::in);
      // initialize the stream:
      ss.write(text.c_str(),text.length());
      // now start reading out values:
      this->initFromStream(ss);
    }

    //------------------------------------
    //! add a new entry to the catalog
    /*!

      \param entry          the entry to be added
      \param updateFPLength (optional) if this is true, our internal
      fingerprint length will also be updated.

    */
    unsigned int addEntry(entryType *entry, bool updateFPLength = true){ 
      PRECONDITION(entry,"bad arguments");
      if (updateFPLength) {
        unsigned int fpl = this->getFPLength();
        entry->setBitId(fpl);
        fpl++;
        this->setFPLength(fpl);
      }
      unsigned int eid = boost::add_vertex(EntryProperty(entry), d_graph);
      orderType etype = entry->getOrder();
      // REVIEW: this initialization is not required: the STL map, in
      // theory, will create a new object when operator[] is called
      // for a new item
      if (d_orderMap.find(etype) == d_orderMap.end()) {
        RDKit::INT_VECT nets;
        d_orderMap[etype] = nets;
      }
      d_orderMap[etype].push_back(eid);
      return eid;
    }

    //------------------------------------
    //! adds an edge between two entries in the catalog
    /*!
      Since we are using a bidirectional graph - the order in 
      which the ids are supplied here makes a difference

      \param id1 index of the edge's beginning
      \param id2 index of the edge's end

    */
    void addEdge(unsigned int id1, unsigned int id2) {
      unsigned int nents = getNumEntries();
      RANGE_CHECK(0, id1, nents-1);
      RANGE_CHECK(0, id2, nents-1);
      // FIX: if we boost::setS for the edgeList BGL will
      // do the checking for duplicity (parallel edges)
      // But for reasons unknown setS results in compile
      // errors while using adjacent_vertices.
      typename CAT_GRAPH_TRAITS::edge_descriptor edge;
      bool found;
      boost::tie(edge,found) = boost::edge(boost::vertex(id1,d_graph),
                                           boost::vertex(id2,d_graph),
                                           d_graph);
      if (!found) {
        boost::add_edge(id1, id2, d_graph);
      }
    }

    //------------------------------------
    //! returns a pointer to our entry with a particular index 
    const entryType *getEntryWithIdx(unsigned int idx) const {
      RANGE_CHECK(0,idx,getNumEntries()-1);
      int vd = boost::vertex(idx, d_graph);
      typename boost::property_map < CatalogGraph, vertex_entry_t>::const_type 
        pMap = boost::get(vertex_entry_t(), d_graph);
      return pMap[vd];
    }

    //------------------------------------
    //! returns a pointer to our entry with a particular bit ID
    const entryType *getEntryWithBitId(unsigned int idx) const {
      RANGE_CHECK(0,idx,this->getFPLength()-1);
      typename boost::property_map < CatalogGraph, vertex_entry_t>::const_type 
        pMap = boost::get(vertex_entry_t(), d_graph);
      const entryType *res=NULL;
      for(unsigned int i=idx;i<this->getNumEntries();i++){
        const entryType *e=pMap[i];
        if(e->getBitId()==static_cast<int>(idx)){
	  res=e;
          break;
        }
      }
      return res;
    }

    //------------------------------------
    //! returns the index of the entry with a particular bit ID
    int getIdOfEntryWithBitId(unsigned int idx) const {
      RANGE_CHECK(0,idx,this->getFPLength()-1);
      typename boost::property_map < CatalogGraph, vertex_entry_t>::const_type 
        pMap = boost::get(vertex_entry_t(), d_graph);
      int res=-1;
      for(unsigned int i=idx;i<this->getNumEntries();i++){
        const entryType *e=pMap[i];
        if(static_cast<unsigned int>(e->getBitId())==idx){
          res=i;
	  break;
        }
      }
      return res;
    }

    //------------------------------------
    //! returns a list of the indices of entries below the one passed in
    RDKit::INT_VECT getDownEntryList(unsigned int idx) const {
      RDKit::INT_VECT res;
      DOWN_ENT_ITER nbrIdx, endIdx;
      boost::tie(nbrIdx, endIdx) = boost::adjacent_vertices(idx, d_graph);
      while (nbrIdx != endIdx) {
        res.push_back(*nbrIdx);
        nbrIdx++;
      }
      //std::cout << res.size() << "\n";
      return res;
    }

    //------------------------------------
    //! returns a list of the indices that have a particular order 
    const RDKit::INT_VECT &getEntriesOfOrder(orderType ord) {
      return d_orderMap[ord];
    }

    //------------------------------------
    //! returns a list of the indices that have a particular order
    /*!
      \overload
    */
    const RDKit::INT_VECT &getEntriesOfOrder(orderType ord) const {
      typename std::map<orderType, RDKit::INT_VECT>::const_iterator elem;
      elem = d_orderMap.find(ord);
      CHECK_INVARIANT(elem!=d_orderMap.end()," catalog does not contain any entries of the order specified");
      return elem->second;
    }

    
  private:
    // graphs that store the entries in the catalog in a hierachical manner
    CatalogGraph d_graph;
    // a  map that maps the order type of entries in the catalog to 
    // a vector of vertex indices in the graphs above
    // e.g.  for a catalog with molecular fragments, the order of a fragment can 
    // simply be the number of bond in it. The list this oder maps to is all the
    // vertex ids of these fragment in the catalog that have this many bonds in them
    std::map<orderType, RDKit::INT_VECT> d_orderMap;

    //------------------------------------
    //! clear any memory that we've used
    void destroy() {
      ENT_ITER_PAIR entItP = boost::vertices(d_graph);
      typename boost::property_map < CatalogGraph, vertex_entry_t>::type 
        pMap = boost::get(vertex_entry_t(), d_graph);
      while (entItP.first != entItP.second) {
        delete pMap[*(entItP.first++)];
      }
    }


    
  };


  //-----------------------------------------------------------------------------
  //! a linear Catalog (analogous to an std::vector)
  /*!
    Here there is no particular hierarchy, simply a
    collection of entries.
  */
  template <class entryType, class orderType> 
  class LinearCatalog : public Catalog <entryType, orderType> {
    // here there is no particular hierarchy of entries
    // we simply model it as a vector of entries
    // FIX: for retrieval purposes a better model map be std::map

  public:
    std::string Serialize();
    
    unsigned int addEntry(entryType *entry, bool updateFPLength = true);
    
    const entryType *getEntryWithIdx(unsigned int idx) const;

  private:
    std::vector<entryType*> d_vector;
  };
}

#endif
