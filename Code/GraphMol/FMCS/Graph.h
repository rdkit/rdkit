// graph topology in terms of indeces in source molecule
#pragma once
#include <boost/graph/adjacency_list.hpp>

namespace RDKit
{
 namespace FMCS
 {
/*
    template<typename T>
    class TMatrix // for scalar value types ! including bool with special STL implementation (no reference to item - bitset used)
    {
        size_t          Size;
        std::vector<T>  Data;
    public:
        inline TMatrix(size_t size=0) : Size(size), Data(cx*cy) {}
        inline size_t getSize()const { return Size;}
        inline bool   empty   ()const { return Data.empty();}
        inline void   resize(size_t size) { Data.resize(size*size); Size = size;}
        inline void set(size_t row, size_t col, T val) { Data[row*Size + col] = val;}
        inline T at(size_t row, size_t col)     { return Data[row*Size + col];}
        inline T at(size_t row, size_t col)const{ return Data[row*Size + col];}
    };
    typedef TMatrix<bool>    AdjacencyMatrix;   // actual graph adjacency matrix, index based
*/

    typedef unsigned AtomIdx_t;
    typedef unsigned BondIdx_t;
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, AtomIdx_t, BondIdx_t> Graph_t; 

    class Graph : public Graph_t
    {
    public:
        void addAtom(unsigned atom)
        {
            Graph::vertex_descriptor which = boost::add_vertex(*this);
            (*this)[which] = atom;
        }
        void addBond(unsigned bond, unsigned beginAtom, unsigned endAtom)
        {
            bool res;
            Graph_t::edge_descriptor which;
            boost::tie(which, res) = boost::add_edge(beginAtom, endAtom, *this);
            (*this)[which] = bond;
        }
    };

 }}