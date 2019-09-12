/*
 * This file is part of the RingDecomposerLib, licensed
 * under BSD New license (see LICENSE in the root directory).
 * Copyright (c) 2016
 * University of Hamburg, ZBH - Center for Bioinformatics
 * Niek Andresen, Florian Flachsenberg, Matthias Rarey
 * 
 * Please cite:
 * 
 * Kolodzik, A.; Urbaczek, S.; Rarey, M.
 * Unique Ring Families: A Chemically Meaningful Description
 * of Molecular Ring Topologies.
 * J. Chem. Inf. Model., 2012, 52 (8), pp 2013-2021
 * 
 * Flachsenberg, F.; Andresen, N.; Rarey, M.
 * RingDecomposerLib: An Open-Source Implementation of
 * Unique Ring Families and Other Cycle Bases.
 * J. Chem. Inf. Model., 2017, 57 (2), pp 122-126
 */

/**
 * @file RingDecomposerLib.h
 * @brief This file contains the API of the RingDecomposerLib
 * library.
 */

#ifndef RING_DECOMPOSER_LIB_H
#define RING_DECOMPOSER_LIB_H

#if (defined( _WIN32 ) && defined( _MSC_VER ) )
    /* Win32 & MS VC ++ */
    #define RDL_API __declspec(dllexport)
#elif __GNUC__ >= 4 || defined(__clang__)
    #define RDL_API __attribute__((visibility("default")))
#else
    #define RDL_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Structure representing a calculation result.
 *
 * This is the central structure, that is used to store
 * the calculation results. Almost all functions
 * in this header work on it.
 */
typedef struct RDL_data RDL_data;

/**
 * @brief Datastructure representing a graph for calculations.
 *
 * Build graph data structure that is used for the calculation.
 * The graph is stored as an adjacency list.
 * The vertices are numbered from 0 to |V|-1. Call the following functions:
 *
 * \ref RDL_initNewGraph(unsigned V) to initialize a new graph with V vertices
 * (returns \ref RDL_graph *)
 *
 * \ref RDL_addUEdge(RDL_graph *, RDL_node from, RDL_node to) to add a new (undirected) edge
 * from the vertex with index "from" to the vertex with index "to".
 *
 * Then \ref RDL_calculate can be called on it.
 *
 */
typedef struct RDL_graph RDL_graph;

/**
 * @brief A node is represented by its index in the graph (0 to |V|-1).
 */
typedef unsigned RDL_node;

/** @brief An edge is represented by an array of size two
 * containing the adjacent nodes.
 */
typedef RDL_node RDL_edge[2];

/**
 * @brief error levels for custom logging functions
*/
typedef enum RDL_ERROR_LEVEL {
  RDL_DEBUG,
  RDL_WARNING,
  RDL_ERROR,
  RDL_INITIALIZE_OUTPUT
} RDL_ERROR_LEVEL;

/*
 * @brief output function for writing errors/warnings
 *
 * @param level `RDL_ERROR_LEVEL` specifying the level of the message
 * @param m message associated with the event
 */
typedef void (*RDL_outputFunction) (RDL_ERROR_LEVEL level, const char* m, ...);

/** @brief the output function for warnings and errors */
RDL_API extern RDL_outputFunction RDL_outputFunc;

/**
 * @brief Set the output function for warning and error messages.
 * @param func: callback function for errors and warnings
 * @warning This function sets a `static` variable. When using multiple
 * threads, only one thread should set the output function and the
 * output function should be thread-safe.
 *
 * Default is no output at all (\ref RDL_writeNothing).
 * Call `RDL_setOutputFunction(RDL_writeToStderr)` to enable output
 * of warning and error messages to `stderr`.
 */
RDL_API
void RDL_setOutputFunction(RDL_outputFunction func);

/**
 * @brief An output function for writing everything to stderr
 */
RDL_API
void RDL_writeToStderr(RDL_ERROR_LEVEL level, const char* fmt, ...);

/**
 * @brief No output function
 */
RDL_API
void RDL_writeNothing(RDL_ERROR_LEVEL level, const char* fmt, ...);

/** @brief Invalid result indicator */
RDL_API extern const unsigned RDL_INVALID_RESULT;
/** @brief Duplicate edge indicator */
RDL_API extern const unsigned RDL_DUPLICATE_EDGE;
/** @brief No ringsystem indicator */
RDL_API extern const unsigned RDL_NO_RINGSYSTEM;
/** @brief Invalid number of RCs */
RDL_API extern const double RDL_INVALID_RC_COUNT;

/**
 * @brief Initializes a new \ref RDL_graph
 * @param nof_nodes the number of nodes in the graph
 * @return pointer to the new \ref RDL_graph structure.
 */
RDL_API
RDL_graph *RDL_initNewGraph(unsigned nof_nodes);

/**
 * @brief Delete \ref RDL_graph.
 * @param graph pointer to \ref RDL_graph to delete
 * @note You don't have to delete the graph after successful calculation!
 */
RDL_API
void RDL_deleteGraph(RDL_graph *graph);

/**
 * @brief Adds an undirected edge to the graph.
 * @param graph pointer to the \ref RDL_graph that the edge will to be added to
 * @param node1, node2 pair of \ref RDL_node which the edge is going to connect
 * @returns internal index of edge if successful,
 * \ref RDL_INVALID_RESULT on failures (invalid index, loop),
 * \ref RDL_DUPLICATE_EDGE if the edge was already present
 */
RDL_API
unsigned RDL_addUEdge(RDL_graph *graph, RDL_node node1, RDL_node node2);

/**
 * @brief Calculates the \ref RDL_data structure of the given graph and returns it.
 * @param input_graph The graph for calculation.
 * @return pointer to the \ref RDL_data containing the calculation result, `NULL` on failure
 * @note The \ref RDL_graph has to be an undirected graph.
 * @note Takes ownership of the \ref RDL_graph, if calculation is successful!
 */
RDL_API
RDL_data *RDL_calculate(RDL_graph *input_graph);

/**
 * @brief Deletes \ref RDL_data from memory, including the \ref RDL_graph.
 * @param data pointer to the \ref RDL_data that is going to be deleted
 * @note also deletes the associated \ref RDL_graph
 */
RDL_API
void RDL_deleteData(RDL_data *data);

/**
 * @brief Returns the number of URFs.
 * @param data Pointer to the \ref RDL_data of which the number of URF is requested
 * @return The number of URFs, \ref RDL_INVALID_RESULT on failure
 */
RDL_API
unsigned RDL_getNofURF(const RDL_data *data);

/**
 * @brief Returns the number of RCFs.
 * @param data Pointer to the \ref RDL_data of which the number of RCF is requested
 * @return The number of RCFs, \ref RDL_INVALID_RESULT on failure
 */
RDL_API
unsigned RDL_getNofRCF(const RDL_data *data);

/**
 * @brief Returns the weight of each cycle in the URF identified by its index.
 * @param data pointer to the \ref RDL_data holding the URFs
 * @param index the index of the URF
 * @return the weight of the URF, \ref RDL_INVALID_RESULT on failure (if `index` is out of range)
 */
RDL_API
unsigned RDL_getWeightForURF(const RDL_data *data, unsigned index);

/**
 * @brief Returns the weight of each cycle in the RCF identified by its index.
 * @param data pointer to the \ref RDL_data holding the URFs
 * @param index the index of the RCF
 * @return the weight of the RCF, \ref RDL_INVALID_RESULT on failure (if `index` is out of range)
 */
RDL_API
unsigned RDL_getWeightForRCF(const RDL_data *data, unsigned index);

/**
 * @brief Gives the nodes of an URF identified with its index in an array of \ref RDL_node.
 * @return the number of nodes in the URF, \ref RDL_INVALID_RESULT on failure
 * @param data pointer to the \ref RDL_data holding the calculation results.
 * @param index the index of the URF
 * @param RDL_node_array_ptr pointer that points to the result array (declare
 * \ref RDL_node * and give address as parameter)
 * @note result has to be deallocated using `free(*RDL_node_array_ptr)`.
 */
RDL_API
unsigned RDL_getNodesForURF(const RDL_data *data, unsigned index,
                            RDL_node **RDL_node_array_ptr);

/**
 * @brief Gives the edges of an URF identified with its index.
 * @return the number of edges in this URF, \ref RDL_INVALID_RESULT on failure
 * @param data pointer to the \ref RDL_data holding the calculation results.
 * @param index the index of the URF
 * @param RDL_edge_array_ptr pointer that points to the result array (declare
 * \ref RDL_edge * and give address as parameter)
 * @note Result has to be deallocated using `free(*RDL_edge_array_ptr)`
 * Gives an array of edges where an edge is represented by two \ref RDL_node that it
 * connects.
 */
RDL_API
unsigned RDL_getEdgesForURF(const RDL_data *data, unsigned index,
                            RDL_edge **RDL_edge_array_ptr);

/**
 * @brief Gives the nodes of an RCF identified with its index in an array of \ref RDL_node.
 * @return the number of nodes in the RCF, \ref RDL_INVALID_RESULT on failure
 * @param data pointer to the \ref RDL_data holding the calculation results.
 * @param index the index of the RCF
 * @param RDL_node_array_ptr pointer that points to the result array (declare
 * \ref RDL_node * and give address as parameter)
 * @note result has to be deallocated using `free(*RDL_node_array_ptr)`.
 */
RDL_API
unsigned RDL_getNodesForRCF(const RDL_data *data, unsigned index,
                            RDL_node **RDL_node_array_ptr);

/**
 * @brief Gives the edges of an RCF identified with its index.
 * @return the number of edges in this RCF, \ref RDL_INVALID_RESULT on failure
 * @param data pointer to the \ref RDL_data holding the calculation results.
 * @param index the index of the RCF
 * @param RDL_edge_array_ptr pointer that points to the result array (declare
 * \ref RDL_edge * and give address as parameter)
 * @note Result has to be deallocated using `free(*RDL_edge_array_ptr)`
 * Gives an array of edges where an edge is represented by two \ref RDL_node that it
 * connects.
 */
RDL_API
unsigned RDL_getEdgesForRCF(const RDL_data *data, unsigned index,
                            RDL_edge **RDL_edge_array_ptr);

/**
 * @brief Returns the number of URFs that contain the given node.
 * @param data pointer to the \ref RDL_data holding the calculation results.
 * @param node the \ref RDL_node to look for in the URFs
 * @return number of URFs containing the node, \ref RDL_INVALID_RESULT on failure
 * @note This functions internally enumerates all nodes for each URF.
 * For repeated checks with different \ref RDL_node it is recommended to
 * enumerate them with \ref RDL_getNodesForURF and check against this array.
 */
RDL_API
unsigned RDL_getNofURFContainingNode(const RDL_data *data, RDL_node node);

/**
 * @brief Returns the number of URFs that contain the \ref RDL_edge defined
 * by the two given nodes
 * @param data Pointer to the \ref RDL_data storing the calculation results.
 * @param node1, node2 pair of \ref RDL_node connected by the edge
 * @return number of URFs containing the edge, \ref RDL_INVALID_RESULT on failure
 * @note This functions internally enumerates all edges for each URF.
 * For repeated checks with different \ref RDL_node it is recommended to
 * enumerate them with \ref RDL_getEdgesForURF and check against this array.
 */
RDL_API
unsigned RDL_getNofURFContainingEdge(const RDL_data *data, RDL_node node1,
                                     RDL_node node2);

/**
 * @brief Returns the number of RCFs that contain the given node.
 * @param data pointer to the \ref RDL_data holding the calculation results.
 * @param node the \ref RDL_node to look for in the RCFs
 * @return number of RCFs containing the node, \ref RDL_INVALID_RESULT on failure
 * @note This functions internally enumerates all nodes for each RCF.
 * For repeated checks with different \ref RDL_node it is recommended to
 * enumerate them with \ref RDL_getNodesForRCF and check against this array.
 */
RDL_API
unsigned RDL_getNofRCFContainingNode(const RDL_data *data, RDL_node node);

/**
 * @brief Returns the number of RCFs that contain the \ref RDL_edge defined
 * by the two given nodes
 * @param data Pointer to the \ref RDL_data storing the calculation results.
 * @param node1, node2 pair of \ref RDL_node connected by the edge
 * @return number of RCFs containing the edge, \ref RDL_INVALID_RESULT on failure
 * @note This functions internally enumerates all edges for each RCF.
 * For repeated checks with different \ref RDL_node it is recommended to
 * enumerate them with \ref RDL_getEdgesForRCF and check against this array.
 */
RDL_API
unsigned RDL_getNofRCFContainingEdge(const RDL_data *data, RDL_node node1,
                                     RDL_node node2);

/**
 * @brief Gives all URFs containing the node.
 * @return the number of URFs containing the node, \ref RDL_INVALID_RESULT on failure
 * @param data pointer to the \ref RDL_data holding the calculation results
 * @param node the \ref RDL_node
 * @param RDL_ids_ptr pointer that points to the result array of integers
 * containing all indices of URFs containing the node. (declare `unsigned *` and give
 * address as parameter)
 * @note The array RDL_ids_ptr has to be to be deallocated with
 * `free(*RDL_ids_ptr)`
 * @note This functions internally enumerates all nodes for each URF.
 * For repeated queries with different \ref RDL_node it is recommended to
 * enumerate them with \ref RDL_getNodesForURF once and use this array.
 */
RDL_API
unsigned RDL_getURFsContainingNode(const RDL_data *data, RDL_node node,
                                   unsigned **RDL_ids_ptr);

/**
 * @brief Gives all URFs containing the edge.
 * @return the number of URFs containing the edge, \ref RDL_INVALID_RESULT on failure
 * @param data pointer to the \ref RDL_data holding the calculation results
 * @param node1 the first \ref RDL_node of the edge
 * @param node2 the the second \ref RDL_node of the edge
 * @param RDL_ids_ptr pointer that points to the result array of integers
 * containing all indices of URFs containing the edge. (declare  `unsigned *` and give
 * address as parameter)
 * @note The array RDL_ids_ptr has to be to be deallocated with
 * `free(*RDL_ids_ptr)`
 * @note This functions internally enumerates all edges for each URF.
 * For repeated queries with different \ref RDL_node it is recommended to
 * enumerate them with \ref RDL_getEdgesForURF once and use this array.
 */
RDL_API
unsigned RDL_getURFsContainingEdge(const RDL_data *data, RDL_node node1,
                                   RDL_node node2, unsigned **RDL_ids_ptr);

/**
 * @brief Gives all RCFs containing the node.
 * @return the number of RCFs containing the node, \ref RDL_INVALID_RESULT on failure
 * @param data pointer to the \ref RDL_data holding the calculation results
 * @param node the \ref RDL_node
 * @param RDL_ids_ptr pointer that points to the result array of integers
 * containing all indices of URFs containing the node. (declare `unsigned *` and give
 * address as parameter)
 * @note The array RDL_ids_ptr has to be to be deallocated with
 * `free(*RDL_ids_ptr)`
 * @note This functions internally enumerates all nodes for each RCF.
 * For repeated queries with different \ref RDL_node it is recommended to
 * enumerate them with \ref RDL_getNodesForRCF once and use this array.
 */
RDL_API
unsigned RDL_getRCFsContainingNode(const RDL_data *data, RDL_node node,
                                   unsigned **RDL_ids_ptr);

/**
 * @brief Gives all RCFs containing the edge.
 * @return the number of RCFs containing the edge, \ref RDL_INVALID_RESULT on failure
 * @param data pointer to the \ref RDL_data holding the calculation results
 * @param node1 the first \ref RDL_node of the edge
 * @param node2 the the second \ref RDL_node of the edge
 * @param RDL_ids_ptr pointer that points to the result array of integers
 * containing all indices of RCFs containing the edge. (declare  `unsigned *` and give
 * address as parameter)
 * @note The array RDL_ids_ptr has to be to be deallocated with
 * `free(*RDL_ids_ptr)`
 * @note This functions internally enumerates all edges for each RCF.
 * For repeated queries with different \ref RDL_node it is recommended to
 * enumerate them with \ref RDL_getEdgesForRCF once and use this array.
 */
RDL_API
unsigned RDL_getRCFsContainingEdge(const RDL_data *data, RDL_node node1,
                                   RDL_node node2, unsigned **RDL_ids_ptr);

/**
 * This structure holds a cycle (or ring) in the graph.
 * It is essentially an array of edges which are represented by pairs of nodes.
 */
typedef struct RDL_cycle {
    /** array of \ref RDL_edge in the cycle */
    RDL_edge *edges;
    /** weight of the cycle (length of the array) */
    unsigned weight;
    /** URF of the cycle */
    unsigned urf;
    /** RCF of the cycle */
    unsigned rcf;
} RDL_cycle;

/**
 * @brief Free memory of \ref RDL_cycle
 * @param cycle The \ref RDL_cycle to delete.
 */
RDL_API
void RDL_deleteCycle(RDL_cycle *cycle);

/**
 * @brief Iterator for relevant cycles.
 *
 * Always check if the iterator is at end BEFORE incrementing
 * or accessing it.
 *
 * Example usage:
 *
 *     RDL_cycle *cycle;
 *     RDL_cycleIterator *it = RDL_getRCyclesIterator(data);
 *     while(!RDL_cycleIteratorAtEnd(it)) {
 *       cycle = RDL_cycleIteratorGetCycle(it);
 *       <do something with the cycle>
 *       RDL_deleteCycle(cycle);
 *       RDL_cycleIteratorNext(it);
 *     }
 *     RDL_deleteCycleIterator(it);
 */
typedef struct RDL_cycleIterator RDL_cycleIterator;

/**
 * @brief Advance the cycle iterator by one.
 * @param it \ref RDL_cycleIterator
 * @returns returns the iterator `it` itself, if successful; `NULL` on error
 * @note DO NOT call without checking \ref RDL_cycleIteratorAtEnd
 */
RDL_API
RDL_cycleIterator* RDL_cycleIteratorNext(RDL_cycleIterator* it);

/**
 * @brief Get the cycle as \ref RDL_cycle
 * @param it \ref RDL_cycleIterator
 * @returns \ref RDL_cycle holding the cycle, `NULL` on error
 * @note The \ref RDL_cycle is constructed and returned as a pointer.
 * You have to delete it using \ref RDL_deleteCycle.
 * @note DO NOT call without checking \ref RDL_cycleIteratorAtEnd
 */
RDL_API
RDL_cycle* RDL_cycleIteratorGetCycle(RDL_cycleIterator* it);

/**
 * @brief Check if iterator is at end (invalid)
 * @param it \ref RDL_cycleIterator
 * @returns 0 if not at end, a non-zero value otherwise
 */
RDL_API
int RDL_cycleIteratorAtEnd(RDL_cycleIterator* it);
/**
 * @brief Free memory of the cycle iterator
 * @param it \ref RDL_cycleIterator
 */
RDL_API
void RDL_deleteCycleIterator(RDL_cycleIterator* it);

/**
 * @brief Gives all relevant cycles of the URF with the given index.
 * @return the number of cycles found, \ref RDL_INVALID_RESULT on failure
 * @param data pointer to the \ref RDL_data holding the calculation results
 * @param index the index of the URF
 * @param RDL_cycle_array_ptr pointer that points to the result array of cycles
 * (declare \ref RDL_cycle ** and give address as parameter)
 * @note RDL_cycle_array_ptr has to be deallocated using
 * \ref RDL_deleteCycles(RDL_cycle **cycles, unsigned number)
 * @note Consider using \ref RDL_getRCyclesForURFIterator instead
 */
RDL_API
unsigned RDL_getRCyclesForURF(const RDL_data *data, unsigned index,
                              RDL_cycle ***RDL_cycle_array_ptr);

/**
 * @brief Get iterator for all relevant cycles of the URF with the given index.
 * @return \ref RDL_cycleIterator, NULL on failure
 * @param data pointer to the \ref RDL_data holding the calculation results
 * @param index the index of the URF
 * @note See \ref RDL_cycleIterator for an example.
 */
RDL_API
RDL_cycleIterator* RDL_getRCyclesForURFIterator(const RDL_data *data, unsigned index);

/**
 * @brief Gives all relevant cycles of the RCF with the given index.
 * @return the number of cycles found, \ref RDL_INVALID_RESULT on failure
 * @param data pointer to the \ref RDL_data holding the calculation results
 * @param index the index of the RCF
 * @param RDL_cycle_array_ptr pointer that points to the result array of cycles
 * (declare \ref RDL_cycle ** and give address as parameter)
 * @note RDL_cycle_array_ptr has to be deallocated using
 * \ref RDL_deleteCycles(RDL_cycle **cycles, unsigned number)
 * @note Consider using \ref RDL_getRCyclesForRCFIterator instead
 */
RDL_API
unsigned RDL_getRCyclesForRCF(const RDL_data *data, unsigned index,
                              RDL_cycle ***RDL_cycle_array_ptr);

/**
 * @brief Get iterator for all relevant cycles of the RCF with the given index.
 * @return \ref RDL_cycleIterator, NULL on failure
 * @param data pointer to the \ref RDL_data holding the calculation results
 * @param index the index of the RCF
 * @note See \ref RDL_cycleIterator for an example.
 */
RDL_API
RDL_cycleIterator* RDL_getRCyclesForRCFIterator(const RDL_data *data, unsigned index);

/**
 * @brief Gives a list of all relevant cycles
 * @return the number of cycles, \ref RDL_INVALID_RESULT on error
 * @param data pointer to the \ref RDL_data holding the calculation results
 * @param RDL_cycle_array_ptr pointer to the result array
 * (declare \ref RDL_cycle ** and give address as parameter)
 * The result is an array of cycles.
 * @note Result has to be deallocated using
 * \ref RDL_deleteCycles(RDL_cycle **cycles, unsigned number)
 * @note Consider using \ref RDL_getRCyclesIterator instead
 */
RDL_API
unsigned RDL_getRCycles(const RDL_data *data, RDL_cycle ***RDL_cycle_array_ptr);

/**
 * @brief Get iterator for all relevant cycles.
 * @return \ref RDL_cycleIterator, NULL on failure
 * @param data pointer to the \ref RDL_data holding the calculation results
 * @note See \ref RDL_cycleIterator for an example.
 */
RDL_API
RDL_cycleIterator* RDL_getRCyclesIterator(const RDL_data *data);

/**
 * @brief Gives the number of relevant cycles in this URF
 * @param data pointer to the \ref RDL_data holding the calculation results
 * @param index the index of the URF
 * @returns the number of relevant cycles in the URF or
 * \ref RDL_INVALID_RC_COUNT if the number is too large to calculate or on failure
 */
RDL_API
double RDL_getNofRCForURF(const RDL_data *data, unsigned index);

/**
 * @brief Gives the number of relevant cycles in this RCF.
 * @param data pointer to the \ref RDL_data holding the calculation results.
 * @param index id of the RCF
 * @returns the number of relevant cycles in the RCF or
 * \ref RDL_INVALID_RC_COUNT if the number is too large to calculate or on failure
 */
RDL_API
double RDL_getNofRCForRCF(const RDL_data *data, unsigned index);

/**
 * @brief Gives the number of relevant cycles
 * @param data pointer to the \ref RDL_data holding the calculation results
 * @returns the number of relevant cycles
 * \ref RDL_INVALID_RC_COUNT if the number is too large to calculate or on failure
 */
RDL_API
double RDL_getNofRC(const RDL_data *data);

/**
 * @brief Gives a set of cycles that forms a Minimal Cycle Basis of the graph.
 * @return the number of cycles returned (`|E|-|V|+1` for connected graphs),
 * RDL_INVALID_RESULT on error
 * @param data pointer to the \ref RDL_data holding the calculation results
 * @param RDL_cycle_array_ptr pointer that points to the result array
 * (declare \ref RDL_cycle ** and give address as parameter)
 * @note Result has to be deallocated using
 * \ref RDL_deleteCycles(RDL_cycle **cycles, unsigned number)
 * The result is an array of cycles.
 */
RDL_API
unsigned RDL_getSSSR(const RDL_data *data, RDL_cycle ***RDL_cycle_array_ptr);

/**
 * @brief Gives a list of relevant cycle prototypes (one for each RCF).
 * @return the number of prototypes, \ref RDL_INVALID_RESULT on error
 * @param data pointer to the \ref RDL_data holding the calculation results
 * @param RDL_cycle_array_ptr pointer to the result array
 * (declare \ref RDL_cycle ** and give address as parameter)
 * The result is an array of cycles.
 * @note Result has to be deallocated using
 * \ref RDL_deleteCycles(RDL_cycle **cycles, unsigned number)
 */
RDL_API
unsigned RDL_getRCPrototypes(const RDL_data *data,
                             RDL_cycle ***RDL_cycle_array_ptr);

/**
 * @brief Deallocates the structure given by \ref RDL_getRCycles(), \ref RDL_getSSSR(),
 * \ref RDL_getRCPrototypes() and \ref RDL_getRCyclesForURF(), if called on its result and
 * return value (the number of cycles)
 * @param cycles The array of cycles.
 * @param number The number of cycles in the array
 */
RDL_API
void RDL_deleteCycles(RDL_cycle **cycles, unsigned number);

/**
 * @brief Translates the results of \ref RDL_getRCycles(), \ref RDL_getSSSR(),
 * \ref RDL_getRCPrototypes() and \ref RDL_getRCyclesForURF() (arrays of \ref RDL_cycle) into an
 * array of cycles represented by arrays of `{0,1}^|E|` (bitsets).
 * @param data pointer to the \ref RDL_data holding the calculation results.
 * @param old_array The resulting structure of the functions named above
 * @param number The return value of the functions named above
 * (the number of cycles given)
 * @param RDL_cycle_array_ptr pointer to the result array
 * (declare char ** and give address as parameter)
 * @return The number of cycles given (same as the parameter 'number'), \ref RDL_INVALID_RESULT on failure
 * @note The initial structure still exists afterwards and still has to be
 * deleted.
 *
 * The resulting array has a 1 at
 * position i if edge i is part of the cycle or 0 otherwise. An edge is identified
 * by the position at which it was added to the graph starting at 0
 * (the return value of \ref RDL_addUEdge()).
 */
RDL_API
unsigned RDL_translateCycArray(const RDL_data *data,
                               RDL_cycle **old_array,
                               unsigned number,
                               char ***RDL_cycle_array_ptr);

/**
 * @brief Deallocates the structure given by \ref RDL_translateCycArray(), if called on
 * its result and return value (the number of cycles).
 * @param cycles The translated array of cycles.
 * @param number The number of cycles in the array.
 */
RDL_API
void RDL_deleteEdgeIdxArray(char **cycles, unsigned number);

/**
 * @brief Gives the edges of the graph.
 * @return the number of edges in the graph, \ref RDL_INVALID_RESULT on error
 * @param data pointer to the \ref RDL_data holding the calculation results.
 * @param RDL_edge_array_ptr pointer that points to the result array
 * (declare \ref RDL_edge * and give address as parameter)
 * @note Result has to be deallocated using `free(*RDL_edge_array_ptr)`
 */
RDL_API
unsigned RDL_getEdgeArray(const RDL_data *data, RDL_edge **RDL_edge_array_ptr);

/**
 * @brief Get the number of edges in the graph.
 * @param data pointer to the \ref RDL_data holding the calculation results.
 * @returns number of edges in the graph, \ref RDL_INVALID_RESULT on error
 */
RDL_API
unsigned RDL_getNofEdges(const RDL_data *data);

/**
 * @brief Get the id of the edge.
 * @param data pointer to the \ref RDL_data holding the calculation results.
 * @param from \ref RDL_node starting at this edge
 * @param to \ref RDL_node starting at this edge
 * @return edge id, \ref RDL_INVALID_RESULT, if not present
 * @note As the graph is stored as an adjacency list, this function
 * has linear runtime in the node degree.
 */
RDL_API
unsigned RDL_getEdgeId(const RDL_data *data, unsigned from, unsigned to);

/**
 * @brief Get the number of ring systems in the graph.
 * @param data pointer to the \ref RDL_data holding the calculation results.
 * @returns number of ring systems in the graph, \ref RDL_INVALID_RESULT on error
 */
RDL_API
unsigned RDL_getNofRingsystems(const RDL_data* data);

/**
 * @brief Get the number of nodes in the ring system.
 * @param data pointer to the \ref RDL_data holding the calculation results.
 * @param idx index of the ring system
 * @returns number of nodes in the ring system, \ref RDL_INVALID_RESULT on error
 * @note A ring system is a 2-connected component of the graph.
 */
RDL_API
unsigned RDL_getNofNodesForRingsystem(const RDL_data *data, unsigned idx);

/**
 * @brief Get the number of edges in the ring system.
 * @param data pointer to the \ref RDL_data holding the calculation results.
 * @param idx index of the ring system
 * @returns number of edges in the ring system, \ref RDL_INVALID_RESULT on error
 * @note A ring system is a 2-connected component of the graph.
 */
RDL_API
unsigned RDL_getNofEdgesForRingsystem(const RDL_data *data, unsigned idx);

/**
 * @brief Get the edges in the ring system.
 * @param data pointer to the \ref RDL_data holding the calculation results.
 * @param idx index of the ring system
 * @param edges pointer that points to the result array
 * (declare \ref RDL_edge * and give address as parameter)
 * @note Result has to be deallocated using `free(*edges)`
 * @returns number of edges in the ring system, \ref RDL_INVALID_RESULT on error
 * @note A ring system is a 2-connected component of the graph.
 */
RDL_API
unsigned RDL_getEdgesForRingsystem(
    const RDL_data *data, unsigned idx, RDL_edge** edges);

/**
 * @brief Get the nodes in the ring system.
 * @param data pointer to the \ref RDL_data holding the calculation results.
 * @param idx index of the ring system
 * @param nodes pointer that points to the result array
 * (declare \ref RDL_edge * and give address as parameter)
 * @note Result has to be deallocated using `free(*nodes)`
 * @returns number of nodes in the ring system, \ref RDL_INVALID_RESULT on error
 * @note A ring system is a 2-connected component of the graph.
 */
RDL_API
unsigned RDL_getNodesForRingsystem(
    const RDL_data *data, unsigned idx, RDL_node** nodes);

/**
 * @brief Get the ring system id for given edge.
 * @param data pointer to the \ref RDL_data holding the calculation results.
 * @param from \ref RDL_node starting at this edge
 * @param to \ref RDL_node starting at this edge
 * @returns number of edges in the ring system, \ref RDL_INVALID_RESULT on error
 *          \ref RDL_NO_RINGSYSTEM if edge is not part of a ring system
 * @note A ring system is a 2-connected component of the graph.
 */
RDL_API
unsigned RDL_getRingsystemForEdge(
    const RDL_data* data, unsigned from, unsigned to);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
