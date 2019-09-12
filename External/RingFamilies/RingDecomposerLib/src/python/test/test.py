# This file is part of the RingDecomposerLib, licensed
# under BSD New license (see LICENSE in the root directory).
# Copyright (c) 2016
# University of Hamburg, ZBH - Center for Bioinformatics
# Niek Andresen, Florian Flachsenberg, Matthias Rarey
# 
# Please cite:
# 
# Kolodzik, A.; Urbaczek, S.; Rarey, M.
# Unique Ring Families: A Chemically Meaningful Description
# of Molecular Ring Topologies.
# J. Chem. Inf. Model., 2012, 52 (8), pp 2013-2021
# 
# Flachsenberg, F.; Andresen, N.; Rarey, M.
# RingDecomposerLib: An Open-Source Implementation of
# Unique Ring Families and Other Cycle Bases.
# J. Chem. Inf. Model., 2017, 57 (2), pp 122-126

import unittest
import operator
import py_rdl
import sys
import os
import re
import argparse


# simple DIMACS parser
def read_dimacs(fname):
    infile = open(fname, 'r')
    nof_nodes = nof_edges = None
    edges = []
    for line in infile:
        if line[0] == 'p':
            m = re.match("p edge (\d+) (\d+)", line)
            assert(m)
            nof_nodes = int(m.group(1))
            nof_edges = int(m.group(2))
        elif line[0] == 'e':
            assert nof_nodes and nof_edges
            m = re.match("e (\d+) (\d+)", line)
            node1 = int(m.group(1))
            node2 = int(m.group(2))
            assert node1 >= 1 and node1 <= nof_nodes \
               and node2 >= 1 and node2 <= nof_nodes
            assert len(edges) < nof_edges
            # use a more complicated representation of edges for testing
            edges.append("N{},N{}".format(node1, node2))
        elif line[0] == 'r':
            break
        elif line[0] == 'c':
            pass
        else:
            assert False
    assert len(edges) == nof_edges
    infile.close()
    return edges, nof_nodes

Testcases = {}

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='DIMACS input', required=True)
    parser.add_argument('unittest_args', nargs='*')
    args = parser.parse_args()

    testname = os.path.basename(args.input).replace('.', '_')
    graph, nof_nodes = read_dimacs(args.input)
    Testcases[testname] = graph, (nof_nodes)
    # we need this trick so we can pass command line arguments to tests
    sys.argv[1:] = args.unittest_args


# meta class trick so we can generate multiple tests for different data
class RDLTestMeta(type):
    def __new__(mcs, name, bases, dict):

        def calculate_urf(g, get_node1, get_node2):
            """
            Calculate URFdata
            """
            return py_rdl.Calculator.get_calculated_result(g, get_node1, get_node2)

        def make_graph(g):
            """
            Calculate with explicit graph
            """
            new_g = []
            for e in g:
                m = re.match("N(\d+),N(\d+)", e)
                new_e = (int(m.group(1)), int(m.group(2)))
                new_g.append(new_e)
            urf_g = py_rdl.Graph.from_edges(new_g)
            return urf_g

        def make_test_build(g, get_node1, get_node2):
            """
            Verify that calculation is sucessful.
            """
            def test(self):
                data1 = calculate_urf(g, get_node1, get_node2)
                self.assertTrue(data1.is_calculated())
                graph = make_graph(g)
                data2 = py_rdl.Calculator(graph)
                self.assertFalse(data2.is_calculated())
                self.assertRaises(py_rdl.RDLError, data2.get_nof_urf)
                data2.calculate()
                self.assertTrue(data2.is_calculated())
                self.assertEqual(data1.get_nof_urf(), data2.get_nof_urf())
                self.assertEqual(data1.get_nof_rcf(), data2.get_nof_rcf())
                with self.assertRaises(IndexError):
                    data1.get_weight_for_urf(1000000000)
                with self.assertRaises(IndexError):
                    data1.get_weight_for_rcf(1000000000)

            return test

        def make_test_urf(g, get_node1, get_node2, z):
            """
            Verify the intergrity of different URF accessors
            """
            def test(self):
                data = calculate_urf(g, get_node1, get_node2)
                self.assertTrue(data.is_calculated())
                if z:
                    self.assertTrue(data.urfs)
                urfs_for_nodes = {}
                urfs_for_edges = {}
                for urf in data:
                    edges = set()
                    nodes = set()
                    for cycle in data.get_relevant_cycles_for_urf(urf):
                        for e in cycle.edges:
                            edges.add(e)
                            if e not in urfs_for_edges:
                                urfs_for_edges[e] = set()
                            urfs_for_edges[e].add(urf.index)

                        for n in cycle.nodes:
                            nodes.add(n)
                            if n not in urfs_for_nodes:
                                urfs_for_nodes[n] = set()
                            urfs_for_nodes[n].add(urf.index)

                        self.assertEqual(cycle.urf.index, urf.index)
                        self.assertEqual(cycle.weight, urf.weight)
                    self.assertEqual(edges, set(urf.edges))
                    self.assertEqual(nodes, set(urf.nodes))
                if z:
                    self.assertTrue(urfs_for_nodes)
                    self.assertTrue(urfs_for_edges)
                for node in urfs_for_nodes:
                    urfs1 = urfs_for_nodes[node]
                    urfs2 = set([urf.index for urf in data.get_urfs_for_node(node)])
                    self.assertEqual(urfs1, urfs2)
                for edge in urfs_for_edges:
                    urfs1 = urfs_for_edges[edge]
                    urfs2 = set([urf.index for urf in data.get_urfs_for_edge(edge)])
                    self.assertEqual(urfs1, urfs2)
            return test

        def make_test_rcf(g, get_node1, get_node2, z):
            """
            Verify the integrity of the RCF results.
            """
            def test(self):
                data = calculate_urf(g, get_node1, get_node2)
                self.assertTrue(data.is_calculated())
                if z:
                    self.assertTrue(data.rcfs)
                rcfs_for_nodes = {}
                rcfs_for_edges = {}
                for rcf in data.rcfs:
                    edges = set()
                    nodes = set()
                    for cycle in data.get_relevant_cycles_for_rcf(rcf):
                        for e in cycle.edges:
                            edges.add(e)
                            if e not in rcfs_for_edges:
                                rcfs_for_edges[e] = set()
                            rcfs_for_edges[e].add(rcf.index)

                        for n in cycle.nodes:
                            nodes.add(n)
                            if n not in rcfs_for_nodes:
                                rcfs_for_nodes[n] = set()
                            rcfs_for_nodes[n].add(rcf.index)

                        self.assertEqual(cycle.rcf.index, rcf.index)
                        self.assertEqual(cycle.weight, rcf.weight)
                    self.assertEqual(edges, set(rcf.edges))
                    self.assertEqual(nodes, set(rcf.nodes))
                if z:
                    self.assertTrue(rcfs_for_nodes)
                    self.assertTrue(rcfs_for_edges)
                for node in rcfs_for_nodes:
                    rcfs1 = rcfs_for_nodes[node]
                    rcfs2 = set([rcf.index for rcf in data.get_rcfs_for_node(node)])
                    self.assertEqual(rcfs1, rcfs2)
                for edge in rcfs_for_edges:
                    rcfs1 = rcfs_for_edges[edge]
                    rcfs2 = set([rcf.index for rcf in data.get_rcfs_for_edge(edge)])
                    self.assertEqual(rcfs1, rcfs2)
            return test

        def make_test_rcp(g, get_node1, get_node2, z):
            """
            Verify integrity of RCPs.
            """
            def test(self):
                data = calculate_urf(g, get_node1, get_node2)
                self.assertTrue(data.is_calculated())
                rcps = data.get_relevant_cycle_prototypes()
                if z:
                    self.assertTrue(rcps)
                for i, rcp in enumerate(rcps):
                    self.assertEqual(i, rcp.rcf.index)
            return test

        def make_test_rc(g, get_node1, get_node2, z):
            """
            Verify integrity of RCs and numbers.
            """
            def test(self):
                data = calculate_urf(g, get_node1, get_node2)
                self.assertTrue(data.is_calculated())
                rcs = list(data.get_relevant_cycles())
                if z:
                    self.assertTrue(rcs)
                self.assertEqual(len(rcs), data.get_nof_relevant_cycles())
                sum = 0
                for urf in data:
                    urf_rc = list(data.get_relevant_cycles_for_urf(urf))
                    self.assertEqual(len(urf_rc), data.get_nof_relevant_cycles_for_urf(urf))
                    sum += data.get_nof_relevant_cycles_for_urf(urf)
                self.assertEqual(sum, data.get_nof_relevant_cycles())
            return test

        def make_test_rs(g, get_node1, get_node2, z):
            """
            Verify integrity of Ringsystems.
            """
            def test(self):
                data = calculate_urf(g, get_node1, get_node2)
                self.assertTrue(data.is_calculated())
                rss = data.get_nof_ringsystems()
                self.assertEqual(rss, len(data.ringsystems))
                if z:
                    self.assertTrue(data.ringsystems)
                for rs in data.ringsystems:
                    self.assertTrue(rs.nodes)
                    self.assertTrue(rs.edges)
            return test

        def make_test_sssr(g, get_node1, get_node2, z):
            """
            Verify integrity of SSSR.
            """
            def test(self):
                data = calculate_urf(g, get_node1, get_node2)
                self.assertTrue(data.is_calculated())
                sssr = data.get_sssr()
                self.assertEqual(z, len(sssr))
                self.assertEqual(z, len(set([c.rcf for c in sssr])))
            return test

        for tname in Testcases:
            g, nof_nodes = Testcases[tname]
            print(tname)

            def get_node1(e):
                return e.split(',')[0]
            def get_node2(e):
                return e.split(',')[1]

            nof_edges = len(g)

            nodes = set()
            for e in g:
                nodes.add(get_node1(e))
                nodes.add(get_node2(e))

            nof_nodes_real = len(nodes)

            def dfs_visit(n, visited):
                visited.add(n)
                for e in g:
                    next = None
                    if get_node1(e) == n:
                        next = get_node2(e)
                    elif get_node2(e) == n:
                        next = get_node1(e)
                    if next is not None and next not in visited:
                        dfs_visit(next, visited)

            nof_cc = 0
            visited = set()
            for n in nodes:
                if n not in visited:
                    nof_cc += 1
                dfs_visit(n, visited)

            z = nof_edges - nof_nodes_real + nof_cc

            test_name = "test_build_{}".format(tname)
            dict[test_name] = make_test_build(g, get_node1, get_node2)
            test_name = "test_urf_{}".format(tname)
            dict[test_name] = make_test_urf(g, get_node1, get_node2, z)
            test_name = "test_rcf_{}".format(tname)
            dict[test_name] = make_test_rcf(g, get_node1, get_node2, z)
            test_name = "test_rcp_{}".format(tname)
            dict[test_name] = make_test_rcp(g, get_node1, get_node2, z)
            test_name = "test_rc_{}".format(tname)
            dict[test_name] = make_test_rc(g, get_node1, get_node2, z)
            test_name = "test_sssr_{}".format(tname)
            dict[test_name] = make_test_sssr(g, get_node1, get_node2, z)
            test_name = "test_rs_{}".format(tname)
            dict[test_name] = make_test_rs(g, get_node1, get_node2, z)
        return type.__new__(mcs, name, bases, dict)


# strange syntax for python2/python3 compatibility
RDLTest = RDLTestMeta(str('RDLTest'), (unittest.TestCase,), {})

if __name__ == '__main__':
    unittest.main()
