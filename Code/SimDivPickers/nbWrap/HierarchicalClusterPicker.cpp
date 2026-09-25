//
//  Copyright (C) 2003-2026 Greg Landrum and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/ndarray.h>

#include <SimDivPickers/DistPicker.h>
#include <SimDivPickers/HierarchicalClusterPicker.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDPickers {

// REVIEW: the poolSize can be pulled from the numeric array
RDKit::INT_VECT HierarchicalPicks(HierarchicalClusterPicker *picker,
                                  nb::ndarray<nb::numpy, const double,
                                              nb::ndim<1>, nb::c_contig> distMat,
                                  int poolSize, int pickSize) {
  if (pickSize >= poolSize) {
    throw nb::value_error("pickSize must be less than poolSize");
  }

  std::vector<double> dMatCopy(distMat.data(), distMat.data() + distMat.shape(0));
  return picker->pick(dMatCopy.data(), poolSize, pickSize);
}

// REVIEW: the poolSize can be pulled from the numeric array
RDKit::VECT_INT_VECT HierarchicalClusters(HierarchicalClusterPicker *picker,
                                          nb::ndarray<nb::numpy, const double,
                                                      nb::ndim<1>, nb::c_contig> distMat,
                                          int poolSize, int pickSize) {
  std::vector<double> dMatCopy(distMat.data(), distMat.data() + distMat.shape(0));
  return picker->cluster(dMatCopy.data(), poolSize, pickSize);
}

}  // namespace RDPickers

void wrap_HierarchCP(nb::module_ &m) {
  nb::class_<RDPickers::HierarchicalClusterPicker>(
      m, "HierarchicalClusterPicker",
      "A class for diversity picking of items using Hierarchical Clustering\n")
      .def(nb::init<RDPickers::HierarchicalClusterPicker::ClusterMethod>(),
           "clusterMethod"_a)
      .def("Pick", RDPickers::HierarchicalPicks,
           "distMat"_a, "poolSize"_a, "pickSize"_a,
           R"DOC(Pick a diverse subset of items from a pool of items using hierarchical clustering

ARGUMENTS:
  - distMat: 1D distance matrix (only the lower triangle elements)
  - poolSize: number of items in the pool
  - pickSize: number of items to pick from the pool
)DOC")
      .def("Cluster", RDPickers::HierarchicalClusters,
           "distMat"_a, "poolSize"_a, "pickSize"_a,
           R"DOC(Return a list of clusters of item from the pool using hierarchical clustering

ARGUMENTS:
  - distMat: 1D distance matrix (only the lower triangle elements)
  - poolSize: number of items in the pool
  - pickSize: number of items to pick from the pool
)DOC");

  nb::enum_<RDPickers::HierarchicalClusterPicker::ClusterMethod>(
      m, "ClusterMethod")
      .value("WARD", RDPickers::HierarchicalClusterPicker::WARD)
      .value("SLINK", RDPickers::HierarchicalClusterPicker::SLINK)
      .value("CLINK", RDPickers::HierarchicalClusterPicker::CLINK)
      .value("UPGMA", RDPickers::HierarchicalClusterPicker::UPGMA)
      .value("MCQUITTY", RDPickers::HierarchicalClusterPicker::MCQUITTY)
      .value("GOWER", RDPickers::HierarchicalClusterPicker::GOWER)
      .value("CENTROID", RDPickers::HierarchicalClusterPicker::CENTROID)
      .export_values();
}
