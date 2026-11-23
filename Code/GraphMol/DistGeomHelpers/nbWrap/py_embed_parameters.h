//
// Created by David Cosgrove on 28/10/2026.
//

#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/optional.h>

#include "RDBoost/Wrap_nb.h"
#include "DistGeom/BoundsMatrix.h"
#include "DistGeom/ChiralSet.h"
#include "GraphMol/ForceFieldHelpers/CrystalFF/TorsionPreferences.h"

namespace nb = nanobind;
using namespace nb::literals;  // Clearly this is bad form in a header, but
                               // neither I nor Gemini could come up with
                               // anything better.

namespace RDKit::DGeomHelpers {
struct PyEmbedParameters : public RDKit::DGeomHelpers::EmbedParameters {
 public:
  PyEmbedParameters() : RDKit::DGeomHelpers::EmbedParameters() {}
  PyEmbedParameters(const RDKit::DGeomHelpers::EmbedParameters &other)
      : RDKit::DGeomHelpers::EmbedParameters(other) {
    if (other.coordMap) {
      d_coordMap.reset(new std::map<int, RDGeom::Point3D>());
      for (const auto &item : *(other.coordMap)) {
        d_coordMap->insert({item.first, item.second});
      }
    }
  }

  void setCoordMap(const nb::dict &cmap) {
    d_coordMap.reset(new std::map<int, RDGeom::Point3D>());
    for (auto item : cmap) {
      (*d_coordMap)[nb::cast<int>(item.first)] =
          nb::cast<RDGeom::Point3D>(item.second);
    }
    coordMap = d_coordMap.get();
  }

  nb::tuple getFailureCounts() const {
    nb::list lst;
    for (auto failure : failures) {
      lst.append(failure);
    }
    return nb::tuple(lst);
  }

  void setCPCI(const nb::dict &CPCIdict) {
    CPCI = std::make_shared<
        std::map<std::pair<unsigned int, unsigned int>, double>>();
    for (auto item : CPCIdict) {
      auto id = nb::cast<nb::tuple>(item.first);
      unsigned int a = nb::cast<unsigned int>(id[0]);
      unsigned int b = nb::cast<unsigned int>(id[1]);
      (*CPCI)[{a, b}] = nb::cast<double>(item.second);
    }
  }

  void setBoundsMatrix(
      nb::ndarray<nb::numpy, double, nb::ndim<2>, nb::c_contig> bm) {
    if (bm.shape(0) != bm.shape(1)) {
      throw nb::value_error("The array has to be square");
    }
    if (bm.shape(0) == 0) {
      throw nb::value_error("The array has to have a nonzero size");
    }
    unsigned int nrows = (unsigned int)bm.shape(0);
    unsigned int dSize = nrows * nrows;
    auto *cData = new double[dSize];
    memcpy(cData, bm.data(), dSize * sizeof(double));
    DistGeom::BoundsMatrix::DATA_SPTR sdata(cData);
    boundsMat = boost::shared_ptr<const DistGeom::BoundsMatrix>(
        new DistGeom::BoundsMatrix(nrows, sdata));
  }

 private:
  std::unique_ptr<std::map<int, RDGeom::Point3D>> d_coordMap;
};
}  // namespace RDKit::DGeomHelpers
