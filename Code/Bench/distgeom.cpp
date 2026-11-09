#include <catch2/catch_all.hpp>

#include <vector>

#include "bench_common.hpp"

#include <DistGeom/BoundsMatrix.h>
#include <DistGeom/TriangleSmooth.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/DistGeomHelpers/BoundsMatrixBuilder.h>

using namespace RDKit;

TEST_CASE("Triangle bounds smoothing", "[distgeom]") {
  // unsmoothed topological bounds, as ETKDG builds them before smoothing
  std::vector<std::vector<double>> bounds;
  std::vector<unsigned int> sizes;
  for (auto &mol : bench_common::load_samples()) {
    RWMol withHs(mol);
    MolOps::addHs(withHs);
    const auto n = withHs.getNumAtoms();
    DistGeom::BoundsMatPtr mat(new DistGeom::BoundsMatrix(n));
    DGeomHelpers::initBoundsMat(mat);
    DGeomHelpers::setTopolBounds(withHs, mat);
    bounds.emplace_back(mat->getData(), mat->getData() + n * n);
    sizes.push_back(n);
  }

  BENCHMARK_ADVANCED("triangleSmoothBounds")
  (Catch::Benchmark::Chronometer meter) {
    std::vector<std::vector<DistGeom::BoundsMatrix>> work(meter.runs());
    for (auto &run : work) {
      for (size_t i = 0; i < bounds.size(); ++i) {
        DistGeom::BoundsMatrix::DATA_SPTR data(new double[bounds[i].size()]);
        std::copy(bounds[i].begin(), bounds[i].end(), data.get());
        run.emplace_back(sizes[i], data);
      }
    }
    meter.measure([&](int run) {
      unsigned int nSmoothed = 0;
      for (auto &mat : work[run]) {
        nSmoothed += DistGeom::triangleSmoothBounds(&mat);
      }
      return nSmoothed;
    });
  };
}
