// benches about benches

#include <catch2/catch_all.hpp>
#include <string>

#include "bench_common.hpp"

TEST_CASE("bench_common::nth_random", "[meta]") {
  BENCHMARK("bench_common::nth_random", i) {
    return bench_common::nth_random(i);
  };
}
