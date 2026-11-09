#include <catch2/catch_all.hpp>
#include <RDGeneral/RDLog.h>

#include <string>

#include "bench_common.hpp"

int main(int argc, char *argv[]) {
  RDLog::InitLogs();

  Catch::Session session;
  std::string smilesPath;
  using namespace Catch::Clara;
  session.cli(session.cli() |
              Opt(smilesPath, "path")["--smiles"](
                  "SMILES file (one molecule per line, first column) to use "
                  "instead of the built-in sample molecules"));

  if (const auto rc = session.applyCommandLine(argc, argv); rc != 0) {
    return rc;
  }
  if (!smilesPath.empty()) {
    bench_common::set_smiles_file(smilesPath);
  }
  return session.run();
}
