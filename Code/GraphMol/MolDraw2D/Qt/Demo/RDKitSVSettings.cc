//
// file RDKitSVSettings.cc
// David Cosgrove
// AstraZeneca
// 4th July 2014
//

#include "RDKitSVSettings.H"

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <iostream>

using namespace std;
namespace po = boost::program_options;

// *****************************************************************************
RDKitSVSettings::RDKitSVSettings(int argc, char **argv) {
  po::options_description desc("Allowed Options");
  build_program_options(desc);

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    cout << desc << endl;
    exit(1);
  }

  ostringstream oss;
  oss << desc;
  usage_text_ = oss.str();
}

// *****************************************************************************
void RDKitSVSettings::build_program_options(po::options_description &desc) {
  desc.add_options()("help", "Produce this help text.")(
      "molecule-file,M", po::value<vector<string>>(&mol_files_),
      "Name of file containing molecules.")(
      "smarts-file,S", po::value<string>(&smarts_file_),
      "Name of file containing SMARTS strings.");
}
