//
// file RDKitSVSettings.H
// David Cosgrove
// AstraZeneca
// 4th July 2014
//
// Command-line parsing for RDKitSV

#ifndef RDKITSVSETTINGS_H
#define RDKITSVSETTINGS_H

#include <string>
#include <vector>

#include <boost/program_options/options_description.hpp>

// ****************************************************************************

class RDKitSVSettings {

public :

  RDKitSVSettings( int argc , char **argv );

  std::string usage_text() const { return usage_text_; }
  std::vector<std::string> mol_files() const { return mol_files_; }
  std::string smarts_file() const { return smarts_file_; }

private :

  std::string usage_text_;
  std::vector<std::string> mol_files_;
  std::string smarts_file_;

  void build_program_options( boost::program_options::options_description &desc );

};

#endif // RDKITSVSETTINGS_H
