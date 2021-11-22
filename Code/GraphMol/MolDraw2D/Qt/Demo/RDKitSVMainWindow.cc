//
// file RDKitSV.cc
// David Cosgrove
// AstraZeneca
//
// 20th June 2014
//

#include "RDKitSVMainWindow.H"
#include "RDKitSVPanel.H"
#include "RDKitSVSettings.H"
#include "QTGet2Strings.H"
#include "QT4SelectItems.H"

#include <QAction>
#include <QApplication>
#include <QFileDialog>
#include <QFileInfo>
#include <QHBoxLayout>
#include <QMenu>
#include <QMenuBar>
#include <QMessageBox>
#include <QStatusBar>
#include <QString>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/lexical_cast.hpp>

#include <fstream>
#include <iostream>
#include <list>

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

using namespace boost;
using namespace std;

namespace RDKitSV {

// ****************************************************************************
RDKitSVMainWindow::RDKitSVMainWindow(int argc, char **argv) : QMainWindow() {
  build_actions();
  build_menubar();
  build_widget();

  parse_args(argc, argv);
}

// ****************************************************************************
void RDKitSVMainWindow::build_actions() {
  build_file_actions();
  build_smarts_actions();
}

// ****************************************************************************
void RDKitSVMainWindow::build_file_actions() {
  file_exit_ = new QAction(tr("Exit"), this);
  file_exit_->setShortcut(QString("Ctrl+Q"));
  connect(file_exit_, &QAction::triggered, this, &RDKitSVMainWindow::slot_exit);

  file_read_mols_ = new QAction(tr("Molecules"), this);
  connect(file_read_mols_, &QAction::triggered, this,
          &RDKitSVMainWindow::slot_read_mols);

  file_read_smarts_ = new QAction(tr("SMARTS"), this);
  connect(file_read_smarts_, &QAction::triggered, this,
          &RDKitSVMainWindow::slot_read_smarts);

  file_write_left_ = new QAction(tr("Molecules (Left Panel)"), this);
  connect(file_write_left_, &QAction::triggered, this,
          &RDKitSVMainWindow::slot_write_left_molecules);

  file_write_right_ = new QAction(tr("Molecules (Right Panel)"), this);
  connect(file_write_right_, &QAction::triggered, this,
          &RDKitSVMainWindow::slot_write_right_molecules);

  file_write_smarts_ = new QAction(tr("SMARTS"), this);
  connect(file_write_smarts_, &QAction::triggered, this,
          &RDKitSVMainWindow::slot_write_smarts);
}

// ****************************************************************************
void RDKitSVMainWindow::build_smarts_actions() {
  smarts_match_ = new QAction(tr("Match"), this);
  smarts_match_->setShortcut(QString("Ctrl+M"));
  connect(smarts_match_, &QAction::triggered, this,
          &RDKitSVMainWindow::slot_match_smarts);

  smarts_edit_ = new QAction(tr("Edit"), this);
  connect(smarts_edit_, &QAction::triggered, this,
          &RDKitSVMainWindow::slot_edit_smarts);

  smarts_new_ = new QAction(tr("New"), this);
  connect(smarts_new_, &QAction::triggered, this,
          &RDKitSVMainWindow::slot_new_smarts);
}

// ****************************************************************************
void RDKitSVMainWindow::build_menubar() {
  QMenu *file_menu = menuBar()->addMenu(tr("File"));
  QMenu *read_menu = file_menu->addMenu(tr("Read"));
  read_menu->addAction(file_read_mols_);
  read_menu->addAction(file_read_smarts_);

  QMenu *write_menu = file_menu->addMenu(tr("Write"));
  write_menu->addAction(file_write_left_);
  write_menu->addAction(file_write_right_);
  write_menu->addAction(file_write_smarts_);

  file_menu->addSeparator();
  file_menu->addAction(file_exit_);

  QMenu *smarts_menu = menuBar()->addMenu(tr("SMARTS"));
  smarts_menu->addAction(smarts_match_);
  smarts_menu->addAction(smarts_edit_);
  smarts_menu->addAction(smarts_new_);
}

// ****************************************************************************
void RDKitSVMainWindow::build_widget() {
  left_panel_ = new RDKitSVPanel;
  right_panel_ = new RDKitSVPanel(false);

  QWidget *cen_wid = new QWidget;
  QHBoxLayout *hbox = new QHBoxLayout;
  hbox->addWidget(left_panel_);
  hbox->addWidget(right_panel_);
  right_panel_->hide();

  cen_wid->setLayout(hbox);
  setCentralWidget(cen_wid);
}

// ****************************************************************************
void RDKitSVMainWindow::parse_args(int argc, char **argv) {
  RDKitSVSettings s(argc, argv);

  vector<string> mol_files = s.mol_files();
  if (!mol_files.empty()) {
    BOOST_FOREACH (string mf, mol_files) { read_mols(mf); }
  }

  string smarts_file = s.smarts_file();
  if (!smarts_file.empty()) {
    read_smarts(smarts_file);
  }
}

// ****************************************************************************
RDKitSVMainWindow::FILE_TYPE RDKitSVMainWindow::get_filetype(
    const string &filename, bool &is_compressed) {
  is_compressed = false;
  if (filename.substr(filename.length() - 3) == string(".gz")) {
    is_compressed = true;
  }

  if (string(".smi") == filename.substr(filename.length() - 4) ||
      string(".smi.gz") == filename.substr(filename.length() - 7) ||
      string(".ism") == filename.substr(filename.length() - 4) ||
      string(".ism.gz") == filename.substr(filename.length() - 7)) {
    return SMILES;
  } else if (string(".sdf") == filename.substr(filename.length() - 4) ||
             string(".sdf.gz") == filename.substr(filename.length() - 7)) {
    return SDF;
  }

  return UNKNOWN;
}

// ****************************************************************************
void RDKitSVMainWindow::read_mols(const string &filename) {
  QApplication::setOverrideCursor(Qt::WaitCursor);
  bool is_compressed;
  FILE_TYPE file_type = get_filetype(filename, is_compressed);

  boost::iostreams::filtering_istream ins;
  if (is_compressed) {
    ins.push(boost::iostreams::gzip_decompressor());
  }
  ins.push(boost::iostreams::file_source(filename.c_str()));
  if (!ins || !ins.good()) {
    QMessageBox::warning(this, tr("File Read Error"),
                         QString("%1 %2 %3 for reading.")
                             .arg(tr("Couldn't open"))
                             .arg(filename.c_str())
                             .arg(tr("for reading.")));
    return;
  }

  vector<RDKit::ROMOL_SPTR> new_mols;
  if (SDF == file_type) {
    RDKit::MolSupplier *mol_supplier =
        new RDKit::ForwardSDMolSupplier(&ins, true);
    while (!mol_supplier->atEnd()) {
      new_mols.push_back(RDKit::ROMOL_SPTR(mol_supplier->next()));
      if (!new_mols.back()) {
        new_mols.pop_back();
        break;
      }
    }
  } else if (SMILES == file_type) {
    // SmilesMolSupplier doesn't seem to work with compressed streams. It might
    // be because of the use of getline.
    // For now, assume the file is a simple one with no header and just SMILES
    // followed by optional name
    istreambuf_iterator<char> isb(ins), eos;

    vector<char> nextline;

    while (isb != eos) {
      nextline.clear();
      while (isb != eos && *isb != '\n') {
        nextline.push_back(*isb);
        ++isb;
      }
      if ('#' != nextline[0]) {
        vector<string> split_line;
        split(split_line, nextline, is_any_of(" \t"));
        string smi(split_line[0]);
        string smi_name;
        if (1 == split_line.size()) {
          smi_name = string("MOL_") + lexical_cast<string>(mols_.size() + 1);
        } else {
          smi_name = string(nextline.begin() + split_line[0].length() + 1,
                            nextline.end());
        }
        boost::trim(smi);
        boost::trim(smi_name);
        RDKit::ROMol *mol = RDKit::SmilesToMol(smi);
        if (mol) {
          mol->setProp("_Name", smi_name);
          mols_.push_back(RDKit::ROMOL_SPTR(mol));
        }
      }
      ++isb;  // get past '\n'
    }
  } else {
    QMessageBox::warning(this, tr("File Read Error"),
                         QString("%1 %2.")
                             .arg(tr("Unrecognised file type for"))
                             .arg(filename.c_str()));
    return;
  }

  mols_.insert(mols_.end(), new_mols.begin(), new_mols.end());
  left_panel_->set_molecules(mols_);

  QApplication::restoreOverrideCursor();
}

// ****************************************************************************
void RDKitSVMainWindow::read_smarts(const string &filename) {
#ifdef NOTYET
  cout << "reading SMARTS from " << filename << endl;
#endif
  ifstream ifs(filename.c_str());
  if (!ifs || !ifs.good()) {
    QMessageBox::warning(this, tr("File Read Error"),
                         QString("%1 %2 %3 for reading.")
                             .arg(tr("Couldn't open"))
                             .arg(filename.c_str())
                             .arg(tr("for reading.")));
    return;
  }

  while (1) {
    string next_line;
    getline(ifs, next_line);
#ifdef NOTYET
    cout << next_line << "XX" << endl;
#endif
    if (!ifs.good() || ifs.eof()) {
      break;
    }
    boost::trim(next_line);
    if (next_line.empty()) {
      continue;
    }
    if ('#' == next_line[0]) {
      continue;  // it's a comment
    }
    list<string> splits;
    boost::split(splits, next_line, boost::is_any_of(" \t"));
    string smt, smt_name;
    if (splits.size() >= 2) {
      smt = splits.front();
      splits.pop_front();
      boost::trim(smt);
      smt_name = boost::join(splits, " ");
      boost::trim(smt_name);
#ifdef NOTYET
      cout << "smt = " << smt << " name = " << smt_name << endl;
#endif
    } else if (1 == splits.size()) {
      smt = splits.front();
      boost::trim(smt);
      smt_name = string("SMARTS_") + lexical_cast<string>(smarts_.size() + 1);
    }
    smarts_.push_back(make_pair(smt_name, smt));
  }
  statusBar()->showMessage(
      QString("%1 %2.").arg(tr("Number of SMARTS now")).arg(smarts_.size()),
      2000);
}

// ****************************************************************************
void RDKitSVMainWindow::write_mols(RDKitSVPanel &panel,
                                   const string &filename) {
  vector<RDKit::ROMOL_SPTR> mols = panel.get_molecules();
  if (mols.empty()) {
    QMessageBox::information(
        this, tr("No molecules."),
        tr("Panel contains no molecules, so nothing to write."));
    return;
  }

  bool is_compressed;
  FILE_TYPE file_type = get_filetype(filename, is_compressed);

  boost::iostreams::filtering_ostream ons;
  if (is_compressed) {
    ons.push(boost::iostreams::gzip_compressor());
  }
  ons.push(boost::iostreams::file_sink(filename.c_str()));
  if (!ons || !ons.good()) {
    QMessageBox::warning(
        this, tr("Can't write to file"),
        QString("%1 %2.").arg(
            tr("Couldn't write to file").arg(filename.c_str())));
    return;
  }

  RDKit::MolWriter *mw;
  if (SMILES == file_type) {
    bool ism = filename.substr(filename.length() - 4) == ".ism" ||
                       filename.substr(filename.length() - 7) == ".ism.gz"
                   ? true
                   : false;
    mw = new RDKit::SmilesWriter(&ons, " ", "Name", true, false, ism);
  } else if (SDF == file_type) {
    mw = new RDKit::SDWriter(&ons, false);
  } else {
    QMessageBox::warning(this, tr("File Write Error"),
                         QString("%1 %2.")
                             .arg(tr("Unrecognised file type for"))
                             .arg(filename.c_str()));
    return;
  }

  BOOST_FOREACH (RDKit::ROMOL_SPTR mol, mols) { mw->write(*mol); }

  delete mw;
}

// ****************************************************************************
void RDKitSVMainWindow::match_smarts(const vector<pair<string, string>> &smts) {
  cout << "Matching " << smts.size() << " SMARTS" << endl;
  vector<vector<int>> smarts_match(mols_.size(), vector<int>());

  QString smts_label;
  for (int i = 0, is = smts.size(); i < is; ++i) {
    cout << "Matching " << i << " : " << smts[i].second << endl;
    smts_label += QString("%1 (%2)")
                      .arg(smts[i].first.c_str())
                      .arg(smts[i].second.c_str());
    if (i < is - 1) {
      smts_label += "|";
    }
    RDKit::RWMol *q = RDKit::SmartsToMol(smts[i].second);
    for (int j = 0, js = mols_.size(); j < js; ++j) {
      vector<RDKit::MatchVectType> hits_vect;
      if (RDKit::SubstructMatch(*mols_[j], *q, hits_vect)) {
        // store the numbers of the atoms that matched
        typedef pair<int, int> INTPAIR;
        BOOST_FOREACH (RDKit::MatchVectType hit, hits_vect) {
          BOOST_FOREACH (INTPAIR ip, hit) {
            smarts_match[j].push_back(ip.second);
          }
        }
        std::sort(smarts_match[j].begin(), smarts_match[j].end());
        smarts_match[j].erase(
            unique(smarts_match[j].begin(), smarts_match[j].end()),
            smarts_match[j].end());
      }
    }
    delete q;
  }

  vector<RDKit::ROMOL_SPTR> match_mols, miss_mols;
  vector<vector<int>> match_atom_matches;
  for (int i = 0, is = mols_.size(); i < is; ++i) {
    if (smarts_match[i].empty()) {
      miss_mols.push_back(mols_[i]);
    } else {
      match_mols.push_back(mols_[i]);
      match_atom_matches.push_back(smarts_match[i]);
    }
  }

  right_panel_->show();
  left_panel_->set_molecules(match_mols, match_atom_matches);
  // Details, details!!
  if (1 == match_mols.size()) {
    left_panel_->set_label(
        QString("%1 %2.").arg(tr("1 molecule matched.").arg(smts_label)));
  } else {
    left_panel_->set_label(QString("%1 %2 %3.")
                               .arg(match_mols.size())
                               .arg(tr("molecules matched"))
                               .arg(smts_label));
  }
  right_panel_->set_molecules(miss_mols);
  if (1 == miss_mols.size()) {
    right_panel_->set_label(tr("1 molecule didn't match"));
  } else {
    right_panel_->set_label(QString("%1 %2")
                                .arg(miss_mols.size())
                                .arg(tr("molecules didn't match.")));
  }
}

// ****************************************************************************
vector<pair<string, string>> RDKitSVMainWindow::select_smarts(
    bool multi_select) {
  vector<pair<string, string>> smts_to_use;

  if (smarts_.empty()) {
    return smts_to_use;
  }

  vector<char> selected_smarts(smarts_.size(), 0);
  vector<QString> smts;
  for (int i = 0, is = smarts_.size(); i < is; ++i) {
#ifdef NOTYET
    cout << i << " : " << smarts_[i].first << " :: " << smarts_[i].second
         << endl;
#endif
    smts.push_back(QString("%1 (%2)")
                       .arg(smarts_[i].first.c_str())
                       .arg(smarts_[i].second.c_str()));
  }

  DACLIB::QT4SelectItems si(string("Select SMARTS to use."), smts,
                            selected_smarts, !multi_select, this);
  if (QDialog::Accepted != si.exec()) {
    return smts_to_use;
  }
  si.get_results(selected_smarts);

  for (int i = 0, is = smarts_.size(); i < is; ++i) {
    if (selected_smarts[i]) {
      smts_to_use.push_back(smarts_[i]);
    }
  }

  return smts_to_use;
}

// ****************************************************************************
void RDKitSVMainWindow::update_smarts(const string &new_name,
                                      const string &new_val) {
  vector<pair<string, string>>::iterator p =
      std::find_if(smarts_.begin(), smarts_.end(),
                   bind(equal_to<string>(), new_name,
                        bind(&pair<string, string>::first, _1)));
  if (p == smarts_.end()) {
    smarts_.push_back(make_pair(new_name, new_val));
  } else {
    p->second = new_val;
  }
}

// ****************************************************************************
void RDKitSVMainWindow::slot_read_mols() {
  QString filename = QFileDialog::getOpenFileName(
      this, "Select molecule file", last_dir_,
      "Mol. files (*.smi *.smi.gz *.ism *.ism.gz *.sdf *.sdf.gz)");
  if (filename.isEmpty() || filename.isNull()) {
    return;
  }

  QFileInfo fi(filename);
  last_dir_ = fi.absolutePath();

  read_mols(filename.toLocal8Bit().data());
}

// ****************************************************************************
void RDKitSVMainWindow::slot_read_smarts() {
  QString filename =
      QFileDialog::getOpenFileName(this, "Select SMARTS file", last_dir_,
                                   "SMARTS files (*.smt);;Any file (*)");
  if (filename.isEmpty() || filename.isNull()) {
    return;
  }

  QFileInfo fi(filename);
  last_dir_ = fi.absolutePath();

  read_smarts(filename.toLocal8Bit().data());
}

// ****************************************************************************
void RDKitSVMainWindow::slot_exit() { exit(0); }

// ****************************************************************************
void RDKitSVMainWindow::slot_match_smarts() {
  vector<pair<string, string>> smts_to_use = select_smarts(true);
  if (smts_to_use.empty()) {
    return;
  }

  match_smarts(smts_to_use);
}

// ****************************************************************************
void RDKitSVMainWindow::slot_edit_smarts() {
  if (smarts_.empty()) {
    slot_new_smarts();
  }

  vector<pair<string, string>> smts_to_use = select_smarts(false);
  if (smts_to_use.empty()) {
    return;
  }

  DACLIB::QTGet2Strings g2s(
      tr("SMARTS name"), smts_to_use.front().first.c_str(), tr("SMARTS value"),
      smts_to_use.front().second.c_str(), this);
  if (QDialog::Accepted != g2s.exec()) {
    return;
  }

  QString new_name, new_smarts;
  g2s.get_values(new_name, new_smarts);

  update_smarts(new_name.toLocal8Bit().data(), new_smarts.toLocal8Bit().data());
}

// ****************************************************************************
void RDKitSVMainWindow::slot_new_smarts() {
  DACLIB::QTGet2Strings g2s(tr("SMARTS name"),
                            QString("SMARTS_%1").arg(smarts_.size() + 1),
                            tr("SMARTS value"), QString(""), this);
  if (QDialog::Accepted != g2s.exec()) {
    return;
  }

  QString new_name, new_smarts;
  g2s.get_values(new_name, new_smarts);

  update_smarts(new_name.toLocal8Bit().data(), new_smarts.toLocal8Bit().data());
}

// ****************************************************************************
void RDKitSVMainWindow::slot_write_left_molecules() {
  QString filename = QFileDialog::getSaveFileName(
      this, "File for left panel molecules", last_dir_,
      "Mol. files (*.smi *.smi.gz *.ism *.ism.gz *.sdf *.sdf.gz)");
  if (filename.isEmpty()) {
    return;
  }

  QFileInfo fi(filename);
  last_dir_ = fi.absolutePath();

  write_mols(*left_panel_, filename.toLocal8Bit().data());
}

// ****************************************************************************
void RDKitSVMainWindow::slot_write_right_molecules() {
  if (right_panel_->isHidden()) {
    // it's empty, so nothing to do.
    QMessageBox::information(
        this, tr("No molecules."),
        tr("Right panel is inactive, so nothing to write."));
    return;
  }

  QString filename = QFileDialog::getSaveFileName(
      this, "File for left panel molecules", last_dir_,
      "Mol. files (*.smi *.smi.gz *.ism *.ism.gz *.sdf *.sdf.gz)");
  if (filename.isEmpty()) {
    return;
  }

  QFileInfo fi(filename);
  last_dir_ = fi.absolutePath();

  write_mols(*right_panel_, filename.toLocal8Bit().data());
}

// ****************************************************************************
void RDKitSVMainWindow::slot_write_smarts() {
  QString filename =
      QFileDialog::getSaveFileName(this, "File for SMARTS.", last_dir_,
                                   "SMARTS files (*.smt);;Any file (*)");
  if (filename.isEmpty()) {
    return;
  }

  QFileInfo fi(filename);
  last_dir_ = fi.absolutePath();

  ofstream ofs(filename.toLocal8Bit().data());
  if (!ofs || !ofs.good()) {
    QMessageBox::warning(
        this, tr("Can't write to file"),
        QString("%1 %2.").arg(tr("Couldn't write to file").arg(filename)));
    return;
  }

  vector<pair<string, string>> smts_to_use = select_smarts(true);
  if (smts_to_use.empty()) {
    return;
  }

  typedef pair<string, string> STRINGPAIR;
  BOOST_FOREACH (STRINGPAIR smt, smts_to_use) {
    ofs << smt.first << " " << smt.second << endl;
  }
}

}  // namespace RDKitSV
