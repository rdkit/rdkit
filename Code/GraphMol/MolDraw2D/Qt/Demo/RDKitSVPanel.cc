//
// file RDKitSVPanel.cc
// David Cosgrove
// AstraZeneca
// 20th June 2014
//

#include "RDKitSVPanel.H"
#include "MolDisplay2DWidget.H"

#include <QLabel>
#include <QLayout>
#include <QSlider>

using namespace std;

namespace RDKitSV {

// ****************************************************************************
RDKitSVPanel::RDKitSVPanel(bool left_slider, QWidget *parent, Qt::WindowFlags f)
    : QWidget(parent, f) {
  build_widget(left_slider);
}

// ****************************************************************************
void RDKitSVPanel::set_molecules(const vector<RDKit::ROMOL_SPTR> &new_mols,
                                 const vector<vector<int>> &highlight_atoms) {
#ifdef NOTYET
  cout << "RDKitSVPanel::set_molecules : " << new_mols.size() << endl;
#endif

  mols_ = new_mols;
  highlight_atoms_ = highlight_atoms;
  if (highlight_atoms_.size() != mols_.size()) {
    highlight_atoms_.clear();
  }

  if (mols_.empty()) {
    mol_slider_->setDisabled(true);
  } else {
    mol_slider_->setEnabled(true);
    mol_slider_->setRange(0, mols_.size() - 1);
    mol_slider_->setValue(0);
  }

  slot_slider_changed();  // to force picture update
}

// ****************************************************************************
void RDKitSVPanel::set_label(const QString &new_label) {
  label_->setText(new_label);
  label_->setWordWrap(true);
  label_->setAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
  if (new_label.isEmpty()) {
    label_->hide();
  } else {
    label_->show();
  }
}

// ****************************************************************************
void RDKitSVPanel::build_widget(bool left_slider) {
  QHBoxLayout *hbox = new QHBoxLayout;
  mol_draw_ = new RDKit::MolDisplay2DWidget;

  mol_slider_ = new QSlider;
  connect(mol_slider_, &QSlider::valueChanged, this,
          &RDKitSVPanel::slot_slider_changed);
  mol_slider_->setPageStep(1);

  if (left_slider) {
    hbox->addWidget(mol_slider_);
    hbox->addWidget(mol_draw_, 1);
  } else {
    hbox->addWidget(mol_draw_, 1);
    hbox->addWidget(mol_slider_);
  }

  label_ = new QLabel;
  QVBoxLayout *vbox = new QVBoxLayout;
  vbox->addLayout(hbox, 1);
  vbox->addWidget(label_);
  setLayout(vbox);

  label_->hide();
}

// ****************************************************************************
void RDKitSVPanel::slot_slider_changed() {
  if (mols_.empty()) {
    mol_draw_->set_display_mol(RDKit::ROMOL_SPTR());
  } else {
    int mol_num = mol_slider_->value();
    mol_draw_->set_display_mol(mols_[mol_num]);
    if (!highlight_atoms_.empty()) {
      mol_draw_->set_selected_atoms(highlight_atoms_[mol_num]);
    }
  }
}

}  // namespace RDKitSV
