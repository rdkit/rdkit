//
// file MolDisplay2DWidget.H
// David Cosgrove
// AstraZeneca
// 6th June 2014
//

#include "MolDisplay2DWidget.H"

#include "stddefs.H"
#include <GraphMol/MolDraw2D/MolDraw2DQt.h>

#include <QColor>
#include <QMouseEvent>
#include <QPainter>
#include <QPaintEvent>
#include <QString>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>

#include <limits>
#include <map>

#include <boost/tuple/tuple_io.hpp>

using namespace boost;
using namespace std;

namespace RDKit {

// ****************************************************************************
MolDisplay2DWidget::MolDisplay2DWidget(QWidget *parent, Qt::WindowFlags f)
    : QWidget(parent, f), pick_circle_rad_(-1), atom_picking_(false) {}

// ****************************************************************************
void MolDisplay2DWidget::set_display_mol(ROMOL_SPTR new_mol) {
  if (!new_mol) {
    disp_mol_.reset();
  } else {
    disp_mol_ = RWMOL_SPTR(new RWMol(*new_mol));
    RDKit::MolOps::Kekulize(*disp_mol_);
    RDDepict::compute2DCoords(*disp_mol_);
    if (!disp_mol_->hasProp("_drawingBondsWedged")) {
      RDKit::Conformer conf = disp_mol_->getConformer();
      RDKit::WedgeMolBonds(*disp_mol_, &conf);
    }
  }
  update();
}

// ****************************************************************************
QSize MolDisplay2DWidget::minimumSize() const { return QSize(100, 100); }

// ****************************************************************************
void MolDisplay2DWidget::set_atom_picking(bool new_val) {
  atom_picking_ = new_val;

  if (atom_picking_) {
    picked_atoms_.clear();
  }
}

// ****************************************************************************
void MolDisplay2DWidget::set_selected_atoms(const vector<int> &sa) {
  picked_atoms_ = sa;
  update();
}

// ****************************************************************************
int MolDisplay2DWidget::pick_circle_rad() const {
  if (-1 == pick_circle_rad_) {
    pick_circle_rad_ = std::min(width() / 100, height() / 100);
  }

  return pick_circle_rad_;
}

// ****************************************************************************
QSize MolDisplay2DWidget::sizeHint() const { return QSize(400, 400); }

// ****************************************************************************
void MolDisplay2DWidget::paintEvent(QPaintEvent *event) {
  pick_circle_rad_ = std::min(width() / 100, height() / 100);

  QPainter qp;
  qp.begin(this);

  QFont font = qp.font();
  font.setPointSize(6);
  qp.setFont(font);

  qp.setRenderHint(QPainter::Antialiasing, true);
  qp.setRenderHint(QPainter::TextAntialiasing, true);

  if (!disp_mol_) {
    return;
  }

  draw_molecule(qp);
}

// ****************************************************************************
void MolDisplay2DWidget::mousePressEvent(QMouseEvent *event) {
  if (atom_picking_) {
    select_atom(event);
  }

  update();
}

// ****************************************************************************
void MolDisplay2DWidget::select_atom(QMouseEvent *event) {
  int na = find_nearest_atom(event->x(), event->y());
  if (-1 != na) {
    vector<int>::iterator p =
        std::find(picked_atoms_.begin(), picked_atoms_.end(), na);
    if (p == picked_atoms_.end()) {
      picked_atoms_.push_back(na);
    } else {
      picked_atoms_.erase(p);
    }
  } else {
    if (!(event->modifiers() & Qt::ControlModifier)) {
      // user missed, doesn't have control button down, so clear
      picked_atoms_.clear();
    }
  }
}

// ****************************************************************************
void MolDisplay2DWidget::draw_molecule(QPainter &qp) {
  string mol_name = disp_mol_->getProp<string>("_Name");
  int h = rect().height();
  if (!mol_name.empty()) {
    h = int(float(rect().height()) * 0.95);
  }

  mol_drawer_.reset(new RDKit::MolDraw2DQt(rect().width(), h, qp));

  vector<int> sa = selected_atoms();
  mol_drawer_->drawMolecule(*disp_mol_, &sa);
  add_molecule_title(qp, mol_name, h);

  identify_selected_atoms(qp);
}

// ****************************************************************************
void MolDisplay2DWidget::add_molecule_title(QPainter &qp,
                                            const string &mol_name,
                                            int label_box_height) {
  qp.setPen("Black");
  if (!mol_name.empty()) {
    int box_height = height() - label_box_height;
    qp.fillRect(0, label_box_height, width(), box_height, qp.background());
    while (true) {
      QRect br = qp.boundingRect(0, label_box_height, width(), box_height,
                                 Qt::AlignHCenter | Qt::AlignVCenter,
                                 mol_name.c_str());
      if (br.height() > box_height) {
        float scale = float(box_height) / float(br.height());
        QFont ft = qp.font();
        float new_fs = scale * ft.pointSizeF();
        ft.setPointSizeF(new_fs);
        qp.setFont(ft);
      } else {
        break;
      }
    }
    qp.drawText(0, label_box_height, width(), height() - label_box_height,
                Qt::AlignHCenter | Qt::AlignVCenter, mol_name.c_str());
  }
}

// ****************************************************************************
void MolDisplay2DWidget::identify_selected_atoms(QPainter &qp) {
  static QPen sel_pen(QColor("Orange"));
  qp.setPen(sel_pen);

  // put an orange square round selected atoms.
  vector<int> sa = selected_atoms();
  for (int i = 0, is = sa.size(); i < is; ++i) {
    Point2D at_cds = mol_drawer()->getDrawCoords(sa[i]);
    qp.drawRect(at_cds.x - pick_circle_rad(), at_cds.y - pick_circle_rad(),
                2 * pick_circle_rad(), 2 * pick_circle_rad());
  }
}

// ****************************************************************************
int MolDisplay2DWidget::find_nearest_atom(int x_screen_pos,
                                          int y_screen_pos) const {
  int nearest_at = -1, nearest_dist = numeric_limits<int>::max();

  for (int i = 0, is = disp_mol_->getNumAtoms(); i < is; ++i) {
    Point2D screen_cds = mol_drawer_->getDrawCoords(i);
    int dist = DACLIB::square(screen_cds.x - x_screen_pos) +
               DACLIB::square(screen_cds.y - y_screen_pos);
    if (dist < nearest_dist) {
      nearest_dist = dist;
      nearest_at = i;
    }
  }

  if (nearest_dist < pick_circle_rad_ * pick_circle_rad_) {
    return nearest_at;
  } else {
    return -1;
  }
}

}  // namespace RDKit
