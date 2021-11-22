//
// file QT4SelectItems.H
// Dave Cosgrove
// AstraZeneca
// 27th June 2006
//
// Puts up a QListWidget in a box allowing selections from a list of strings.
// Similar to QTSelectItems.H but in Qt4 speak.
// This must always be used as a modal dialog, since the address of the vector
// holding the selections is kept, and this must be in scope when it is filled.
// It is filled when either the ok or cancel buttons are pressed.

#include <QLabel>
#include <QLayout>
#include <QListWidget>
#include <QListWidgetItem>
#include <QPushButton>

#include "QT4SelectItems.H"

#include <iostream>

using namespace std;

namespace DACLIB {

// ***************************************************************************
QT4SelectItems::QT4SelectItems(const string &label,
                               vector<QString> &item_labels,
                               vector<char> &selected_items, bool radio_box,
                               QWidget *parent)
    : QDialog(parent) {
  setWindowTitle("Select Items");

  QWidget *vbox = new QWidget(this);

  list_widget_ = new QListWidget();

  if (radio_box) {
    list_widget_->setSelectionMode(QAbstractItemView::SingleSelection);
  } else {
    list_widget_->setSelectionMode(QAbstractItemView::ExtendedSelection);
  }

  vector<QString>::iterator p, ps;
  int i = 0;
  for (p = item_labels.begin(), ps = item_labels.end(); p != ps; ++p, ++i) {
    QListWidgetItem *item = new QListWidgetItem(*p, list_widget_);
    list_widget_->setItemSelected(item, selected_items[i]);
  }

  connect(list_widget_, SIGNAL(itemDoubleClicked(QListWidgetItem *)), this,
          SLOT(slot_list_double_clicked(QListWidgetItem *)));

  build_action_box();

  vlayout_ = new QVBoxLayout();
  vlayout_->setDirection(QBoxLayout::BottomToTop);
  vlayout_->addWidget(action_box_);
  vlayout_->addWidget(list_widget_);
  vlayout_->addWidget(new QLabel(label.c_str()));

  vbox->setLayout(vlayout_);
  vbox->show();
}

// *****************************************************************************
void QT4SelectItems::get_results(vector<char> &selected_items) const {
  selected_items = vector<char>(list_widget_->count(), 0);
  for (int i = 0, is = list_widget_->count(); i < is; ++i) {
    selected_items[i] = list_widget_->isItemSelected(list_widget_->item(i));
  }
}

// *****************************************************************************
void QT4SelectItems::build_action_box() {
  action_box_ = new QWidget(this);
  QHBoxLayout *hlayout = new QHBoxLayout();

  QPushButton *button = new QPushButton("Ok");
  button->setDefault(true);
  hlayout->addWidget(button);
  connect(button, SIGNAL(clicked()), this, SLOT(accept()));

  button = new QPushButton("Cancel");
  hlayout->addWidget(button);
  connect(button, SIGNAL(clicked()), this, SLOT(reject()));

  action_box_->setLayout(hlayout);
}

// *****************************************************************************
// select item and out
void QT4SelectItems::slot_list_double_clicked(QListWidgetItem *item) {
  item->setSelected(true);
  accept();
}

}  // namespace DACLIB
