//
// file QTGet2Strings.cc
// David Cosgrove
// AstraZeneca
// 24th April 2014
//

#include "QTGet2Strings.H"

#include <QBoxLayout>
#include <QFormLayout>
#include <QFrame>
#include <QLineEdit>
#include <QPushButton>
#include <QString>

namespace DACLIB {

// ****************************************************************************
QTGet2Strings::QTGet2Strings(QString prompt1, QString initval1, QString prompt2,
                             QString initval2, QWidget *parent,
                             Qt::WindowFlags f)
    : QDialog(parent, f) {
  build_widget(prompt1, initval1, prompt2, initval2);
}

// ****************************************************************************
void QTGet2Strings::get_values(QString &string1, QString &string2) {
  string1 = le1_->text();
  string2 = le2_->text();
}

// ****************************************************************************
void QTGet2Strings::build_widget(QString prompt1, QString initval1,
                                 QString prompt2, QString initval2) {
  QVBoxLayout *vbox = new QVBoxLayout;

  QFormLayout *fl = new QFormLayout;
  le1_ = new QLineEdit;
  fl->addRow(prompt1, le1_);
  le1_->setText(initval1);
  le2_ = new QLineEdit;
  fl->addRow(prompt2, le2_);
  le2_->setText(initval2);

  vbox->addLayout(fl);
  vbox->addWidget(build_action_box());
  setLayout(vbox);
}

// ****************************************************************************
QWidget *QTGet2Strings::build_action_box() {
  QFrame *action_frame = new QFrame;
  action_frame->setFrameStyle(QFrame::Box);

  QHBoxLayout *hlayout = new QHBoxLayout;

  QPushButton *button = new QPushButton("Ok");
  hlayout->addWidget(button);
  connect(button, SIGNAL(clicked()), this, SLOT(accept()));
  button->setDefault(true);
  button->setAutoDefault(true);

  button = new QPushButton("Cancel");
  hlayout->addWidget(button);
  connect(button, SIGNAL(clicked()), this, SLOT(reject()));
  button->setDefault(false);
  button->setAutoDefault(false);

  action_frame->setLayout(hlayout);

  return action_frame;
}

}  // namespace DACLIB
