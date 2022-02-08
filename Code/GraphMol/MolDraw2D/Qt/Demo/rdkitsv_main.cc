//
// file rdkitsv_main.cc
// David Cosgrove
// AstraZeneca
// 20th June 2014
//
// This is the main function for the program rdkitsv, a simple program
// demonstrating the use of the DrawMol2D classes for drawing into a Qt
// widget.

#include "RDKitSVMainWindow.H"

#include <RDGeneral/versions.h>

#include <iostream>

#include <QApplication>

using namespace std;

// ****************************************************************************

int main(int argc, char **argv) {
  QApplication a(argc, argv);
  cout << "Built with Qt version " << QT_VERSION_STR << endl
       << "Running with Qt version " << qVersion() << endl
       << "Using RDKit version " << RDKit::rdkitVersion << endl
       << endl;

  RDKitSV::RDKitSVMainWindow *mw = new RDKitSV::RDKitSVMainWindow(argc, argv);
  mw->setWindowTitle(QObject::tr("RDKit SV"));
  mw->setGeometry(0, 0, 1000, 1000);
  mw->show();

  return a.exec();
}
