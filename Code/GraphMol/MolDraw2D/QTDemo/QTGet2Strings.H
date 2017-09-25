//
// file QTGet2Strings.H
// David Cosgrove
// AstraZeneca
// 24th April 2014
//
// Dialog that gets 2 strings from the user.

#ifndef QTGET2STRINGS_H
#define QTGET2STRINGS_H

#include <QDialog>

// ****************************************************************************

class QLineEdit;
class QString;

// ****************************************************************************

namespace DACLIB {

class QTGet2Strings : public QDialog {

public :

  QTGet2Strings( QString prompt1 , QString initval1 ,
                 QString prompt2 , QString initval2 ,
                 QWidget *parent = 0 , Qt::WindowFlags f = 0 );

  void get_values( QString &string1 , QString &string2 );

private :

  QLineEdit *le1_ , *le2_;

  void build_widget( QString prompt1 , QString initval1 ,
                     QString prompt2 , QString initval2 );
  QWidget *build_action_box();

};

} // EO namespace DACLIB

#endif // QTGET2STRINGS_H
