//
// file QT4SelectItems.H
// Dave Cosgrove
// AstraZeneca
// 27th June 2006
//
// Puts up a QListWidget in a box allowing selections from a list of strings.

#ifndef DAC_QT4_SELECT_ITEMS__
#define DAC_QT4_SELECT_ITEMS__

#include <string>
#include <vector>

#include <QDialog>

// ***************************************************************************

class QListWidget;
class QListWidgetItem;
class QVBoxLayout;

namespace DACLIB {

class QT4SelectItems : public QDialog {

  Q_OBJECT

public :

  // selected_items needs to be big enough to hold the result, and can hold
  // pre-selections if required. I.e. any item in item_labels that has a
  // corresponding 1 in selected_items will be set as selected when the dialog
  // pops up.
  QT4SelectItems( const std::string &label , std::vector<QString> &item_labels ,
                  std::vector<char> &selected_items , bool radio_box = false ,
                  QWidget *parent = 0 );
  void get_results( std::vector<char> &selected_items ) const;

protected slots :

  void slot_list_double_clicked( QListWidgetItem *item ); // select item and out

protected :

  QWidget       *action_box_;
  QVBoxLayout   *vlayout_;
  QListWidget   *list_widget_;

  void build_action_box();

};

} // EO namespace DACLIB

#endif
