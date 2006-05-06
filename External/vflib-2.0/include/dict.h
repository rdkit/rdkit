/*----------------------------------------------------
 * dict.h
 * A simple dictionary, not very efficient,
 * implemented with linked lists.
 * Author: P. Foggia
 * $Id: dict.h,v 1.1 2001/10/24 01:20:11 glandrum Exp $
 ----------------------------------------------------*/

/*----------------------------------------------------
 * REVISION HISTORY
 *   $Log: dict.h,v $
 *   Revision 1.1  2001/10/24 01:20:11  glandrum
 *   added
 *
 *   Revision 1.2  1998/12/19 11:35:45  foggia
 *   Minor changes
 *
 *   Revision 1.1  1998/12/13 13:37:24  foggia
 *   Initial revision
 *
 ---------------------------------------------------*/

#ifndef DICT_H
#define DICT_H

#include <stddef.h>
#include <iostream.h>

#include "error.h"

#ifdef NEED_BOOL
#ifndef BOOL_DEFINED
#define BOOL_DEFINED
typedef int bool;
const int false=0;
const int true=1
#endif
#endif


template <class Key, class Value>
struct DictionaryNode
  { Key key;
    Value value;
    DictionaryNode<Key, Value> *next;
  };


template <class Key, class Value>
class DictionaryIterator;

template <class Key, class Value>
class Dictionary
  { private:
      typedef DictionaryNode<Key, Value> Node;
      Node *head;

      friend class DictionaryIterator<Key, Value>;

    public:
      typedef DictionaryIterator<Key, Value> iterator_type;

      Dictionary(); 
      ~Dictionary();
      void clear();
      Value *get(const Key &k);
      void put(const Key &k, const Value &v);
      iterator_type *iterator();
  };

template <class Key, class Value>
class DictionaryIterator
  { private:
      DictionaryNode<Key, Value>* node;

    public:
      DictionaryIterator(Dictionary<Key, Value> &dict)
                  { node=dict.head; }
      bool more() { return node!=NULL; }
      void next() { if (node!=NULL) node=node->next; }
      const Key &key() { return node->key; }
      Value &value()   { return node->value; }
  };

template <class Key, class Value>
istream& operator>>(istream& in, Dictionary<Key, Value> &dict); 
template <class Key, class Value>
ostream& operator<<(ostream& out, Dictionary<Key, Value> &dict); 

/*-----------------------------------------------------------
 * Methods of class Dictionary
 ---------------------------------------------------------*/
template <class Key, class Value>
Dictionary<Key,Value>::Dictionary()
  { head=NULL; 
  }

template <class Key, class Value>
Dictionary<Key,Value>::~Dictionary()
  { clear();
  }

template <class Key, class Value>
void Dictionary<Key,Value>::clear()
  { while (head!=NULL)
      { Node *alias=head;
        head=head->next;
        delete alias;
       }
  }

template <class Key, class Value>
Value *Dictionary<Key,Value>::get(const Key &k)
  { Node *p;
    for(p=head; p!=NULL && p->key!=k; p=p->next)
      ; 
    if (p==NULL)
      return NULL;
    else
      return &p->value;
  }

template <class Key, class Value>
void Dictionary<Key,Value>::put(const Key &k, const Value &v)
  { Node *p;
    Value *vp=get(k);

    if (vp!=NULL)
      *vp=v;
    else
      { p=new Node;
        if (p==NULL)
          error("Out of memory");
        p->key=k;
        p->value=v;
        p->next=head;
        head=p;
      }
  }

template <class Key, class Value>
DictionaryIterator<Key,Value> *
Dictionary<Key,Value>::iterator()
  { DictionaryIterator<Key,Value> *p=new DictionaryIterator<Key,Value>(*this);
    if (p==NULL)
      error("Out of memory");
    return p;
  }


template <class Key, class Value>
istream& operator>>(istream& in, Dictionary<Key, Value> &dict)
  { int n, i;
    Key k;
    Value v;

    in>>n;
    dict.clear();

    for(i=0; i<n; i++)
      { in >> k >> v;
        dict.put(k, v);
      }

    return in;
  }

template <class Key, class Value>
ostream& operator<<(ostream& out, Dictionary<Key, Value> &dict)
  { DictionaryIterator<Key, Value> *iter;

    int ct=0;
    iter=dict.iterator();
    while (iter->more())
      { ct++;
        iter->next();
      }
    delete iter;

    out << ct ;

    iter=dict.iterator();
    while (iter->more())
      { out << ' ' << iter->key() << ' ' << iter->value();
        iter->next();
      }
    delete iter;

    return out;

  }

#endif
