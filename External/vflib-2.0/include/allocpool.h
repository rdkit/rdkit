/*-----------------------------------------------------
 * allocpool.h
 * Definition of an allocation pool class.
 *
 * Author: P. Foggia
 * $Id: allocpool.h,v 1.1 2001/10/24 01:20:11 glandrum Exp $
 ----------------------------------------------------*/

/*----------------------------------------------------
 * REVISION HISTORY
 *   $Log: allocpool.h,v $
 *   Revision 1.1  2001/10/24 01:20:11  glandrum
 *   added
 *
 *   Revision 1.4  1998/12/19 11:35:24  foggia
 *   Minor changes
 *
 *   Revision 1.3  1998/12/12 12:16:33  foggia
 *   Name changes
 *
 *   Revision 1.2  1998/09/26 09:02:32  foggia
 *   Minor changes
 *
 *   Revision 1.1  1998/09/19 14:40:35  foggia
 *   Initial revision
 *
 ---------------------------------------------------*/

/*----------------------------------------------------------------------
 *   CLASS DESCRIPTION
 * An allocation pool is an efficient way to allocate on heap
 * a set of objects of the same type that will be deleted
 * at the same time. Instead of allocating each object separately,
 * thus incurring in both space and time overhead, "chunks" of 
 * objects are allocated. The allocation pool mantains a list of
 * the chunks, that will be deleted by the pool distructor.
 *
 * Two classes are provided: an Allocator, which is an
 * interface independent of actual allocation strategy, and 
 * AllocationPool which is the real class. 
 ---------------------------------------------------------------------*/

#ifndef ALLOCPOOL_H
#define ALLOCPOOL_H


#include <stddef.h>
#include "error.h"

/*-----------------------------------------------------------
 * Declaration of class Allocator,
 * which represents an abstract way of allocating objects
 *---------------------------------------------------------*/

template <class T>
class Allocator
  { protected:
      virtual T *alloc() = 0;
    public:
      virtual ~Allocator() {}
      virtual T *Allocate() 
          { T*p=alloc();
            if (p==NULL)
              error("Out of memory");
            return p;
          }
  };


/*-----------------------------------------------------------
 * Declaration of class NewAllocator,
 * which simply use new to allocate objects
 ----------------------------------------------------------*/
template <class T>
class NewAllocator: public Allocator<T>
  { protected:
      virtual T* alloc() { return new T; }
  };

/*-----------------------------------------------------------
 * Declaration of class NullAllocator, which always
 * return a null pointer
 ----------------------------------------------------------*/
template <class T>
class NullAllocator: public Allocator<T>
  { public:
      virtual T *Allocate() { return NULL; }
	protected:
	  virtual T *alloc() { return NULL; }
  };


/*-----------------------------------------------------------
 * Declaration of class AllocationPool
 *---------------------------------------------------------*/

template <class T, int CHUNK_SIZE>
class AllocationPool: public Allocator<T>
  { private:
      struct chunk
        { chunk *link;
          T content[CHUNK_SIZE];
        };
      chunk *chunkList;
      int rest;
      void grow();
	protected:
      virtual T *alloc();
    public:
      AllocationPool();
      ~AllocationPool();
  };



/*----------------------------------------------------
 * inline methods of class AllocationPool
 ---------------------------------------------------*/
template <class T, int CHUNK_SIZE>
AllocationPool<T, CHUNK_SIZE>::AllocationPool() 
  { chunkList=0; 
    rest=0; 
  }

template <class T, int CHUNK_SIZE>
AllocationPool<T, CHUNK_SIZE>::~AllocationPool() 
  { chunk *p=chunkList;
    while (p!=0)
      { chunkList=p->link;
        delete p;
        p=chunkList;
      }
  }

template <class T, int CHUNK_SIZE>
T * AllocationPool<T, CHUNK_SIZE>::alloc() 
  { if (rest==0)
      grow();
    if (rest>0)
      return chunkList->content + (--rest);
    else
      return 0;
  }

template <class T, int CHUNK_SIZE>
void AllocationPool<T, CHUNK_SIZE>::grow() 
  { chunk *p=new chunk;
    if (p!=0)
      { p->link=chunkList;
        chunkList=p;
        rest=CHUNK_SIZE;
      }
  }

#endif
