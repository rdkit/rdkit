//
//  Copyright (C) 2003-2007 Greg Landrum and Rational Discovery LLC
//  Copyright (C) 2017-2019 Greg Landrum and NextMove Software
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RD_LEADERPICKER_H
#define RD_LEADERPICKER_H

#include <RDGeneral/types.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Exceptions.h>
#include <cstdlib>
#include "DistPicker.h"

namespace RDPickers {

/*! \brief Implements the Leader algorithm for picking a subset of item from a
 *pool
 *
 *  This class inherits from the DistPicker and implements a specific picking
 *strategy aimed at diversity. See documentation for "pick()" member function
 *for the algorithm details
 */
class LeaderPicker : public DistPicker {
 public:
  double default_threshold{0.0};
  int default_nthreads{1};

  /*! \brief Default Constructor
   *
   */
  LeaderPicker()  {}
  LeaderPicker(double threshold)
      : default_threshold(threshold), default_nthreads(1) {}
  LeaderPicker(double threshold, int nthreads)
      : default_threshold(threshold), default_nthreads(nthreads) {}

  /*! \brief Contains the implementation for a lazy Leader diversity picker
   *
   * See the documentation for the pick() method for details about the algorithm
   *
   *   \param func - a function (or functor) taking two unsigned ints as
   *arguments and returning the distance (as a double) between those two
   *elements.
   *
   *   \param poolSize - the size of the pool to pick the items from. It is
   *assumed that the distance matrix above contains the right number of
   *elements; i.e. poolSize*(poolSize-1)
   *
   *   \param pickSize - the number items to pick from pool (<= poolSize)
   *
   *   \param firstPicks - (optional)the first items in the pick list
   *
   *   \param seed - (optional) seed for the random number generator. If this is
   *<0 the generator will be seeded with a random number.
   */
  template <typename T>
  RDKit::INT_VECT lazyPick(T &func, unsigned int poolSize,
                           unsigned int pickSize) const;

  template <typename T>
  RDKit::INT_VECT lazyPick(T &func, unsigned int poolSize,
                           unsigned int pickSize, double threshold) const;

  template <typename T>
  RDKit::INT_VECT lazyPick(T &func, unsigned int poolSize,
                           unsigned int pickSize,
                           const RDKit::INT_VECT &firstPicks,
                           double threshold) const;

  template <typename T>
  RDKit::INT_VECT lazyPick(T &func, unsigned int poolSize,
                           unsigned int pickSize,
                           const RDKit::INT_VECT &firstPicks, double threshold,
                           int nthreads) const;

  /*! \brief Contains the implementation for the Leader diversity picker
   *
   *   \param distMat - distance matrix - a vector of double. It is assumed that
   *only the lower triangle element of the matrix are supplied in a 1D array\n
   *
   *   \param poolSize - the size of the pool to pick the items from. It is
   *assumed that the distance matrix above contains the right number of
   *elements; i.e. poolSize*(poolSize-1) \n
   *
   *   \param pickSize - maximum number items to pick from pool (<= poolSize)
   *
   *   \param firstPicks - indices of the items used to seed the pick set.
   */
  RDKit::INT_VECT pick(const double *distMat, unsigned int poolSize,
                       unsigned int pickSize, const RDKit::INT_VECT &firstPicks,
                       double threshold, int nthreads) const {
    CHECK_INVARIANT(distMat, "Invalid Distance Matrix");
    if (!poolSize) throw ValueErrorException("empty pool to pick from");
    if (poolSize < pickSize)
      throw ValueErrorException("pickSize cannot be larger than the poolSize");
    distmatFunctor functor(distMat);
    return this->lazyPick(functor, poolSize, pickSize, firstPicks, threshold,
                          nthreads);
  }

  /*! \overload */
  RDKit::INT_VECT pick(const double *distMat, unsigned int poolSize,
                       unsigned int pickSize) const {
    RDKit::INT_VECT iv;
    return pick(distMat, poolSize, pickSize, iv, default_threshold,
                default_nthreads);
  }
};

#ifdef USE_THREADED_LEADERPICKER
// Note that this block of code currently only works on linux (which is why it's
// disabled by default) We will revisit this during the 2020.03 release cycle in
// order to get a multi-threaded version of the LeaderPicker that works on all
// supported platforms
template <typename T>
void *LeaderPickerWork(void *arg);

template <typename T>
struct LeaderPickerState {
  typedef struct {
    int *ptr;
    unsigned int capacity;
    unsigned int len;
    unsigned int next[2];
  } LeaderPickerBlock;
  typedef struct {
    LeaderPickerState<T> *stat;
    pthread_t tid;
    unsigned int id;
  } LeaderPickerThread;

  std::vector<LeaderPickerThread> threads;
  std::vector<LeaderPickerBlock> blocks;
  pthread_barrier_t wait;
  pthread_barrier_t done;
  std::vector<int> v;
  LeaderPickerBlock *head_block;
  unsigned int thread_op;
  unsigned int nthreads;
  unsigned int tick;
  double threshold;
  int query;
  T *func;

  LeaderPickerState(unsigned int count, int nt) {
    v.resize(count);
    for (unsigned int i = 0; i < count; i++) v[i] = i;

    // InitializeBlocks
    unsigned int bcount;
    unsigned int bsize;
    if (nt > 1) {
      bsize = 4096;
      bcount = (count + (bsize - 1)) / bsize;
      unsigned int tasks = (bcount + 1) / 2;
      // limit number of threads to available work
      if (nt > (int)tasks) nt = tasks;
    } else {
      bsize = 32768;
      bcount = (count + (bsize - 1)) / bsize;
    }
    blocks.resize(bcount);
    head_block = &blocks[0];
    tick = 0;

    if (bcount > 1) {
      int *ptr = &v[0];
      unsigned int len = count;
      for (unsigned int i = 0; i < bcount; i++) {
        LeaderPickerBlock *block = &blocks[i];
        block->ptr = ptr;
        if (len > bsize) {
          block->capacity = bsize;
          block->len = bsize;
          block->next[0] = i + 1;
        } else {
          block->capacity = len;
          block->len = len;
          block->next[0] = 0;
          break;
        }
        ptr += bsize;
        len -= bsize;
      }
    } else {
      head_block->capacity = count;
      head_block->len = count;
      head_block->next[0] = 0;
      head_block->next[1] = 0;
      head_block->ptr = &v[0];
    }

    // InitializeThreads
    if (nt > 1) {
      nthreads = nt;
      pthread_barrier_init(&wait, NULL, nthreads + 1);
      pthread_barrier_init(&done, NULL, nthreads + 1);

      threads.resize(nt);
      for (unsigned int i = 0; i < nthreads; i++) {
        threads[i].id = i;
        threads[i].stat = this;
        pthread_create(&threads[i].tid, NULL, LeaderPickerWork<T>,
                       (void *)&threads[i]);
      }
    } else
      nthreads = 1;
  }

  ~LeaderPickerState() {
    if (nthreads > 1) {
      thread_op = 1;
      pthread_barrier_wait(&wait);
      for (unsigned int i = 0; i < nthreads; i++)
        pthread_join(threads[i].tid, 0);
      pthread_barrier_destroy(&wait);
      pthread_barrier_destroy(&done);
    }
  }

  bool empty() {
    while (head_block) {
      if (head_block->len) return false;
      unsigned int next_tick = head_block->next[tick];
      if (!next_tick) return true;
      head_block = &blocks[next_tick];
    }
    return true;
  }

  unsigned int compact(int *dst, int *src, unsigned int len) {
    unsigned int count = 0;
    for (unsigned int i = 0; i < len; i++) {
      if ((*func)(query, src[i]) > threshold) dst[count++] = src[i];
    }
    return count;
  }

  void compact_job(unsigned int cycle) {
    // On entry, next[tick] for each block is the current linked list.
    // On exit, next[tock] is the linked list for the next iteration.
    unsigned int tock = tick ^ 1;

    LeaderPickerBlock *list = head_block;
    for (;;) {
      unsigned int next_tick = list->next[tick];
      if (next_tick) {
        LeaderPickerBlock *next = &blocks[next_tick];
        unsigned int next_next_tick = next->next[tick];
        if (cycle == 0) {
          list->len = compact(list->ptr, list->ptr, list->len);
          if (list->len + next->len <= list->capacity) {
            list->len += compact(list->ptr + list->len, next->ptr, next->len);
            list->next[tock] = next_next_tick;
          } else {
            next->len = compact(next->ptr, next->ptr, next->len);
            if (next->len) {
              list->next[tock] = next_tick;
              next->next[tock] = next_next_tick;
            } else
              list->next[tock] = next_next_tick;
          }
          cycle = nthreads - 1;
        } else
          cycle--;
        if (next_next_tick) {
          list = &blocks[next_next_tick];
        } else
          break;
      } else {
        if (cycle == 0) {
          list->len = compact(list->ptr, list->ptr, list->len);
          list->next[tock] = 0;
        }
        break;
      }
    }
  }

  void compact(int pick) {
    query = pick;
    if (nthreads > 1) {
      thread_op = 0;
      pthread_barrier_wait(&wait);
      pthread_barrier_wait(&done);
    } else
      compact_job(0);
    tick ^= 1;
  }

  int compact_next() {
    compact(head_block->ptr[0]);
    return query;
  }
};

// This is the loop the worker threads run
template <typename T>
void *LeaderPickerWork(void *arg) {
  typename LeaderPickerState<T>::LeaderPickerThread *thread;
  thread = (typename LeaderPickerState<T>::LeaderPickerThread *)arg;
  LeaderPickerState<T> *stat = thread->stat;

  for (;;) {
    pthread_barrier_wait(&stat->wait);
    if (stat->thread_op) return (void *)0;
    stat->compact_job(thread->id);
    pthread_barrier_wait(&stat->done);
  }
}
#else

template <typename T>
struct LeaderPickerState {
  std::vector<int> v;
  unsigned int left;
  double threshold;
  int query;
  T *func;

  LeaderPickerState(unsigned int count, int) : left(count), threshold(0.0), query(0), func(nullptr) {
    v.resize(count);
    for (unsigned int i = 0; i < count; i++) v[i] = i;
  }

  bool empty() { return left == 0; }

  unsigned int compact(int *dst, int *src, unsigned int len) {
    unsigned int count = 0;
    for (unsigned int i = 0; i < len; i++) {
      double ld = (*func)(query, src[i]);
      // std::cerr << query << "-" << src[i] << " " << ld << std::endl;
      if (ld > threshold) dst[count++] = src[i];
    }
    return count;
  }

  void compact(int pick) {
    query = pick;
    left = compact(&v[0], &v[0], left);
  }

  int compact_next() {
    query = v[0];
    left = compact(&v[0], &v[1], left - 1);
    return query;
  }
};

#endif
// we implement this here in order to allow arbitrary functors without link
// errors
template <typename T>
RDKit::INT_VECT LeaderPicker::lazyPick(T &func, unsigned int poolSize,
                                       unsigned int pickSize,
                                       const RDKit::INT_VECT &firstPicks,
                                       double threshold, int nthreads) const {
  if (!poolSize) throw ValueErrorException("empty pool to pick from");

  if (poolSize < pickSize)
    throw ValueErrorException("pickSize cannot be larger than the poolSize");

  if (!pickSize) pickSize = poolSize;
  RDKit::INT_VECT picks;

  LeaderPickerState<T> stat(poolSize, nthreads);
  stat.threshold = threshold;
  stat.func = &func;

  unsigned int picked = 0;  // picks.size()
  unsigned int pick = 0;

  if (!firstPicks.empty()) {
    for (RDKit::INT_VECT::const_iterator pIdx = firstPicks.begin();
         pIdx != firstPicks.end(); ++pIdx) {
      pick = static_cast<unsigned int>(*pIdx);
      if (pick >= poolSize) {
        throw ValueErrorException("pick index was larger than the poolSize");
      }
      picks.push_back(pick);
      stat.compact(pick);
      picked++;
    }
  }

  while (picked < pickSize && !stat.empty()) {
    pick = stat.compact_next();
    picks.push_back(pick);
    picked++;
  }
  return picks;
}

template <typename T>
RDKit::INT_VECT LeaderPicker::lazyPick(T &func, unsigned int poolSize,
                                       unsigned int pickSize) const {
  RDKit::INT_VECT firstPicks;
  return LeaderPicker::lazyPick(func, poolSize, pickSize, firstPicks,
                                default_threshold, default_nthreads);
}

template <typename T>
RDKit::INT_VECT LeaderPicker::lazyPick(T &func, unsigned int poolSize,
                                       unsigned int pickSize,
                                       double threshold) const {
  RDKit::INT_VECT firstPicks;
  return LeaderPicker::lazyPick(func, poolSize, pickSize, firstPicks, threshold,
                                default_nthreads);
}
template <typename T>
RDKit::INT_VECT LeaderPicker::lazyPick(T &func, unsigned int poolSize,
                                       unsigned int pickSize,
                                       const RDKit::INT_VECT &firstPicks,
                                       double threshold) const {
  return LeaderPicker::lazyPick(func, poolSize, pickSize, firstPicks, threshold,
                                default_nthreads);
}

};  // namespace RDPickers

#endif
