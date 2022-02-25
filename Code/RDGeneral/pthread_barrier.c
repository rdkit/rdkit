/*
 * Copyright (c) 2015, Aleksey Demakov
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// Code obtained from https://github.com/ademakov/DarwinPthreadBarrier.git
#include "pthread_barrier.h"

#include <errno.h>

#ifdef __APPLE__

#define __unused __attribute__((unused))

int
pthread_barrierattr_init(pthread_barrierattr_t *attr __unused)
{
	return 0;
}

int
pthread_barrierattr_destroy(pthread_barrierattr_t *attr __unused)
{
	return 0;
}

int
pthread_barrierattr_getpshared(const pthread_barrierattr_t *__restrict attr __unused,
			       int *__restrict pshared)
{
	*pshared = PTHREAD_PROCESS_PRIVATE;
	return 0;
}

int
pthread_barrierattr_setpshared(pthread_barrierattr_t *attr __unused,
			       int pshared)
{
	if (pshared != PTHREAD_PROCESS_PRIVATE) {
		errno = EINVAL;
		return -1;
	}
	return 0;
}

int
pthread_barrier_init(pthread_barrier_t *__restrict barrier,
		     const pthread_barrierattr_t *__restrict attr __unused,
		     unsigned count)
{
	if (count == 0) {
		errno = EINVAL;
		return -1;
	}

	if (pthread_mutex_init(&barrier->mutex, 0) < 0) {
		return -1;
	}
	if (pthread_cond_init(&barrier->cond, 0) < 0) {
		int errno_save = errno;
		pthread_mutex_destroy(&barrier->mutex);
		errno = errno_save;
		return -1;
	}

	barrier->limit = count;
	barrier->count = 0;
	barrier->phase = 0;

	return 0;
}

int
pthread_barrier_destroy(pthread_barrier_t *barrier)
{
    pthread_mutex_destroy(&barrier->mutex);
    pthread_cond_destroy(&barrier->cond);
    return 0;
}

int
pthread_barrier_wait(pthread_barrier_t *barrier)
{
	pthread_mutex_lock(&barrier->mutex);
	barrier->count++;
	if (barrier->count >= barrier->limit) {
		barrier->phase++;
		barrier->count = 0;
		pthread_cond_broadcast(&barrier->cond);
		pthread_mutex_unlock(&barrier->mutex);
		return PTHREAD_BARRIER_SERIAL_THREAD;
	} else {
		unsigned phase = barrier->phase;
		do
			pthread_cond_wait(&barrier->cond, &barrier->mutex);
		while (phase == barrier->phase);
		pthread_mutex_unlock(&barrier->mutex);
		return 0;
	}
}

#endif /* __APPLE__ */
