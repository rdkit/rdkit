// Added by Glydade: Reimplementation of cs_swapbytes

// BSD 3-Clause License
// 
// Copyright (c) 2025, Glysade Inc
// 
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
#ifndef CS_SWAPBYTES_H
#define CS_SWAPBYTES_H

#if defined(PLATFORM_LITTLEENDIAN)
inline INT16 SwapBytes(INT16 m) {
  return m;
}

inline UINT16 SwapBytes(UINT16 m) {
  return m;
}

inline INT32 SwapBytes(INT32 m) {
  return m;
}

inline UINT32 SwapBytes(UINT32 m) {
  return m;
}


inline double SwapBytes(double m) {
  return m;
}

inline float SwapBytes(float m) {
  return m;
}

inline void SwapBytes(const float *src, char * dst) {
  memcpy(dst, src, sizeof(float));
}

inline void SwapBytes(const double *src, char * dst) {
  memcpy(dst, src, sizeof(double));
}
#elif defined(PLATFORM_BIGENDIAN)
// this function swap the bytes of values given it's size as a template
// parameter (could sizeof be used?).
template <class T, unsigned int size>
inline T SwapBytes(T value) {
  if (size < 2) {
    return value;
  }

  union {
    T value;
    char bytes[size];
  } in, out;

  in.value = value;

  for (unsigned int i = 0; i < size; ++i) {
    out.bytes[i] = in.bytes[size - 1 - i];
  }

  return out.value;
}
inline INT16 SwapBytes(INT16 m) {
  return SwapBytes<INT16, sizeof(INT16)>(m);
}

inline UINT16 SwapBytes(UINT16 m) {
  return SwapBytes<UINT16, sizeof(UINT16)>(m);
}

inline INT32 SwapBytes(INT32 m) {
  return SwapBytes<INT32, sizeof(INT32)>(m);
}

inline UINT32 SwapBytes(UINT32 m) {
  return SwapBytes<UINT32, sizeof(UINT32)>(m);
}


inline double SwapBytes(double m) {
  return SwapBytes<double,sizeof(double))>(m);
}

inline float SwapBytes(float m) {
  return SwapBytes<foud,sizeof(float))>(m);
}

inline void SwapBytes(const float *src, char * dst) {
  memcpy(dst, src, sizeof(float));
  SwapBytes(*((float*)dst));
}

inline void SwapBytes(const double *src, char * dst) {
  memcpy(dst, src, sizeof(double));
  SwapBytes(*((double*)dst));
}
#else
#error neither PLATFORM_LITTLEENDIAN nor PLATFORM_BIGENDIAN defined
#endif

#endif
