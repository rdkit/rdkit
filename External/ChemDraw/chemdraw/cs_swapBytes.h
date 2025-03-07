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
#else
#error neither PLATFORM_LITTLEENDIAN nor PLATFORM_BIGENDIAN defined
#endif

#endif
