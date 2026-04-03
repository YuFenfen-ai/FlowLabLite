// now_clock.c
// Cross-compiler minimal implementation of now_seconds()
// Windows: uses QueryPerformanceCounter for high-resolution time
// Fallback: returns 0.0

#ifdef _WIN32
#include <windows.h>

__declspec(dllexport) double now_seconds(void) {
  static LARGE_INTEGER freq = {0};
  LARGE_INTEGER cnt;
  if (freq.QuadPart == 0) QueryPerformanceFrequency(&freq);
  QueryPerformanceCounter(&cnt);
  return (double)cnt.QuadPart / (double)freq.QuadPart;
}
#else
#include <time.h>

double now_seconds(void) {
  struct timespec ts;
#if defined(CLOCK_REALTIME)
  if (clock_gettime(CLOCK_REALTIME, &ts) == 0) {
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
  }
#endif
  return 0.0;
}
#endif
