#include <cstdlib>
#include <sys/time.h>

#include "Timer.hpp"

double Timer::update_time(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  prevtime = curtime;
  curtime = tv.tv_sec + 1e-6 * tv.tv_usec;
  return curtime - prevtime;
}

Timer::Timer(void) {
  update_time();
}
