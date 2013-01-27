#ifndef TIMER_HPP
#define TIMER_HPP

class Timer {
public:
  double prevtime, curtime;
  double update_time(void);
  Timer(void);
};

#endif
