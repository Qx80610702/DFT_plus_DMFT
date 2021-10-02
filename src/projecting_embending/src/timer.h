#pragma once

#include <sys/time.h>

namespace timer
{
  void timestamp(double& time_start);
  void convert(const struct timeval T, double& time);

  void get_time(const double time_strat, double& seconds);

  void get_time(const double time_start, double& seconds, int& minutes);

  void get_time(const double time_start, double& seconds, int& minutes, int& hours);

  void get_time(const double time_start, double& seconds, int& minutes, int& hours,  int& days);
  
}
