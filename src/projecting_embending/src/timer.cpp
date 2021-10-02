
#include "timer.h"

namespace timer
{

void timestamp(double& time_start)
{
  struct timeval T;

  gettimeofday(&T,nullptr);
  convert(T, time_start);

  return;
}

void convert(const struct timeval T, double& time)
{
  time = (double)T.tv_sec+ (double)T.tv_usec/1.0e6;  //in second
  return;
}

void get_time(const double time_start, double& seconds)
{
  struct timeval T;
  double time_now;

  gettimeofday(&T,nullptr);
  convert(T, time_now);

  seconds=time_now-time_start;

  return;
}


void get_time(const double time_start, double& seconds, int& minutes)
{
  struct timeval T;
  double time_now, diff;

  gettimeofday(&T,nullptr);

  convert(T, time_now);

  diff = time_now-time_start;

  minutes=(int)diff/60;
  seconds=diff-minutes*60;

  return;
}

void get_time(const double time_start, double& seconds, int& minutes, int& hours)
{
  struct timeval T;
  double time_now, diff;

  gettimeofday(&T,nullptr);

  convert(T, time_now);

  diff = time_now-time_start;

  hours=(int)(diff/3600);
  minutes=(int)((diff-hours*3600)/60);
  seconds=diff-hours*3600-minutes*60;

  return;
}


void get_time(const double time_start, double& seconds, int& minutes, int& hours,  int& days)
{
  struct timeval T;
  double time_now, diff;

  gettimeofday(&T,nullptr);

  convert(T, time_now);

  diff = time_now-time_start;

  days=(int)(diff/(3600*24));
  hours=(int)((diff-days*3600*24)/3600);
  minutes=(int)((diff-days*3600*24-hours*3600)/60);
  seconds=diff-days*3600*24-hours*3600-minutes*60;

  return;
}


}