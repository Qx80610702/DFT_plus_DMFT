
#include "timer.h"

#include <ctime>
#include <sstream>
#include <iomanip>

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

void get_date_time(std::string& date)
{
  time_t nowtime;
  struct tm* p;
  time(&nowtime);
  p = localtime(&nowtime);
   
  std::stringstream ss;
  int year = p->tm_year + 1900;
  int month = p->tm_mon + 1;
  ss << year << "." << std::setfill('0') << std::setw(2) << month 
     << "." << std::setfill('0') << std::setw(2) << p->tm_mday 
     << "--" << p->tm_hour 
     << ":" <<p->tm_min
     << ":" <<p->tm_sec;
  
  date = ss.str();
  return;
}

}