#ifndef TIMER_H
#define TIMER_H
#include <sys/resource.h>
class t_timer
{

    struct rusage usage2,usage1;
    double time;
    public:
    t_timer() : time(0) {} ;
    void start() { getrusage(RUSAGE_SELF,&usage1);};
    void stop() 
    { 
	getrusage(RUSAGE_SELF,&usage2); 
	time += (double)(usage2.ru_utime.tv_sec-usage1.ru_utime.tv_sec);
        time += (double)(usage2.ru_utime.tv_usec-usage1.ru_utime.tv_usec)*1e-6;
    }
    double gettime() { return time; };
    void reset() { time=0;};

};
#endif
