#ifndef TIMER_H
#define TIMER_H
#include <sys/resource.h>
#include <sys/time.h>

class t_timer
{

    struct rusage usage2,usage1;
    double time;
    double real_time;
    timeval tval1, tval2;
    public:
    t_timer() : time(0), real_time(0) {} ;
    void start()
    {
        getrusage(RUSAGE_SELF,&usage1);
        gettimeofday(&tval1, NULL);
    };
    void stop() 
    { 
	getrusage(RUSAGE_SELF,&usage2); 
	time += (double)(usage2.ru_utime.tv_sec-usage1.ru_utime.tv_sec);
        time += (double)(usage2.ru_utime.tv_usec-usage1.ru_utime.tv_usec)*1e-6;

        gettimeofday(&tval2, NULL);
        real_time += (double)(tval2.tv_sec-tval1.tv_sec);
        real_time += (double)(tval2.tv_usec-tval1.tv_usec)*1e-6;
    }
    double get_cpu_time() { return time; };
    double get_real_time() { return real_time; };
    void reset() { time=0;};

};
#endif
