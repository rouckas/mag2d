#include <iostream>
#include "histogram.hpp"
#include "random.cpp"
#include "timer.hpp"

int main()
{
    Histogram histogram(10000, -1e-8, 1e-7);
    t_random rng(1234);
    t_timer timer;

    //speed tests
    int niter = 20000000;
    float volatile dum;
    for(int j=0; j<3; j++)
    {
        timer.reset();
        timer.start();
        for(int i=0; i<niter; i++)
        {
            dum = rng.uni();
        }
        timer.stop();
        cout << "uni() RNG speed = " << niter/timer.gettime() << " samples/sec\n";

        timer.reset();
        timer.start();
        for(int i=0; i<niter; i++)
        {
            dum = rng.uni_small();
        }
        timer.stop();
        cout << "uni_small() RNG speed = " << niter/timer.gettime() << " samples/sec\n";

        cout << endl;
    }


    exit(0);
    for(long int i=0; i<1e12; i++)
    {
        histogram.add(rng.uni());
        if(i%100000000 == 0)
            histogram.print("test_hist.txt");
    }
    return 0;
}

