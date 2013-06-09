#include <iostream>
#include "gas.hpp"
using namespace std;

double RAND()
{
    return (2*(rand()/(double)RAND_MAX)-1);
}

double RANDP()
{
    srand(rand());
    return (rand()/(double)RAND_MAX);
}

int main()
{
    world w;
    time_t tt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    srand(tt);
    //w.addParticle(particle(1,vect3d(0,0,0),vect3d(0,0,0),1,0,0,1));
    //w.addParticle(particle(1,vect3d(10,0.1,0),vect3d(-2,0,0),1,0,1,1));
    for(int i = 0;i<1000;++i)
    {
        i += w.addParticle(particle(0.15,vect3d(10*RAND(),10*RAND(),0),vect3d(9*RAND(),9*RAND(),0),1.0,1,0,0));
    }
    w.animate(0.002,_FOREVER);
    cout << "Hello world!" << endl;
    return 0;
}
