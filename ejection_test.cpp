#include <iostream>
#include "Ejection.h"
#include "MemAlloc.hpp"


int main(int argc, char** argv)
{
    SOENS::Ejection eject(40.0, 12.0, 6.0, 5.0, 1.0);
    eject.setEinc(10.0);
   
    double vmag, vx, vy, vz;
    int* pdf = memAllocInit<int>(100.0, 0);
    double de = 10.0/100.0;

    for(int i=0; i<10000; i++) {
	eject.velocitySD(vmag, vx, vy, vz);

	int ii = vmag/de;

	++pdf[ii];
    }

    for(int i=0; i<100; i++) {
	std::cout << (i+0.5)*de << "\t" << pdf[i] << std::endl;
    }

    return 0;
}
