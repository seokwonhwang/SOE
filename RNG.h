#ifndef RNG_H
#define RNG_H

// NameSpace of Simulation Of Everything //
namespace SOENS
{

class RNG
{
private:
    long int seed_;

public:
    RNG() {
	seed_ = 1;
    }
    RNG(long int seed) {
	seed_= seed;
    }
    void resetSeed(long int seed) {
	seed_ = seed;
    }

    // get the random number
    double uniform() {
	long int a = 16807;
	long int m = 2147483647;
	long int q = 127773;
	long int r = 2836;
	long int hi, lo;

	hi = seed_/q;
	lo = seed_ - q*hi;
	seed_ = a*lo - r*hi;

	/* "seed_" will always be a legal integer of 32 bits (including sign). */
	if(seed_ <= 0) seed_ = seed_ + m;

	return (seed_/2147483646.0);
    }
     
};

} // End of namespace, SOENS


#endif
