#ifndef EJECTION_H
#define EJECTION_H

#include "RNG.h"
#include "Constants.h"
#include <cmath>
#include <iostream>

// Ejected particles have Stepanova and Dew's distribution,
// Only normally incident particle is valid
// Reference : K. Nanbu, T. Ohshita, Vacuum 87 (2013) 103-108

// NameSpace of Simulation Of Everything //
namespace SOENS {

class Ejection
{
    friend std::ostream& operator<<(std::ostream& os,
				    const Ejection& ejection);
protected:
    // Projectile mass, [u]
    double mi_;

    // Target mass, [u]
    double mt_;

    // Binding energy of the target solid, [eV]
    double U_;
    
    // Threshold energy for sputtering, [eV]
    double Eth_;

    // Parameters for S-D distribution
    double A_;
    double q1_;
    double q2_;
    double m_;
    double emax2einc_;

    // a parameter shaping the distribution to be over or under-cosine, respectilvely
    // g(theta) ~ cos(theta)^(alpha^2), default value of 1
    double alpha2_;

    // the energy of the incident particle, [eV]
    double einc_;

    // the maximum value of S-D distribution function
    double fmax_;

    // uniform number generator
    RNG R1, R2, R3, R4, R5;

private:
    inline double funcSD(double einc, double eps, double theta) const;
    inline double getMaximumFuncSD(double einc) const;

public:
    // Constructor : define incident mass, the properties of the target
    // and parameters for distribution function
    Ejection(double mi, double mt, double Eth, double U, double alpha);

    // Set the energy of incident particles, [eV]
    void setEinc(double einc);

    // Get the velocity of the distribution, S-D
    void velocitySD(double& eps, double& vx, double& vy, double& vz);

};

} // End of namespace, SOENS

#endif
