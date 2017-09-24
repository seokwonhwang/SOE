#include "Ejection.h"

using namespace SOENS;

Ejection::Ejection(double mi, double mt, double Eth, double U, double alpha)
: R1(1), R2(2), R3(3), R4(4), R5(5)
{
    mi_ = mi;
    mt_ = mt;
    U_ = U;
    Eth_ = Eth;

    // default values for S-D Distribution
    // the ratio of the maximum energy of ejected particles to incident energy
    // the maximum energy is related to the threshold energy of sputtering, Eth
    // Emax/E0 = U/Eth, Emax/E0 = 0.35 for Cu atoms under Ar+ bombardment
    emax2einc_ = U_/Eth_;

    A_ = 13.0;
    q1_ = 2.0 - 0.25*mt_/mi_;
    q2_ = 0.55;
    m_ = 0.2;

    // cosine parameter for theta distribution
    alpha2_ = alpha*alpha;
}

void Ejection::setEinc(double einc)
{
    einc_ = einc;
    fmax_ = getMaximumFuncSD(einc);
}

double Ejection::funcSD(double einc, double eps, double theta) const
// Calculation of the S-D distribution function 
// which does not include the normalization constant, C(theta)
// einc : the energy of incident particle
// eps : the energy of ejected particles
// theta : the theta angle of ejected particle
{
    double emax = emax2einc_*einc;

    return (eps*pow((eps+U_), -3.0+2.0*m_)*(1.0 - (eps+U_)/(emax+U_))
	*exp( -A_*pow( (mi_*(eps*(pow(cos(theta),q1_))+U_))/(mt_*einc) , q2_) ));
}

double Ejection::getMaximumFuncSD(double einc) const
// The function finds the maximum value of funcSD(eps, pi/2)
// einc : the energy of incident particle
{

    int ne = 100;
    double emax = emax2einc_*einc;
    double de = emax/ne;

    double eps = 0.0;
    double fmax = funcSD(einc, eps, 0.5*M_PI);

    for(int i=1; i<=ne; i++) {
	eps += de;
	double fvalue = funcSD(einc, eps, 0.5*M_PI);
	
	if(fvalue > fmax) {
	    fmax = fvalue;
	}
    }

    return fmax;
}

void Ejection::velocitySD(double& eps, double& vx, double& vy, double& vz)
// The energy of ejected particle, eps
{
    // To get a polar angle, Acceptance-Rejection Method is used.
    // The pdf of the function, g(theta) = A*cos(theta)^(alpha^2)*sin(theta)
    // g(theta) has the maximum value at theta = atan(1/alpha)
    double theta;
    double thetaFmax = atan(1/sqrt(alpha2_));
    double gmax = pow(cos(thetaFmax), alpha2_)*sin(thetaFmax);

    while(1) {
	// 1.5707963267949 = 0.5*PI
	theta = 1.5707963267949*R1.uniform();

	if( (gmax*R2.uniform()) < 
		(pow(cos(theta), alpha2_)*sin(theta)) ) {
	    break;
	}
    }

    // Acceptance-Rejection Method is also used to get the energy of ejected particle.
    // Instead of calculating the constant function, C(theta),
    // the ratio, f(epsilon|theta) = F(epsilon|theta)/F(epsilon|pi/2)
    // f(epsilion|theta) has the maximum value at theta = pi/2
    // If einc doesn't equals to the previous value of einc, 
    // the maximum value at theta = pi/2, f(epsilon|pi/2, einc) is recalculated.

    // double fmax = getMaximumFuncSD(einc);
    // maximum energy of ejected particles
    double emax = emax2einc_*einc_;

    while(1) {
	eps = R3.uniform()*emax;
	double fvalue = funcSD(einc_, eps, theta);

	if( (fmax_*R4.uniform()) < fvalue ) {
	    break;
	}
    }

    // Azimuthal angle
    double phi = 2.0*M_PI*R5.uniform(); 

    double vmag = sqrt(2*eps*ECHARGE/(mt_*UAMU));
    double sintheta = sin(theta);
    double costheta = cos(theta);
    double sinphi = sin(phi);
    double cosphi = cos(phi);

    vx = vmag*sintheta*cosphi;
    vy = vmag*sintheta*sinphi;
    vz = vmag*costheta;
}



std::ostream& operator<<(std::ostream& os, const SOENS::Ejection& ejection)
{
    os << "!------------------------ Ejection --------------------------!" << std::endl;
    //os << "M_i = " << ejection.mi_ << std::endl;
    //os << "M_t = " << ejection.mt_ << std::endl;
    os << "!-------------------- End of Ejection -----------------------!" << std::endl;

    return os;
}
