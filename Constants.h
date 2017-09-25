#ifndef CONSTANTS_H
#define CONSTANTS_H

// NameSpace of Simulation Of Everything //
namespace SOENS {

/* PI */
#ifndef M_PI
#define M_PI	    3.14159265358979
#endif

/* BOLTZMANN CONSTANT [J/K] */
#ifndef KB
#define KB          1.3806488E-23
#endif

/* VACUUM PERMITTIVITY [F/m] */
#ifndef EPS0
#define EPS0        8.8541E-12
#endif

/* ELEMENTARY CHARGE [C] */
#ifndef ECHARGE
#define ECHARGE     1.602176565E-19
#endif

/* UNIFIED ATOMIC MASS UNIT [kg] */
#ifndef UAMU
#define UAMU 1.660539040E-27
#endif

/* ELECTRON MASS [kg] */
#ifndef EMASS
#define EMASS       9.10938291E-31
#endif

/* MASS OF NEUTRON [kg] */
#ifndef NMASS
#define NMASS 1.674927351E-27
#endif

/* MASS OF POSITRON [kg] */
#ifndef PMASS
#define PMASS 1.672621777e-27
#endif

/* 1ATM [735.5592 TORR, 1.01325E5 Pa] */
#ifndef ATM
#define ATM 1.01325E5
#endif

/* STANDARD TEMPERATURE, IUPAC VERSION */
#ifndef ST_IUPAC
#define ST_IUPAC 273.15
#endif

/* 1TORR = 133.3224 [Pa = N/m^2] */
#ifndef TORR_MKS
#define TORR_MKS    133.3224
#endif

/* LOSCHMIDIT CONSTANT [1/m^3] */
#ifndef LOSCHMIDT
#define LOSCHMIDT   2.6867774E25 
#endif

/* AVOGADRO CONSTANT [#/mol] */
#ifndef NA
#define NA	6.02214129E23 
#endif

#ifndef CLIGHT
#define CLIGHT 2.99792458E8
#endif

} // End of namespace, SOENS

#endif
