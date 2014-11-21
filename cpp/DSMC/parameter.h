#if !defined(__PARAMETER_H)
#define __PARAMETER_H

// Boltzmann constant.
#define BOLTZ			1.380622e-23
#define PI			3.141593


// Argon gas properties.
#define DIAMETER		4.17e-10
#define REF_TEMP		273.
#define VIS_TEMP_INDEX		0.81
#define MASS			6.63e-26


// CROSS_SECTION
#define MIX_CROSS_SECTION 0.25*PI*((2*DIAMETER)*(2*DIAMETER))
// REFERENCE_TEMPERATURE
#define MIX_RT  0.5*(2*REF_TEMP)
// VISCOSITY_TEMPERATURE_INDEX
#define MIX_VTI 0.5*(2*VIS_TEMP_INDEX)
// REDUCE_MASS
#define MIX_REDUCE_MASS (MASS*MASS)/(MASS+MASS)


#endif
