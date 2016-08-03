/***************************************************************************************/
/*                                    Struct.h

    This code creates the many data storage structures used in other pieces of code.
    
    Based on code written by Coire Cadeau and modified by Sharon Morsink and 
    Abigail Stevens.
    
    Does Coire still get to claim copyright? I don't think we personally own our research 
    at the University of Alberta. Anyways, here it is:
    (C) Coire Cadeau, 2007; Source (C) Coire Cadeau 2007, all rights reserved.
*/
/***************************************************************************************/

#ifndef STRUCT_H
#define STRUCT_H

#include <exception>
#include <vector>
#include <float.h>

#define NN 20            // lookup table for bending angle (deflection angle) calculation
#define MAX_NUMBINS 512  // REMEMBER TO CHANGE THIS IN CHI.H AS WELL!! how many time bins the light curve is cut up into
#define NCURVES 100        // REMEMBER TO CHANGE THIS IN CHI.H AS WELL!! number of different light curves that it will calculate


struct Parameters {      // local bit of spot information
  double theta;          // Angle between the NS spin axis and the centre of the spot; in radians
  double phi_0;          // Azimuthal angular location of the centre of the spot; in radians
  double dS;             // Area of grid bit
  double incl;           // Inclination angle (between NS spin axis and observer's line of sight); in radians
  double aniso;          // Anisotropy
  double Gamma;          // Angle between true normal to surface and radial vector
  double Gamma1;         // Not related to above gamma; spectral index
  double Gamma2;         // Not related to above gamma; spectral index
  double Gamma3;         // Not related to above gamma; spectral index
  double temperature;    // Temperature of the spot, in the spot's frame; in keV
  double mass;           // Mass of the star; unitless in here
  double radius;         // Radius of the spot; unitless in here
  double req;            // Radius at the equator
  double rspot;          // Radius at the spot
  double mass_over_r;    // Dimensionless mass divided by radius
  double omega;          // Spin frequency of the star; unitless in here
  double cosgamma;       // Gamma is angle between true normal to surface and radial vector; probably in radians
  double bbrat;          // Ratio of blackbody-to-comptonization of spectrum
  double ts;             // Time shift between 0 and 1 (beginning of period, end of period)
  double E_band_lower_1; // Lower bound of energy band for flux calculation; in keV
  double E_band_upper_1; // Upper bound of energy band for flux calculation; in keV
  double E_band_lower_2; // Lower bound of energy band for flux calculation; in keV
  double E_band_upper_2; // Upper bound of energy band for flux calculation; in keV
  double distance;       // Distance from earth to NS; in meters
  double rsc;            // Scattering radius; adds scattering junk; in meters
  double Isc;            // Scattering intensity; adds scattering junk
  double bmodel;         // who knows?
  double E0; // NICER
  double E1; // NICER
  double E2; // NICER
  double DeltaE; // NICER 
};

struct Flags {
	double shift_t;               // for shifting of time, to match light curve phases
	bool infile_is_set;           // if an input file has been set
	bool ignore_time_delays;      // if we should ignore time delays
	unsigned int spectral_model;  // stating which model we're using -- definitions of models given elsewhere
	unsigned int beaming_model;   // stating which model we're using -- definitions of models given elsewhere
};


class Defl {
	public:                       // allocates the memory for the lookup table -- not evenly spaced
	double psi_b[3*NN+1];         // table where given psi, look up b
	double b_psi[3*NN+1];         // table where given b, look up psi
	double psi_max;               // largest possible value of psi
	double b_max;                 // largest possible value of b
	//class OblDeflectionTOA defltoa;
};

class LightCurve {                     // Stores all the data about the light curve!
	public:
	double t[MAX_NUMBINS];                 // one-dimensional array that hold the value of time of emission; normalized between 0 and 1
	double f[NCURVES][MAX_NUMBINS];        // two-dimensional array of fluxes (one one-dimensional array for each energy curve)
	bool visible[MAX_NUMBINS];             // is the spot visible at that point
	double t_o[MAX_NUMBINS];               // the time in the observer's frame; takes into account the light travel time
	double cosbeta[MAX_NUMBINS];           // as seen in MLCB17 (cos of zenith angle, between the norm vector and initial photon direction)
	double eta[MAX_NUMBINS];               // doppler shift factor; MLCB33
	double psi[MAX_NUMBINS];               // bending angle; MLCB15
	double R_dpsi_db[MAX_NUMBINS];         // derivative with respect to b of MLCB20 times the radius
	double b[MAX_NUMBINS];                 // impact parameter; defined in dimensionless units
	double dcosalpha_dcospsi[MAX_NUMBINS]; // appears in MLCB30
	double dOmega_s[MAX_NUMBINS];          // solid angle weight factor; MLCB30
	struct Parameters para;                // parameters from above; para is like i, Parameters is like Integer
	struct Flags flags;                    // flags from above
	class Defl defl;                       // deflection from above
	unsigned int numbins;                  // Number of time or phase bins for one spin period; Also the number of flux data points
	unsigned int numbands;
	bool eclipse;                          // True if an eclipse occurs
	bool ingoing;                          // True if one or more photons are ingoing
	bool problem;                          // True if a problem occurs
	double maxFlux[NCURVES];               // true (continuous) maximum flux values for each light curve
	double minFlux[NCURVES];               // true (continuous) minimum flux values for each light curve
	double pulseFraction[NCURVES];         // Pulse fraction of the light curve
	double norm[NCURVES];                  // The average flux value of a light curve, used to normalize a light curve to 1
	double asym[NCURVES];                  // Asymmetry between the rise and fall times for the light curve. =0 is rise=fall
	unsigned int count;                    // for outputting command line args in Chisquare, chi.cpp
};

// Need to have t, f, and err defined as pointers to arrays in this way:
/*for (unsigned int y(0); y < NCURVES; y++) {
 	obsdata.t = new double[numbins];
 	obsdata.f[y] = new double[numbins];
 	obsdata.err[y] = new double[numbins];
 }*/

class DataStruct {             // if reading in data, this would be the experimental data
	public:
	double *t;                       // time
	double *f[NCURVES];              // flux; an array of double pointers
	double *err[NCURVES];            // error bars
	double chisquare;                // chi squared
	double shift;                    // if we need it; unused
	unsigned int numbins;            // Number of time or phase bins for one spin period; Also the number of flux data points
};


// Shouldn't be needed anymore, but leaving in just in case
/*
class DataStruct {             // if reading in data, this would be the experimental data; needed for reading in in Spot
	public:
	double t[MAX_NUMBINS];             // time
	double f[NCURVES][MAX_NUMBINS];    // flux
	double err[NCURVES][MAX_NUMBINS];  // error bars
	double chisquare;                  // chi squared
	double shift;                      // if we need it; unused
	unsigned int numbins;              // Number of time or phase bins for one spin period; Also the number of flux data points
};
*/

#endif // STRUCT_H
