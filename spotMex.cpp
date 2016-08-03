/***************************************************************************************/
/*                                   SpotMex.cpp

    This code produces a pulse profile once a set of parameter describing the star, 
    spectrum, and hot spot have been inputed.
    
    This is the matlab executable version of Spot.cpp, by Abigail Stevens (2012/2013)
    
*/
/***************************************************************************************/

// INCLUDE ALL THE THINGS! 
// If you do not get this reference, see http://hyperboleandahalf.blogspot.ca/2010/06/this-is-why-ill-never-be-adult.html
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <exception>
#include <vector>
#include <string>
#include <string.h>
#include <stdio.h>
#include "OblDeflectionTOA.h"
#include "Chi.h"
#include "PolyOblModelNHQS.h"
#include "PolyOblModelCFLQS.h"
#include "SphericalOblModel.h"
#include "OblModelBase.h"
#include "Units.h"
#include "Exception.h"
#include "Struct.h"
#include "time.h"
#include "io64.h"
#include "mex.h"

// MAIN

void mexFunction ( int numOutputs, mxArray *theOutput[], int numInputs, const mxArray *theInput[] ) {
//theOutput[] is an array containing Null pointers.  We need to allocate the memory when returning the results
//numOutputs is the expected number of output pointers, i.e., the number of variables to be returned (unlike C, MatLab can have more than one)
//theInput[] is an array containing pointers to data of mxArray type
//numInputs is the number of input pointers, i.e., the number of variables passed to the function

	/*********************************************/
    /* VARIABLE DECLARATIONS AND INITIALIZATIONS */
    /*********************************************/

    double incl_1(90.0),               // Inclination angle of the observer, in degrees
    	   incl_2(90.0),               // PI - incl_2; needed for computing flux from second hot spot, since cannot have a theta greater than 
           theta_1(90.0),              // Emission angle (latitude) of the first upper spot, in degrees, down from spin pole
           theta_2(90.0),              // Emission angle (latitude) of the second lower spot, in degrees, up from spin pole (180 across from first spot)
           mass,                       // Mass of the star, in M_sun
           rspot,                      // Radius of the star at the spot, in km
           omega,                      // Frequency of the spin of the star, in Hz
           req,                        // Radius of the star at the equator, in km
           bbrat(1.0),                 // Ratio of blackbody to Compton scattering effects, unitless
           ts(0.0),                    // Phase shift or time off-set from data; Used in chi^2 calculation
           spot_temperature(0.0),      // Inner temperature of the spot, in the star's frame, in keV
           phi_0_1,                    // Equatorial azimuth at the center of the piece of the spot that we're looking at
           phi_0_2,                    // Equatorial azimuth at the center of the piece of the second hot spot that we're looking at
           theta_0_1,                  // Latitude at the center of the piece of the spot we're looking at
           theta_0_2,                  // Latitude at the center of the piece of the second hot spot that we're looking at
           rho(0.0),                   // Angular radius of the inner bullseye part of the spot, in degrees (converted to radians)
           dphi(1.0),                  // Each chunk of azimuthal angle projected onto equator, when broken up into the bins (see numphi)
           dtheta(1.0),                // Each chunk of latitudinal angle, when broken up into the bins (see numtheta)
           phi_edge_1(0.0),            // Equatorial azimuth at the edge of the spot at some latitude theta_0_1
           phi_edge_2(0.0),            // Equatorial azimuth at the edge of the second spot at some latitude theta_0_2
           theta_edge_1(0.0),          // Latitudinal edge of upper hot spot at -rho
           theta_edge_2(0.0),          // Latitudinal edge of lower (second) hot spot at -rho
           aniso(0.586),               // Anisotropy parameter
           Gamma1(2.0),                // currently unused variable; for spectral model #2 I think?
           Gamma2(2.0),                // ask Sharon
           Gamma3(2.0),                // ask Sharon
           mu_1(1.0),                  // = cos(theta_1), unitless
           mu_2(1.0),                  // = cos(theta_2), unitless
           cosgamma,                   // Cos of the angle between the radial vector and the vector normal to the surface; defined in equation 13, MLCB
           Flux[NCURVES][MAX_NUMBINS], // Array of fluxes; Each curve gets its own vector of fluxes based on photons in bins.
           trueSurfArea(0.0),          // The true geometric surface area of the spot
           E_band_lower_1(2.0),        // Lower bound of first energy band to calculate flux over, in keV.
           E_band_upper_1(3.0),        // Upper bound of first energy band to calculate flux over, in keV.
           E_band_lower_2(5.0),        // Lower bound of second energy band to calculate flux over, in keV.
           E_band_upper_2(6.0),        // Upper bound of second energy band to calculate flux over, in keV.
           //T_mesh[30][30],           // Temperature mesh over the spot; same mesh as theta and phi bins, assuming square mesh
           distance(3.0857e22),        // Distance from earth to the NS, in meters; default is 10kpc
           chisquared(0.0),            // chi^2 value for how well the computed light curve fits the data set
           *curveOut,
           *chiOut;
    
    unsigned int NS_model(3),          // Specifies oblateness (option 3 is spherical)
                 spectral_model(0),    // Spectral model choice (initialized to blackbody)
                 beaming_model(0),     // Beaming model choice (initialized to isotropic)
                 numbins(MAX_NUMBINS), // Number of time or phase bins for one spin period; Also the number of flux data points
                 numphi(1),            // Number of azimuthal (projected) angular bins per spot
                 numtheta(1),          // Number of latitudinal angular bins per spot
                 gray_flag(0);         // graybody (1) or isotropic (0)
	
	bool ignore_time_delays(false); // true if we're ignoring time delays, false if we are not ignoring time delays. default false.
    bool negative_theta(false);

/***************************************************************************************/
	/* Need to change filenameheader here below: */
  //  char filenameheader[]="RunAstar";
    // also need to change this in init.m and outputFerret.m
/***************************************************************************************/
	
    // Create LightCurve data structure
    class LightCurve curve, normcurve;  // variables curve and normalized curve, of type LightCurve
    class DataStruct obsdata;           // observational data as read in from a file
    
    // Setting up the output parameters
    int dimSize[2];
    dimSize[0] = 1;
    dimSize[1] = 1;

    // output needs to be an array so set it up as a [1,1] array holding one value
    theOutput[0] = mxCreateNumericArray(2, dimSize, mxDOUBLE_CLASS, mxREAL);

    chiOut = mxGetPr(theOutput[0]);
    
    // the following makes the current date a string called 'todaysdate'
   /* time_t tim;
	struct tm *ptr;
	char todaysdate[80];
	time( &tim );
	ptr = localtime( &tim );
	strftime(todaysdate,80,"%d-%b-%Y",ptr);*/
	


	/********************************************************/
    /* READ INFORMATION PASSED BY MEX FUNCTION              */
    /* MAKE SURE CORRECT NUMBER OF PARAMETERS ARE PASSED IN */
    /********************************************************/
   
	if ( numInputs != 27 ) { // need to go through and figure out how many inputs we should have
		//mexPrintf("numInputs = %d, expecting 26", numInputs);
		mexErrMsgTxt("Wrong number of inputs.");
	}
	
	numbins = mxGetScalar(theInput[0]); // int
	obsdata.t = mxGetPr(theInput[1]); // array of double
	NS_model = mxGetScalar(theInput[2]); // int
	mass = mxGetScalar(theInput[3]); // double
	rspot = mxGetScalar(theInput[4]); // rspot; double
	omega = mxGetScalar(theInput[5]); // double
	incl_1 = mxGetScalar(theInput[6]); // double
	theta_1 = mxGetScalar(theInput[7]); // double
	rho = mxGetScalar(theInput[8]); // double 
    numtheta = mxGetScalar(theInput[9]); //int
    spot_temperature = mxGetScalar(theInput[10]); // double
	E_band_lower_1 = mxGetScalar(theInput[11]); // double 
	E_band_upper_1 = mxGetScalar(theInput[12]); // double
	E_band_lower_2 = mxGetScalar(theInput[13]); // double 
	E_band_upper_2 = mxGetScalar(theInput[14]); // double 
	aniso = mxGetScalar(theInput[15]); // double
	ts = mxGetScalar(theInput[16]); // double
	bbrat = mxGetScalar(theInput[17]); // double
	gray_flag = mxGetScalar(theInput[18]); // int
	Gamma1 = mxGetScalar(theInput[19]); // double
	Gamma2 = mxGetScalar(theInput[20]); // double
	Gamma3 = mxGetScalar(theInput[21]); // double
	distance = mxGetScalar(theInput[22]); // double
	obsdata.f[1] = mxGetPr(theInput[23]); // array of double
    obsdata.f[2] = mxGetPr(theInput[24]); // array of double
    obsdata.err[1] = mxGetPr(theInput[25]); // array of double
    obsdata.err[2] = mxGetPr(theInput[26]); // array of double
  	//numtheta = 1;
  

	
    /******************************************/
    /* SENSIBILITY CHECKS ON INPUT PARAMETERS */
    /******************************************/
    
  //  if ( numtheta > 30 || numtheta < 1 ) {
  //      mexErrMsgTxt("Illegal number of theta bins. Must be between 1 and 30, inclusive. Exiting.");
  //  }
    if ( spot_temperature < 0.0 || rho < 0.0 || mass < 0.0 || rspot < 0.0 || omega < 0.0 ) {
        mexErrMsgTxt("Cannot have a negative spot temperature, spot angular radius, NS mass, NS radius, or NS spin frequency. Exiting." );
    }
  //  if ( rho == 0.0 && numtheta != 1 ) {
  //      std::cout << "\nWarning: Setting theta bin to 1 for a trivially sized spot.\n" << std::endl;
  //      numtheta = 1;
  //  }
    if ( E_band_lower_1 >= E_band_upper_1 ) {
    	mexErrMsgTxt("Must have E_band_lower_1 < E_band_upper_1 (both in keV). Exiting." );
    }
    if ( E_band_lower_2 >= E_band_upper_2 ) {
    	mexErrMsgTxt("Must have E_band_lower_2 < E_band_upper_2 (both in keV). Exiting." );
    }
    if ( numbins > MAX_NUMBINS || numbins <= 0 ) {
    	mexErrMsgTxt("Illegal number of phase bins. Must be between 1 and MAX_NUMBINS, inclusive. Exiting.");
    }
   	//std::cout << "Pre unit conversions." << std::endl;
 
    /*****************************************************/
    /* UNIT CONVERSIONS -- MAKE EVERYTHING DIMENSIONLESS */
    /*****************************************************/
    
    incl_1 *= (Units::PI / 180.0);  // radians
    incl_2 = incl_1; // radians
    theta_1 *= (Units::PI / 180.0); // radians
    theta_2 = theta_1; // radians
    rho *= (Units::PI / 180.0);  // radians
    mu_1 = cos( theta_1 ); // unitless
    mass = Units::cgs_to_nounits( mass*Units::MSUN, Units::MASS ); // from Msun to unitless
    rspot = Units::cgs_to_nounits( rspot*1.0e5, Units::LENGTH ); // from km to unitless
    omega = Units::cgs_to_nounits( 2.0*Units::PI*omega, Units::INVTIME ); // from Hz to unitless
    distance = Units::cgs_to_nounits( distance*100, Units::LENGTH ); // from m to unitless

	//std::cout << "Pre loading values into structure." << std::endl;

    /************************************************************/
    /* PASSING VALUES INTO THE STRUCTURE, ASSIGNING SOME VALUES */
    /************************************************************/   
	
	obsdata.numbins = numbins;
	curve.numbins = numbins;
    curve.para.mass = mass;
    curve.para.radius = rspot;
    curve.para.omega = omega;
    curve.para.theta = theta_1;
    curve.para.incl = incl_1;
    curve.para.aniso = aniso;
    curve.para.bbrat = bbrat;
    curve.para.Gamma1 = Gamma1;
    curve.para.Gamma2 = Gamma2;
    curve.para.Gamma3 = Gamma3;
    curve.para.temperature = spot_temperature;
    curve.para.ts = ts;
    curve.para.E_band_lower_1 = E_band_lower_1;
    curve.para.E_band_upper_1 = E_band_upper_1;
    curve.para.E_band_lower_2 = E_band_lower_2;
    curve.para.E_band_upper_2 = E_band_upper_2;
    curve.para.distance = distance;
  
    curve.flags.beaming_model = gray_flag;
	curve.flags.ignore_time_delays = ignore_time_delays;
	//curve.flags.spectral_model = spectral_model;
    
  	/***************************/
	/* START SETTING THINGS UP */
	/***************************/ 
	
	numphi = numtheta; // code currently only handles a square mesh over the hotspot
	
	//std::cout << "Pre setting up model." << std::endl;

	
    // Calculate the Equatorial Radius of the star.
    req = calcreq( omega, mass, theta_1, rspot );  // implementation of MLCB11

    /*********************************************************************************/
	/* Set up model describing the shape of the NS; oblate, funky quark, & spherical */
	/*********************************************************************************/
	OblModelBase* model;

    if ( NS_model == 1 ) {
        // Default model for oblateness
        model = new PolyOblModelNHQS( rspot, req,
		   		    PolyOblModelBase::zetaparam(mass,req),
				    PolyOblModelBase::epsparam(omega, mass, req) );
    }
    else if ( NS_model == 2 ) {
        // Alternative model for quark stars (not very different)
        model = new PolyOblModelCFLQS( rspot, req,
				     PolyOblModelBase::zetaparam(mass,rspot),
				     PolyOblModelBase::epsparam(omega, mass, rspot) );
    }
    else if ( NS_model == 3 ) {
        // Use a spherical star
        model = new SphericalOblModel( rspot );
    }
    else {
        mexErrMsgTxt("\nInvalid NS_model parameter. Exiting.\n");
    }

	//ts = obsdata.t[0]; // same thing is happening here. commenting this out so that ts does something.
	//obsdata.shift = obsdata.t[0];
	obsdata.shift = ts; // comment this out if you want to have ts = obsdata.t[0]
	
    // defltoa is a structure that "points" to routines in the file "OblDeflectionTOA.cpp"
    // used to compute deflection angles and times of arrivals 
    // defltoa is the deflection time of arrival
    OblDeflectionTOA* defltoa = new OblDeflectionTOA(model, mass); // defltoa is a pointer (of type OblDeclectionTOA) to a group of functions    
	
	cosgamma = model->cos_gamma(mu_1);
    curve.para.cosgamma = cosgamma;

	//std::cout << "Pre b lookup table." << std::endl;

    /**********************************************************/
	/* Compute maximum deflection for purely outgoing photons */
	/**********************************************************/
	
    double  b_mid;  // the value of b, the impact parameter, at 90% of b_max
    curve.defl.b_max =  defltoa->bmax_outgoing(rspot); // telling us the largest value of b
    curve.defl.psi_max = defltoa->psi_max_outgoing(curve.defl.b_max, rspot, &curve.problem); // telling us the largest value of psi

    /********************************************************************/
	/* COMPUTE b VS psi LOOKUP TABLE, GOOD FOR THE SPECIFIED M/R AND mu */
	/********************************************************************/
	
    b_mid = curve.defl.b_max * 0.9; // since we want to split up the table between coarse and fine spacing. 0 - 90% is large spacing, 90% - 100% is small spacing. b_mid is this value at 90%.
    curve.defl.b_psi[0] = 0.0; // definitions of the actual look-up table for b when psi = 0
    curve.defl.psi_b[0] = 0.0; // definitions of the actual look-up table for psi when b = 0
    
    for ( unsigned int i(1); i < NN+1; i++ ) { // compute table of b vs psi points
        curve.defl.b_psi[i] = b_mid * i / (NN * 1.0);
        curve.defl.psi_b[i] = defltoa->psi_outgoing(curve.defl.b_psi[i], rspot, curve.defl.b_max, curve.defl.psi_max, &curve.problem); // calculates the integral
    }
    
    // For arcane reasons, the table is not evenly spaced.
    for ( unsigned int i(NN+1); i < 3*NN; i++ ) {  // compute table of b vs psi points
        curve.defl.b_psi[i] = b_mid + (curve.defl.b_max - b_mid) / 2.0 * (i - NN) / (NN * 1.0); // spacing for the part where the points are closer together
        curve.defl.psi_b[i] = defltoa->psi_outgoing(curve.defl.b_psi[i], rspot, curve.defl.b_max, curve.defl.psi_max, &curve.problem); // referenced same as above
    }
    
    curve.defl.b_psi[3*NN] = curve.defl.b_max;   // maximums
    curve.defl.psi_b[3*NN] = curve.defl.psi_max;
    // Finished computing lookup table


	//std::cout << "Pre initializing time and flux." << std::endl;

	
    /****************************/
	/* Initialize time and flux */
	/****************************/
	
    for ( unsigned int i(0); i < numbins; i++ ) {
        curve.t[i] = i / (1.0 * numbins);// + ts;   // defining the time used in the lightcurves
        for ( unsigned int p(0); p < NCURVES; p++ ) {
            curve.f[p][i] = 0.0;                       // initializing flux to 0
            Flux[p][i] = 0.0;
        }
    } 

    /****************************************/
	/* Set up the star and calculate fluxes */
	/****************************************/
	
    if ( rho == 0.0 ) { // Default is an infinitesmal spot, defined to be 0 degrees in radius.
      	phi_0_1 = 0.0;
		theta_0_1 = theta_1;
	    curve.para.dS = trueSurfArea = 0.0; // dS is the area of the particular bin of spot we're looking at right now
   	}   
    else {   // for a nontrivially-sized spot
	    dtheta = 2.0 * rho / (numtheta * 1.0);    // center of spot at theta, rho is angular radius of spot; starts at theta_1-rho to theta_1+rho
       	theta_edge_1 = theta_1 - rho;
       	trueSurfArea = 2 * Units::PI * pow(rspot,2) * (1 - cos(rho));
    }
    
       	//std::cout << "Pre calculating spot." << std::endl;

    /************************************/
	/* LOCATION OF THE SPOT ON THE STAR */
	/************************************/
    	
    /*********************************************/
	/* Spot is TRIVIALLY SIZED ON GEOMETRIC POLE */
	/*********************************************/
		
	if ( theta_1 == 0.0 && rho == 0.0 ) {
		for ( unsigned int i(0); i < numbins; i++ ) {
	    	for ( unsigned int p(0); p < NCURVES; p++ )
		       	curve.f[p][i] = 0.0;
	    }
	    
		// Add curves
	    for (unsigned int i(0); i < numbins; i++) {
		    for (unsigned int p(0); p < NCURVES; p++) {
		        //std::cout << Flux[p][i] << ", " << curve.f[p][i] << std::endl;
		        Flux[p][i] += curve.f[p][i];
		        //if (p == NCURVES-1) std::cout << "Time = " << curve.t[i] << ", i = " << i << ", Flux = " << Flux[p][i] << std::endl;
		        //if( std::isnan(Flux[p][i]) || Flux[p][i] == 0) std::cout << "Flux is NAN or 0 at p="<<p<<", i="<<i << std::endl;
		    }
	    } 
	} // ending trivial spot on pole
	
	/*****************************************/
	/* Spot is SYMMETRIC OVER GEOMETRIC POLE */
	/*****************************************/
		
	else if ( theta_1 == 0.0 && rho != 0.0 ) {
   	 	// Looping through the mesh of the spot
		for ( unsigned int k(0); k < numtheta; k++ ) {
			theta_0_1 = theta_1 - rho + k * dtheta + 0.5 * dtheta;   // note: theta_0_1 changes depending on how we've cut up the spot
			curve.para.theta = theta_0_1;
			// don't need a phi_edge the way i'm doing the phi_0_1 calculation
			dphi = Units::PI / numphi;
			for ( unsigned int j(0); j < numphi; j++ ) {
				phi_0_1 = -Units::PI/2 + dphi * j + 0.5 * dphi;

				//std::cout << "\nk = " << k << ", j = " << j << std::endl;
	          /*std::cout << "\ntheta_0_1 = " << theta_0_1*180.0/Units::PI << " degrees" << std::endl;
	            std::cout << "dtheta = " << dtheta*180.0/Units::PI << " degrees" << std::endl;
	            std::cout << "phi_0_1 = " << phi_0_1*180.0/Units::PI << " degrees" << std::endl;
	            std::cout << "dphi = " << dphi*180.0/Units::PI << " degrees" << std::endl;
	            std::cout << "fabs(theta_0_1) = " << fabs(theta_0_1)*180.0/Units::PI << " degrees" << std::endl;
	            std::cout << "sin(fabs(theta_0_1)) = " << sin(fabs(theta_0_1)) << std::endl;
	            //std::cout << "r spot = " << rspot << std::endl;
				std::cout << "theta_edge_1 = " << theta_edge_1*180.0/Units::PI << " degrees" << std::endl;
			 	std::cout << std::endl;*/
				
				curve.para.dS = pow(rspot,2) * sin(fabs(theta_0_1)) * dtheta * dphi; // assigning partial dS here
	            if ( numtheta == 1 ) 
	            	curve.para.dS = trueSurfArea;
	            if ( NS_model == 1 || NS_model == 2)
	            	curve.para.dS /= curve.para.cosgamma;
	            curve.para.phi_0 = phi_0_1;
	            
	            curve = ComputeAngles(&curve, defltoa);  // Computing the parameters it needs to compute the light curve; defltoa has the routines for what we want to do
	        	curve = ComputeCurve(&curve);
			
				if (curve.para.temperature == 0.0) { 
	        		for (unsigned int i(0); i < numbins; i++) {
		        		for (unsigned int p(0); p < NCURVES; p++)
		           			curve.f[p][i] = 0.0;  // if temperature is 0, then there is no flux!
	        		}
	        	}
	        		
				// Add curves
       			for (unsigned int i(0); i < numbins; i++) {
        			for (unsigned int p(0); p < NCURVES; p++) {
        				//std::cout << Flux[p][i] << ", " << curve.f[p][i] << std::endl;
            			Flux[p][i] += curve.f[p][i];
	           			//if (p == NCURVES-1) std::cout << "Time = " << curve.t[i] << ", i = " << i << ", Flux = " << Flux[p][i] << std::endl;
		       			//if( std::isnan(Flux[p][i]) || Flux[p][i] == 0) std::cout << "Flux is NAN or 0 at p="<<p<<", i="<<i << std::endl;
		      		}
	        	} // end 'add curves'
			}
		}
	} // ending symmetric spot over pole
		
	// new
/*********************************************/
	/* SPOT IS ASYMMETRIC OVER GEOMETRIC POLE */
	/*********************************************/
	
	else if ( (theta_1 - rho) <= 0 ) {
	
		if ( numtheta == 1 ) {
		
			curve.para.theta = theta_1;
			curve.para.phi_0 = 0.0;
			curve.para.dS = trueSurfArea;
			if ( !T_mesh_in ) {
	           		T_mesh[0][0] = spot_temperature;
			}
			curve.para.temperature = T_mesh[0][0];

			curve = ComputeAngles(&curve, defltoa);  // Computing the parameters it needs to compute the light curve; defltoa has the routines for what we want to do
			curve = ComputeCurve(&curve);
	        	
			if ( curve.para.temperature == 0.0 ) { 
			  for ( unsigned int i(0); i < numbins; i++ ) {
			    for ( unsigned int p(0); p < NCURVES; p++ )
			      curve.f[p][i] = 0.0;
			  }
			}
			// Add curves, load into Flux array
			for ( unsigned int i(0); i < numbins; i++ ) {
			  for ( unsigned int p(0); p < NCURVES; p++ ) {
			    Flux[p][i] += curve.f[p][i];
			  }
			} // end Add curves
		} // ending if numtheta == 1
			
		else { // if numtheta != 1
		
		  if ( T_mesh_in ) {
		    std::cout << "WARNING: code can't handle a spot asymmetric over the pole with a temperature mesh." << std::endl;
		    spot_temperature = 2;
		  }
		  curve.para.temperature = spot_temperature;

		  // Looping through the mesh of the spot
		  for (unsigned int k(0); k < numtheta; k++) { 
		    negative_theta = false;
		    dtheta = 2.0 * rho / (1.0*numtheta);
		    theta_0_1 = theta_1 - rho + 0.5 * dtheta + k * dtheta;   // note: theta_0_1 changes depending on how we've cut up the spot
		    if (theta_0_1 < 0.0) {
		      negative_theta = true;
		      theta_0_1 = fabs(theta_0_1);
		    }
		    if (theta_0_1 == 0.0) 
		      theta_0_1 = 0.0 + DBL_EPSILON;

		    curve.para.theta = theta_0_1;   // passing this value into curve theta, so that we can use it in another routine; somewhat confusing names

		    double cos_phi_edge = (cos(rho) - cos(theta_1)*cos(theta_0_1))/(sin(theta_1)*sin(theta_0_1));
		    if ( cos_phi_edge > 1.0 ) { 
		      cos_phi_edge = 1.0;
		    }

		    if ( fabs( sin(theta_1) * sin(theta_0_1) ) > 0.0) { // checking for a divide by 0
		      phi_edge_1 = acos( cos_phi_edge );  
		      dphi = 2.0 * phi_edge_1 / ( numphi * 1.0 );
		    }
		    else {  // trying to divide by zero
		      mexErrMsgTxt("Tried to divide by zero. Likely sin(theta_1) or sin(theta_0_1) is 0. Exiting.");
		    }

		    for ( unsigned int j(0); j < numphi; j++ ) {   // looping through the phi divisions
		      phi_0_1 = -phi_edge_1 + 0.5 * dphi + j * dphi;   // phi_0_1 is at the center of the piece of the spot that we're looking at; 
		      //defined as distance from the y-axis as projected onto the equator

		      if ( negative_theta )
			phi_0_1 = phi_0_1 + Units::PI;
				
		      curve.para.phi_0 = phi_0_1;


			 	
		      curve.para.dS = pow(rspot,2) * sin(fabs(theta_0_1)) * dtheta * dphi;
		      if ( NS_model == 1 || NS_model == 2 )
			curve.para.dS /= curve.para.cosgamma;
		      // Need to multiply by R^2 here because of my if numtheta=1 statement, 
		      // which sets dS = true surface area
		      // alternative method for computing dS: dS=truesurfarea/pow(numtheta,2);

		     


		      curve = ComputeAngles(&curve, defltoa); 
		      curve = ComputeCurve(&curve);
	        	
		      if ( curve.para.temperature == 0.0 ) { 
			for ( unsigned int i(0); i < numbins; i++ ) {
			  for ( unsigned int p(0); p < NCURVES; p++ )
			    curve.f[p][i] = 0.0;
			}
		      }
	       		// Add curves, load into Flux array
		      for ( unsigned int i(0); i < numbins; i++ ) {
			for ( unsigned int p(0); p < NCURVES; p++ ) {
			  Flux[p][i] += curve.f[p][i];
			}
		      } // ending Add curves
		    } // end for-j-loop
		  } // closing for loop through theta divisions
		} // closing if numtheta != 1
	} // ending antisymmetric spot over pole





    	  




	/****************************************/
	/* Spot DOES NOT GO OVER GEOMETRIC POLE */
	/****************************************/
		
	else {
		//std::cout << "Within spot not over pole, pre computations and calls." << std::endl;

   	 	// Looping through the mesh of the spot
    	for ( unsigned int k(0); k < numtheta; k++ ) { // looping through the theta divisions
    		//std::cout << "Within spot not over pole, within k loop, pre computations and calls." << std::endl;

			theta_0_1 = theta_1 - rho + 0.5 * dtheta + k * dtheta;   // note: theta_0_1 changes depending on how we've cut up the spot
			curve.para.theta = theta_0_1;   // passing this value into curve theta, so that we can use it in another routine; somewhat confusing names
        
        	double cos_phi_edge = (cos(rho) - cos(theta_1)*cos(theta_0_1))/(sin(theta_1)*sin(theta_0_1));
	        if ( cos_phi_edge > 1.0 || cos_phi_edge < -1.0 ) 
	        	cos_phi_edge = 1.0;
        
	        if ( fabs( sin(theta_1) * sin(theta_0_1) ) > 0.0 ) { // checking for a divide by 0
	           	phi_edge_1 = acos( cos_phi_edge );   // value of phi (a.k.a. azimuth projected onto equatorial plane) at the edge of the circular spot at some latitude theta_0_1
		       	dphi = 2.0 * phi_edge_1 / ( numphi * 1.0 );
		    }
	        else {  // trying to divide by zero
    	       	mexErrMsgTxt("Tried to divide by zero in calculation of phi_edge_1. Likely, theta_0_1 = 0. Exiting.");
        	}

    		for ( unsigned int j(0); j < numphi; j++ ) {   // looping through the phi divisions
    			//std::cout << "Within spot not over pole, within j loop, pre computations and calls." << std::endl;

	        	phi_0_1 = -phi_edge_1 + 0.5 * dphi + j * dphi;   // phi_0_1 is at the center of the piece of the spot that we're looking at; defined as distance from the y-axis as projected onto the equator
	        	/*if ( only_second_spot ) 
	        		phi_0_1 = Units::PI - phi_edge_1 + 0.5 * dphi + j * dphi;
	            */		
	            /* std::cout << "\nk = " << k << ", j = " << j << std::endl;
            	std::cout << "\ntheta_0_1 = " << theta_0_1*180.0/Units::PI << " degrees" << std::endl;
            	std::cout << "inclination = " << incl_1*180.0/Units::PI << " degrees" << std::endl;
            	std::cout << "dtheta = " << dtheta*180.0/Units::PI << " degrees" << std::endl;
            	std::cout << "phi_0_1 = " << phi_0_1*180.0/Units::PI << " degrees" << std::endl;
            	std::cout << "dphi = " << dphi*180.0/Units::PI << " degrees" << std::endl;
            	std::cout << "fabs(theta_0_1) = " << fabs(theta_0_1)*180.0/Units::PI << " degrees" << std::endl;
            	std::cout << "sin(fabs(theta_0_1)) = " << sin(fabs(theta_0_1)) << std::endl;
            	//std::cout << "r spot = " << rspot << std::endl;
				std::cout << "theta_edge_1 = " << theta_edge_1*180.0/Units::PI << " degrees" << std::endl;
				std::cout << "phi_edge_1 = " << phi_edge_1*180.0/Units::PI << " degrees" << std::endl; 
		 		std::cout << std::endl;*/
          
            	curve.para.dS = pow(rspot,2) * sin(fabs(theta_0_1)) * dtheta * dphi;
            	if( numtheta == 1 ) 
            		curve.para.dS = trueSurfArea;
            	if ( NS_model == 1 || NS_model == 2 )
            		curve.para.dS /= curve.para.cosgamma;
            	//std::cout << "dS = " << curve.para.dS << std::endl;
            	curve.para.phi_0 = phi_0_1;
            	//std::cout << "Within spot not over pole, pre calls to Chi.cpp." << std::endl;

            	curve = ComputeAngles(&curve, defltoa);  // Computing the parameters it needs to compute the light curve; defltoa has the routines for what we want to do
        		curve = ComputeCurve(&curve);            // Compute Light Curve, for each separate mesh bit
    			//std::cout << "dOmega = " << curve.dOmega_s[0] << std::endl;
                //std::cout << "curve.f[3][2] = " << curve.f[3][2] << std::endl;
	    		
                /*for (unsigned int i(0); i < numbins; i++) {
		    		for (unsigned int p(0); p < NCURVES; p++)
		        		if( std::isnan(curve.f[p][i]) || curve.f[p][i] == 0) std::cout << "curve.f[p="<<p<<"][i="<<i<<"] = " << curve.f[p][i] << std::endl;
	        	}*/
	        	if ( curve.para.temperature == 0.0 ) { 
	        		for ( unsigned int i(0); i < numbins; i++ ) {
		        		for ( unsigned int p(0); p < NCURVES; p++ )
		           			curve.f[p][i] = 0.0;
	        		}
	        	}
	       		// Add curves
	       		//std::cout << "Within spot not over pole, pre 'add curves'." << std::endl;

	        	for ( unsigned int i(0); i < numbins; i++ ) {
		        	for ( unsigned int p(0); p < NCURVES; p++ ) {
		        		//std::cout << Flux[p][i] << ", " << curve.f[p][i] << std::endl;
		            	Flux[p][i] += curve.f[p][i];
		            	//if (p == NCURVES-1) std::cout << "Time = " << curve.t[i] << ", i = " << i << ", Flux = " << Flux[p][i] << std::endl;
		            	//if( std::isnan(Flux[p][i]) || Flux[p][i] == 0) std::cout << "Flux is NAN or 0 at p="<<p<<", i="<<i << std::endl;
		            }
	        	} // close add curves   
	       	} // closing for loop through phi divisions
	    } // closing for loop through theta divisions
    
    } // closing spot does not go over the pole or is antisymmetric over the pole
    
  /*  for ( unsigned int i(0); i < numbins; i++ ) {
     	for ( unsigned int p(0); p < NCURVES; p++ ) {
     		normcurve.f[p][i] = 0.0;
	        //if( std::isnan(normcurve.f[p][i]) || normcurve.f[p][i] == 0) std::cout << "normcurve.f is NAN or 0 at p="<<p<<", i="<<i << std::endl;
		}
    }  */
    
	/***********************************************************************************/
	/* SECOND HOT SPOT -- Note that this cannot handle going over the geometric pole,  */
	/*                    but doing so would be quite easy. This has precisely the     */
	/*                    same parameters as the first hot spot, just on the antipodal */
	/*                    magnetic pole.                                               */
	/***********************************************************************************/
/*
    if ( two_spots ) {
       	incl_2 = Units::PI - incl_1;
    	curve.para.incl = incl_2;
    	
    	if ( rho == 0.0 ) { // Default is an infinitesmal spot, defined to be 0 degrees in radius.
        	phi_0_2 = Units::PI;
        	theta_0_2 = theta_2;
        	curve.para.dS = trueSurfArea = 0.0; // dS is the area of the particular bin of spot we're looking at right now
        }  
   		else {   // for a nontrivially-sized spot
       		dtheta = 2.0 * rho / (numtheta * 1.0);    // center of spot at theta, rho is angular radius of spot; starts at theta_1-rho to theta_1+rho
       		theta_edge_2 = theta_2 - rho;
       		trueSurfArea = 2 * Units::PI * pow(rspot,2) * (1 - cos(rho));
   		}
    
       	// Looping through the mesh of the second spot
   		for ( unsigned int k(0); k < numtheta; k++ ) { // looping through the theta divisions
			theta_0_2 = theta_2 - rho + 0.5 * dtheta + k * dtheta;   // note: theta_0_2 changes depending on how we've cut up the spot
       		//std::cout << "Theta 2 = " << theta_0_2 * 180 / Units::PI << std::endl; // so that it prints to the screen in degrees for us!
       		curve.para.theta = theta_0_2;   // passing this value into curve theta, so that we can use it in another routine; somewhat confusing names
       		//std::cout //<< "\n " k = " << k << ", theta_0_1 = " << theta_0_1 * 180.0 / Units::PI << std::endl;
        
       		double cos_phi_edge = (cos(rho) - cos(theta_2)*cos(theta_0_2))/(sin(theta_2)*sin(theta_0_2));
    		if ( cos_phi_edge > 1.0 ) 
    			cos_phi_edge = 1.0;
        
        	if ( fabs( sin(theta_2) * sin(theta_0_2) ) > 0.0 ) { // checking for a divide by 0
           		phi_edge_2 = acos ( cos_phi_edge );   // value of phi (a.k.a. azimuth projected onto equatorial plane) at the edge of the circular spot at some latitude theta_0_2
	            dphi = 2.0 * phi_edge_2 / ( numphi * 1.0 ); // want to keep the same dphi as before
	        }
        	else {  // trying to divide by zero
            	mexErrMsgTxt("Tried to divide by zero in calculation of phi_edge_2. \n Check values of sin(theta_2) and sin(theta_0_2). Exiting.");
        	}     
        
       		for ( unsigned int j(0); j < numphi; j++ ) {   // looping through the phi divisions
           		phi_0_2 = Units::PI - phi_edge_2 + 0.5 * dphi + j * dphi;
 */          		
           		/*std::cout << "\nk = " << k << ", j = " << j << std::endl;
           		std::cout << "inclination = " << incl_2 * 180 / Units::PI << std::endl;
           		std::cout << "theta_0_2 = " << theta_0_2*180/Units::PI << std::endl;
           		std::cout << "dtheta = " << dtheta*180/Units::PI << std::endl;
				std::cout << "phi_0_2 = " << phi_0_2*180/Units::PI << std::endl;
           		std::cout << "dphi = " << dphi*180/Units::PI << std::endl;
           		//std::cout << "fabs(theta_0_2) = " << fabs(theta_0_2) *180./Units::PI << std::endl;
           		//std::cout << "sin(fabs(theta_0_2)) = " << sin(fabs(theta_0_2)) << std::endl;
        		//std::cout << "r spot = " << rspot << std::endl;
				std::cout << "theta_edge_2 = " << theta_edge_2*180/Units::PI << std::endl;
				std::cout << "phi_edge_2 = " << phi_edge_2*180/Units::PI << std::endl;
				std::cout << std::endl;*/
  /*          
           		//curve.para.temperature = T_mesh[k][j];
           		curve.para.temperature = spot_temperature;
        		curve.para.phi_0 = phi_0_2;
        		curve.para.dS = pow(rspot,2) * sin(fabs(theta_0_2)) * dtheta * dphi;
        		if( numtheta == 1 ) 
            		curve.para.dS = trueSurfArea;
            	if ( NS_model == 1 || NS_model == 2 )
            		curve.para.dS /= curve.para.cosgamma;
        		//std::cout << "dS = " << curve.para.dS << std::endl;
        		//std::cout << sin(theta_1) * dtheta * dzeta << std::endl;
				curve = ComputeAngles(&curve, defltoa);  // Computing the parameters it needs to compute the light curve; defltoa has the routines for what we want to do
        		curve = ComputeCurve(&curve);            // Compute Light Curve, for each separate mesh bit
        	
       			if ( curve.para.temperature == 0.0 ) { 
       				for ( unsigned int i(0); i < numbins; i++ ) {
        				for ( unsigned int p(0); p < NCURVES; p++ )
	        				curve.f[p][i] += 0.0;
        			}
        		}
        		// Add curves
				for ( unsigned int i(0); i < numbins; i++ ) {
	    			for ( unsigned int p(0); p < NCURVES; p++ ) {
	        			//std::cout << Flux[p][i] << ", " << curve.f[p][i] << std::endl;
	            		Flux[p][i] += curve.f[p][i];
	            		//if (p == NCURVES-1) std::cout << "Time = " << curve.t[i] << ", i = " << i << ", Flux = " << Flux[p][i] << std::endl;
	            		//if( std::isnan(Flux[p][i]) || Flux[p][i] == 0) std::cout << "Flux is NAN or 0 at p="<<p<<", i="<<i << std::endl;
	            	}
       			}
       		} // closing for loop through phi divisions
		} // closing for loop through theta divisions
    } // closing if two spots

*/
 	//std::cout << "mex1 curve.f[3][28] = " << curve.f[3][28] << std::endl;

 	/*******************************/
	/* NORMALIZING THE FLUXES TO 1 */
	/*******************************/
	//std::cout << "Pre normalizing." << std::endl;
    //std::cout << "Flux[p=3][i=2] = " << Flux[3][2] << std::endl;
 	//if ( normalize_flux ) {   // Normalize to 1
    	normcurve = Normalize( Flux, numbins );
     	// fun little way around declaring Normalize as returning a matrix! so that we can still print out Flux
     	for ( unsigned int i(0); i < numbins; i++ ) {
     		for ( unsigned int p(0); p < NCURVES; p++ ) {
     			Flux[p][i] = normcurve.f[p][i];
	            //if( std::isnan(normcurve.f[p][i]) || normcurve.f[p][i] == 0) std::cout << "normcurve.f[p="<<p<<"][i="<<i<<"] = " << normcurve.f[p][i] << std::endl;
     		}
     	}  
    //}
    //std::cout << "Pre loading Flux back into curve struct." << std::endl;

    // Loading the total Flux back into the structure curve, so the chisquared is done properly
    for ( unsigned int i(0); i < numbins; i++ ) {
	    for ( unsigned int p(0); p < NCURVES; p++ ) {
	    	//std::cout << Flux[p][i] << ", " << curve.f[p][i] << std::endl;
	        curve.f[p][i] = Flux[p][i];
            //if ( p == 3 ) std::cout << "Flux[p="<<p<<"][i="<<i<<"] = " << Flux[p][i] << std::endl;
	        //if (p == NCURVES-1) std::cout << "Time = " << curve.t[i] << ", i = " << i << ", Flux = " << Flux[p][i] << std::endl;
	    	//if( std::isnan(Flux[p][i]) || Flux[p][i] == 0) std::cout << "Flux is NAN or 0 at p="<<p<<", i="<<i << std::endl;
	    }
    }
 	//std::cout << "mex_norm curve.f[3][28] = " << curve.f[3][28] << std::endl;
	/*******************************/
	/* Calculating pulse fractions */
	/*******************************/
	
//    double avgPulseFraction(0.0); // average of pulse fractions across all energy bands
 /*   double sumMaxFlux(0.0), sumMinFlux(0.0);
    double overallPulseFraction(0.0);
    for ( unsigned int j(0); j < NCURVES; j++ ) {
    	sumMaxFlux += curve.maxFlux[j];
    	sumMinFlux += curve.minFlux[j];
        avgPulseFraction += curve.pulseFraction[j];
        //std::cout << "Maximum = " << curve.maxFlux[j] << std::endl;
        //std::cout << "Minimum = " << curve.minFlux[j] << std::endl; 
        std::cout << "Pulse Fraction ["<<j<<"] = " << curve.pulseFraction[j] << std::endl; 
    }
	overallPulseFraction = (sumMaxFlux - sumMinFlux) / (sumMaxFlux + sumMinFlux);
    avgPulseFraction /= (NCURVES);
    //std::cout << "\nAveraged Pulse Fraction = " << avgPulseFraction << std::endl;
    //std::cout << "\nPulse Fraction = " << overallPulseFraction << std::endl;
	*/
	//std::cout << "Pre chisquared computation." << std::endl;

	chisquared = ChiSquare( &obsdata, &curve );
	
	if (std::isnan(chisquared)) {
		//mexPrintf("Chisquared is nan! Setting it to 1,000,000.\n");
		chisquared = 10000000.0;
	}
	else {
		//mexPrintf("X^2 = %f\n\n", chisquared);
	}
	
	if (ts > 1.0) ts -= 1.0; // to make sure ts is in [0:1], for printing out
    //std::cout << "Pre printing to everything." << std::endl;

	/*****************************************************************/
    /* DUMPING DATA INTO MATLAB OUTPUT                               */
    /* For filewriting, need to use fopen, fprintf, and fclose (etc) */
    /*****************************************************************/
    //std::cout << "Pre output." << std::endl;

    chiOut[0] = chisquared; // saved this to theOutput[0] at top just after declarations
    //std::cout << "Output 1." << std::endl;
    dimSize[1] = 3; // number of columns (1: time, 2: flux in 1st energy band, 3: flux in 2nd energy band)
    //std::cout << "Output 2." << std::endl;

    dimSize[0] = (int)numbins; // number of rows ( = numbins)
    //std::cout << "Output 3." << std::endl;

    theOutput[1] = mxCreateNumericArray(2, dimSize, mxDOUBLE_CLASS, mxREAL); // formatting/setting up the output
    //std::cout << "Output 4." << std::endl;

    curveOut = mxGetPr(theOutput[1]); // matlab and mex love pointers
    //std::cout << "Output 5." << std::endl;

    /*FILE * printCurveOut;  // printing curveOut values to a file
    printCurveOut = fopen("run_data/curveOut_flux.txt", "w+"); // writing to curveOut_flux.txt
    
    if ( printCurveOut == NULL) {
        mexErrMsgTxt("\n printCurveOut failed. Exiting.\n");
    }
    
    FILE * printMexFlux;  // printing Flux values to a file
    printMexFlux = fopen("run_data/mex_flux.txt", "w+"); // writing to mex_flux.txt
    
    if ( printMexFlux == NULL ) {
        mexErrMsgTxt("\n printMexFlux failed. Exiting.\n");
    }
    
    fprintf(printCurveOut,"# What gets saved in curveOut. Should match what is put into mex_flux.txt \n");
	fprintf(printCurveOut,"# Column 1: time \n# Column 2: Flux, first energy band \n# Column 3: Flux, second energy band \n");
    fprintf(printCurveOut,"# \n");
    
    fprintf(printMexFlux,"# The flux calculated in spotMex.cpp. Should match what is put into curveOut_flux.txt \n");
    fprintf(printMexFlux,"# Column 1: time \n# Column 2: Flux, first energy band \n# Column 3: Flux, second energy band \n");
    fprintf(printMexFlux,"# \n");
    */
    for ( unsigned int i(0); i < numbins; i++ ) {

    	for ( int column(0); column < dimSize[1]; column++ ) {
    		    //std::cout << "Output 6: i="<<i<<", column="<< column << std::endl;

    		if ( column == 0 ) { // time goes in first column
    			curveOut[i + (int)numbins*column] = curve.t[i];
    			//fprintf(printCurveOut, "%f \t",curveOut[i + (int)numbins*column]);
    			//fprintf(printMexFlux, "%f \t",curve.t[i]);
    		}
    		else if ( column == 1 ) { // first energy band flux goes in second column
    			if ( std::isnan(Flux[NCURVES-2][i]) ){
    				curveOut[i + (int)numbins*column] = 0.0;
    				//mexPrintf("ERROR! FLUX BAND 1 = NAN. Setting it to 0.\n");
    			}
    			else {
    				curveOut[i + (int)numbins*column] = Flux[NCURVES-2][i]; // first energy band
    			}
    			//fprintf(printCurveOut, "%f \t",curveOut[i + (int)numbins*column]);
    			//fprintf(printMexFlux, "%f \t",curve.f[NCURVES-2][i]);
    		}
    		else { // second energy band flux goes in third column
    			//std::cout << "i + (int)numbins*column = " << i + (int)numbins*column << std::endl;
    			if ( std::isnan(Flux[NCURVES-1][i]) ){
    				curveOut[i + (int)numbins*column] = 0.0;
    				//mexPrintf("ERROR! FLUX BAND 2 = NAN. Setting it to 0.\n");
    			}
    			else {
    				curveOut[i + (int)numbins*column] = Flux[NCURVES-1][i]; // second energy band
    			}
    			//std::cout << "Printed to curveOut correctly." << std::endl;
    			//fprintf(printCurveOut,"%f \n",curveOut[i + (int)numbins*column]);
    			//fprintf(printMexFlux,"%f \n",curve.f[NCURVES-1][i]);
    		}
    	}
    }
    //std::cout << "Output done." << std::endl;
    //fclose(printCurveOut);
	//fclose(printMexFlux);

   /* out << "# Photon Flux for isotropic blackbody spot on a spherical star, D = 10 kpc. \n"
        << "# R = " << Units::nounits_to_cgs(rspot, Units::LENGTH )*1.0e-5 << " km; "
        << "# M = " << Units::nounits_to_cgs(mass, Units::MASS)/Units::MSUN << " Msun; "
        << "# Spin = " << Units::nounits_to_cgs(omega, Units::INVTIME)/(2.0*Units::PI) << " Hz \n"
        << "# Gravitational Redshift 1+z = " << 1.0/sqrt(1.0 - 2*mass/rspot) << "\n"
        << "# Inclination Angle = " << incl_1 * 180.0/Units::PI << " degrees \n"
        << "# Spot 1 at Angle = " << theta_1 * 180.0/Units::PI << " degrees \n"
        //<< "# Spot 2 at Angle = " << theta_2 * 180.0/Units::PI << " degrees \n"
        << "# Angular Radius of Spot = " << rho * 180.0/Units::PI << " degrees \n"
        << "# Spot Temperature, (star's frame) kT = " << spot_temperature << " keV \n" 
        //<< "# Pulse Fraction = " << avgPulseFraction << "\n"
        << std::endl;

    out << "# Column 1: phase bins (0 to 1)\n"
        << "# Column 2: Bolometric Number flux, photons/(cm^2 s) \n"
        << "# Columns 3, 4, 5: Monochromatic Number flux (photons/(cm^2 s keV) measured at energies (at infinity) of 2, 6, 12 keV\n" 
        << "# Column 6: Number flux (photons/(cm^2 s)) in the energy band " << E_band_lower_1 << " eV to " << E_band_upper_1 << " eV \n"
        << "#"
        << std::endl; 

    for ( unsigned int i(0); i < numbins; i++ ) {
        out << curve.t[i]<< "\t" ;
        for ( unsigned int p(0); p < NCURVES; p++ ) { // 3 monochromatic energies + 1 band + bolometric
        //tw	if( std::isnan(Flux[p][i]) || Flux[p][i] == 0) std::cout << "Flux is NAN or 0 at p="<<p<<", i="<<i << std::endl;
            out << Flux[p][i] << "\t" ;
        }*/
    
    //std::cout << "Pre deleting models." << std::endl;

    delete defltoa;
    delete model;
    //std::cout << "Everything is done." << std::endl;
}
