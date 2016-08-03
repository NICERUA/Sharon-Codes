/***************************************************************************************/
/*                                     Chi.cpp

    This holds the computational methods/functions called in Spot.cpp and probably
    elsewhere.
    
    And now it computes chi^2 too!
    
    PGxx means equation xx in Poutanen & Gierlinski 2003, MNRAS 343
    MLCBxx means equation xx in Morsink, Leahy, Cadeau & Braga 2007, ApJ 663 
*/
/***************************************************************************************/

#include "matpack.h"
#include <exception>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Chi.h"
#include "OblDeflectionTOA.h"
#include "OblModelBase.h"
#include "PolyOblModelNHQS.h"
#include "Exception.h"
#include "Units.h"
#include "Struct.h"
#include "time.h"
#include <stdio.h>
using namespace std;

/**************************************************************************************/
/* ChiSquare:                                                                         */
/*           computes the chi^2 fit of data vs a simulation                           */
/*																					  */
/* pass: obsdata = light curve fluxes from observational data                         */
/*       ts = time shift or phase shift (normalized)                                  */
/**************************************************************************************/
double ChiSquare ( class DataStruct* obsdata, class LightCurve* curve) {
        
    /***************************************/
   	/* VARIABLE DECLARATIONS FOR ChiSquare */
    /***************************************/
    
    unsigned int numbins;  // Number of phase (or time) bins the light curve is cut into
    
    int k,      // Array index variable
    	new_b,  // 
    	n;      // Array index variable
    
    double ts,                              // time shift, so the phase of the simulation matches the phase of the data
           new_shift,                       // A time shift (for rebinning?)
    	   ja_check,                        // 
    	   chisquare(0.0),                  // Computed chi^2
    	   tempflux[NCURVES][MAX_NUMBINS],  // Temporary array to store the flux
    	   min_location;                    // 
    
    numbins = obsdata->numbins;
    ts = curve->para.ts;
    
    for ( unsigned int z(1); z<=1 ; z++ ) { // for different epochs
        
        while ( ts < 0.0 ) {
            ts += 1.0;
        }
        while ( ts > 1.0 ) {
            ts -= 1.0;
        }
        
        min_location = ts * numbins; // real version of the bin with the num
        new_shift = modf(min_location, &ja_check) / (numbins*1.0); // makes min_location an int
        new_b = ja_check;
        
        /*if ( ts < 1.0 ) { // ?? looks like it's for cycling around
            if ( ( ts+1 / (numbins*1.0) ) > 1.0 ) { // changed to numbins from 64 // ???
                new_b = 0; // ??
                new_shift = ts - 1.0; // ??
            }
        }*/
        
        // Rebinning the data and store shifted data back in Flux
        
        for ( unsigned int i(0); i < numbins; i++ ) {
            k = i - new_b; //May changed
            if (k > static_cast<int>(numbins)-1) k -= numbins;
            if (k < 0) k += numbins;
            
            unsigned int p = 1;
        	unsigned int q = 2;
            tempflux[p][i] = curve->f[q][k]; // putting things from curve->f's k bin into tempflux's i bin
            p = 2;
        	q = 3;
        	tempflux[p][i] = curve->f[q][k]; // putting things from curve->f's k bin into tempflux's i bin
        }
        
        for ( unsigned int i(0); i < numbins; i++ ) {
            n = i - 1;
            if ( n < 0 ) n += numbins;
            // n = i+1;
            if ( n > static_cast<int>(numbins) - 1 ) n -= numbins;
            
            unsigned int p = 1;
        	unsigned int q = 2;
            curve->f[q][i] = tempflux[p][i] + (tempflux[p][n]-tempflux[p][i]) * new_shift * numbins;
            p = 2;
        	q = 3;
        	curve->f[q][i] = tempflux[p][i] + (tempflux[p][n]-tempflux[p][i]) * new_shift * numbins;
        }

        
        // Compute chisquare for shifted data
        
        // energy band 1
        unsigned int p = 1;
        unsigned int q = 2;
        for ( unsigned int i(0); i < numbins; i++ ) {
            chisquare += pow( (obsdata->f[p][i] - curve->f[q][i])/obsdata->err[p][i], 2);
            // using different 'p' and 'q' because we're not comparing the exact same columns of obsdata.f and Flux
        }

		// energy band 2
        p = 2;
        q = 3;
        for ( unsigned int i(0); i < numbins; i++ ) {
            chisquare += pow( (obsdata->f[p][i] - curve->f[q][i])/obsdata->err[p][i], 2);
            // using different 'p' and 'q' because we're not comparing the exact same columns of obsdata.f and Flux
        }
    }
    
    //chisquare += pow( (bbrat - 0.32)/(0.064) , 2);
    
    //std::cout << "M = " << Units::nounits_to_cgs(curve->para.mass, Units::MASS)/Units::MSUN << ", R = " << Units::nounits_to_cgs(curve->para.radius, Units::LENGTH )*1.0e-5 << ", i = " 
    //<< curve->para.incl * 180.0 / Units::PI <<", e = " << curve->para.theta * 180.0 / Units::PI <<", ts = " << ts << std::endl;

    //std::cout << "\nChisquare = " << chisquare << std::endl;

    return chisquare;
    
}

/**************************************************************************************/
/* ComputeAngles:                                                                     */
/*              computes all angles necessary to create the x-ray light curve         */
/*																					  */
/* pass: incurve = contains needed values like mass, radius, inclination, theta_0     */
/*       defltoa =                                                                    */
/**************************************************************************************/
class LightCurve ComputeAngles ( class LightCurve* incurve,
				                 class OblDeflectionTOA* defltoa) {

	/*******************************************/
	/* VARIABLE DECLARATIONS FOR ComputeAngles */
	/*******************************************/
	
    class LightCurve curve;
    curve = *incurve;

    double theta_0,                   // Emission angle (latitude) of the spot, in radians          
           incl,                      // inclination angle of the observer, in radians
           mass,                      // Mass of the star, in M_sun
           radius,                    // Radius of the star (at whatever little bit we're evaluating at)
      mass_over_r,
           omega,                     // Spin frequency of the neutron star, in Hz
           cosgamma,                  // Cos of the angle between the radial vector and the surface normal vector
           shift_t,                   // Time off-set from data
           b_guess(0.0),              // Impact parameter; starting off with reasonable guess then refining it
           mu(1.0),                   // = cos(theta_0), unitless
           speed(0.0),                // Velocity of the spot, defined in MLCB34
           alpha(0.0),                // Zenith angle, in radians
           sinalpha(0.0),             // Sin of zenith angle, defined in MLCB19
           cosalpha(1.0),             // Cos of zenith angle, used in MLCB30
           b(0.0),                    // Photon's impact parameter
           toa_val(0.0),              // Time of arrival, MLCB38
           dpsi_db_val(0.0),          // Derivative of MLCB20 with respect to b
           phi_0,                     // Azimuthal location of the spot, in radians
           dS,                        // Surface area of the spot, defined in MLCB2; computed in Spot.cpp
           distance;                  // Distance from Earth to NS, inputted in meters

    double eps, epspsi,dcosa_dcosp;
           
    unsigned int numbins(MAX_NUMBINS);// Number of phase bins the light curve is split into; same as number of flux data points
    numbins = curve.numbins;


    bool ingoing(false);
    bool infile_is_set(false);

    std::vector< double > phi_em(numbins, 0.0);   // Azimuth as measured from the star's emission frame; one for each phase bin
    std::vector< double > psi(numbins, 0.0);      // Bending angle, as defined in MLCB20

    // These are calculated in the second loop.
    std::vector< double > dcosalpha_dcospsi(numbins, 0.0);    // Used in MLCB30
    std::vector< double > cosdelta(numbins, 0.0);             // 
    std::vector< double > cosxi(numbins, 0.0);                // Used in Doppler boost factor, defined in MLCB35

    // vectors for 4-point interpolation
    std::vector< double > psi_k(4, 0.0);       // Bending angle, sub k?
    std::vector< double > b_k(4, 0.0);         // Impact parameter, sub k?
    
    /************************************************************************************/
    /* SETTING THINGS UP - keep in mind that these variables are in dimensionless units */
    /************************************************************************************/
    
    dS = curve.para.dS;
    theta_0 = curve.para.theta;
    incl = curve.para.incl;
    phi_0 = curve.para.phi_0;
    mass = curve.para.mass;
    radius = curve.para.radius;
    mass_over_r = curve.para.mass_over_r;
    omega = curve.para.omega;
    cosgamma = curve.para.cosgamma;  // for spherical, cosgamma = 0
    distance = curve.para.distance;
    shift_t = curve.para.ts;
    infile_is_set = curve.flags.infile_is_set;
    speed = omega*radius*sin(theta_0) / sqrt( 1.0 - 2.0*mass_over_r ); // MLCB34
    //std::cout << "Speed = " << speed << std::endl;
    mu = cos(theta_0);

    std::cout << "ComputeAngles: radius = " << radius << std::endl;


    //initial assumptions
    curve.eclipse = false;
    curve.ingoing = false;
    curve.problem = false;
	
    /********************************************************/
    /* COMPUTE EMISSION TIME, PHASE BINS, AND LIGHT BENDING */
    /********************************************************/



    for ( unsigned int i(0); i < numbins; i++ ) { // opening For-Loop-1
        // SMM: Added an offset of phi_0
        // SMM: Time is normalized to the spin period so 0 < t_e < 1 
        // curve.t[i] = i/(1.0*numbins) + shift_t; (This is computed in Spot.cpp)
        phi_em.at(i) = phi_0 + (2.0*Units::PI) * curve.t[i]; // phi_em is the same thing as "phi" in PG; changes with t
        psi.at(i) = acos(cos(incl)*cos(theta_0) + sin(incl)*sin(theta_0)*cos(phi_em.at(i))); // PG1; this theta_0 is the theta_0 from spot.cpp
       
    } // closing For-Loop-1

    int j(0);
	
    for ( unsigned int i(0); i < numbins; i++ ) { // opening For-Loop-2
        int sign(0);
        double bval(0.0);
        bool result(false);
        double b1(0.0), b2(0.0), psi1(0.0), psi2(0.0);
        double xb(0.0);
        int k(0);
        //j=0;
        /**************************************************************************/
	/* TEST FOR VISIBILITY FOR EACH VALUE OF b, THE PHOTON'S IMPACT PARAMETER */
	/**************************************************************************/

        if ( psi.at(i) < curve.defl.psi_max ) {
            if ( psi.at(i) > curve.defl.psi_b[j] )
	            while ( psi.at(i) > curve.defl.psi_b[j] ) {
	                j++;      
                }
            else {
	            while ( psi.at(i) < curve.defl.psi_b[j] ) {
	              j--;
	            }
	           j++;
            }
      
            b1 = curve.defl.b_psi[j-1];
            b2 = curve.defl.b_psi[j];
            psi1 = curve.defl.psi_b[j-1];
            psi2 = curve.defl.psi_b[j];
            k = j - 2;
      
            if ( j == 1 ) k = 0;
            if ( j == 3 * NN ) k = 3 * NN - 3;

            for ( j = 0; j < 4; j++ ) {
	            b_k.at(j) = curve.defl.b_psi[k+j];
	            psi_k.at(j) = curve.defl.psi_b[k+j];
            }
  
            // 4-pt interpolation to find the correct value of b given psi.
            xb = psi.at(i);
            b_guess = (xb-psi_k.at(1))*(xb-psi_k.at(2))*(xb-psi_k.at(3))*b_k.at(0)/
	                  ((psi_k.at(0)-psi_k.at(1))*(psi_k.at(0)-psi_k.at(2))*(psi_k.at(0)-psi_k.at(3)))
	                  +(xb-psi_k.at(0))*(xb-psi_k.at(2))*(xb-psi_k.at(3))*b_k.at(1)/
	                  ((psi_k.at(1)-psi_k.at(0))*(psi_k.at(1)-psi_k.at(2))*(psi_k.at(1)-psi_k.at(3)))
	                  +(xb-psi_k.at(0))*(xb-psi_k.at(1))*(xb-psi_k.at(3))*b_k.at(2)/
	                  ((psi_k.at(2)-psi_k.at(0))*(psi_k.at(2)-psi_k.at(1))*(psi_k.at(2)-psi_k.at(3)))
	                  +(xb-psi_k.at(0))*(xb-psi_k.at(1))*(xb-psi_k.at(2))*b_k.at(3)/
                	  ((psi_k.at(3)-psi_k.at(0))*(psi_k.at(3)-psi_k.at(1))*(psi_k.at(3)-psi_k.at(2)));
        } // ending psi.at(i) < curve.defl.psi_max
        
        /***********************************************/
	/* FINDING IF A SOLUTION EXISTS, SETTING FLAGS */
	/***********************************************/
		
        result = defltoa->b_from_psi( fabs(psi.at(i)), radius, mu, bval, sign, curve.defl.b_max, 
        		 curve.defl.psi_max, b_guess, fabs(psi.at(i)), b2, fabs(psi.at(i))-psi2, 
        		 &curve.problem );
        if ( result == false ) { 
            curve.visible[i] = false;
            curve.t_o[i] = curve.t[i] + curve.t_o[i-1] - curve.t[i-1] ;
            curve.dOmega_s[i] = 0.0;
            curve.eclipse = true;
        }
		
        else { // there is a solution
            b = bval;
            if ( sign < 0 ) { // if the photon is initially ingoing (only a problem in oblate models)
	            ingoing = true;
	            curve.ingoing = true;
	            //std::cout << "ingoing!"<< std::endl;
            }
            else if ( sign > 0 ) {
	            ingoing = false;
            }
            else {
	            throw( Exception("Chi.cpp: sign not returned as + or - with success.") ); // used to say "ObFluxApp.cpp"
            }
            
			double b_maximum = radius/sqrt(1.0 - 2.0*mass_over_r);
			if ( (fabs(b-b_maximum) < 1e-7) && (b > 0.0) && (b > b_maximum) ) { 
			// this corrects for b being ever so slightly over bmax, which yields all kinds of errors in OblDeflectionTOA
				std::cout << "Setting b = b_max." << std::endl;
				b = b_maximum - DBL_EPSILON;
			}
            curve.b[i] = b/radius;
            
            /*******************************************************/
	    /* IF A SOLUTION EXISTS, CALCULATING ANGLES:  sinalpha */
	    /*                    						  cosalpha */
	    /*                    						  alpha    */
	    /*                   						  cosdelta */
	    /*                							  cosbeta  */
	    /*******************************************************/

            sinalpha =  b * sqrt( 1.0 - 2.0 * mass_over_r ) / radius;  // PG4, reconfigured
            cosalpha = sqrt( 1.0 - sinalpha * sinalpha ); 
            alpha    = asin( sinalpha );
            

            if ( sign < 0 ) { // alpha is greater than pi/2
	            alpha = Units::PI - alpha;
            }

            cosdelta.at(i) =  (cos(incl) - cos(theta_0)*cos(psi.at(i))) / (sin(theta_0)*sin(psi.at(i))) ;
            if ( theta_0 == 0.0 )  // cosdelta needs to be redone if theta_0 = 0
            	cosdelta.at(i) = sqrt( 1 - pow( sin(incl) * sin(phi_em.at(i)) / sin(psi.at(i)) ,2) ); // law of sines from MLCB Fig 1 and in MLCB17 and just above
     
            if ( (cos(theta_0) < 0) && (cos(incl) < 0 ) ) {
	            cosdelta.at(i) *= -1.0;
            }
            if ( sin(psi.at(i)) != 0.0 ) {
	            curve.cosbeta[i] = cosalpha * cosgamma + sinalpha * sqrt( 1.0 - pow( cosgamma, 2.0 )) * cosdelta.at(i);
	            //if( std::isnan(curve.cosbeta[i]) ) std::cout << "cosdelta.at(i="<<i<<") = " << cosdelta.at(i) << std::endl;
            }
            else {
	            curve.cosbeta[i] = cosalpha * cosgamma;
            }

            if ( cosalpha < 0.01 ) {
	            curve.cosbeta[i] = (Units::PI/2.0 - alpha + sqrt(2) * sqrt(1.0-cosgamma) * cosdelta.at(i));
  			}
            if ( curve.cosbeta[i] < 0.0 && curve.cosbeta[i] > 1.0 ) { // cosbeta > 1.0 added by Abbie, Feb 22 2013
	            std::cerr << "i = " << i
	        	          << ", Not visible at phi_em = " << 180.0 * phi_em.at(i) / Units::PI 
	                   	  << ", cos(beta) = " << curve.cosbeta[i] 
		                  << ", cos(alpha) = " << cosalpha
		                  << ", cos(gamma) = " << cosgamma
		                  << ", cos(delta) = " << cosdelta.at(i)
		                  << " (visibility condition)." << std::endl << std::endl;
	            curve.visible[i] = false;
            }
            else {
	            curve.visible[i] = true;
            }
			
			/********************************************************/
			/* IF THE PHOTON IS VISIBLE, COMPUTE: dpsi_db           */     
			/*                                    toa               */
			/*									  cosxi             */
			/*									  eta               */
			/*									  t_o               */
			/*									  psi               */
			/*									  R_dpsi_db         */
			/*									  dcosalpha_dcospsi */
			/*									  dOmega_s          */
			/********************************************************/
			
            if ( curve.visible[i] ) { // visible 
            	if (alpha == 0.0 && psi.at(i) == 0.0 & phi_em.at(i) == 0.0) 
            		cosxi.at(i) = 0.0; // to avoid NAN errors from dividing by 0; appears when incl = theta at i=0
	            else 
	            	cosxi.at(i) = - sinalpha * sin(incl) * sin(phi_em.at(i)) / sin(fabs(psi.at(i)));  // PG11
	            curve.eta[i] = sqrt( 1.0 - speed*speed ) / (1.0 - speed*cosxi.at(i) ); // Doppler boost factor, MLCB33
	            
	            if ((1.0 - speed*cosxi.at(i)) == 0.0) 
	            	std::cout << "dividing by zero" << std::endl;
	            if((std::isnan(speed) || speed == 0.0) && theta_0 != 0.0 )
	            	std::cout << "speed = " << speed << " at i = " << i << std::endl;
	            //if(std::isnan(cosxi.at(i)) || cosxi.at(i) == 0.0) 
	            //	std::cout << "cosxi(i="<<i<<") = " << cosxi.at(i) << std::endl;
	            

	            if ( ingoing ) {
	             //  std::cout << "Ingoing b = " << b << std::endl;
		      dpsi_db_val = defltoa->dpsi_db_ingoing( b, radius, mu, &curve.problem );
		      toa_val = defltoa->toa_ingoing( b, radius, mu, &curve.problem );
	            }
                else {

		  if (b != curve.defl.b_max ){
	                dpsi_db_val = defltoa->dpsi_db_outgoing_u( b, radius, &curve.problem );
	                //if (i == 0) std::cout << "b = " << b <<", dpsi_db = " << dpsi_db_val << std::endl;

	                toa_val = defltoa->toa_outgoing_u( b, radius, &curve.problem );
	                //std::cout << "dpsi_db_val = " << dpsi_db_val << ", toa_val = " << toa_val << std::endl;
		  }
		  else{
		    toa_val = defltoa->toa_outgoing_u( b, radius, &curve.problem );
		    	eps = 1e-6;
			b = curve.defl.b_max * sqrt(1.0 - eps);
			epspsi = defltoa->psi_outgoing( b, radius, curve.defl.b_max, curve.defl.psi_max, &curve.problem);
			dpsi_db_val = defltoa->dpsi_db_outgoing( b, radius, &curve.problem );
			dcosa_dcosp = sqrt(1.0-2*mass_over_r) / (sqrt(eps) * sin(fabs(epspsi)) * radius * dpsi_db_val) * sqrt(1.0-eps);

		  }


		}

	            //std::cout << "Done computing TOA " << std::endl;
	            curve.t_o[i] = curve.t[i] + (omega * toa_val) / (2.0 * Units::PI);
	            curve.psi[i] = psi.at(i);
	            curve.R_dpsi_db[i] = dpsi_db_val * radius;

		    if (b != curve.defl.b_max){
		    if ( psi.at(i) == 0 && alpha == 0 ) 
		      curve.dcosalpha_dcospsi[i] = fabs( (1.0 - 2.0 * mass_over_r) / curve.R_dpsi_db[i]);
		    //if (psi.at(i) == 0 && alpha == 0 ) curve.dcosalpha_dcospsi[i] = 0.0;
	            else 
		      curve.dcosalpha_dcospsi[i] = fabs( sinalpha/cosalpha * sqrt(1.0 - 2.0*mass_over_r) / (sin(fabs(psi.at(i))) * curve.R_dpsi_db[i]) );
		    }
		    else
		      curve.dcosalpha_dcospsi[i] = dcosa_dcosp; 
		   

	            curve.dOmega_s[i] = (dS / (distance * distance)) 
	                               * (1.0 / (1.0 - 2.0 * mass_over_r)) 
	                               * curve.cosbeta[i] 
	                               * curve.dcosalpha_dcospsi[i];  // PG8
		    //std::cout << "t = " << curve.t[i] << " dOmega = " << curve.dOmega_s[i] << std::endl;

	            /***********************************/
				/* FLAGS IF A VALUE IS NAN OR ZERO */
				/***********************************/
	            if (std::isnan(dpsi_db_val) || dpsi_db_val == 0) std::cout << "dpsi_db_val = " << dpsi_db_val << "at i = " << i << std::endl;
				if (std::isnan(psi.at(i))) std::cout << "psi.at(i="<<i<<") = " << psi.at(i) << std::endl;
				if (std::isnan(curve.dOmega_s[i])) std::cout << "dOmega is NAN at i = " << i << std::endl;
				if (std::isnan(curve.cosbeta[i])) std::cout << "cosbeta is NAN at i = " << i << std::endl;
				if (std::isnan(curve.dcosalpha_dcospsi[i])) std::cout << "dcosalpha_dcospsi is NAN at i = " << i<< std::endl;
				if (std::isnan(sinalpha)) std::cout << "sinalpha is NAN at i = " << i << std::endl;
				if (std::isnan(cosalpha)) std::cout << "cosalpha is NAN at i = " << i << std::endl;
				if (std::isnan(psi.at(i))) std::cout << "psi is NAN at i = " << i << std::endl;
				if (std::isnan(mass)) std::cout << "mass is NAN at i = " << i << std::endl;
				if (std::isnan(radius)) std::cout << "radius is NAN at i = " << i << std::endl;
				

            } // end visible
            
            /**************************************************************/
			/* IF THE PHOTON IS NOT VISIBLE, POLITELY CRASH OR GO TO ZERO */
			/**************************************************************/
      
            else { // not visible; we think that it shouldn't matter if it's not visible at i=0
	      curve.t_o[i] = curve.t[i] + (curve.t_o[i-1] - curve.t[i-1]) ; // t_o is not defined properly, so we'll set it to emission time
	      //curve.t_o[i] = curve.t[i] ;
	      curve.dOmega_s[i] = 0.0;    // don't see the spot, so dOmega = 0
	      curve.cosbeta[i] = 0.0;     // doesn't matter, doesn't enter into calculation
	      curve.eta[i] = 1.0;	        // doesn't matter, doesn't enter into calculation
	      std::cout << "toa = " << (curve.t_o[i-1] - curve.t[i-1]) << std::endl;
	            //}
            } // end not visible
        } // end "there is a solution"
    }  // closing For-Loop-2

    return curve;

} // End ComputeAngles

/**************************************************************************************/
/* Bend:                                                                              */
/*              Creates a Table of alpha, b, psi, dpsi_db                             */
/*																					  */
/* pass: incurve = contains needed values like mass, radius, inclination, theta_0     */
/*       defltoa =                                                                    */
/**************************************************************************************/
void Bend ( class LightCurve* incurve,
				                 class OblDeflectionTOA* defltoa) {

	/*******************************************/
	/* VARIABLE DECLARATIONS FOR ComputeAngles */
	/*******************************************/
	
    class LightCurve curve;
    curve = *incurve;

    std::ofstream bend;      // output stream; printing information to the output file

    double 
      eps,
           mass,                      // Mass of the star, in M_sun
           radius,                    // Radius of the star (at whatever little bit we're evaluating at)
      mass_over_r,
           alpha(0.0),                // Zenith angle, in radians
           sinalpha(0.0),             // Sin of zenith angle, defined in MLCB19
           cosalpha(1.0),             // Cos of zenith angle, used in MLCB30
           b(0.0),                    // Photon's impact parameter
      psi(0.0),
      b_max,
      psi_max,
      toa_val(0.0),              // Time of arrival, MLCB38
      dpsi_db_val(0.0),          // Derivative of MLCB20 with respect to b
      dcosa_dcosp;
          
    unsigned int numbins(MAX_NUMBINS);// Number of phase bins the light curve is split into; same as number of flux data points
    numbins = 1000;


    std::cout << "Entering Bend" << std::endl;

    /************************************************************************************/
    /* SETTING THINGS UP - keep in mind that these variables are in dimensionless units */
    /************************************************************************************/
    
    mass = curve.para.mass;
    radius = curve.para.radius;
    mass_over_r = curve.para.mass_over_r;

    b_max = curve.defl.b_max;
    // psi_max = curve.defl.psi_max;

    psi_max = defltoa->psi_max_outgoing_u(b_max,radius, &curve.problem);
    std::cout << "BEND: psi_max = " << psi_max << std::endl;

    //initial assumptions
    curve.problem = false;
	
    // Open output file and print out header

    bend.open("bend.txt", std::ios_base::trunc);
    bend.precision(10);
    bend << "# Mass = " <<  Units::nounits_to_cgs(mass, Units::MASS)/Units::MSUN << " Msun "    << std::endl;
    bend << "# R_sp = " << Units::nounits_to_cgs(radius, Units::LENGTH )*1.0e-5 << " km " << std::endl;
    bend << "# M/R = " << curve.para.mass_over_r << std::endl;
    bend << "# R/M = " << 1.0/curve.para.mass_over_r << std::endl;
    bend << "#alpha #b/R #psi #dcosalpha/dcospsi #toa*C/R" << std::endl;

    for (unsigned int i(0); i < numbins; i++ ) { 

      alpha = i/(numbins-1.0) * Units::PI/2.0; 
      sinalpha = sin(alpha);
      cosalpha = sqrt( 1.0 - sinalpha * sinalpha ); 

      b = sinalpha * radius / sqrt( 1.0 - 2.0 * mass_over_r );
 
      psi = defltoa->psi_outgoing_u( b, radius, b_max, psi_max, &curve.problem);

      dpsi_db_val = defltoa->dpsi_db_outgoing_u( b, radius, &curve.problem );
	               
      toa_val = defltoa->toa_outgoing_u( b, radius, &curve.problem );

      std::cout << "bend: alpha = " << alpha << " b/r = " << b/radius << " psi = " << psi << " dpsi = " << dpsi_db_val << std::endl; 

      if (psi==0 && alpha == 0)
	dcosa_dcosp =  fabs( (1.0 - 2.0 * mass_over_r) / (radius * dpsi_db_val));
      	//dcosa_dcosp =  fabs( (1.0 - 2.0 * mass_over_r));
      else
	dcosa_dcosp = fabs( sinalpha/cosalpha * sqrt(1.0 - 2.0*mass_over_r) / (sin(fabs(psi)) * radius*dpsi_db_val));

       if (cosalpha == 0){

	eps = 5e-7;
	b = b_max * sqrt(1.0 - eps);
	psi = defltoa->psi_outgoing_u( b, radius, b_max, psi_max, &curve.problem);
	dpsi_db_val = defltoa->dpsi_db_outgoing_u( b, radius, &curve.problem );
	dcosa_dcosp = sqrt(1.0-2*mass_over_r) / (sqrt(eps) * sin(fabs(psi)) * radius * dpsi_db_val) * sqrt(1.0-eps);

	b = b_max;
	psi = defltoa->psi_outgoing_u( b, radius, b_max, psi_max, &curve.problem);

	}
		
      bend 
	//<< i << " " 
	<< alpha << " " 
	<< b/radius << " " 
	//	<< cosalpha << " "
	<< psi << " " 
	//<< cos(psi) << " "
	//<< 1.0/dcosa_dcosp << " " 
	<< dcosa_dcosp << " " 
	//<< Units::nounits_to_cgs(toa_val,Units::TIME)  << " " 
	<< toa_val / radius 
	<< std::endl;
	       			   
    }


} // End Bend


/**************************************************************************************/
/* ComputeCurve:                                                                      */
/*              computes the flux of each light curve                                 */
/*																					  */
/* pass: angles = all the angles necessary to compute the flux correctly;             */
/*                computed in the routine/method/function above [radians or unitless] */
/**************************************************************************************/
class LightCurve ComputeCurve( class LightCurve* angles ) {
	
	/******************************************/
	/* VARIABLE DECLARATIONS FOR ComputeCurve */
	/******************************************/
	
    class LightCurve curve;

    std::ifstream input;
    std::ofstream output;
    std::ofstream ttt;


    bool infile_is_set;   // If the user has specified an input file
    double 
           mass,               // Mass of the star, in M_sun
           radius,             // Radius of the star at the spot, in km
      mass_over_r,
           temperature,        // Temperature of the spot, in keV
      //           E_mono,             // Energy observed, in keV, for each monochromatic energy
           E_band_lower_1,     // Lower bound of energy band for flux integration, in keV
           E_band_upper_1,     // Upper bound of energy band for flux integration, in keV
           E_band_lower_2,     // Lower bound of energy band for flux integration, in keV
           E_band_upper_2,     // Upper bound of energy band for flux integration, in keV
           redshift,           // Gravitational redshift = 1 + z = (1-2M/R)^{-1/2}
           bolo,               // Bolometric flux; bolo = sigma T^4/pi
           gray(1.0);          // Graybody factor (when = 1, not effective)
        
    double E0, E1, E2, DeltaE, E_obs;

    unsigned int numbins(MAX_NUMBINS);  // Time bins of light curve (usually 128)
    unsigned int numbands(NCURVES);  // Number of Energy Bands
        
    std::vector< double > totflux(MAX_NUMBINS, 0.0); // integrated flux

    std::vector< bool > nullcurve(NCURVES,true); // true means that the curve is zero everywhere

    //std::vector< double > softbb(MAX_NUMBINS, 0.0);  // blackbody soft flux
    //std::vector< double > softcm(MAX_NUMBINS, 0.0);  // compton soft flux

    // double softbbave(0.0), softcmave(0.0), hardave(0.0);  // average value of softbb and softcm and high energy band

    /*********************/
    /* SETTING THINGS UP */
    /*********************/
    
    // One monochromatic energy, hardwired value, in keV
    //    E_mono = 1.0;

    curve = (*angles);
    mass = curve.para.mass;                     // unitless
    radius = curve.para.radius;                 // unitless
    mass_over_r = curve.para.mass_over_r;
    temperature = curve.para.temperature;       // T in keV 
    infile_is_set = curve.flags.infile_is_set;
    numbins = curve.numbins;
    numbands = curve.numbands;
    E_band_lower_1 = curve.para.E_band_lower_1;     // in keV
    E_band_upper_1 = curve.para.E_band_upper_1;     // in keV
    E_band_lower_2 = curve.para.E_band_lower_2;     // in keV
    E_band_upper_2 = curve.para.E_band_upper_2;     // in keV
    E0 = curve.para.E0;
    E1 = curve.para.E1;
    E2 = curve.para.E2;
    DeltaE = curve.para.DeltaE;
   
    redshift = 1.0 / sqrt( 1 - 2.0 * mass_over_r);

    //bolo = 2.0e12/15.0 * pow(Units::H_PLANCK,-3) * pow(Units::C,-2) * pow( Units::PI * Units::EV * temperature, 4);
    //bolo *= 1.0e-3/Units::EV;

    bolo = 2.0e9 * 2.404 * Units::C * pow(temperature * Units::EV/(Units::H_PLANCK * Units::C) , 3); // use this one! probably!
    // the e9 in the beginning is for changing T^3 from keV to eV
    // 2.404 comes from evaluating Bradt equation 6.17 (modified, for photon number count units), using the Riemann zeta function for z=3

    for ( unsigned int i(0); i < numbins; i++ ) { // Compute flux for each phase bin

      if ( curve.dOmega_s[i] != 0.0 ) {
	if ( curve.flags.beaming_model == 1 ) {
	  gray = Gray(curve.cosbeta[i]*curve.eta[i]); // gets the graybody factor
	}
	else
	  gray = 1.0; // Isotropic, the same at all angles

	if (std::isnan(gray) || gray == 0) std::cout << "gray = " << gray << std::endl;
	if (std::isnan(curve.eta[i]) || curve.eta[i] == 0) std::cout << "eta[i="<<i<<"] = " << curve.eta[i] << std::endl;
	if (std::isnan(redshift) || redshift == 0) std::cout << "redshift = " << redshift << std::endl;


	if (curve.flags.spectral_model == 0){ // Monochromatic Observation of Blackbody
	/*******************************************************************/
	/* COMPUTING BLACKBODY LIGHT CURVE FOR MONOCHROMATIC ENERGY, p = 0 */
	/*      First computes in [erg/(s cm^2 Hz), converts to            */
	/*		photons/(s cm^2 keV)                               */
	/*******************************************************************/
	  //std::cout << "E_obs = " << E0 << std::endl;

	  // Moonochromatic light curve in energy flux erg/(s cm^2 Hz)
	  curve.f[0][i] = gray * curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * BlackBody(temperature,E0*redshift/curve.eta[i]); 
	  // Units: erg/(s cm^2 Hz)
	  //Convert to photons/(s cm^2 keV)
	  curve.f[0][i] *= (1.0 / ( E0 * Units::H_PLANCK )); // Units: photons/(s cm^2 keV)

	  if (curve.f[0][i] != 0.0) nullcurve[0] = false;
	}

	if (curve.flags.spectral_model == 1){ // Funny Line Emission for NICER
	  for (unsigned int p=0; p<numbands; p++){
	    
	    E_obs = E0 + p*DeltaE;
	    //curve.f[p][i] = gray * curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * Line(temperature,E_obs*redshift/curve.eta[i],E1,E2); 
	    // Units: erg/(s cm^2 Hz)
	    //Convert to photons/(s cm^2 keV)
	    //curve.f[p][i] *= (1.0 / ( E_obs * Units::H_PLANCK )); // Units: photons/(s cm^2 keV)
	    // Multiply by the width of the energy band;
	    //curve.f[p][i] *= DeltaE; // Units: photons/(s cm^2)

	    //curve.eta[i] = 1.0;

	    curve.f[p][i] = gray * curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * LineBandFlux(temperature, (E_obs-0.5*DeltaE)*redshift/curve.eta[i], (E_obs+0.5*DeltaE)*redshift/curve.eta[i], E1, E2); // Units: photon/(s cm^2)

	    if (curve.f[p][i] != 0.0) nullcurve[0] = false;

	  }
	}

	if (curve.flags.spectral_model == 7) { // Not Used Right Now.
			
	  /*******************************************/
	  /* COMPUTING BOLOMETRIC LIGHT CURVE, p = 0 */
	  /* Units: photons/(cm^2 s)                 */
	  /*******************************************/
    		
	  // Bolometric Light Curve for energy flux erg/(cm^2 s)
	  curve.f[0][i] =  bolo * gray * curve.dOmega_s[i] * pow(curve.eta[i],5) * pow(redshift,-4); // Units: erg/(cm^2 s)

	  // Bolometric Light Curve for photon number flux photons/(cm^2 s)
	  curve.f[0][i] = bolo * gray * curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3); // Units: photons/(cm^2 s)
                                     
	  /***************************************************/
	  /* COMPUTING PHOTON NUMBER FLUX FOR AN ENERGY BAND */
	  /*   	1 < p < NCURVES-1                          */
	  /*      Units: photons/(s cm^2)                    */
	  /***************************************************/
    		
	  // First energy band
	  curve.f[NCURVES-2][i] = gray * curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * EnergyBandFlux(temperature, E_band_lower_1*redshift/curve.eta[i], E_band_upper_1*redshift/curve.eta[i]); // Units: photon/(s cm^2)
	  // Second energy band
	  curve.f[NCURVES-1][i] = gray * curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * EnergyBandFlux(temperature, E_band_lower_2*redshift/curve.eta[i], E_band_upper_2*redshift/curve.eta[i]); // Units: photon/(s cm^2)
	}		
      }
      else { // if curve.dOmega_s[i] == 0.0
	for ( unsigned int p(0); p < numbands; p++) {
	  curve.f[p][i] = 0.0;
	}
      }

    } // ending the for(i) loop
    

	
    /***********************************************************/
    /* DEALING WITH THE TIME DELAYS, REBINNING THE LIGHT CURVE */
    /* This is where the jumpy problems tend to be.            */
    /***********************************************************/

    curve.flags.ignore_time_delays = false;
    		
    if ( !curve.flags.ignore_time_delays ) { // if we are not ignoring the time delays        
      int k(0), j(0); // index placeholders; approximately, k is i+1 and j is i-1
      // but time delays mean that j isn't always i-1
      // used in the linear interpolation
        
      /********************************/
      /* LOOP THROUGH THE LIGHTCURVES */
      /********************************/
    	
      for (unsigned int p(0); p < numbands; p++) {
      // for (unsigned int p(1); p < 2; p++) {

	if (!nullcurve[p]) {
	

	/*********************/
	/* MORE DECLARATIONS */
	/*********************/
    		
	  std::vector< double > newflux(MAX_NUMBINS, 0.0);                  // rebinned flux (to account for photon travel time)
	  unsigned int imax(0), imin(0);                                    // index of element with maximum flux value, minimum flux value
	  double max_discrete_flux(0.0), min_discrete_flux(curve.f[p][0]);  // maximum flux value and minimum flux value assigned to a grid point (discrete)
	  double tx(0), ta(0), tb(0), tc(0);                                // a,b,c: three-point interpolation on a parabola for min and max areas (time a, time b, time c)
	  double fa(0), fb(0), fc(0);                                       // corresponding fluxes for three-point interpolation
	  int ia(0), ic(0);                                                 // indices in array for where flux is fa, fc
	  double temporary1(0.0), temporary2(0.0);                          // makes math easier
	  double maximum(0.0), minimum(100000.0);                           // true (continuous) maximum and minimum flux values
	  double tmin(0.0);                                                 // value of t at the true minimum

	  int ec1(0), ec2(0),j1,j2;
	  double te1, te2, tt1, tt2, slope1, slope2;

	  /**********************************************************/
	  /* START WITH FINDING THE MAXIMUM FLUX OF THE LIGHT CURVE */
	  /**********************************************************/
			
	  /********************************************************/
	  /* FINDING THE DISCRETE MAXIMUM FLUX OF THE LIGHT CURVE */
	  /********************************************************/
    		
	  for ( unsigned int i(0); i < numbins; i++ ) {
	    if ( curve.f[p][i] > max_discrete_flux ) { // tells you where the maximum is
	      imax = i;  
	      max_discrete_flux = curve.f[p][i];
	    }
	  }
	        
	  /******************************************/
	  /* FINDING THE TRUE MAXIMUM VALUE OF FLUX */
	  /******************************************/
    		
	  if ( imax == 0 ) {
	    ia = numbins - 1;
	    ta = curve.t_o[ia] - 1.0;
	  }
	  else { // else imax == -1
	    ia = imax - 1;
	    ta = curve.t_o[ia];
	  }
	  tb = curve.t_o[imax];
	  if ( imax == numbins - 1 ) {
	    ic = 0;
	    tc = curve.t_o[ic] + 1.0;
	  }
	  else { //else imax == 0
	    ic = imax + 1;
	    tc = curve.t_o[ic];
	  }
	  fa = curve.f[p][ia];
	  fb = curve.f[p][imax];
	  fc = curve.f[p][ic];
            
	  /**********************************************/
	  /* NUMERICAL RECIPES, PARABOLIC INTERPOLATION */
	  /* Equation 10.3.1 (our t is their x)         */
	  /**********************************************/
    		
	  if ( ( (tb-ta)*(fb-fc) - (tb-tc)*(fb-fa) ) != 0.0 ) { // otherwise you get a big fat NAN for a flux
	    tx = tb - 0.5 * (pow(tb-ta,2)*(fb-fc) - pow(tb-tc,2)*(fb-fa)) / ((tb-ta)*(fb-fc) - (tb-tc)*(fb-fa));
	    temporary1 =  (fa-fc)/( pow(tc-tx,2) - pow(ta-tx,2)) ;
	  }
	  if ( ( (tb-ta)*(fb-fc) - (tb-tc)*(fb-fa) ) == 0.0 ) { // to avoid dividing by zero
	    /********/      tx = 0.0; // is this what it should be?
	  }
	  if ( ( pow(tc-tx,2) - pow(ta-tx,2) ) != 0.0 ) {
	    temporary1 =  (fa-fc)/( pow(tc-tx,2) - pow(ta-tx,2)) ;
	  }
	  if ( ( pow(tc-tx,2) - pow(ta-tx,2) ) == 0.0 ) {
	    temporary1 = 0.0;
	  }
	  maximum = fb - temporary1*pow(tb-tx,2);
	  // maximum is the maximum flux for the specific "p" light curve
	  // tx is the time of maximum 

	  /************************************************/
	  /* NOW FIND THE MINIMUM FLUX OF THE LIGHT CURVE */
	  /* Note: real light curves won't have eclipses  */
	  /* If eclipsed, then min = 0                    */
	  /************************************************/
    		
	  if ( !curve.eclipse ){// && !curve.ingoing ) {   
            
	    /********************************************************/
	    /* FINDING THE DISCRETE MINIMUM FLUX OF THE LIGHT CURVE */
	    /* If not eclipsed                                      */
	    /********************************************************/
    			
	    for ( unsigned int i(0); i < numbins; i++) {
	      if (curve.f[p][i] < min_discrete_flux) {
		imin = i;
		min_discrete_flux = curve.f[p][i];
	      }
	    }
	        
	    //std::cout << "imin = " << imin << "minflux discrete = " << min_discrete_flux << std::endl;
    
	    /******************************************/
	    /* FINDING THE TRUE MINIMUM VALUE OF FLUX */
	    /* If not eclipsed                        */
	    /******************************************/
    			
	    if ( imin == 0 ) {
	      ia = numbins - 1;
	      ta = curve.t_o[ia] - 1.0;
	    }
	    else {  //else imin == 1
	      ia = imin - 1;
	      ta = curve.t_o[ia];
	    }
	    tb = curve.t_o[imin];
	    if ( imin == numbins-1 ) {
	      ic = 0;
	      tc = curve.t_o[ic]+1.0;
	    }
	    else {  //else imin == 2
	      ic = imin+1;
	      tc = curve.t_o[ic];
	    }
	    fa = curve.f[p][ia];
	    fb = curve.f[p][imin];
	    fc = curve.f[p][ic];
                	
	    /**********************************************/
	    /* NUMERICAL RECIPES, PARABOLIC INTERPOLATION */
	    /* Equation 10.3.1 (our t is their x)         */
	    /**********************************************/

	    if ( ( (tb-ta)*(fb-fc) - (tb-tc)*(fb-fa) ) != 0.0 ) { // otherwise you get a big fat NAN for a flux
	      tmin = tb - 0.5*(pow(tb-ta,2)*(fb-fc) - pow(tb-tc,2)*(fb-fa)) / ((tb-ta)*(fb-fc) - (tb-tc)*(fb-fa));  
	    }
	    if ( ( pow(tc-tmin,2) - pow(ta-tmin,2) ) != 0.0 ) {
	      temporary2 =  (fa-fc)/( pow(tc-tmin,2) - pow(ta-tmin,2)) ;
	    }
	    if ( ( (tb-ta)*(fb-fc) - (tb-tc)*(fb-fa) ) == 0.0) {  // what should happen for dividing by zero?
	      tmin = 0.0;
	    }
	    if ( ( pow(tc-tmin,2) - pow(ta-tmin,2) ) == 0.0 ) { // what should happen for dividing by zero?
	      temporary2 = 0.0;
	    }
           	 	
	    minimum = fb - temporary2 * pow(tb-tmin,2);

	    /*	    std::cout << "imin = " << imin 
		      << "tmin = " << tmin
		      << "True minflux = " << minimum << std::endl;*/
           	 
            } // ending "not eclipsed, not ingoing" section

	  // minimum is the minimum value of flux
	  // tmin is the time of the minimum
	  // only true of not eclipsed

            
	  /*******************************************/
	  /* FOR ECLIPSING LIGHT CURVES, MINIMUM = 0 */
	  /*******************************************/
    		
	  else if ( curve.eclipse ) {
	    minimum = 0.0;

	    // Find the first eclipsed bin; set ec1=i;
	    for (unsigned int i(0);i<numbins;i++){
	      if (curve.f[p][i]==0.0){
		ec1=i;
		break;
	      }
	    }
	    
	    //	    std::cout << "eclipse begins at bin" << ec1 << std::endl;
	    // Extrapolate to find beginning of eclipse (linear extrapolation)
	    j1 = ec1-2;
	    j2 = ec1-1;

	    if (ec1==0){
	      j1+=numbins;
	      j2+=numbins;
	    }
	    if (ec1==1){
	      j1+=numbins;
	    }
	    tt1 = curve.t_o[j1];
	    if (ec1<2)
	      tt1 +=-1.0;
	    tt2 = curve.t_o[j2];
	    if (ec1==0)
	      tt2+=-1.0;
	    


	    slope1 = (curve.f[p][j1]-curve.f[p][j2])/(tt1-tt2);
	
	    if (slope1 != 0.0)
	      te1  = - (curve.f[p][j2])/slope1 + tt2 ;
	    else
	      te1 = tt1;


	    // After the eclipse, find the first bin with non-zero flux
	    for (unsigned int i(ec1);i<numbins;i++){
	      if (curve.f[p][i]>0.0){
		ec2=i;
		break;
	      }
	    }
	    //std::cout << "eclipse ends at bin" << ec2 << std::endl;
	    j1=ec2;
	    j2=ec2+1;
	    if (ec2==numbins-1){
	      j2+=-numbins;
	    }
	    tt1=curve.t_o[j1];
	    tt2=curve.t_o[j2];
	    if (ec2==numbins-1)
	      tt2+=1.0;

	    slope2 = (curve.f[p][j1]-curve.f[p][j2])/(tt1-tt2);
	    te2 = tt1 - curve.f[p][j1]/slope2;
	    //std::cout << "Eclipse ends at time = " << te2 << std::endl;


	  } // ending "yes eclipsed
            
            /****************************************************************************/
	  /* FOR NOT ECLIPSED AND YES INGOING PHOTONS, POLITELY SET TO ZERO AND CRASH */
	  /****************************************************************************/
    		
	  else {    //if there is no eclipse and the photon is intially ingoing
	    throw( Exception(" Photons are ingoing and not eclipsing. This is an issue that the code cannot handle.") );
	    minimum = 0.0;
	    maximum = 0.0;
	    break;
	  } // ending "not eclipsed, yes ingoing" section
            
	  /****************************************************/
	  /* COMPUTING THE PULSE FRACTION FOR THE LIGHT CURVE */
	  /****************************************************/
   			
	  curve.minFlux[p] = minimum;
	  curve.maxFlux[p] = maximum;
	  curve.pulseFraction[p] = (curve.maxFlux[p] - curve.minFlux[p]) / (curve.maxFlux[p] + curve.minFlux[p]);

	  curve.asym[p] = (tmin - tx) - 0.5;
	  if (curve.asym[p] < 0.0) curve.asym[p]+=1.0;


	  // Initializing totflux
	  for ( unsigned int i(0); i < numbins; i++ )
	    totflux.at(i) = 0.0;
			
		           		
	  /**************************************************************/
	  /* ADDING FLUXES FROM ALL THE PHASE BINS AND OTHER FUN THINGS */
	  /**************************************************************/

	  for ( unsigned int i(0); i < numbins; i++ ) {  // for-i-loop, looping through the phase bins
            
	    k = i + 1;
	    j = i;		 
	    if ( k == static_cast<int>(numbins) ) 
	      k = 0;
	 		
	    /* if (curve.eclipse) {  // begin eclipse section
	      if (curve.f[p][k] == 0.0) {
		if (curve.f[p][j] != 0.0) {
		  k = i;
		  j = i-1;
		}
		else if (curve.f[p][j-1] != 0.0) {
		  k = j - 1;
		  j = j - 2;
		}
	      }
	      else {
		if ( curve.f[p][j] == 0.0 ){
		  j = k;
		  k = k + 1;
		}
	      } 
	      } // finished eclipse section*/

	      // Check to see if we're near the maximum
	    if ( (i <= imax + 1 && i >= imax -1) || 
		 (imax == 0 && (i <= 1 || i == numbins-1)) || 
		 (imax == numbins-1 && (i==0 || i >= numbins-2))) { // parabolic interpolation near the maximum
	      if ( imax == numbins-1 && i == 0 ) {
		newflux.at(i) = maximum - temporary1 * pow(curve.t[i] + 1.0 - tx,2);  // intermediate value of flux
	      }
	      else {
		if ( imax == 0 && i == numbins-1 )
		  newflux.at(i) = maximum - temporary1 * pow(curve.t[i] - 1.0 - tx,2);
		else
		  newflux.at(i) = maximum - temporary1 * pow(curve.t[i] - tx,2);
	      }
	      //std::cout << "MAX i = " << i << " time = " << curve.t[i] << " flux  = " << newflux.at(i) << std::endl;
	    }
	            
	    // Not near the maximum
	    else {	            
	      	     
		double t1, t2;
		// Find point to the left of the "ith" point 
		// Time delays mean that j isn't always i-1

		// Check to see if t_o[0] > t_o[i]

		if ( curve.t_o[0] > curve.t[i]){ // if thing
		  if ( i==0){
		    if (curve.t_o[numbins-1] > 1.0){
		      j = numbins-2;
		      t1 = curve.t_o[j]-1.0;

		      k = numbins -1;
		      t2 = curve.t_o[k]-1.0;
		    }
		    else{
		      j = numbins-1;
		      t1 = curve.t_o[j]-1.0;
		      k = 0;
		      t2 = curve.t_o[k];
		    }
		  } // end if i==0
		  else{ //i>0
		    j = i-2;
		    if (j<0) j += numbins;
		    t1 = curve.t_o[j];
		    if (t1 > 1.0) t1 += -1.0;
		    k = j+1;
		    if (k>=numbins) k += -numbins;
		    t2 = curve.t_o[k];
		  } 
		}
		else { //Default
		  j = 0;
		  while ( (curve.t_o[j] <= curve.t[i]) && (j < static_cast<int>(numbins)) ) 
		    j++;
		  j--;
		  if ( j < 0 ) // because otherwise the computer gives us garbage and the flux looks ridiculous
		    j = numbins - abs(j);
		  t1 = curve.t_o[j]; // time to the left of the time we're interested in

		  if ( j == static_cast<int>(numbins) - 1 ) {
		    k = 0;
		  }
		  else { // else5; effectively, if i != 0 because that would make j != numbins-1
		    k = j + 1;
		    while (curve.t_o[k] <= curve.t[i] && k < static_cast<int>(numbins)) 
		      k++;
		  }
		  t2 = curve.t_o[k]; // time to the right of the point we're interested in
		  if (k==0) t2 += 1.0;

		}

		/*if (i==0){
		  std::cout << "bin 0 " << std::endl;
		  std::cout << "j=" << j << " t_o[j}=" << curve.t_o[j] << " t1 = " << t1 << std::endl;
		  std::cout << "k=" << k << " t_o[k}=" << curve.t_o[k] << " t2 = " << t2 << std::endl;
		  }*/

		newflux.at(i) = curve.f[p][j] + (curve.f[p][k]-curve.f[p][j])/(t2-t1) * (curve.t[i]-t1); // linear interpolation!

		if (curve.eclipse){
		  if (curve.t[i] < te1 && curve.t[i] > te1 - 1.0/numbins){
		    //std::cout << "time just before eclipse at te1=" << te1 <<"! i=" << i << " t = " << curve.t[i] << std::endl;
		    newflux.at(i) = slope1 * (curve.t[i] - te1);
		  }
		  if (curve.t[i] < te1 + 1.0/numbins && curve.t[i] > te1)
		    newflux.at(i) = 0.0;

		  if (curve.t[i] > te2 && curve.t[i] < te2 + 1.0/numbins){
		    // std::cout << "eclipse ends at te2=" << te2 <<"! i=" << i << " t = " << curve.t[i] << std::endl;
		    newflux.at(i) = slope2 * (curve.t[i] - te2);
		  }

		  if (curve.t[i] > te1 && curve.t[i] < te2)
		    newflux.at(i) = 0.0;
		}

		if (minimum != 0.0)
		  if ( fabs(newflux.at(i) - minimum)/minimum <= 0.000001) {

		    //std::cout << "near minimum!!!! " << std::endl;
		    if ( i == numbins-1 ) 
		      newflux.at(i) = minimum - temporary2 * pow(curve.t[i] - 1 - tmin,2); 
		    else {
		      newflux.at(i) = minimum - temporary2 * pow(curve.t[i] - tmin,2);
		    }
		    //std::cout << "MIN i = " << i << " time = " << curve.t[i] << " flux  = " << newflux.at(i) << std::endl;



		  }

		//std::cout << "i = " << i << " j = " << j << " k = " << k 
		//	  << "t=" << curve.t[i] << " f=" << newflux.at(i) 
		//	  <<std::endl;
		//else
		//newflux.at(i) = curve.f[p][i];

		// }// end else-not-near the minimum
	    }
	
	    /* newflux vs t_e corresponds to the new re-binned light curve.
	       It corresponds to the same light curve as bolflux vs t_o */
	    
	    if ( newflux.at(i) < 0.0 )
	      newflux.at(i) = 0.0;
	  
	  } // closes the for(i) loop, going through the curves.
    	
	  /****************************************************************/
	  /* SUMMING THE FLUXES AT EACH PHASE BIN ACROSS ALL LIGHT CURVES */
	  /****************************************************************/
    	
	  for ( unsigned int i(0); i < numbins; i++ ) {
            totflux.at(i) += newflux.at(i);   // setting totflux = newflux
	  }
	  //if (p==0)
	  //ttt.open("time.txt", std::ios_base::trunc);

	  for ( unsigned int i(0); i < numbins; i++ ) {
	    /*if ( p==0 ) 
	      ttt << curve.t_o[i] << " " << curve.f[p][i] << " " << curve.t[i] << " " << totflux.at(i) 
	    	  << " " << i
	    	  << std::endl;*/
	    
            curve.f[p][i] = totflux.at(i);
	  }
	  // only plotting versus evenly spaced time -- not using t_o
	}
      }// end for-p-loop
    } // end time delay section
    
    return curve;

} // end ComputeCurve


/**************************************************************************************/
/* ShiftCurve:                                                                      */
/*              Shifts curve by azimuthal angle phi                                 */
/*																					  */
/* pass: angles = all the angles necessary to compute the flux correctly;             */
/*                computed in the routine/method/function above [radians or unitless] */
/**************************************************************************************/
class LightCurve ShiftCurve( class LightCurve* angles, double phishift) {
		
    class LightCurve curve;

    std::ifstream input;
    std::ofstream output;
    std::ofstream ttt;

    unsigned int numbins(MAX_NUMBINS);  // Time bins of light curve (usually 128)
    unsigned int numbands(NCURVES);  // Number of Energy Bands
        
    std::vector< double > totflux(MAX_NUMBINS, 0.0); // integrated flux
    std::vector< bool > nullcurve(NCURVES,false); // true means that the curve is zero everywhere

    double timeshift;

    /*********************/
    /* SETTING THINGS UP */
    /*********************/

    curve = (*angles);
    numbins = curve.numbins;
    numbands = curve.numbands;

    timeshift = phishift/(2.0*Units::PI);

    //std::cout << "TimeShift = " << timeshift <<std::endl;

    // Shift light curve by angle phishift<(2pi)/numbins
  
    int k(0), j(0); // index placeholders; approximately, k is i+1 and j is i-1

        
      /********************************/
      /* LOOP THROUGH THE LIGHTCURVES */
      /********************************/
    	
    for (unsigned int p(0); p < numbands; p++) {
  
	if (!nullcurve[p]) {
	

	/*********************/
	/* MORE DECLARATIONS */
	/*********************/
    		
	  std::vector< double > newflux(MAX_NUMBINS, 0.0);                  // rebinned flux (to account for photon travel time)
	  double minimum(100000.0);                           // true (continuous) maximum and minimum flux values
	 
	  int ec1(0), ec2(0),j1,j2;
	  double te1, te2, tt1, tt2, slope1, slope2;

  	                   
	  /*******************************************/
	  /* FOR ECLIPSING LIGHT CURVES, MINIMUM = 0 */
	  /*******************************************/
    		
	  if ( curve.eclipse ) {
	    minimum = 0.0;

	    // Find the first eclipsed bin; set ec1=i;
	    for (unsigned int i(0);i<numbins;i++){
	      if (curve.f[p][i]==0.0){
		ec1=i;
		break;
	      }
	    }
	    
	    //	    std::cout << "eclipse begins at bin" << ec1 << std::endl;
	    // Extrapolate to find beginning of eclipse (linear extrapolation)
	    j1 = ec1-2;
	    j2 = ec1-1;

	    if (ec1==0){
	      j1+=numbins;
	      j2+=numbins;
	    }
	    if (ec1==1){
	      j1+=numbins;
	    }
	    tt1 = curve.t[j1];
	    if (ec1<2)
	      tt1 +=-1.0;
	    tt2 = curve.t[j2];
	    if (ec1==0)
	      tt2+=-1.0;
	    
	    slope1 = (curve.f[p][j1]-curve.f[p][j2])/(tt1-tt2);
	
	    if (slope1 != 0.0)
	      te1  = - (curve.f[p][j2])/slope1 + tt2 ;
	    else
	      te1 = tt1;


	    // After the eclipse, find the first bin with non-zero flux
	    for (unsigned int i(ec1);i<numbins;i++){
	      if (curve.f[p][i]>0.0){
		ec2=i;
		break;
	      }
	    }
	    //std::cout << "eclipse ends at bin" << ec2 << std::endl;
	    j1=ec2;
	    j2=ec2+1;
	    if (ec2==numbins-1){
	      j2+=-numbins;
	    }
	    tt1=curve.t[j1];
	    tt2=curve.t[j2];
	    if (ec2==numbins-1)
	      tt2+=1.0;

	    slope2 = (curve.f[p][j1]-curve.f[p][j2])/(tt1-tt2);
	    te2 = tt1 - curve.f[p][j1]/slope2;
	    //std::cout << "Eclipse ends at time = " << te2 << std::endl;
	  } // ending "yes eclipsed
            
		           	   
	  for ( unsigned int i(0); i < numbins; i++ ) {  // for-i-loop, looping through the phase bins
            
		double t1, t2;
		// Find point to the left of the "ith" point; call it j	       
		if ( i==0){		     
		  j = numbins-1;
		  t1 = curve.t[j]-1.0;
		  k = 0;
		  t2 = curve.t[k];		    
		} // end if i==0		 
		else { //Default
		  j = i-1;
		  t1 = curve.t[j]; // time to the left of the time we're interested in
		  k = i;
		  t2 = curve.t[k]; // time to the right of the point we're interested in
		}
		newflux.at(i) = curve.f[p][k] - (curve.f[p][k]-curve.f[p][j])/(t2-t1) * (timeshift); // linear interpolation!

		if (curve.eclipse){
		  if (curve.t[i] < te1 && curve.t[i] > te1 - 1.0/numbins){
		    //std::cout << "time just before eclipse at te1=" << te1 <<"! i=" << i << " t = " << curve.t[i] << std::endl;
		    newflux.at(i) = -slope1 * (timeshift) + curve.f[p][k];
		  }
		  if (curve.t[i] < te1 + 1.0/numbins && curve.t[i] > te1)
		    newflux.at(i) = 0.0;

		  if (curve.t[i] > te2 && curve.t[i] < te2 + 1.0/numbins){
		    // std::cout << "eclipse ends at te2=" << te2 <<"! i=" << i << " t = " << curve.t[i] << std::endl;
		    newflux.at(i) = -slope2 * (timeshift) + curve.f[p][k];
		  }

		  if (curve.t[i] > te1 && curve.t[i] < te2)
		    newflux.at(i) = 0.0;
		}

		
	  
	
		/* newflux vs t corresponds to the new shifted light curve.
		   It corresponds to the old curve is curve.f */
	    
		if ( newflux.at(i) < 0.0 )
		  newflux.at(i) = 0.0;
	  
	  } // closes the for(i) loop, going through the curves.
    	
	 
	  // if (p==0)
	  //ttt.open("shift.txt", std::ios_base::trunc);

	  for ( unsigned int i(0); i < numbins; i++ ) {
	    /*if ( p==0 ) 
	      ttt << curve.t[i] << " " << curve.f[p][i] << " " << curve.t[i] << " " << newflux.at(i) 
	    	  << " " << i
	    	  << std::endl;
	    */
            curve.f[p][i] = newflux.at(i);
	  }
	 

	} // end NOT NULL
    }// end for-p-loop
   
    
    return curve;

} // end ShiftCurve


/**************************************************************************************/
/* Normalize1:                                                                         */
/*           normalizes the fluxes in each energy band to 1 by dividing by 
	     the average flux of each light curve													  */
/*																					  */
/* pass: Flux = the flux for each part of the light curve [units don't matter]        */
/*       numbins = number of phase bins the light curve is divided into [Unitless]    */
/**************************************************************************************/
class LightCurve Normalize1( double Flux[NCURVES][MAX_NUMBINS], unsigned int numbins ) {

  /***************************************/
  /* VARIABLE DECLARATIONS FOR Normalize */
  /***************************************/
	
    class LightCurve newcurve;  // new curve structure to hold the normalized lightcurve flux values
    unsigned int i, p;          // loop variables
    double norm[NCURVES];       // normalization factor, one for each curve

    for ( p = 0; p < NCURVES; p++ ) {
        norm[p] = 0.0;
        for ( i = 0; i < numbins; i++ ) {
            norm[p] += Flux[p][i]; 
            newcurve.f[p][i] = Flux[p][i];
        }
        if ( norm[p] != 0.0) norm[p] /= (numbins*1.0); // makes norm the average value for each curve
    }


    for ( p = 0; p < NCURVES; p++ ) {
        for ( i = 0; i < numbins; i++ ) {
            if ( norm[p] != 0.0)  newcurve.f[p][i] /= norm[p]; 
            else newcurve.f[p][i] = 1.0;
        }
    }
    
    return newcurve;
} // end Normalize1

/**************************************************************************************/
/* Normalize2:                                                                         */
/*           normalizes the fluxes in low energy band to 1 by dividing by 
	     the average flux of the low energy light curve		       	  */
/*									       	  */
/* pass: Flux = the flux for each part of the light curve [units don't matter]        */
/*       numbins = number of phase bins the light curve is divided into [Unitless]    */
/**************************************************************************************/
class LightCurve Normalize2( double Flux[NCURVES][MAX_NUMBINS], unsigned int numbins ) {

  /***************************************/
  /* VARIABLE DECLARATIONS FOR Normalize */
  /***************************************/
	
    class LightCurve newcurve;  // new curve structure to hold the normalized lightcurve flux values
    unsigned int i, p;          // loop variables
    double norm[NCURVES];       // normalization factor, one for each curve

    for ( p = 0; p < NCURVES; p++ ) {
        norm[p] = 0.0;
        for ( i = 0; i < numbins; i++ ) {
            norm[p] += Flux[p][i]; 
            newcurve.f[p][i] = Flux[p][i];
        }
        if ( norm[p] != 0.0) norm[p] /= (numbins*1.0); // makes norm the average value for each curve
    }


    for ( p = 0; p < NCURVES; p++ ) {
        for ( i = 0; i < numbins; i++ ) {
            if ( norm[2] != 0.0)  newcurve.f[p][i] /= norm[2]; 
            else newcurve.f[p][i] = 1.0;
        }
    }
    
    return newcurve;
} // end Normalize2



#define TINY 1.0e-10
#define NMAX 300
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
#define NDIM 5
#define MPTS 6

/**************************************************************************************/
/* recalc: 																			  */
/*        this does not get called in Spot.cpp										  */
/*																					  */
/* pass: curve = flux in the light curve?                                             */
/*       omega = spin frequency of the NS [unitless]                                  */
/*       mass = mass of the NS [unitless]                                             */
/*       theta_0 = angle from spin axis to center of spot bit [radians]               */
/*       rspot = radius of the NS at the spot [unitless]                              */
/**************************************************************************************/
class OblDeflectionTOA* recalc( class LightCurve* curve, double omega, double mass, 
                                double theta_0, double rspot ) {

    /************************************/
   	/* VARIABLE DECLARATIONS FOR recalc */
    /************************************/
    
    double mu,        // = cos(theta_0), unitless
	       cosgamma,  // Cos of the angle between the radial vector and surface normal vector, defined in MLCB13
	       req,       // Radius of the neutron star at the equator, in km (but maybe unitless here)
	       zeta,      // Defined in MLCB9
	       epsilon,   // Doppler boost factor, defined in MLCB33
	       P2,        // Legendre polynomial
	       P4,        // Legendre polynomial
	       b0,        // Coefficient in MLCB Table 1
	       b2,        // Coefficient in MLCB Table 1
	       b4,        // Coefficient in MLCB Table 1
	       c2,        // Coefficient in MLCB Table 1
	       c4;        // Coefficient in MLCB Table 1
    
    double mass_over_r;
	           
	//G = 1.32746E+11;
	mu = cos( theta_0 );
	//    pi = 4.0 * atan(1.0);
	zeta = mass / rspot; // Since G = c = 1, with mass in gravitational units
    epsilon = pow( rspot * omega, 2.0 ) / zeta;

    b0 = 0.18*epsilon - 0.23*zeta*epsilon + 0.18*pow( epsilon, 2 );
    b2 = 0.39*epsilon - 0.29*zeta*epsilon + 0.42*pow( epsilon, 2 );
    b4 = -0.04*epsilon + 0.15*zeta*epsilon - 0.13*pow(epsilon, 2);
    //b4 = -8.0/3.0 * (b0 - 0.5*b2);

    c2 = 0.60 * pow( epsilon, 2 );
    c4 = -0.12 * pow( epsilon, 2 );

    P2 = LegP2( mu );
    P4 = LegP4( mu );

    req = rspot * (1.0 + b0 + b2*P2 + b4*P4 + P2*(c2*P2 + c4*P4));
	/*
	std::cout << "recalc: theta_0 = " << theta_0*180.0/Units::PI;
    std::cout << "rspot_recalc = " << rspot <<std::endl;
    std::cout << "req_recalc = " << req << std::endl;
    */
    OblModelBase* model;
    model = new PolyOblModelNHQS( rspot, req,
				PolyOblModelBase::zetaparam(mass,req),
				PolyOblModelBase::epsparam(omega, mass, req)
				);

    cosgamma = model->cos_gamma( mu );
    curve->para.cosgamma = cosgamma;
    // defltoa is a structure that "points" to routines in the file "OblDeflectionTOA.cpp"
    // used to compute deflection angles and times of arrivals 
    
    class OblDeflectionTOA* dt;

    mass_over_r = mass/rspot;

    dt = new OblDeflectionTOA( model, mass, mass_over_r, rspot );

    return (dt);

} // end recalc

/**************************************************************************************/
/* Blackbody:                                                                         */
/*           computes the monochromatic blackbody flux in units of erg/cm^2			  */
/*																					  */
/* pass: T = the temperature of the hot spot, in keV                                  */
/*       E = monochromatic energy in keV * redshift / eta                             */
/**************************************************************************************/
double BlackBody( double T, double E ) {   // Blackbody flux in units of erg/cm^2
    return ( 2.0e9 / pow(Units::C * Units::H_PLANCK, 2) * pow(E * Units::EV, 3) / (exp(E/T) - 1) ); // shouldn't it have a pi?
    // the e9 is to switch E from keV to eV; Units::EV gets it from eV to erg, since it's first computed in erg units.
    // the switch from erg units to photon count units happens above just after this is called.
} // end Blackbody


/**************************************************************************************/
/* Line:                                                                         */
/*           computes the monochromatic blackbody flux in units of erg/cm^2*/
/*           for a small range of emission frame energies. */
/*																					  */
/* pass: T = the temperature of the hot spot, in keV                                  */
/*       E = monochromatic energy in keV * redshift / eta                             */
/**************************************************************************************/
double Line( double T, double E , double E1, double E2) {   // Blackbody flux in units of erg/cm^2
  // If  E1 <= E <= E2 then blackbody flux. Zero otherwise  
  // E, E1, E2 are all as locally measured in the rest frame of the star
  if ( E <= E2 && E >= E1)
    return ( 2.0e9 / pow(Units::C * Units::H_PLANCK, 2) * pow(E * Units::EV, 3) / (exp(E/T) - 1) ); 

  else
    return 0.0;
   
} // end line

/**************************************************************************************/
/* EnergyBandFlux:                                                                    */
/*                computes the blackbody flux in units of erg/cm^2 using trapezoidal  */
/*                rule for approximating an integral                                  */
/*                variant of Bradt equation 6.6                                       */
/*                T, E1, E2 put into eV in this routine                               */
/*																					  */
/* pass: T = the temperature of the hot spot, in keV                                  */
/*       E1 = lower bound of energy band in keV * redshift / eta                      */
/*       E2 = upper bound of energy band in keV * redshift / eta                      */
/*       L1 = lower limit of emitted energy in star's frame                           */
/*       L2 = lower limit of emitted energy in star's frame                           */
/**************************************************************************************/
double LineBandFlux( double T, double E1, double E2, double L1, double L2 ) {
  T *= 1e3; // from keV to eV
  // x = E / T
  E1 *= 1e3; // from keV to eV
  E2 *= 1e3; // from keV to eV
  L1 *= 1e3;
  L2 *= 1e3;

  /********************************************/
  /* VARIABLE DECLARATIONS FOR LINEBandFlux */
  /********************************************/
	
	// a, b, x, n, h as defined by Mathematical Handbook eqn 15.16 (Trapezoidal rule to approximate definite integrals)
	double a = E1 / T;          // lower bound of integration
	double b = E2 / T;          // upper bound of integration
	double current_x(0.0);      // current value of x, at which we are evaluating the integrand; x = E / T; unitless
	unsigned int current_n(0);  // current step
	unsigned int n_steps(400); // total number of steps
	double h = (b - a) / n_steps;     // step amount for numerical integration; the size of each step
	double integral_constants = 2.0 * pow(T*Units::EV,3) / pow(Units::C,2) / pow(Units::H_PLANCK,3); // what comes before the integral when calculating flux using Bradt eqn 6.6 (in units of photons/cm^2/s)
	double flux(0.0);           // the resultant energy flux density; Bradt eqn 6.17
	
	// begin trapezoidal rule
	current_x = a + h * current_n;
	if ( L1/T <= current_x && current_x <= L2/T)
	  flux = Bradt_flux_integrand(current_x);

	for ( current_n = 1; current_n < n_steps-1; current_n++ ) {
		current_x = a + h * current_n;
		if ( L1/T <= current_x && current_x <= L2/T)
		  flux += 2.0 * Bradt_flux_integrand(current_x);
	}
	
	current_x = a + h * current_n;
	if ( L1/T <= current_x && current_x <= L2/T)
	  flux += Bradt_flux_integrand(current_x);
	
	flux *= h/2.0;	
	// end trapezoidal rule; numerical integration complete!
	
	flux *= integral_constants;

	return flux;
} // end LINEBandFlux

/**************************************************************************************/
/* EnergyBandFlux:                                                                    */
/*                computes the blackbody flux in units of erg/cm^2 using trapezoidal  */
/*                rule for approximating an integral                                  */
/*                variant of Bradt equation 6.6                                       */
/*                T, E1, E2 put into eV in this routine                               */
/*																					  */
/* pass: T = the temperature of the hot spot, in keV                                  */
/*       E1 = lower bound of energy band in keV * redshift / eta                      */
/*       E2 = upper bound of energy band in keV * redshift / eta                      */
/**************************************************************************************/
double EnergyBandFlux( double T, double E1, double E2 ) {
	T *= 1e3; // from keV to eV
	// x = E / T
	E1 *= 1e3; // from keV to eV
	E2 *= 1e3; // from keV to eV
	
	/********************************************/
   	/* VARIABLE DECLARATIONS FOR EnergyBandFlux */
    /********************************************/
	
	// a, b, x, n, h as defined by Mathematical Handbook eqn 15.16 (Trapezoidal rule to approximate definite integrals)
	double a = E1 / T;          // lower bound of integration
	double b = E2 / T;          // upper bound of integration
	double current_x(0.0);      // current value of x, at which we are evaluating the integrand; x = E / T; unitless
	unsigned int current_n(0);  // current step
	unsigned int n_steps(3000); // total number of steps
	double h = (b - a) / n_steps;     // step amount for numerical integration; the size of each step
	double integral_constants = 2.0 * pow(T*Units::EV,3) / pow(Units::C,2) / pow(Units::H_PLANCK,3); // what comes before the integral when calculating flux using Bradt eqn 6.6 (in units of photons/cm^2/s)
	double flux(0.0);           // the resultant energy flux density; Bradt eqn 6.17
	
	// begin trapezoidal rule
	current_x = a + h * current_n;
	flux = Bradt_flux_integrand(current_x);

	for ( current_n = 1; current_n < n_steps-1; current_n++ ) {
		current_x = a + h * current_n;
		flux += 2.0 * Bradt_flux_integrand(current_x);
	}
	
	current_x = a + h * current_n;
	flux += Bradt_flux_integrand(current_x);
	
	flux *= h/2.0;	
	// end trapezoidal rule; numerical integration complete!
	
	flux *= integral_constants;

	return flux;
} // end EnergyBandFlux

/**************************************************************************************/
/* Bradt_flux_integrand:                                                              */
/*                      integrand of Bradt eqn 6.6 when integrating over nu, modified */
/*                      so the exponent is 2 not 3, so that it comes out as photon    */
/*                      number flux, instead of erg flux                              */
/*																					  */
/* pass: x = current_x from above routine                                             */
/**************************************************************************************/
double Bradt_flux_integrand( double x ) {
	return ( pow(x,2) / (exp(x) - 1) );  // 2 (not 3) for photon number flux
} // end Bradt_flux_integrand

/**************************************************************************************/
/* Gray:																			  */
/*		computes and returns the limb darkening factors for a Gray electron -         */
/*		scattering atmosphere (Hopf function) using values from Chandrasekhar's       */
/*		"Radiative Transfer", table 13												  */
/*																					  */
/* pass: cosine = curve.cosbeta[i] * curve.eta[i]                                     */
/**************************************************************************************/
double Gray( double cosine ) {
	
	/**********************************/
   	/* VARIABLE DECLARATIONS FOR Gray */
    /**********************************/
	
    double F[11],   // limb darkening values for electron-scattering atmosphere
           mu[11];  // a table of cos(theta) values
    int i;          // loop variable

    for ( i = 0; i <= 10; i++ ) {
        mu[i] = 0.0 + i*0.1;
    }
    /*
    // Value from Mihalas's "Stellar Atmospheres"
    F[0] = 0.4330;
    F[1] = 0.5401;
    F[2] = 0.6280;
    F[3] = 0.7112;
    F[4] = 0.7921;
    F[5] = 0.8716;
    F[6] = 0.9501;
    F[7] = 1.0280;
    F[8] = 1.1053;
    F[9] = 1.1824;
    F[10] = 1.2591;
    */
  
    // Values from Chandrasekhar's "Radiative Transfer", Table 13
    F[0] = 1.0;
    F[1] = 1.2647;
    F[2] = 1.48009;
    F[3] = 1.68355;
    F[4] = 1.88105;
    F[5] = 2.07496;
    F[6] = 2.2665;
    F[7] = 2.45639;
    F[8] = 2.64503;
    F[9] = 2.83274;
    F[10] = 3.01973;
  
    for ( i = 0; i <= 10; i++ ) {
        F[i] *= (0.5/1.194);
    }

    int index(0);
    //i=1;

    if (cosine == 1.0)
    	index = 10;
    else {
        while ( cosine >= mu[index]) 
        	index++;
    }
    // cosine lies in the range   mu[index-1] <= cosine <= mu[index]

    if (cosine == mu[index-1]) 
    	return F[index-1];
    else {
        if (cosine == mu[index]) 
        	return F[index];
        else {   
            return ( F[index-1] + (cosine - mu[index-1]) * (F[index] - F[index-1])/(mu[index]-mu[index-1]) );
        }
    }
} // end Gray


/**************************************************************************************/
/* LegP2:																			  */
/*       calculating a Legendre polynomial P2										  */
/*       called in the routine/method/function recalc								  */
/*																					  */
/* pass: costheta = cos of angle of the hotspot down from geom. north pole [unitless] */
/**************************************************************************************/
double LegP2( double costheta ) {
    double P2;
    P2 = 0.5 * ( 3.0 * pow( costheta, 2 ) - 1.0 );
    return P2;
} // end LegP2

/**************************************************************************************/
/* LegP4:																			  */
/*       calculating a Legendre polynomial P4										  */
/*       called in the routine/method/function recalc								  */
/*																					  */
/* pass: costheta = cos of angle of the hotspot down from geom. north pole [unitless] */
/**************************************************************************************/
double LegP4( double costheta ) {
    double P4;
    P4 = 0.125 * ( 35 * pow( costheta, 4 ) - 30 * pow( costheta, 2 ) + 3);
    return P4;
} // end LegP4

/**************************************************************************************/
/* calcreq: 																		  */
/*		   calculating the radius of the NS at its equator                            */
/*		   for oblate stars, useful; for spherical stars, redundant                   */
/*																					  */		   
/* pass: omega = spin frequency of the NS [unitless]                                  */
/*       mass = mass of the NS [unitless]											  */
/*       theta_0 = angle of the hotspot down from the geometric north pole [radians]  */
/*       rspot = radius of the NS at the hotspot [unitless]							  */
/**************************************************************************************/
double calcreq( double omega, double mass, double theta_0, double rspot ) {
   
    /*************************************/
   	/* VARIABLE DECLARATIONS FOR calcreq */
    /*************************************/
   
	double mu,       // = cos(theta_0)
           req,      // Radius of the neutron star at the equator
           pi,       // PI, the constant (3.1415926...)
           zeta,     // Defined in MLCB9
           epsilon,  // Doppler boost factor, defined in MLCB33
           G,        // Gravitational constant
           P2,       // Legendre polynomial
           P4,       // Legendre polynomial
           b0,       // Coefficient from MLCB Table 1
           b2,       // Coefficient from MLCB Table 1
           b4,       // Coefficient from MLCB Table 1
           c2,       // Coefficient from MLCB Table 1
           c4;       // Coefficient from MLCB Table 1
	    
	G = 1.32746E+11;
    mu = cos(theta_0);
    pi = 4.0 * atan(1.0);
    zeta = mass / rspot;
    epsilon = pow( rspot * omega, 2.0 ) / zeta;

    b0 = 0.18*epsilon - 0.23*zeta*epsilon + 0.18*pow(epsilon,2);
    b2 = 0.39*epsilon - 0.29*zeta*epsilon + 0.42*pow(epsilon,2);
    	    b4 = -0.04*epsilon + 0.15*zeta*epsilon - 0.13*pow(epsilon,2);
	    //b4 = -8.0/3.0 * (b0 - 0.5*b2);

    c2 = 0.60*pow(epsilon,2);
    c4 = -0.12*pow(epsilon,2);
    
    P2 = LegP2( mu );
    P4 = LegP4( mu );

    req = rspot * ( 1.0 + b0 + b2*P2 + b4*P4 + P2*(c2*P2 + c4*P4) );
    	    
    //std::cout << "rspot_calcreq = " << rspot <<std::endl;
    //std::cout << "req_calcreq = " << req <<std::endl;
	  
	return req;
} // end calcreq

// end Chi.cpp

