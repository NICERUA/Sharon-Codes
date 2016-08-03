/***************************************************************************************/
/*                               OblDeflectionTOA.cpp

    This code calculates the relationships between impact parameter (b), deflection angle,
    and ray time of arrival (TOA), given for an oblate NS model.
    
    Based on code written by Coire Cadeau and modified by Sharon Morsink and 
    Abigail Stevens.
    
    (C) Coire Cadeau, 2007; Source (C) Coire Cadeau 2007, all rights reserved.
    Permission is granted for private use only, and not distribution, either verbatim or
    of derivative works, in whole or in part.
    This code is not thoroughly tested or guaranteed for any particular use.
*/
/***************************************************************************************/

#include "matpack.h"

#include <exception>
#include <cmath>
#include <iostream>
#include "OblDeflectionTOA.h"
#include "OblModelBase.h"
#include "Exception.h"
#include "Units.h"
#include "Struct.h" //Jan 21 (year? prior to 2012)
// Globals are bad, but there's no easy way to get around it
// for this particular case (need to pass a member function
// to a Matpack routine which does not have a signature to accomodate
// passing in the object pointer).

class LightCurve curve; //Jan 21 (year? prior to 2012)

const OblDeflectionTOA* OblDeflectionTOA_object;
double OblDeflectionTOA_b_value;
double OblDeflectionTOA_b_over_r;
double OblDeflectionTOA_costheta_value;
double OblDeflectionTOA_psi_value;
double OblDeflectionTOA_b_max_value;
double OblDeflectionTOA_psi_max_value;
double OblDeflectionTOA_b_guess;
double OblDeflectionTOA_psi_guess;

/**********************************************************/
/* OblDeflectionTOA_psi_integrand_wrapper:                */
/*                                                        */
/* Passed r (NS radius, unitless), *prob (problem toggle) */
/* Returns integrand if it is not nan.                    */
/**********************************************************/
double OblDeflectionTOA_psi_integrand_wrapper ( double r,  bool *prob ) {
  double integrand( OblDeflectionTOA_object->psi_integrand( OblDeflectionTOA_b_value, r)); 
  //OblDeflectionTOA_object is looking up psi_integrand by feeding it OblDeflectionTOA_b_value, which is defined in psi_outgoing
  	if( std::isnan(integrand) ) {
    	std::cout << "integrand is nan!" << std::endl;
    	std::cout << "r (km) = " << r  << std::endl;
    	std::cerr << "ERROR in OblDeflectionTOA_psi_integrand_wrapper(): returned NaN." << std::endl;
    	*prob = true;
    	integrand = -7888.0;
    	return integrand;
  	}
  	return integrand;
}

double OblDeflectionTOA_psi_integrand_wrapper_u ( double u,  bool *prob ) {
  double integrand( OblDeflectionTOA_object->psi_integrand_u( OblDeflectionTOA_b_over_r, u)); 
  //OblDeflectionTOA_object is looking up psi_integrand by feeding it OblDeflectionTOA_b_over_r, which is defined in psi_outgoing
  	if( std::isnan(integrand) ) {
    	std::cout << "integrand is nan!" << std::endl;
    	std::cout << "u = " << u  << std::endl;
    	std::cerr << "ERROR in OblDeflectionTOA_psi_integrand_wrapper_u(): returned NaN." << std::endl;
    	*prob = true;
    	integrand = -7888.0;
    	return integrand;
  	}
  	return integrand;
}


/************************************************************/
/* OblDeflectionTOA_dpsi_db_integrand_wrapper:              */
/*                                                          */
/* Passed r (NS radius, unitless), *prob (problem toggle)   */
/* Returns integrand if it is not nan.                      */
/************************************************************/
double OblDeflectionTOA_dpsi_db_integrand_wrapper ( double r, bool *prob ) {
	double integrand( OblDeflectionTOA_object->dpsi_db_integrand( OblDeflectionTOA_b_value, r ) );
  
  	if ( std::isnan(integrand) ) {
    	std::cout << "dpsi_db integrand is nan!" << std::endl;
    	std::cout << "r (km) = " << Units::nounits_to_cgs(r, Units::LENGTH)/1.0e5 << std::endl;
    	std::cerr << "ERROR in OblDeflectionTOA_dpsi_db_integrand_wrapper(): returned NaN." << std::endl;
    	*prob = true;
    	integrand = -7888.0;
    	return integrand;
  	}
  
  	return integrand;
}

/************************************************************/
/* OblDeflectionTOA_dpsi_db_integrand_wrapper:              */
/*                                                          */
/* Passed r (NS radius, unitless), *prob (problem toggle)   */
/* Returns integrand if it is not nan.                      */
/************************************************************/
double OblDeflectionTOA_dpsi_db_integrand_wrapper_u ( double u, bool *prob ) {
	double integrand( OblDeflectionTOA_object->dpsi_db_integrand_u( OblDeflectionTOA_b_over_r, u ) );
  
  	if ( std::isnan(integrand) ) {
    	std::cout << "dpsi_db_u integrand is nan!" << std::endl;
    	std::cout << "u (km) = " << Units::nounits_to_cgs(u, Units::LENGTH)/1.0e5 << std::endl;
    	std::cerr << "ERROR in OblDeflectionTOA_dpsi_db_integrand_wrapper(): returned NaN." << std::endl;
    	*prob = true;
    	integrand = -7888.0;
    	return integrand;
  	}
  
  	return integrand;
}

/***********************************************************/
/* OblDeflectionTOA_toa_integrand_wrapper:                 */
/*                                                         */
/* Passed r (NS radius, unitless), *prob (problem toggle)  */
/* Returns integrand if it is not nan.                     */
/***********************************************************/
double OblDeflectionTOA_toa_integrand_wrapper ( double r, bool *prob ) {
  	double integrand( OblDeflectionTOA_object->toa_integrand( OblDeflectionTOA_b_value, r ) );
  	if ( std::isnan(integrand) ) {
    	std::cout << "toa_minus_toa integrand is nan!" << std::endl;
    	std::cout << "r (km) = " << Units::nounits_to_cgs(r, Units::LENGTH)/1.0e5 << std::endl;
    	std::cerr << "ERROR in OblDeflectionTOA_toa_integrand_wrapper(): returned NaN." << std::endl;
    	*prob = true;
    	integrand = -7888.0;
    	return integrand;
  	}
  	return integrand;
}

/***********************************************************/
/* OblDeflectionTOA_toa_integrand_minus_b0_wrapper:        */
/*                                                         */
/* Passed r (NS radius, unitless), *prob (problem toggle)  */
/* Returns integrand if it is not nan.                     */
/***********************************************************/
double OblDeflectionTOA_toa_integrand_minus_b0_wrapper ( double r, bool *prob ) {
  	double integrand( OblDeflectionTOA_object->toa_integrand_minus_b0( OblDeflectionTOA_b_value, r ) );
  	if ( std::isnan(integrand) ) {
    	std::cout << "toa_minus_b0 integrand is nan!" << std::endl;
    	std::cout << "r (km) = " << Units::nounits_to_cgs(r, Units::LENGTH)/1.0e5 << std::endl;
    	std::cerr << "ERROR in OblDeflectionTOA_toa_integrand_minus_b0_wrapper(): returned NaN." << std::endl;
    	*prob = true;
    	integrand = -7888.0;
    	return integrand;
  	}
  	return integrand;
}

/***********************************************************/
/* OblDeflectionTOA_toa_integrand_minus_b0_wrapper:        */
/*                                                         */
/* Passed r (NS radius, unitless), *prob (problem toggle)  */
/* Returns integrand if it is not nan.                     */
/***********************************************************/
double OblDeflectionTOA_toa_integrand_minus_b0_wrapper_u ( double u, bool *prob ) {
  	double integrand( OblDeflectionTOA_object->toa_integrand_minus_b0_u( OblDeflectionTOA_b_over_r, u ) );
  	if ( std::isnan(integrand) ) {
    	std::cout << "toa_minus_b0_u integrand is nan!" << std::endl;
    	std::cout << "r (km) = " << Units::nounits_to_cgs(u, Units::LENGTH)/1.0e5 << std::endl;
    	std::cerr << "ERROR in OblDeflectionTOA_toa_integrand_minus_b0_wrapper(): returned NaN." << std::endl;
    	*prob = true;
    	integrand = -7888.0;
    	return integrand;
  	}
  	return integrand;
}

/*****************************************************/
/*****************************************************/
double OblDeflectionTOA_rcrit_zero_func_wrapper ( double rc ) {
  	return double ( OblDeflectionTOA_object->rcrit_zero_func( rc, 
  	                	OblDeflectionTOA_b_value ) );
}

/***********************************************************/
/* OblDeflectionTOA_b_from_psi_ingoing_zero_func_wrapper:  */
/*                                                         */
/* Passed b (photon impact parameter)                      */
/* Returns                      */
/***********************************************************/
double OblDeflectionTOA_b_from_psi_ingoing_zero_func_wrapper ( double b ) {
  	return double ( OblDeflectionTOA_object->b_from_psi_ingoing_zero_func( b, 
						OblDeflectionTOA_costheta_value,
						OblDeflectionTOA_psi_value ) );
}

/*****************************************************/
/*****************************************************/
double OblDeflectionTOA_b_from_psi_outgoing_zero_func_wrapper ( double b ) {
  // std::cerr << "b = " << Units::nounits_to_cgs(b, Units::LENGTH)/1.0e5 <<std::endl;
	return double ( OblDeflectionTOA_object->b_from_psi_outgoing_zero_func( b, 
						OblDeflectionTOA_costheta_value,
						OblDeflectionTOA_psi_value,
						OblDeflectionTOA_b_max_value,
						OblDeflectionTOA_psi_max_value,
						OblDeflectionTOA_b_guess,
						OblDeflectionTOA_psi_guess ) );
}

// End of global pollution.


// Defining constants
//const double OblDeflectionTOA::INTEGRAL_EPS = 1.0e-7;
const double OblDeflectionTOA::FINDZERO_EPS = 1.0e-6;
const double OblDeflectionTOA::RFINAL_MASS_MULTIPLE = 1.0e7;
//const double OblDeflectionTOA::DIVERGENCE_GUARD = 2.0e-2; // set to 0 to turn off
const double OblDeflectionTOA::DIVERGENCE_GUARD = 0.0;
const long int OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_N_MAX = 100000;
const long int OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_N_MAX_1 = 1000;
const long int OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_N = 1000;
const long int OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_N_1 = 1000;
const long int OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_INGOING_N = 400;
const double OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_POWER = 4.0; 
//const double OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_POWER = 8.0;

OblDeflectionTOA::OblDeflectionTOA ( OblModelBase* modptr, const double& mass_nounits, const double& mass_over_r_nounits, const double& radius_nounits ) 
  : r_final( RFINAL_MASS_MULTIPLE * mass_nounits ) {
  this->modptr = modptr;
  mass = mass_nounits;
  mass_over_r = mass_over_r_nounits;
  rspot = radius_nounits;
}

/************************************************************/
/* OblDeflectionTOA::bmax_outgoing                          */
/*														    */
/* Computes the maximal value of b for outgoing light rays  */
/* Passed rspot (radius of NS at spot, unitless)            */
/* Returns maximal value of b (in funny unitless units!)    */
/************************************************************/
double OblDeflectionTOA::bmax_outgoing ( const double& rspot ) const {

  //double bmax_outgoing ( rspot / sqrt( 1.0 - 2.0 * get_mass() / rspot ) );
  double bmax_outgoing ( rspot / sqrt( 1.0 - 2.0 * get_mass_over_r() ) );
  	return bmax_outgoing;
}

/*******************************************************************************/
/* OblDeflectionTOA::bmin_ingoing                                              */
/*														                       */
/* Computes value of b corresponding to the most ingoing light ray allowed,    */
/*   i.e. the ray that just grazes the surface and moves toward axis with b=0  */
/* Passed cos_theta (cosine of theta (what is theta here?))                    */
/* Returns minimum value of b allowed                                          */
/*******************************************************************************/
double OblDeflectionTOA::bmin_ingoing ( const double& rspot, const double& cos_theta ) const { 
  //	costheta_check( cos_theta );

  	// Need to return the value of b corresponding to dR/dtheta of the surface.

  	double drdth ( modptr->Dtheta_R( cos_theta ) );
	double b ( rspot / sqrt( (1.0 - 2.0 * get_mass_over_r()) + pow( drdth / rspot , 2.0 ) ) );
 
  	// NaN check:
  	if ( std::isnan(b) ) {
    	std::cerr << "ERROR in OblDeflectionTOA::bmin_ingoing(): returned NaN." << std::endl;
    	b = -7888.0;
    	return b;
  	}

  	return b;    // returning b!!!
}

/*************************************************************/
/* OblDeflectionTOA::ingoing_allowed                         */
/*														     */
/* Determines if an ingoing light ray is allowed             */
/*   for locations where dR/dtheta \neq 0                    */
/* Passed cos_theta (cosine of theta (what is theta here?))  */
/* Returns boolean (true if allowed)                         */
/*************************************************************/
bool OblDeflectionTOA::ingoing_allowed ( const double& cos_theta ) {
  	// ingoing rays can be a consideration for locations where dR/dtheta \neq 0
  	//costheta_check( cos_theta );
  	return bool ( modptr->Dtheta_R( cos_theta ) != 0.0 );
}

/************************************************************************/
/* OblDeflectionTOA::rcrit                                              */
/*														                */
/* Determines if an ingoing light ray is allowed                        */
/*   for locations where dR/dtheta \neq 0                               */
/* Passed b (photon impact parameter), cos_theta (cos of spot angle?),  */
/*        prob (problem toggle)                                         */
/* Returns                                    */
/************************************************************************/
double OblDeflectionTOA::rcrit ( const double& b, const double& cos_theta, bool *prob ) const {
  	double candidate;
       
  	double r ( modptr->R_at_costheta( cos_theta ) );
 
  	if ( b == bmax_outgoing(r) ) {
    	return double(r);
  	}

  	OblDeflectionTOA_object = this;
  	OblDeflectionTOA_b_value = b;

 	candidate = MATPACK::FindZero(modptr->R_at_costheta(1.0),
				modptr->R_at_costheta(0.0),
				OblDeflectionTOA_rcrit_zero_func_wrapper); // rcrit_guess

 
  	return candidate;
}

/*****************************************************/
/*****************************************************/
double OblDeflectionTOA::psi_outgoing ( const double& b, const double& rspot,
				       					const double& b_max, const double& psi_max, 
				       					bool *prob ) const{

  	double dummy;

  	if ( b > b_max || b < 0.0 ) {
    	std::cerr << "ERROR in OblDeflectionTOA::psi_outgoing(): b out-of-range." << std::endl;
    	*prob = true;
    	dummy = -7888.0;
    	return dummy;
  	}

  	if ( fabs( b - b_max ) < 1e-7 ) { // essentially the same as b=b_max
    	return psi_max;
	}
  	else {  // set globals
    	OblDeflectionTOA_object = this;
    	OblDeflectionTOA_b_value = b;     // this will be fed to another function

    	double psi(0.0);

    	if ( b != 0.0 ) {
      		psi = Integration( rspot, get_rfinal(), OblDeflectionTOA_psi_integrand_wrapper ); // integrating from r_surf to ~infinity
    	}				
    	return psi;
  	}
}

/*****************************************************/
double OblDeflectionTOA::psi_outgoing_u ( const double& b, const double& rspot,
				       					const double& b_max, const double& psi_max, 
				       					bool *prob ) const{

  	double dummy;

  	if ( b > b_max || b < 0.0 ) {
    	std::cerr << "ERROR in OblDeflectionTOA::psi_outgoing_u(): b out-of-range." << std::endl;
    	*prob = true;
    	dummy = -7888.0;
    	return dummy;
  	}

  	if ( fabs( b - b_max ) < 1e-7 ) { // essentially the same as b=b_max
    	return psi_max;
	}
  	else {  // set globals
    	OblDeflectionTOA_object = this;
    	OblDeflectionTOA_b_over_r = b/rspot;     // this will be fed to another function

    	double psi(0.0);

    	if ( b != 0.0 ) {
	  //psi = Integration( 0.0, 1.0, OblDeflectionTOA_psi_integrand_wrapper_u ); // integrating from r_surf to ~infinity
	  if ( b > 0.995*b_max){
	    double split(0.95);  	
	    psi = Integration( 0.0, split, OblDeflectionTOA_psi_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N );
	    psi += Integration( split, 1.0, OblDeflectionTOA_psi_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N_1 );
	  }
	  else
	    psi = Integration( 0.0, 1.0, OblDeflectionTOA_psi_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N );
    	}				
    	return psi;
  	}
}

/*****************************************************/
// Compute the largest outgoing deflection angle
/*****************************************************/
double OblDeflectionTOA::psi_max_outgoing ( const double& b, const double& rspot, bool *prob ) {
  	double dummy;
  	if ( b > bmax_outgoing(rspot) || b < 0.0 ) {
    	std::cerr << "ERROR in OblDeflectionTOA::psi_outgoing(): b out-of-range." << std::endl;
    	*prob = true;
    	dummy = -7888.0;
    	return dummy;
  	}

  	// set globals
  	OblDeflectionTOA_object = this;
  	OblDeflectionTOA_b_value = b;
  	
  	double psi = Integration( rspot, get_rfinal(), OblDeflectionTOA_psi_integrand_wrapper );
					
	std::cout << "Psi_max: b/r = " << b/rspot << " rspot = " << rspot << " r_final = " << get_rfinal() << std::endl;
	std::cout << "psi = " << psi << std::endl;

  	return psi;
}

double OblDeflectionTOA::psi_max_outgoing_u ( const double& b, const double& rspot, bool *prob ) const{
  	double dummy;
  	if ( b > bmax_outgoing(rspot) || b < 0.0 ) {
    	std::cerr << "ERROR in OblDeflectionTOA::psi_outgoing(): b out-of-range." << std::endl;
    	*prob = true;
    	dummy = -7888.0;
    	return dummy;
  	}

  	// set globals
  	OblDeflectionTOA_object = this;
  	OblDeflectionTOA_b_over_r = b/rspot;

	double split(0.9);  	
  	double psi = Integration( 0.0, split, OblDeflectionTOA_psi_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N_MAX_1 );
	psi += Integration( split, 1.0, OblDeflectionTOA_psi_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N_MAX );
					
	std::cout << "Psi_max_u: b/r = " << b/rspot << " rspot = " << rspot << " r_final = " << get_rfinal() << std::endl;
	std::cout << "psi = " << psi << std::endl;

  	return psi;
}

/*****************************************************/
/*****************************************************/
double OblDeflectionTOA::psi_ingoing ( const double& b, const double& cos_theta, bool *prob ) const {
  //costheta_check( cos_theta );
  
  	double rspot ( modptr->R_at_costheta( cos_theta ) );
  	//double dummy;
  	// std::cout << "psi_ingoing: cos_theta = " << cos_theta << " b = " << b << " r = " << rspot << std::endl;
	/*
  	if ( b > bmax_outgoing(rspot) || b < bmin_ingoing( rspot, cos_theta ) ) {
    	std::cerr << "ERROR in OblDeflectionTOA::psi_ingoing(): b out-of-range." << std::endl;
    	*prob = true;
    	dummy = -7888.0;
    	return dummy;
 	}
	*/
  	// set globals
  	OblDeflectionTOA_object = this;
  	OblDeflectionTOA_b_value = b;

 
  	// See psi_outgoing. Use an approximate formula for the integral near rcrit, and the real formula elsewhere
  
  	double rcrit = this->rcrit( b, cos_theta, &curve.problem );

  	double psi_in = 2.0 * sqrt( 2.0 * (rspot - rcrit) / (rcrit - 3.0 * get_mass()) );

  	return double( psi_in + this->psi_outgoing( b, rspot, 100.0, 100.0, &curve.problem ) );
}

/*****************************************************/
// Major Change to b_from_psi
// a guess for b is given to this routine
// if the photon is purely outgoing nothing is changed
// if the photon is initially ingoing a calculation is done
/*****************************************************/
bool OblDeflectionTOA::b_from_psi ( const double& psi, const double& rspot, const double& cos_theta, double& b,
									int& rdot, const double& bmax_out, 
									const double& psi_out_max, const double& b_guess, 
									const double& psi_guess, const double& b2, 
									const double& psi2, bool *prob ) {
  //costheta_check( cos_theta );

  // return bool indicating whether a solution was found.
  // store result in b
  // store -1 in rdot if its an ingoing solution, +1 otherwise.

  	double bcand;    //
  	double dummyb;   //
  	double dummypsi; //

    OblDeflectionTOA_b_max_value = bmax_out;
    OblDeflectionTOA_psi_max_value = psi_out_max;
    OblDeflectionTOA_b_guess = b_guess;
    OblDeflectionTOA_psi_guess = psi_guess;

  	if ( std::isinf(psi_out_max) ) {
    	std::cerr << "ERROR in OblDeflectionTOA::b_from_psi(): psi_out_max = infinity" << std::endl;
    	*prob = true;
    	dummyb = -7888.0;
    	return dummyb;
  	}

  	if ( psi < 0 ) {
    	std::cerr << "ERROR in OblDeflectionTOA::b_from_psi: need psi >= 0" << std::endl;
    	*prob = true;
    	dummypsi = -7888.0;
    	return dummypsi;
  	}
  	else if ( psi == 0.0 ) {
    	b = 0.0;
    	rdot = 1;
    	return true;
  	}
  	else if ( psi == psi_out_max ) {
    	b = bmax_out;
    	rdot = 1;
    	return true;
  	}
  	else if ( psi < psi_out_max ) { // begin normal outgoing case, psi < psi_out_max
    	OblDeflectionTOA_object = this;
    	OblDeflectionTOA_costheta_value = cos_theta;
    	OblDeflectionTOA_psi_value = psi;

    	bcand = b_guess;

 

    	if ( fabs(bcand) <= std::numeric_limits<double>::epsilon() || fabs(bmax_out - bcand) <= std::numeric_limits<double>::epsilon() ) { // this indicates no soln
      		std::cerr << "ERROR in OblDeflectionTOA::b_from_psi(): outgoing returned no solution?" << std::endl;
      		*prob = true;
      		b = -7888.0;
      		return b;
    	}
    	else {
      		b = bcand;
      		rdot = 1;
      		return true;
    	}
  	} // end normal outgoing case, psi < psi_out_max
  	
  	else { // psi > psi_out_max 
    	// Test to see if ingoing photons are allowed
    	bool ingoing_allowed( this->ingoing_allowed(cos_theta) );
    	double bmin_in,    //
    	       psi_in_max; //

    	//std::cout << "b_from_psi: psi > psi_out_max" << std::endl;

    	if ( ingoing_allowed ) {
	  bmin_in = bmin_ingoing( rspot, cos_theta );
       
	  psi_in_max = psi_ingoing( bmin_in, cos_theta, &curve.problem );

	  if ( psi > psi_in_max ) {
	    return false;
	  }
	  else if ( psi == psi_in_max ) {
	    b = bmin_in;
	    rdot = -1;
	    return true;
	  }
	  else { // psi_out_max < psi < psi_in_max
	    OblDeflectionTOA_object = this;
	    OblDeflectionTOA_costheta_value = cos_theta;
	    OblDeflectionTOA_psi_value = psi;

	    bcand = MATPACK::FindZero(bmin_in, bmax_out,
				      OblDeflectionTOA_b_from_psi_ingoing_zero_func_wrapper,
				      OblDeflectionTOA::FINDZERO_EPS);
 
	    if( fabs(bmin_in - bcand) <= std::numeric_limits<double>::epsilon() || fabs(bmax_out - bcand) <= std::numeric_limits<double>::epsilon() ) { // this indicates no soln
	      std::cerr << "ERROR in OblDeflectionTOA::b_from_psi(): ingoing returned no solution?" << std::endl;
	      *prob = true;
	      b = -7888.0;
	      return b;
	    }
	    else {
	      b = bcand;
	      rdot = -1;
	      return true;
	    }
	  }
    	} // end ingoing photons are allowed
    	else { // ingoing photons not allowed.
      		//std::cerr << "Something is really wrong in b-from-psi!!" << std::endl;
      		return false;
    	}
  	}
  	std::cerr << "ERROR in OblDeflectionTOA::b_from_psi(): reached end of function?" << std::endl;
  	*prob = true;
  	return false;
}

/*****************************************************/
// warning: this integral will diverge near bmax_out. You can't remove
  	// the divergence in the same manner as the one for the deflection.
  	// check the output!
/*****************************************************/
double OblDeflectionTOA::dpsi_db_outgoing( const double& b, const double& rspot, bool *prob ) {
  	//costheta_check(cos_theta);
  
  	double dummy; //

  	if ( (b > bmax_outgoing(rspot) || b < 0.0) ) { 
    	std::cerr << "ERROR inOblDeflectionTOA::dpsi_db_outgoing(): b out-of-range." << std::endl;
    	std::cerr << "b = " << b << ", bmax = " << bmax_outgoing(rspot) << ", |b - bmax| = " << fabs(b - bmax_outgoing(rspot)) << std::endl;
    	*prob = true;
    	dummy = -7888.0;
    	return dummy;
  	}
  
  	// set globals
  	OblDeflectionTOA_object = this;
  	OblDeflectionTOA_b_value = b;

  	//double rsurf = modptr->R_at_costheta(cos_theta);

  	//rsurf = rspot;

  	double dpsidb = Integration( rspot, get_rfinal(), OblDeflectionTOA_dpsi_db_integrand_wrapper );
  	
  	return dpsidb;
}

/*****************************************************/
// warning: this integral will diverge near bmax_out. You can't remove
  	// the divergence in the same manner as the one for the deflection.
  	// check the output!
/*****************************************************/
double OblDeflectionTOA::dpsi_db_outgoing_u( const double& b, const double& rspot, bool *prob ) const {
  	//costheta_check(cos_theta);
  
  	double dummy; //

	double b_max = bmax_outgoing(rspot);

  	if ( (b > b_max || b < 0.0) ) { 
    	std::cerr << "ERROR inOblDeflectionTOA::dpsi_db_outgoing(): b out-of-range." << std::endl;
    	std::cerr << "b = " << b << ", bmax = " << bmax_outgoing(rspot) << ", |b - bmax| = " << fabs(b - bmax_outgoing(rspot)) << std::endl;
    	*prob = true;
    	dummy = -7888.0;
    	return dummy;
  	}
  
  	// set globals
  	OblDeflectionTOA_object = this;
  	OblDeflectionTOA_b_over_r = b/rspot;

  	//double rsurf = modptr->R_at_costheta(cos_theta);

  	//rsurf = rspot;

	double dpsidb(0.0);
	double split(0.9);

	if (b < 0.95*b_max)
	  dpsidb = Integration( 0.0, 1.0, OblDeflectionTOA_dpsi_db_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N );
	else
	  if ( b < 0.99995*b_max){
	      	
	    dpsidb = Integration( 0.0, split, OblDeflectionTOA_dpsi_db_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N );
	    dpsidb += Integration( split, 1.0, OblDeflectionTOA_dpsi_db_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N_1 );
	    //std::cout << "b/r = " << b/rspot << " b_max/r = " << b_max/rspot 
	    //	      << " b/b_max = " << b/b_max << " dpsidb = " << dpsidb << std::endl; 

	  }
	  else{
	    dpsidb = Integration( 0.0, split, OblDeflectionTOA_dpsi_db_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N );
	    dpsidb += Integration( split, 1.0, OblDeflectionTOA_dpsi_db_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N_1*100 );
	    //std::cout << "*******b/r = " << b/rspot << " b_max/r = " << b_max/rspot 
	    //	      << " b/b_max = " << b/b_max << " dpsidb = " << dpsidb << std::endl; 


	  }
	  
	    

	  //double dpsidb = Integration( 0.0, 1.0, OblDeflectionTOA_dpsi_db_integrand_wrapper_u );
  	
  	return dpsidb;
}

/*****************************************************/
// warning: i would have preferred to do this directly,
// but the integrand diverges at rcrit, and so it is difficult
// to see how to get the boundary term of this derivative to
// come out properly. however, we have a closed-form analytical
// approximation of the integral between rcrit and Rsurf, 
// so we use it instead for the derivative.
/*****************************************************/
double OblDeflectionTOA::dpsi_db_ingoing( const double& b, const double& rspot, const double& cos_theta, bool *prob ) {

  	double rcrit = this->rcrit( b, cos_theta, &curve.problem );
  	double m = this->get_mass();
  	double drcrit_db = sqrt( 1.0 - 2.0 * m / rcrit ) / (1.0 - ( (b / rcrit) * (m / rcrit)
	        		   / sqrt( 1.0 - 2.0 * m / rcrit )) );

	double dpsidb_in = -sqrt(2.0) * (drcrit_db) * (rspot - 3.0 * m)
					   / (sqrt( (rspot - rcrit) / (rcrit - 3.0 * m) ) 
					   * pow( rcrit - 3.0 * m , 2.0 ) );

  	return double ( dpsidb_in + this->dpsi_db_outgoing( b, rspot, &curve.problem ) );
}

/********************************************************************************/
/* OblDeflectionTOA::toa_outgoing                                               */
/*														                        */
/* Computes the time-of-arrival for initially outgoing photons                  */
/* Passed b (photon impact parameter), rspot (radius of NS at spot, unitless),  */
/*        prob (problem toggle)                                                 */
/* Returns the time of arrival of the photon                                    */
/********************************************************************************/
double OblDeflectionTOA::toa_outgoing ( const double& b, const double& rspot, bool *prob ) {
	// costheta_check(cos_theta);
  	//double dummy;

  	// set globals
  	OblDeflectionTOA_object = this;
  	OblDeflectionTOA_b_value = b;

  	// Note: Use an approximation to the integral near the surface
  	// to avoid divergence problems, and then use the real
  	// integral afterwards. See astro-ph/07030123 for the formula.

  	// The idea is to compute the time for a ray emitted from the
  	// surface, and subtract the time for a b=0 ray from r=rpole.
  	// To avoid large numbers, part of the subtraction is accomplished in the integrand.

  	// Note that for a b=0 ray, Delta T = int(1/(1-2*M/r), r),
  	// which is r + 2 M ln (r - 2M)

  	//double rsurf = modptr->R_at_costheta(cos_theta);
 
  	double rsurf = rspot;
  	double rpole = modptr->R_at_costheta(1.0);

  	//std::cout << "toa_outgoing: rpole = " << rpole << std::endl;

  	double toa = Integration( rsurf, get_rfinal(), OblDeflectionTOA_toa_integrand_minus_b0_wrapper );

  	double toa_b0_polesurf = (rsurf - rpole) + 2.0 * get_mass() * ( log( rsurf - 2.0 * get_mass() )
		       				 - log( rpole - 2.0 * get_mass() ) );

  	return double( toa - toa_b0_polesurf);
}

double OblDeflectionTOA::toa_outgoing_u ( const double& b, const double& rspot, bool *prob ) {
	// costheta_check(cos_theta);
  	//double dummy;

  	// set globals

  double b_max=bmax_outgoing(rspot);

  OblDeflectionTOA_object = this;
  OblDeflectionTOA_b_over_r = b/rspot;

  	// Note: Use an approximation to the integral near the surface
  	// to avoid divergence problems, and then use the real
  	// integral afterwards. See astro-ph/07030123 for the formula.

  	// The idea is to compute the time for a ray emitted from the
  	// surface, and subtract the time for a b=0 ray from r=rpole.
  	// To avoid large numbers, part of the subtraction is accomplished in the integrand.

  	// Note that for a b=0 ray, Delta T = int(1/(1-2*M/r), r),
  	// which is r + 2 M ln (r - 2M)

  	//double rsurf = modptr->R_at_costheta(cos_theta);
 
  	double rsurf = rspot;
  	double rpole = modptr->R_at_costheta(1.0);

  	//std::cout << "toa_outgoing: rpole = " << rpole << std::endl;
	double toa(0.0);
	double split(0.99);

	if (b < 0.95*b_max)
	  toa = Integration( 0.0, 1.0, OblDeflectionTOA_toa_integrand_minus_b0_wrapper_u,TRAPEZOIDAL_INTEGRAL_N/10 );
	else
	  if ( b < 0.9999999*b_max){
	      	
	    toa = Integration( 0.0, split, OblDeflectionTOA_toa_integrand_minus_b0_wrapper_u,TRAPEZOIDAL_INTEGRAL_N );
	    toa += Integration( split, 1.0, OblDeflectionTOA_toa_integrand_minus_b0_wrapper_u,TRAPEZOIDAL_INTEGRAL_N_1/10 );
	    //std::cout << "b/r = " << b/rspot << " b_max/r = " << b_max/rspot 
	    //	      << " b/b_max = " << b/b_max << " toa = " << toa << std::endl; 

	  }
	  else{
	    toa = Integration( 0.0, split, OblDeflectionTOA_toa_integrand_minus_b0_wrapper_u,TRAPEZOIDAL_INTEGRAL_N );
	    toa += Integration( split, 1.0, OblDeflectionTOA_toa_integrand_minus_b0_wrapper_u,TRAPEZOIDAL_INTEGRAL_N_1*10 );
	    //std::cout << "*******b/r = " << b/rspot << " b_max/r = " << b_max/rspot 
	    //	      << " b/b_max = " << b/b_max << " toa = " << toa << std::endl; 


	  }
	  
	    
  	//double toa = Integration( 0.0, 1.0, OblDeflectionTOA_toa_integrand_minus_b0_wrapper_u );

  	double toa_b0_polesurf = (rsurf - rpole) + 2.0 * get_mass() * ( log( rsurf - 2.0 * get_mass() )
		       				 - log( rpole - 2.0 * get_mass() ) );

  	return double( toa - toa_b0_polesurf);
}


/*****************************************************/
/*****************************************************/
double OblDeflectionTOA::toa_ingoing ( const double& b, const double& rspot, const double& cos_theta, bool *prob ) {
  //double dummy;
  	//costheta_check(cos_theta);
	/*
  	if ( b > bmax_outgoing(rspot) || b < bmin_ingoing(rspot, cos_theta) ) {
    	std::cerr << "ERROR in OblDeflectionTOA::toa_ingoing(): b out-of-range." << std::endl;
    	*prob = true;
    	dummy = -7888.0;
    	return dummy;
  	}
	*/
  	// std::cerr << "DEBUG: rcrit (km) = " << Units::nounits_to_cgs(this->rcrit(b,cos_theta), Units::LENGTH)/1.0e5 << std::endl;

  	OblDeflectionTOA_object = this;
  	OblDeflectionTOA_b_value = b;

  	// See psi_ingoing. Use an approximate formula for the integral near rcrit, and the real formula elsewhere
  
  	double rcrit = this->rcrit( b, cos_theta, &curve.problem );
 
  	double toa_in = rcrit / sqrt( 1 - 2.0 * get_mass() / rcrit ) * 2.0 * 
  					sqrt( 2.0 * (modptr->R_at_costheta(cos_theta) - rcrit) / 
  					(rcrit - 3.0 * get_mass()) );

	//  	double rspot( modptr->R_at_costheta(cos_theta) );

  	return double( toa_in + this->toa_outgoing( b, rspot, &curve.problem ) );

}

/*****************************************************/
// gives the integrand from equation 20, MLCB
/*****************************************************/
double OblDeflectionTOA::psi_integrand ( const double& b, const double& r) const { 
  double integrand( b / (r * r * sqrt( 1.0 - pow( b / r, 2.0 ) * (1.0 - 2.0 * get_mass_over_r() * get_rspot()/r ) )) );
  	return integrand;
}

/*****************************************************/
// Integrand with respect to u=R/r. b_R = b/R
/*****************************************************/
double OblDeflectionTOA::psi_integrand_u ( const double& b_R, const double& u) const { 
  // double integrand( b / (r * r * sqrt( 1.0 - pow( b / r, 2.0 ) * (1.0 - 2.0 * get_mass_over_r() * get_rspot()/r ) )) );
  double integrand ( b_R / sqrt( 1.0 - pow(u*b_R,2)*(1.0-2.0*get_mass_over_r()*u) ));
  	return integrand;
}

/*****************************************************/
// integrand for the derivative of equation 20, MLCB
/*****************************************************/
double OblDeflectionTOA::dpsi_db_integrand ( const double& b, const double& r ) const { 
  // double integrand( (1.0 / ( r * r * sqrt( 1.0 - pow( b / r , 2.0 ) * (1.0 - 2.0 * get_mass_over_r() * get_rspot()/r) )) ) 
  //		    + (b * b * (1.0 - 2.0 * get_mass_over_r() * get_rspot()/r) 
  //		       / (pow( r , 4.0 ) * pow( 1.0 - pow( b / r , 2.0 ) * (1.0 - 2.0 * get_mass_over_r() * get_rspot()/r) , 1.5) ) ) );

  double integrand ( 1.0/(r*r) *  
		     pow( 1.0 - pow( b / r , 2.0 ) * (1.0 - 2.0 * get_mass_over_r() * get_rspot()/r) , -1.5));

  	return integrand;
}

/*****************************************************/
// integrand for the derivative of equation 20, MLCB
/*****************************************************/
double OblDeflectionTOA::dpsi_db_integrand_u ( const double& b_R, const double& u ) const { 
  // double integrand( (1.0 / ( r * r * sqrt( 1.0 - pow( b / r , 2.0 ) * (1.0 - 2.0 * get_mass_over_r() * get_rspot()/r) )) ) 
  //		    + (b * b * (1.0 - 2.0 * get_mass_over_r() * get_rspot()/r) 
  //		       / (pow( r , 4.0 ) * pow( 1.0 - pow( b / r , 2.0 ) * (1.0 - 2.0 * get_mass_over_r() * get_rspot()/r) , 1.5) ) ) );

  double integrand ( 1.0/get_rspot() *  pow( 1.0 - pow( b_R*u , 2.0 ) * (1.0 - 2.0 * get_mass_over_r() * u) , -1.5));

  	return integrand;
}


/*****************************************************/
// integrand from equation 38, MLCB, without the -1
/*****************************************************/
double OblDeflectionTOA::toa_integrand ( const double& b, const double& r ) const { 
  double integrand( (1.0 / (1.0 - 2.0 * get_mass_over_r() * get_rspot()/r)) * ((1.0 / sqrt( 1.0 - pow( b / r , 2.0 )
											    * (1.0 - 2.0 * get_mass_over_r()*get_rspot()/r) )) ) );
  	return integrand;
}

/*****************************************************/
// gives the integrand from equation 38, MLCB
/*****************************************************/
double OblDeflectionTOA::toa_integrand_minus_b0 ( const double& b, const double& r ) const {   
  double integrand( (1.0 / (1.0 - 2.0 * get_mass_over_r()*get_rspot()/r)) 
		    * ( (1.0 / sqrt( 1.0 - pow( b / r , 2.0 ) * (1.0 - 2.0 * get_mass_over_r()*get_rspot()/r) )) - 1.0) );
  	return integrand;
}

/*****************************************************/
// gives the integrand from equation 38, MLCB
/*****************************************************/
double OblDeflectionTOA::toa_integrand_minus_b0_u ( const double& b_R, const double& u ) const {   
  double integrand(  get_rspot()/(u*u) * (1.0 / (1.0 - 2.0 * get_mass_over_r()*u)) 
		    * ( (1.0 / sqrt( 1.0 - pow( b_R * u , 2.0 ) * (1.0 - 2.0 * get_mass_over_r()*u) )) - 1.0) );
  	return integrand;
}

/*****************************************************/
// used in oblate shape of star
/*****************************************************/
double OblDeflectionTOA::rcrit_zero_func ( const double& rc, const double& b ) const { 
  	return double( rc - b * sqrt( 1.0 - 2.0 * get_mass() / rc ) );
}

/*****************************************************/
// used in oblate shape of star
/*****************************************************/
double OblDeflectionTOA::b_from_psi_ingoing_zero_func ( const double& b, 
														const double& cos_theta, 
														const double& psi ) const { 
  	return double( psi - this->psi_ingoing(b, cos_theta, &curve.problem) );
}

/*****************************************************/
/*****************************************************/
double OblDeflectionTOA::b_from_psi_outgoing_zero_func ( const double& b, 
														 const double& cos_theta, 
														 const double& psi,
														 const double& b_max, 
														 const double& psi_max, 
														 const double& b_guess, 
														 const double& psi_guess ) const {
  	return double( psi - this->psi_outgoing( b, cos_theta, b_max, psi_max, &curve.problem ) );
}

/*****************************************************/
/* OblDeflectionTOA::TrapezoidalInteg                */
/*                                                   */
/* Numerically integrates using the trapezoid rule.  */
/*****************************************************/
double OblDeflectionTOA::TrapezoidalInteg ( const double& a, const double& b, 
                                            double (*func)(double x, bool *prob),
					   						const long int& N ) {
  	if ( a == b ) return double(0.0);

  	double integral(0.0);
  	//const long int N( OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_N );
  	const double power( OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_POWER );

  	for ( int i(1); i <= N; i++ ) {
    	double left, right, mid, width, fmid;
    	left = TrapezoidalInteg_pt( a, b, power, N, i-1, &curve.problem );
    	right = TrapezoidalInteg_pt( a, b, power, N, i, &curve.problem );
    	mid = (left + right) / 2.0;
    	width = right - left;
    	fmid = (*func)(mid, &curve.problem); //Double Check
    
    	integral += (width * fmid);
  	}

  	return integral;
}

/*****************************************************/
/* OblDeflectionTOA::TrapezoidalInteg                */
/*                                                   */
/* Numerically integrates using the trapezoid rule.  */


/*******************************************************/
/* OblDeflectionTOA::TrapezoidalInteg_pt               */
/*                                                     */
/* Called in TrapezoidalInteg                          */
/* Power spacing of subdivisions for TrapezoidalInteg  */
/* i ranges from 0 to N                                */
/*******************************************************/
double OblDeflectionTOA::TrapezoidalInteg_pt ( const double& a, const double& b, 
					      					   const double& power, const long int& N,
					      					   const long int& i, bool *prob ) {
  	double dummy;
  	if ( N <= 0 ) {
    	std::cerr << "ERROR in OblDeflectionTOA::TrapezoidalInteg_pt: N <= 0." << std::endl;
    	*prob = true;
    	dummy = -7888.0;
    	return dummy;
  	}
  	if ( i < 0 || i > N ) {
    	std::cerr << "ERROR in OblDeflectionTOA::TrapezoidalInteg_pt: i out of range." << std::endl;
    	*prob = true;
    	dummy = -7888.0;
    	return dummy;
  	}
  	return double( a + (b - a) * pow( 1.0 * i / N , power ) );
}
/*******************************************************/
/* OblDeflectionTOA::Integration                       */
/*                                                     */
/* Integrates! By calling TrapezoidalInteg             */
/* Here, N = TRAPEZOIDAL_INTEGRAL_N                    */
/*******************************************************/
inline double OblDeflectionTOA::Integration ( const double& a, const double& b, 
                                              double (*func)(double x, bool *prob),
					     					  const long int& N ) {
  	return TrapezoidalInteg( a, b, func, N );
  	//return MATPACK::AdaptiveSimpson( a, b, func, OblDeflectionTOA::INTEGRAL_EPS );
}

