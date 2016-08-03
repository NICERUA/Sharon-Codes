/***************************************************************************************/
/*                               OblDeflectionTOA.h

    This is the header file for OblDeflectionTOA.cpp.
    
    Based on code written by Coire Cadeau and modified by Sharon Morsink and 
    Abigail Stevens.
    
    (C) Coire Cadeau, 2007; Source (C) Coire Cadeau 2007, all rights reserved.
    Permission is granted for private use only, and not distribution, either verbatim or
    of derivative works, in whole or in part.
    This code is not thoroughly tested or guaranteed for any particular use.
*/
/***************************************************************************************/

#ifndef OBLDEFLECTIONTOA_H
#define OBLDEFLECTIONTOA_H

#include "OblModelBase.h"
#include "Exception.h"

// Globals are bad, but there's no easy way to get around it
// for this particular case (need to pass member functions
// to Matpack routines which do not have a signature to accomodate
// passing in the object pointer).

double OblDeflectionTOA_psi_integrand_wrapper( double r, bool *prob );
double OblDeflectionTOA_dpsi_db_integrand_wrapper( double r, bool *prob );
double OblDeflectionTOA_toa_integrand_wrapper( double r, bool *prob );
double OblDeflectionTOA_toa_integrand_minus_b0_wrapper( double r, bool *prob );
double OblDeflectionTOA_rcrit_zero_func_wrapper( double rc );
double OblDeflectionTOA_b_from_psi_ingoing_zero_func_wrapper( double b );
double OblDeflectionTOA_b_from_psi_outgoing_zero_func_wrapper( double b );

double OblDeflectionTOA_psi_integrand_wrapper_u ( double u,  bool *prob );
double OblDeflectionTOA_dpsi_db_integrand_wrapper_u ( double u, bool *prob );
double OblDeflectionTOA_toa_integrand_minus_b0_wrapper_u ( double u, bool *prob );


// End global pollution.

class OblDeflectionTOA {
	static const double INTEGRAL_EPS;                     //
	static const double FINDZERO_EPS;                     //
	static const double RFINAL_MASS_MULTIPLE;             //
	static const double DIVERGENCE_GUARD;                 //
  	static const long int TRAPEZOIDAL_INTEGRAL_N;         //
  	static const long int TRAPEZOIDAL_INTEGRAL_N_MAX;         //
  	static const long int TRAPEZOIDAL_INTEGRAL_N_1;         //
  	static const long int TRAPEZOIDAL_INTEGRAL_N_MAX_1;         //
	static const long int TRAPEZOIDAL_INTEGRAL_INGOING_N; //
  	static const double TRAPEZOIDAL_INTEGRAL_POWER;       //

 	private:
  		OblModelBase* modptr; // pointer to the model?
  		double mass;          // NS mass should be stored in internal (unitless?) units.
		double mass_over_r;
		double rspot;
  		const double r_final; //

  		void costheta_check( const double& cos_theta ) const throw(std::exception){
    		if( fabs(cos_theta) > 1 )
      		throw(Exception("cos_theta out of range."));
  		}

 	public: // all of the member functions for which we provide hooks for MATPACK
 		double psi_integrand ( const double& b, const double& r) const;
  		double dpsi_db_integrand ( const double& b, const double& r ) const;
  		double toa_integrand ( const double& b, const double& r ) const;
  		double toa_integrand_minus_b0 ( const double& b, const double& r ) const;
  		double rcrit_zero_func( const double& rc, const double& b ) const;
  		double b_from_psi_ingoing_zero_func ( const double& b, const double& cos_theta, 
  								              const double& psi) const;
  		double b_from_psi_outgoing_zero_func ( const double& b, const double& cos_theta, 
  		                                       const double& psi, const double& b_max, 
  		                                       const double& psi_max, 
  		                                       const double &b_guess, 
  		                                       const double& psi_guess ) const;
		double psi_integrand_u ( const double& b_R, const double& u) const;
		double dpsi_db_integrand_u ( const double& b_R, const double& u ) const ;
		double toa_integrand_minus_b0_u ( const double& b_R, const double& u ) const;

	public:
 		OblDeflectionTOA ( OblModelBase* modptr, const double& mass_nounits ,const double& mass_over_r_nounits, const double& radius_nounits);
  		double bmax_outgoing ( const double& rspot ) const;
  		double bmin_ingoing ( const double& rspot, const double& cos_theta ) const;
  		bool ingoing_allowed ( const double& cos_theta );
  		double get_mass() const { return mass; }
		double get_mass_over_r() const { return mass_over_r; }
		double get_rspot() const { return rspot;}
  		double get_rfinal() const { return r_final; }

  		double rcrit ( const double& b, const double& cos_theta, bool *prob ) const;
  
  		double psi_outgoing ( const double& b, const double& rspot, const double& b_max, 
  		                      const double& psi_max, bool *prob ) const;
  		double psi_max_outgoing ( const double& b, const double& rspot, bool *prob );
  		double psi_ingoing ( const double& b, const double& cos_theta, bool *prob ) const;
  		bool b_from_psi ( const double& psi, const double& rspot, const double& cos_theta, double& b, 
  		                  int& rdot, const double& bmax_out, const double& psi_out_max, 
		   				  const double& b_guess, const double& psi_guess,
		   				  const double& b2, const double&psi2, bool *prob );
  		double dpsi_db_outgoing ( const double& b, const double& rspot, bool *prob );
  		double dpsi_db_ingoing ( const double& b, const double& rspot, const double& cos_theta, bool *prob ); //changed GC

		double psi_outgoing_u ( const double& b, const double& rspot,
					const double& b_max, const double& psi_max, bool *prob ) const;
		double psi_max_outgoing_u ( const double& b, const double& rspot, bool *prob ) const;
		double dpsi_db_outgoing_u( const double& b, const double& rspot, bool *prob ) const;
		double toa_outgoing_u ( const double& b, const double& rspot, bool *prob );

  		double toa_outgoing ( const double& b, const double& cos_theta, bool *prob );
  		double toa_ingoing ( const double& b, const double& rspot, const double& cos_theta, bool *prob );

 		// a trapezoidal integrator which will not evaluate func at the first endpoint, a
  		static double TrapezoidalInteg ( const double& a, const double& b, 
  										 double (*func)(double x, bool *prob),
				  						 const long int& N );
  	
  		// power spacing of subdivisions for above integrator, i ranges from 0 to N.
  		static double TrapezoidalInteg_pt ( const double& a, const double& b, 
				     					    const double& power, const long int& N,
				     					    const long int& i, bool *prob);
  		inline static double Integration ( const double& a, const double& b, 
  										   double (*func)(double x, bool *prob),
				    					   const long int& N = TRAPEZOIDAL_INTEGRAL_N );

  	
};

#endif // OBLDEFLECTIONTOA_H
