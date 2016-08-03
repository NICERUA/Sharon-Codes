/***************************************************************************************/
/*                              PolyOblModelBase.cpp

    This code is the interface for a polynomial oblate NS model.
    
    Based on code written by Coire Cadeau and modified by Sharon Morsink and 
    Abigail Stevens.
    
    (C) Coire Cadeau, 2007; Source (C) Coire Cadeau 2007, all rights reserved.
    Permission is granted for private use only, and not distribution, either verbatim or
    of derivative works, in whole or in part.
    This code is not thoroughly tested or guaranteed for any particular use.
*/
/***************************************************************************************/

#include "PolyOblModelBase.h"
#include <cmath>
#include "Exception.h"

PolyOblModelBase::PolyOblModelBase( const double& Rspot_nounits_value, 
				    const double& Req_nounits_value, const double& zetaval, const double& epsval )
  : Rspot_nounits(Rspot_nounits_value), Req_nounits(Req_nounits_value), zeta(zetaval), eps(epsval) { }

double PolyOblModelBase::R_at_costheta( const double& costheta ) const throw(std::exception) {
  	// Return R(theta) in "nounits".
  	// note that the user supplies cos(theta) and not theta.

  	if (costheta == 0.0)
    	return double( get_Req_nounits()*( 1.0 + a0() + a2()*(-0.5) + a4()*(3.0/8.0) ) ); 
  	else{
    	if (costheta == 1.0)
      		return double( get_Req_nounits()*( 1.0 + a0() + a2() + a4() ) );  
    	else
      		return double( get_Rspot_nounits() );
  	}

    //return double( get_Req_nounits()*( 1.0 + a0()*P0(costheta) + a2()*P2(costheta) + a4()*P4(costheta) ) ); 
}

double PolyOblModelBase::Dtheta_R( const double& costheta ) const throw(std::exception) {
  	// Return dR(theta) / dtheta in "nounits".
  	// note that the user supplies cos(theta) and not theta.
  	if(costheta == 0.0 || fabs(costheta) == 1.0) return double(0.0);
  	else {
   		// -sin(theta) * Dcostheta(R)
    	return double( -sqrt(1.0-costheta*costheta) * get_Req_nounits()
		   *( a0()*Dmu_P0(costheta) + a2()*Dmu_P2(costheta) + a4()*Dmu_P4(costheta) ) );
	}
}

double PolyOblModelBase::z( const double& costheta ) const {
  	return double(1.0/sqrt(1.0 - 2.0*get_zeta()*get_Req_nounits()/R_at_costheta(costheta) ) - 1.0);
}

double PolyOblModelBase::f( const double& costheta ) const throw(std::exception) {
  	return double( (1.0 + z(costheta))*Dtheta_R(costheta)/R_at_costheta(costheta) );
}

double PolyOblModelBase::cos_gamma(const double& costheta) const throw(std::exception) {
  	return double( 1.0/sqrt(1.0 + pow(f(costheta),2.0)) );
}

double PolyOblModelBase::zetaparam( const double& Mass_nounits, const double& Req_nounits ) {
  	return double( Mass_nounits / Req_nounits );
}

double PolyOblModelBase::epsparam ( const double& Omega_nounits, const double& Mass_nounits, const double& Req_nounits ) {
  	return double( pow(Omega_nounits,2.0) * pow(Req_nounits,3.0) / Mass_nounits );
}

double PolyOblModelBase::P0( const double& mu ) {
  	return double( 1.0 );
}

double PolyOblModelBase::P2( const double& mu ) {
  	return double( (3.0 * mu * mu - 1.0) / 2.0 );
}

double PolyOblModelBase::P4( const double& mu ) {
  	return double( (35.0 * pow( mu , 4.0 ) - 30.0 * mu * mu + 3.0) / 8.0 );
}

double PolyOblModelBase::Dmu_P0( const double& mu ) {
  	return double( 0.0 );
}

double PolyOblModelBase::Dmu_P2( const double& mu ) {
  	return double( 3.0*mu ); 
}

double PolyOblModelBase::Dmu_P4( const double& mu ) {
  	return double( mu*(35.0*mu*mu - 15.0)/2.0 );
}

double PolyOblModelBase::get_Req_nounits() const { 
  	return double( Req_nounits );
}

double PolyOblModelBase::get_Rspot_nounits() const { 
  	return double( Rspot_nounits );
}

double PolyOblModelBase::get_zeta() const {
  	return double( zeta );
}

double PolyOblModelBase::get_eps() const {
  	return double( eps );
}
