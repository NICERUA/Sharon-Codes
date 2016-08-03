/***************************************************************************************/
/*                               PolyOblModelBase.h

    This is the header file for PolyOblModelBase.cpp.
    
    Based on code written by Coire Cadeau and modified by Sharon Morsink and 
    Abigail Stevens.
    
    (C) Coire Cadeau, 2007; Source (C) Coire Cadeau 2007, all rights reserved.
    Permission is granted for private use only, and not distribution, either verbatim or
    of derivative works, in whole or in part.
    This code is not thoroughly tested or guaranteed for any particular use.
*/
/***************************************************************************************/

// SMM: Nov 30, 2009
// Added Rspot_nounits as a parameter and function.

#ifndef POLYOBLMODELBASE_H
#define POLYOBLMODELBASE_H

#include "OblModelBase.h"
#include <exception>

class PolyOblModelBase : public OblModelBase {
 	private:
  		double Rspot_nounits;
  		double Req_nounits, zeta, eps;
  		double z(const double& costheta) const;

 	public:
  		PolyOblModelBase( const double& Rspot_nounits, const double& Req_nounits, const double& zeta, const double& eps );
  		double R_at_costheta( const double& costheta ) const throw(std::exception);
  		double Dtheta_R( const double& costheta ) const throw(std::exception);
  		double f(const double& costheta)  const throw(std::exception);
  		double cos_gamma(const double& costheta) const throw(std::exception);
  		static double zetaparam( const double& Mass_nounits, const double& Req_nounits );
  		static double epsparam( const double& Omega_nounits, const double& Mass_nounits, const double& Req_nounits );
  		virtual ~PolyOblModelBase() { }

 	protected:
  		static double P0(const double& mu);
  		static double P2(const double& mu);
  		static double P4(const double& mu);
  		static double Dmu_P0(const double& mu);
  		static double Dmu_P2(const double& mu);
  		static double Dmu_P4(const double& mu);

  	virtual double a0() const = 0;
  	virtual double a2() const = 0;
  	virtual double a4() const = 0;
  
  	double get_Req_nounits() const;
  	double get_zeta() const;
  	double get_eps() const;
  	double get_Rspot_nounits() const;
};

#endif // POLYOBLMODELBASE_H
