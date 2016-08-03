/***************************************************************************************/
/*                               PolyOblModelNHQS.cpp

    This code sets up the model for an oblate NS (neutron and hybrid quark star model).
    
    Based on code written by Coire Cadeau and modified by Sharon Morsink and 
    Abigail Stevens.
    
    (C) Coire Cadeau, 2007; Source (C) Coire Cadeau 2007, all rights reserved.
    Permission is granted for private use only, and not distribution, either verbatim or
    of derivative works, in whole or in part.
    This code is not thoroughly tested or guaranteed for any particular use.
*/
/***************************************************************************************/

#include "PolyOblModelNHQS.h"

PolyOblModelNHQS::PolyOblModelNHQS( const double& Rspot_nounits, const double& Req_nounits, const double& zeta, const double& eps )
  : PolyOblModelBase(Rspot_nounits, Req_nounits, zeta, eps) { }

double PolyOblModelNHQS::a0() const {
  double eps(this->get_eps());
  double zeta(this->get_zeta());
  return double(-0.18*eps + 0.23*zeta*eps - 0.05*eps*eps);
}

double PolyOblModelNHQS::a2() const {
  double eps(this->get_eps());
  double zeta(this->get_zeta());
  return double(-0.39*eps + 0.29*zeta*eps + 0.13*eps*eps);
}

double PolyOblModelNHQS::a4() const {
  double eps(this->get_eps());
  double zeta(this->get_zeta());
    return double(0.04*eps - 0.15*zeta*eps + 0.07*eps*eps );

  //return double( -8.0/3.0*(-0.18*eps + 0.23*zeta*eps - 0.05*eps*eps - 0.5*(-0.39*eps + 0.29*zeta*eps + 0.13*eps*eps)));

}
