/***************************************************************************************/
/*                               PolyOblModelNHQS.h

    This is the header file for PolyOblModelNHQS.cpp, which sets up the model for an 
    oblate NS (neutron and hybrid quark star model).
    
    Based on code written by Coire Cadeau and modified by Sharon Morsink and 
    Abigail Stevens.
    
    (C) Coire Cadeau, 2007; Source (C) Coire Cadeau 2007, all rights reserved.
    Permission is granted for private use only, and not distribution, either verbatim or
    of derivative works, in whole or in part.
    This code is not thoroughly tested or guaranteed for any particular use.
*/
/***************************************************************************************/

#ifndef POLYOBLMODELNHQS_H
#define POLYOBLMODELNHQS_H

#include "PolyOblModelBase.h"

class PolyOblModelNHQS : public PolyOblModelBase {
 	public:
  		PolyOblModelNHQS( const double& Rspot_nounits, const double& Req_nounits, 
  						  const double& zeta, const double& eps );
 	protected:
  	double a0() const;
  	double a2() const;
  	double a4() const;
};

#endif // POLYOBLMODELNHQS_H
