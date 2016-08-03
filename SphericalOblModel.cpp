/***************************************************************************************/
/*                               SphericalOblModel.cpp

    This code sets up the model for a spherical NS.
    
    Based on code written by Coire Cadeau and modified by Sharon Morsink and 
    Abigail Stevens.
    
    (C) Coire Cadeau, 2007; Source (C) Coire Cadeau 2007, all rights reserved.
    Permission is granted for private use only, and not distribution, either verbatim or
    of derivative works, in whole or in part.
    This code is not thoroughly tested or guaranteed for any particular use.
*/
/***************************************************************************************/

#include "SphericalOblModel.h"

SphericalOblModel::SphericalOblModel(const double& Req_val ) : Req(Req_val) { }

// For spherical model, radius is same at all parts on NS surface
double SphericalOblModel::R_at_costheta(const double& costheta ) const throw(std::exception) { 
	return Req; 
}

// For spherical model, 
double SphericalOblModel::Dtheta_R(const double& costheta ) const throw(std::exception) { 
	return double(0.0); 
}

// For spherical model,
double SphericalOblModel::f(const double& costheta) const throw(std::exception) {  // this f is the one from MLCB3
	return double(0.0); 
}

// For spherical model, cos_gamma = 0
double SphericalOblModel::cos_gamma(const double& costheta) const throw(std::exception) { 
	return double(1.0); 
}
