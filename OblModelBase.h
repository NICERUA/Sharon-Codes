// OblModelBase.h
//
// Interface describing an Oblateness model
// (C) Coire Cadeau, 2007

// Source (C) Coire Cadeau 2007, all rights reserved.
//
// Permission is granted for private use only, and not
// distribution, either verbatim or of derivative works,
// in whole or in part.
//
// The code is not thoroughly tested or guaranteed for
// any particular use.

#ifndef OBLMODELBASE_H
#define OBLMODELBASE_H

#include <exception>

class OblModelBase {
 public:
  virtual double R_at_costheta( const double& costheta ) const throw(std::exception) = 0;
  virtual double Dtheta_R( const double& costheta ) const throw(std::exception) = 0;
  virtual double f(const double& costheta) const throw(std::exception) =0; // return (1+z) (dR/dtheta) / R
  virtual double cos_gamma(const double& costheta) const throw(std::exception) = 0; // return cosine of angle between radial and normal
  virtual ~OblModelBase() { }
};

#endif // OBLMODELBASE_H
