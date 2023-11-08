// variance.cpp
// author: Cristian Castiglione
// creation: 08/11/2023
// last change: 08/11/2023

#include "variance.h"

using namespace glm;

// Variance functions
arma::mat Constant::varfun (const arma::mat & mu) {return arma::ones(arma::size(mu));}
arma::mat Linear::varfun (const arma::mat & mu) {return mu;}
arma::mat Squared::varfun (const arma::mat & mu) {return mu % mu;}
arma::mat Cubic::varfun (const arma::mat & mu) {return mu % mu % mu;}
arma::mat cSquared::varfun (const arma::mat & mu) {return mu % (1 - mu);}
