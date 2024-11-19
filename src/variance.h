// variance.h
// author: Cristian Castiglione
// creation: 08/11/2023
// last change: 19/11/2024

#ifndef VARIANCE_H
#define VARIANCE_H

#include <RcppArmadillo.h>
#include "utils.h"

namespace glm {

class Variance {
    public:
        std::string varf = "Variance";
        virtual bool validmu (const arma::mat & mu) = 0;
        virtual arma::mat varfun (const arma::mat & mu, const double & phi) = 0;
        virtual arma::mat initfun (const arma::mat & y) = 0;
        virtual arma::mat devfun (const arma::mat & y, const arma::mat & mu, const double & phi) = 0;
        virtual ~Variance () {}
};

class Constant : public Variance {
    public:
        bool validmu (const arma::mat & mu);
        arma::mat initfun (const arma::mat & y);
        arma::mat varfun (const arma::mat & mu, const double & phi);
        arma::mat devfun (const arma::mat & y, const arma::mat & mu, const double & phi);
        Constant () {this->varf = "const";}
};

class Linear : public Variance {
    public:
        bool validmu (const arma::mat & mu);
        arma::mat initfun (const arma::mat & y);
        arma::mat varfun (const arma::mat & mu, const double & phi);
        arma::mat devfun (const arma::mat & y, const arma::mat & mu, const double & phi);
        Linear () {this->varf = "mu";}
};

class Squared : public Variance {
    public:
        bool validmu (const arma::mat & mu);
        arma::mat initfun (const arma::mat & y);
        arma::mat varfun (const arma::mat & mu, const double & phi);
        arma::mat devfun (const arma::mat & y, const arma::mat & mu, const double & phi);
        Squared () {this->varf = "mu^2";}
};

class Cubic : public Variance {
    public:
        bool validmu (const arma::mat & mu);
        arma::mat initfun (const arma::mat & y);
        arma::mat varfun (const arma::mat & mu, const double & phi);
        arma::mat devfun (const arma::mat & y, const arma::mat & mu, const double & phi);
        Cubic () {this->varf = "mu^3";}
};

class cSquared : public Variance {
    public:
        bool validmu (const arma::mat & mu);
        arma::mat initfun (const arma::mat & y);
        arma::mat varfun (const arma::mat & mu, const double & phi);
        arma::mat devfun (const arma::mat & y, const arma::mat & mu, const double & phi);
        cSquared () {this->varf = "mu(1-mu)";}
};

class NBVariance : public Variance {
    public:
        bool validmu (const arma::mat & mu);
        arma::mat initfun (const arma::mat & y);
        arma::mat varfun (const arma::mat & mu, const double & phi);
        arma::mat devfun (const arma::mat & y, const arma::mat & mu, const double & phi);
        NBVariance () {this->varf = "mu(1+t*mu)";}
};

}

#endif