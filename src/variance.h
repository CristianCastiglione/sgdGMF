// variance.h
// author: Cristian Castiglione
// creation: 08/11/2023
// last change: 08/11/2023

#ifndef VARIANCE_H
#define VARIANCE_H

#include <RcppArmadillo.h>

namespace glm {

class Variance {
    public:
        std::string var = "Variance";
        virtual arma::mat varfun (const arma::mat & mu) = 0;
        virtual ~Variance () {}
};

class Constant : public Variance {
    public:
        arma::mat varfun (const arma::mat & mu);
        Constant () {this->var = "const";}
};

class Linear : public Variance {
    public:
        arma::mat varfun (const arma::mat & mu);
        Linear () {this->var = "mu";}
};

class Squared : public Variance {
    public:
        arma::mat varfun (const arma::mat & mu);
        Squared () {this->var = "mu^2";}
};

class Cubic : public Variance {
    public:
        arma::mat varfun (const arma::mat & mu);
        Cubic () {this->var = "mu^3";}
};

class cSquared : public Variance {
    public:
        arma::mat varfun (const arma::mat & mu);
        cSquared () {this->var = "mu(1-mu)";}
};

}

#endif