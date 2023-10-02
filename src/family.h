// family.h
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 01/10/2023

#ifndef FAMILY_H
#define FAMILY_H

#include <RcppArmadillo.h>

namespace Family {

class Family {
    public:
        std::string family;
        double dispersion;
        virtual arma::mat variance (const arma::mat & mu) = 0;
        virtual arma::mat initialize (const arma::mat & y) = 0;
        virtual arma::mat devresid (const arma::mat & y, const arma::mat & mu) = 0;
        virtual bool validmu (const arma::mat & mu) = 0;
        virtual bool valideta (const arma::mat & eta) = 0;
        Family () {this->family = "Family"; this->dispersion = 1;}
        virtual ~Family () {}
};

class Gaussian : public Family {
    public:
        arma::mat variance (const arma::mat & mu);
        arma::mat initialize (const arma::mat & y);
        arma::mat devresid (const arma::mat & y, const arma::mat & mu);
        bool validmu (const arma::mat & mu);
        bool valideta (const arma::mat & eta);
        Gaussian () {this->family = "Gaussian";}
};

class Binomial : public Family {
    public:
        arma::mat variance (const arma::mat & mu);
        arma::mat initialize (const arma::mat & y);
        arma::mat devresid (const arma::mat & y, const arma::mat & mu);
        bool validmu (const arma::mat & mu);
        bool valideta (const arma::mat & eta);
        Binomial () {this->family = "Binomial";}
};

class Poisson : public Family {
    public:
        arma::mat variance (const arma::mat & mu);
        arma::mat initialize (const arma::mat & y);
        arma::mat devresid (const arma::mat & y, const arma::mat & mu);
        bool validmu (const arma::mat & mu);
        bool valideta (const arma::mat & eta);
        Poisson () {this->family = "Poisson";}
};

class Gamma : public Family {
    public:
        arma::mat variance (const arma::mat & mu);
        arma::mat initialize (const arma::mat & y);
        arma::mat devresid (const arma::mat & y, const arma::mat & mu);
        bool validmu (const arma::mat & mu);
        bool valideta (const arma::mat & eta);
        Gamma () {this->family = "Gamma";}
};

}

#endif