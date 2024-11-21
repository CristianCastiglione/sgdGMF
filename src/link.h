// link.h
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 21/10/2024

#ifndef LINK_H
#define LINK_H

#include <RcppArmadillo.h>
#include "utils.h"

namespace glm {

class Link {
    public:
        std::string link = "Link";
        virtual bool valideta (const arma::mat & eta) = 0;
        virtual arma::mat linkfun (const arma::mat & mu) = 0;
        virtual arma::mat linkinv (const arma::mat & eta) = 0;
        virtual arma::mat mueta (const arma::mat & eta) = 0;
        virtual ~Link () {}
};

class Identity : public Link {
    public:
        bool valideta (const arma::mat & eta);
        arma::mat linkfun (const arma::mat & mu);
        arma::mat linkinv (const arma::mat & eta);
        arma::mat mueta (const arma::mat & eta);
        Identity () {this->link = "Identity";}
};

class Logit : public Link {
    public:
        bool valideta (const arma::mat & eta);
        arma::mat linkfun (const arma::mat & mu);
        arma::mat linkinv (const arma::mat & eta);
        arma::mat mueta (const arma::mat & eta);
        Logit () {this->link = "Logit";}
};

class Probit : public Link {
    public:
        bool valideta (const arma::mat & eta);
        arma::mat linkfun (const arma::mat & mu);
        arma::mat linkinv (const arma::mat & eta);
        arma::mat mueta (const arma::mat & eta);
        Probit () {this->link = "Probit";}
};

class Cauchit : public Link {
    public:
        bool valideta (const arma::mat & eta);
        arma::mat linkfun (const arma::mat & mu);
        arma::mat linkinv (const arma::mat & eta);
        arma::mat mueta (const arma::mat & eta);
        Cauchit () {this->link = "Cauchit";}
};

class cLogLog : public Link {
    public:
        bool valideta (const arma::mat & eta);
        arma::mat linkfun (const arma::mat & mu);
        arma::mat linkinv (const arma::mat & eta);
        arma::mat mueta (const arma::mat & eta);
        cLogLog () {this->link = "cLogLog";}
};

class Log : public Link {
    public:
        bool valideta (const arma::mat & eta);
        arma::mat linkfun (const arma::mat & mu);
        arma::mat linkinv (const arma::mat & eta);
        arma::mat mueta (const arma::mat & eta);
        Log () {this->link = "Log";}
};

class Inverse : public Link {
    public:
        bool valideta (const arma::mat & eta);
        arma::mat linkfun (const arma::mat & mu);
        arma::mat linkinv (const arma::mat & eta);
        arma::mat mueta (const arma::mat & eta);
        Inverse () {this->link = "Inverse";}
};

class SquaredInverse : public Link {
    public:
        bool valideta (const arma::mat & eta);
        arma::mat linkfun (const arma::mat & mu);
        arma::mat linkinv (const arma::mat & eta);
        arma::mat mueta (const arma::mat & eta);
        SquaredInverse () {this->link = "1/mu^2";}
};

class Sqrt : public Link {
    public:
        bool valideta (const arma::mat & eta);
        arma::mat linkfun (const arma::mat & mu);
        arma::mat linkinv (const arma::mat & eta);
        arma::mat mueta (const arma::mat & eta);
        Sqrt () {this->link = "Sqrt";}
};

}

#endif