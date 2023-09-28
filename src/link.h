// link.h
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 28/09/2023

#ifndef LINK_H
#define LINK_H

#include <RcppArmadillo.h>
#include "utils.h"

namespace Link {

class Link {
    public:
        std::string link;
        virtual arma::mat linkfun (const arma::mat & mu) const = 0;
        virtual arma::mat linkinv (const arma::mat & eta) const = 0;
        virtual arma::mat mueta (const arma::mat & eta) const = 0;
        virtual ~Link () {}
};

class Identity : public Link {
    public:
        std::string link = "Identity";
        arma::mat linkfun (const arma::mat & mu) const;
        arma::mat linkinv (const arma::mat & eta) const;
        arma::mat mueta (const arma::mat & eta) const;
};

class Logit : public Link {
    public:
        std::string link = "Logit";
        arma::mat linkfun (const arma::mat & mu) const;
        arma::mat linkinv (const arma::mat & eta) const;
        arma::mat mueta (const arma::mat & eta) const;
};

class Probit : public Link {
    public:
        std::string link = "Probit";
        arma::mat linkfun (const arma::mat & mu) const;
        arma::mat linkinv (const arma::mat & eta) const;
        arma::mat mueta (const arma::mat & eta) const;
};

class Cauchy : public Link {
    public:
        std::string link = "Cauchy";
        arma::mat linkfun (const arma::mat & mu) const;
        arma::mat linkinv (const arma::mat & eta) const;
        arma::mat mueta (const arma::mat & eta) const;
};

class cLogLog : public Link {
    public:
        std::string link = "cLogLog";
        arma::mat linkfun (const arma::mat & mu) const;
        arma::mat linkinv (const arma::mat & eta) const;
        arma::mat mueta (const arma::mat & eta) const;
};

class Log : public Link {
    public:
        std::string link = "Log";
        arma::mat linkfun (const arma::mat & mu) const;
        arma::mat linkinv (const arma::mat & eta) const;
        arma::mat mueta (const arma::mat & eta) const;
};

class Inverse : public Link {
    public:
        std::string link = "Inverse";
        arma::mat linkfun (const arma::mat & mu) const;
        arma::mat linkinv (const arma::mat & eta) const;
        arma::mat mueta (const arma::mat & eta) const;
};

class Sqrt : public Link {
    public:
        std::string link = "Sqrt";
        arma::mat linkfun (const arma::mat & mu) const;
        arma::mat linkinv (const arma::mat & eta) const;
        arma::mat mueta (const arma::mat & eta) const;
};

}

#endif