// family.h
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 28/09/2023

#ifndef FAMILY_H
#define FAMILY_H

#include <RcppArmadillo.h>
#include "link.h"
#include <memory>

namespace Family {

template<class L>
class Family {
    public:
        std::string family;
        // std::string link = linkobj->link;
        arma::mat linkfun (const arma::mat & mu) const {return linkobj.linkfun(mu);}
        arma::mat linkinv (const arma::mat & eta) const {return linkobj.linkinv(eta);}
        arma::mat mueta (const arma::mat & eta) const {return linkobj.mueta(eta);}
        virtual arma::mat variance (const arma::mat & mu) const = 0;
        virtual arma::mat devresid (const arma::mat & y, const arma::mat & mu) const = 0;
        virtual arma::mat devresid (const arma::mat & y, const arma::mat & mu, const arma::mat & wt) const = 0;
        virtual arma::mat initialize (const arma::mat & y) const = 0;
        virtual bool validmu (const arma::mat & mu) const = 0;
        virtual bool valideta (const arma::mat & eta) const = 0;
        Family (const L & link) : linkobj(link) {}
        virtual ~Family () {}
    private:
        L linkobj;
};

template<class L>
class Gaussian : public Family<L> {
    public:
        std::string family = "Gaussian";
        arma::mat variance (const arma::mat & mu) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu, const arma::mat & wt) const;
        arma::mat initialize (const arma::mat & y) const;
        bool validmu (const arma::mat & mu) const;
        bool valideta (const arma::mat & eta) const;
        Gaussian (const L & link) : Family<L>(link) {}
};

template<class L>
class Binomial : public Family<L> {
    public:
        std::string family = "Binomial";
        arma::mat variance (const arma::mat & mu) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu, const arma::mat & wt) const;
        arma::mat initialize (const arma::mat & y) const;
        bool validmu (const arma::mat & mu) const;
        bool valideta (const arma::mat & eta) const;
        Binomial (const L & link) : Family<L>(link) {}
};

template<class L>
class Poisson : public Family<L> {
    public:
        std::string family = "Poisson";
        arma::mat variance (const arma::mat & mu) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu, const arma::mat & wt) const;
        arma::mat initialize (const arma::mat & y) const;
        bool validmu (const arma::mat & mu) const;
        bool valideta (const arma::mat & eta) const;
        Poisson (const L & link) : Family<L>(link) {}
};

template<class L>
class Gamma : public Family<L> {
    public:
        std::string family = "Gamma";
        arma::mat variance (const arma::mat & mu) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu, const arma::mat & wt) const;
        arma::mat initialize (const arma::mat & y) const;
        bool validmu (const arma::mat & mu) const;
        bool valideta (const arma::mat & eta) const;
        Gamma (const L & link) : Family<L>(link) {}
};

// Yet to implement: InverseGaussian, Quasi, QuasiBinomial, QuasiPoisson, NegativeBinomial, ...

}

#endif