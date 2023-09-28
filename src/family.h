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

class Family {
    public:
        const std::string family;
        const std::string link = linkclass->link;
        arma::mat linkfun (const arma::mat & mu) const {return linkclass->linkfun(mu);}
        arma::mat linkinv (const arma::mat & eta) const {return linkclass->linkinv(eta);}
        arma::mat mueta (const arma::mat & eta) const {return linkclass->mueta(eta);}
        virtual arma::mat variance (const arma::mat & mu) const = 0;
        virtual arma::mat devresid (const arma::mat & y, const arma::mat & mu) const = 0;
        virtual arma::mat devresid (const arma::mat & y, const arma::mat & mu, const arma::mat & wt) const = 0;
        virtual arma::mat initialize (const arma::mat & y) const = 0;
        virtual bool validmu (const arma::mat & mu) const = 0;
        virtual bool valideta (const arma::mat & eta) const = 0;
        Family (std::unique_ptr<Link::Link> & linkobj) : linkclass(std::move(linkobj)) {}
        virtual ~Family () {}
    private:
        std::unique_ptr<Link::Link> linkclass;
};

class Gaussian : public Family {
    public:
        const std::string family = "Gaussian";
        arma::mat variance (const arma::mat & mu) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu, const arma::mat & wt) const;
        arma::mat initialize (const arma::mat & y) const;
        bool validmu (const arma::mat & mu) const;
        bool valideta (const arma::mat & eta) const;
        Gaussian (std::unique_ptr<Link::Link> & linkobj) : Family(linkobj) {}
};

class Binomial : public Family {
    public:
        const std::string family = "Binomial";
        arma::mat variance (const arma::mat & mu) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu, const arma::mat & wt) const;
        arma::mat initialize (const arma::mat & y) const;
        bool validmu (const arma::mat & mu) const;
        bool valideta (const arma::mat & eta) const;
        Binomial (std::unique_ptr<Link::Link> & linkobj) : Family(linkobj) {}
};

class Poisson : public Family {
    public:
        const std::string family = "Poisson";
        arma::mat variance (const arma::mat & mu) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu, const arma::mat & wt) const;
        arma::mat initialize (const arma::mat & y) const;
        bool validmu (const arma::mat & mu) const;
        bool valideta (const arma::mat & eta) const;
        Poisson (std::unique_ptr<Link::Link> & linkobj) : Family(linkobj) {}
};

class Gamma : public Family {
    public:
        const std::string family = "Gamma";
        arma::mat variance (const arma::mat & mu) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu, const arma::mat & wt) const;
        arma::mat initialize (const arma::mat & y) const;
        bool validmu (const arma::mat & mu) const;
        bool valideta (const arma::mat & eta) const;
        Gamma (std::unique_ptr<Link::Link> & linkobj) : Family(linkobj) {}
};

// Yet to implement: InverseGaussian, Quasi, QuasiBinomial, QuasiPoisson, NegativeBinomial, ...

}

#endif