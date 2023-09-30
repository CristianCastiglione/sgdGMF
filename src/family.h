// family.h
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 30/09/2023

#ifndef FAMILY_H
#define FAMILY_H

#include <RcppArmadillo.h>
#include <memory>
// memory permits to manage pointer objects of the form std::unique_ptr<typename>

#include "link.h"


namespace Family {

class Family {
    private:
        std::unique_ptr<Link::Link> linkobj;
    public:
        std::string family;
        std::string link = linkobj->link;
        double dispersion;
        arma::mat linkfun (const arma::mat & mu) const {return linkobj->linkfun(mu);};
        arma::mat linkinv (const arma::mat & eta) const {return linkobj->linkinv(eta);};
        arma::mat mueta (const arma::mat & eta) const {return linkobj->mueta(eta);};
        virtual arma::mat variance (const arma::mat & mu) const = 0;
        virtual arma::mat initialize (const arma::mat & y) const = 0;
        virtual arma::mat devresid (const arma::mat & y, const arma::mat & mu) const = 0;
        virtual bool validmu (const arma::mat & mu) const = 0;
        virtual bool valideta (const arma::mat & eta) const = 0;
        Family (std::unique_ptr<Link::Link> & link) : linkobj(std::move(link)) {}
        virtual ~Family () {}
};

class Gaussian : public Family {
    public:
        std::string family = "Gaussian";
        double dispersion = 1;
        arma::mat variance (const arma::mat & mu) const;
        arma::mat initialize (const arma::mat & y) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        bool validmu (const arma::mat & mu) const;
        bool valideta (const arma::mat & eta) const;
        Gaussian (std::unique_ptr<Link::Link> & link) : Family(link) {}
};

class Binomial : public Family {
    public:
        std::string family = "Binomial";
        double dispersion = 1;
        arma::mat variance (const arma::mat & mu) const;
        arma::mat initialize (const arma::mat & y) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        bool validmu (const arma::mat & mu) const;
        bool valideta (const arma::mat & eta) const;
        Binomial (std::unique_ptr<Link::Link> & link) : Family(link) {}
};

class Poisson : public Family {
    public:
        std::string family = "Poisson";
        double dispersion = 1;
        arma::mat variance (const arma::mat & mu) const;
        arma::mat initialize (const arma::mat & y) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        bool validmu (const arma::mat & mu) const;
        bool valideta (const arma::mat & eta) const;
        Poisson (std::unique_ptr<Link::Link> & link) : Family(link) {}
};

class Gamma : public Family {
    public:
        std::string family = "Gamma";
        double dispersion = 1;
        arma::mat variance (const arma::mat & mu) const;
        arma::mat initialize (const arma::mat & y) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        bool validmu (const arma::mat & mu) const;
        bool valideta (const arma::mat & eta) const;
        Gamma (std::unique_ptr<Link::Link> & link) : Family(link) {}
};

}

#endif