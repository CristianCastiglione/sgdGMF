// family.h
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 13/11/2023

#ifndef FAMILY_H
#define FAMILY_H

#include <RcppArmadillo.h>
#include <memory>
#include "link.h"

namespace glm {

class Family {
    private:
        std::unique_ptr<Link> linkobj;

    protected:
        std::string family; // glm-family names
        std::string link; // glm-link name
        bool estdispp; // should the dispersion be estimated?
        double dispersion; // dispersion parameter

    public:
        std::string getfamily() const {return this->family;}
        std::string getlink() const {return this->link;}

        bool estdisp() const {return this->estdispp;}
        double getdisp() const {return this->dispersion;}
        void setdisp(const double & disp) {this->dispersion = disp;}

        arma::mat linkfun (const arma::mat & mu) const {return linkobj->linkfun(mu);};
        arma::mat linkinv (const arma::mat & eta) const {return linkobj->linkinv(eta);};
        arma::mat mueta (const arma::mat & eta) const {return linkobj->mueta(eta);};
        
        virtual arma::mat variance (const arma::mat & mu) const = 0;
        virtual arma::mat initialize (const arma::mat & y) const = 0;
        virtual arma::mat devresid (const arma::mat & y, const arma::mat & mu) const = 0;
        
        virtual bool validmu (const arma::mat & mu) const = 0;
        virtual bool valideta (const arma::mat & eta) const = 0;
        
        Family (std::unique_ptr<Link> & link) : linkobj(std::move(link)) {
            this->family = "Family";
            this->link = linkobj->link;
            this->estdispp = false;
            this->dispersion = 1;
        }
        
        virtual ~Family () {}
};

class Gaussian : public Family {
    public:
        arma::mat variance (const arma::mat & mu) const;
        arma::mat initialize (const arma::mat & y) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        bool validmu (const arma::mat & mu) const;
        bool valideta (const arma::mat & eta) const;
        
        Gaussian (std::unique_ptr<Link> & link) : Family(link) {
            this->family = "Gaussian";
            this->estdispp = true;
        }
};

class Binomial : public Family {
    public:
        arma::mat variance (const arma::mat & mu) const;
        arma::mat initialize (const arma::mat & y) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        bool validmu (const arma::mat & mu) const;
        bool valideta (const arma::mat & eta) const;
        
        Binomial (std::unique_ptr<Link> & link) : Family(link) {
            this->family = "Binomial";
        }
};

class Poisson : public Family {
    public:
        arma::mat variance (const arma::mat & mu) const;
        arma::mat initialize (const arma::mat & y) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        bool validmu (const arma::mat & mu) const;
        bool valideta (const arma::mat & eta) const;
        
        Poisson (std::unique_ptr<Link> & link) : Family(link) {
            this->family = "Poisson";
        }
};

class Gamma : public Family {
    public:
        arma::mat variance (const arma::mat & mu) const;
        arma::mat initialize (const arma::mat & y) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        bool validmu (const arma::mat & mu) const;
        bool valideta (const arma::mat & eta) const;
        
        Gamma (std::unique_ptr<Link> & link) : Family(link) {
            this->family = "Gamma";
            this->estdispp = true;
        }
};

class NegativeBinomial : public Family {
    public:
        arma::mat variance (const arma::mat & mu) const;
        arma::mat initialize (const arma::mat & y) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        bool validmu (const arma::mat & mu) const;
        bool valideta (const arma::mat & eta) const;
        
        NegativeBinomial (std::unique_ptr<Link> & link) : Family(link) {
            this->family = "NegativeBinomial";
            this->estdispp = true;
            this->dispersion = 10;
        }
};

class QuasiBinomial : public Family {
    public:
        arma::mat variance (const arma::mat & mu) const;
        arma::mat initialize (const arma::mat & y) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        bool validmu (const arma::mat & mu) const;
        bool valideta (const arma::mat & eta) const;
        
        QuasiBinomial (std::unique_ptr<Link> & link) : Family(link) {
            this->family = "QuasiBinomial";
            this->estdispp = true;
        }
};

class QuasiPoisson : public Family {
    public:
        arma::mat variance (const arma::mat & mu) const;
        arma::mat initialize (const arma::mat & y) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        bool validmu (const arma::mat & mu) const;
        bool valideta (const arma::mat & eta) const;
        
        QuasiPoisson (std::unique_ptr<Link> & link) : Family(link) {
            this->family = "QuasiPoisson";
            this->estdispp = true;
        }
};

}

#endif