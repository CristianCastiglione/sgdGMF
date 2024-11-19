// family.h
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 19/11/2024

#ifndef FAMILY_H
#define FAMILY_H

#include <RcppArmadillo.h>
#include <memory>
#include "link.h"
#include "variance.h"

// void check_link (const std::string & familyname, const std::string & linkname);
void check_varf (const std::string & familyname, const std::string & varfname);

namespace glm {

class Family {
    private:
        std::unique_ptr<Link> linkobj;
        std::unique_ptr<Variance> varfobj;

    protected:
        std::string family; // glm-family name
        std::string link; // glm-link name
        std::string varf; // glm-variance function name

        bool estdispp; // should the dispersion be estimated?
        double dispersion; // dispersion parameter

    public:
        // Get family, link and variance names
        std::string getfamily() const {return this->family;}
        std::string getlink() const {return this->link;}
        std::string getvarf() const {return this->varf;}

        // Check, get and set the dispersion parameter
        bool estdisp() const {return this->estdispp;}
        double getdisp() const {return this->dispersion;}
        void setdisp(const double & disp) {this->dispersion = disp;}

        // Check if mu and eta are valid wrt the link and variance functions
        bool validmu (const arma::mat & mu) const {return varfobj->validmu(mu);}
        bool valideta (const arma::mat & eta) const {return linkobj->valideta(eta);}

        // Compute variance-related functions
        arma::mat initfun (const arma::mat & y) const {return varfobj->initfun(y);}
        arma::mat varfun (const arma::mat & mu) const {return varfobj->varfun(mu, this->dispersion);}
        arma::mat devfun (const arma::mat & y, const arma::mat & mu) const {
            return varfobj->devfun(y, mu, this->dispersion);
        }

        // Compute link-related functions
        arma::mat linkfun (const arma::mat & mu) const {return linkobj->linkfun(mu);}
        arma::mat linkinv (const arma::mat & eta) const {return linkobj->linkinv(eta);}
        arma::mat mueta (const arma::mat & eta) const {return linkobj->mueta(eta);}
        
        // Compute variance, initial and deviance values
        virtual arma::mat initialize (const arma::mat & y) const = 0;
        virtual arma::mat variance (const arma::mat & mu) const = 0;
        virtual arma::mat devresid (const arma::mat & y, const arma::mat & mu) const = 0;
        
        // Initialize the class using a link and a variance function
        Family (std::unique_ptr<Link> & link, std::unique_ptr<Variance> & varf) 
            : linkobj(std::move(link)), varfobj(std::move(varf)) {
            this->family = "Family";
            this->link = linkobj->link;
            this->varf = varfobj->varf;
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
        
        Gaussian (
            std::unique_ptr<Link> & link, 
            std::unique_ptr<Variance> & varf
        ) : Family(link, varf) {
            this->family = "Gaussian";
            this->varf = "const";
            this->estdispp = true;
        }
};

class Binomial : public Family {
    public:
        arma::mat variance (const arma::mat & mu) const;
        arma::mat initialize (const arma::mat & y) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        
        Binomial (
            std::unique_ptr<Link> & link, 
            std::unique_ptr<Variance> & varf
        ) : Family(link, varf) {
            this->family = "Binomial";
            this->varf = "mu(1-mu)";
        }
};

class Poisson : public Family {
    public:
        arma::mat variance (const arma::mat & mu) const;
        arma::mat initialize (const arma::mat & y) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        
        Poisson (
            std::unique_ptr<Link> & link, 
            std::unique_ptr<Variance> & varf
        ) : Family(link, varf) {
            this->family = "Poisson";
            this->varf = "mu";
        }
};

class Gamma : public Family {
    public:
        arma::mat variance (const arma::mat & mu) const;
        arma::mat initialize (const arma::mat & y) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        
        Gamma (
            std::unique_ptr<Link> & link, 
            std::unique_ptr<Variance> & varf
        ) : Family(link, varf) {
            this->family = "Gamma";
            this->varf = "mu^2";
            this->estdispp = true;
        }
};

class NegativeBinomial : public Family {
    public:
        arma::mat variance (const arma::mat & mu) const;
        arma::mat initialize (const arma::mat & y) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        
        NegativeBinomial (
            std::unique_ptr<Link> & link, 
            std::unique_ptr<Variance> & varf
        ) : Family(link, varf) {
            this->family = "NegativeBinomial";
            this->varf = "mu(1+t*mu)";
            this->estdispp = true;
            this->dispersion = 10;
        }
};

class QuasiBinomial : public Family {
    public:
        arma::mat variance (const arma::mat & mu) const;
        arma::mat initialize (const arma::mat & y) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        
        QuasiBinomial (
            std::unique_ptr<Link> & link, 
            std::unique_ptr<Variance> & varf
        ) : Family(link, varf) {
            this->family = "QuasiBinomial";
            this->varf = "mu(1-mu)";
            this->estdispp = true;
        }
};

class QuasiPoisson : public Family {
    public:
        arma::mat variance (const arma::mat & mu) const;
        arma::mat initialize (const arma::mat & y) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        
        QuasiPoisson (
            std::unique_ptr<Link> & link, 
            std::unique_ptr<Variance> & varf
        ) : Family(link, varf) {
            this->family = "QuasiPoisson";
            this->varf = "mu";
            this->estdispp = true;
        }
};

class Quasi : public Family {
    public:
        arma::mat variance (const arma::mat & mu) const;
        arma::mat initialize (const arma::mat & y) const;
        arma::mat devresid (const arma::mat & y, const arma::mat & mu) const;
        
        Quasi (
            std::unique_ptr<Link> & link, 
            std::unique_ptr<Variance> & varf
        ) : Family(link, varf) {
            this->family = "Quasi";
            this->estdispp = true;
        }
};

}

#endif