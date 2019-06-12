//
//  TMohrAdjust.h
//  plastic_adjust
//
//  Created by Manouchehr Sanei on 4/23/19.
//
//

#ifndef TMohrAdjust_h
#define TMohrAdjust_h

#include <stdio.h>
#include "pzfmatrix.h"
#include "StubFunctions.h"
#include "TPZYCMohrCoulombPV.h"


class TMohrAdjust
{
    
private:
    
   
    TPZYCMohrCoulombPV m_mohr;
    TPZFMatrix<STATE> m_sign_tau;

    REAL m_coh;
    REAL m_phi;
    
public:
    
    /// Default constructor
    TMohrAdjust();
    
    /// Copy constructor
    TMohrAdjust(const TMohrAdjust &other);
    
    /// Assignement constructor
    const TMohrAdjust &operator=(const TMohrAdjust &other);
    
    /// Destructor
    ~TMohrAdjust();
    
    /// Get the cohesion
    REAL Cohval(){return m_coh;}
    
    /// Set the cohesion
    void SetCohval(REAL cohesion){m_coh = cohesion;}
    
    /// Get the cohesion
    REAL Phival(){return m_phi;}
    
    /// Set the friction
    void SetPhival(REAL friction){m_phi = friction;}
    
    /// Method to initialize the class with possible data
    void Populate();
    
    /// Function to represent tan(phi)
    REAL compute_phi_initial();
    
    /// Function to represent cohesion
    REAL compute_coh_initial();
    
    /// Method to represent the error of cost function for NLopt
    STATE errorfunctionMC(const std::vector<STATE> &input);
    
    /// Method to represent the gradient of object function for NLopt
    void gradientfunctionMC(const std::vector<STATE> &input, std::vector<double> &grad);
    
    /// Method to adjust the parameters phi and coh using NLopt libarary
    void Adjust();
    
    
    /// Second derivative (Hessian) of objective function
    void Hessian(TPZFMatrix<REAL> &Residual, REAL signVal, REAL tauVal);
    
    /// First derivative (Gradient) of objective function
    void Residual(TPZFMatrix<REAL> &Residual, REAL signVal, REAL tauVal);
    
    /// Method to assemble Hessian and Residual
    STATE Assemble(TPZFMatrix<REAL> &sign_tau, TPZFMatrix<REAL> &hessian, TPZFMatrix<REAL> &res);
    
    /// Method to adjust the parameters using Quasi-Newton method
    void Adjust2();
    
    /// Method to modifiy the parameters
    void LoadCorrection(TPZFMatrix<REAL> &delx);
    
    /// Method to compute shear stress
    STATE ComputeShearStress(REAL &sig_n, REAL &phi, REAL &coh);
    
    
    
};


#endif /* TMohrAdjust_h */
