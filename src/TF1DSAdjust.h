//
//  TF1DSAdjust.h
//  plastic_adjust
//
//  Created by Manouchehr Sanei on 4/8/19.
//
//

#ifndef TF1DSAdjust_h
#define TF1DSAdjust_h

#include <stdio.h>
#include "pzfmatrix.h"
#include "TPZSandlerExtended.h"
#include "StubFunctions.h"


class TF1DSAdjust
{
    
private:
    
    TPZSandlerExtended m_Sandler;
    TPZFMatrix<STATE> m_I1_SqJ2;
    
    REAL m_Aval;
    REAL m_Bval;
    REAL m_Cval;

    
public:
    
    /// Default constructor
    TF1DSAdjust();
    
    /// Copy constructor
    TF1DSAdjust(const TF1DSAdjust &other);
    
    /// Assignement constructor
    const TF1DSAdjust &operator=(const TF1DSAdjust &other);
    
    /// Destructor
    ~TF1DSAdjust();
    
    /// Method to initialize the class with possible data
    void Populate();
    
    
    /// Get the A value
    REAL Aval(){return m_Aval;}
    
    /// Set the A value
    void SetAval(REAL aval){m_Aval = aval;}
    
    /// Get the B value
    REAL Bval(){return m_Bval;}
    
    /// Set the B value
    void SetBval(REAL bval){m_Bval = bval;}
    
    /// Get the C value
    REAL Cval(){return m_Cval;}
    
    /// Set the C value
    void SetCval(REAL cval){m_Cval = cval;}
    
    /// Function to represent mean value of B
    REAL computeB_F1();
    
    /// Function to represent mean value of C
    REAL computeC_F1();
    
    /// Function to represent mean value of A
    REAL computeA_F1();
    
    /// Method to represent the error of cost function
    STATE errorfunction(const std::vector<STATE> &input);
    
    /// Method to adjust the parameters using NEWUOA unconstrained optimization (NLopt)
    void Adjust();
    
    
    /// Second derivative (Hessian) of objective function
    void Hessian(TPZFMatrix<REAL> & Hessian, REAL I1val, REAL SqJ2Val);
    
    /// First derivative (Gradient) of objective function
    void Residual(TPZFMatrix<REAL> &Residual, REAL I1val, REAL SqJ2Val);
    
    /// Method to assemble Hessian and Residual
    STATE Assemble(TPZFMatrix<REAL> &I1_SqJ2, TPZFMatrix<REAL> &hessian, TPZFMatrix<REAL> &res);
    
    /// Method to adjust the parameters using Quasi-Newton method
    void Adjust2();
    
    /// Method to modifiy the parameters
    void LoadCorrection(TPZFMatrix<REAL> &delx);
    
    

 

};


#endif /* TF1DSAdjust_h */
