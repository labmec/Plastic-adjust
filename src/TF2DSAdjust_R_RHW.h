//
//  TF2DSAdjust_R_RHW.h
//  plastic_adjust
//
//  Created by Manouchehr Sanei on 6/25/19.
//
//

#ifndef TF2DSAdjust_R_RHW_h
#define TF2DSAdjust_R_RHW_h

#include <stdio.h>
#include "pzfmatrix.h"
#include "TPZSandlerExtended.h"
#include "StubFunctions.h"


class TF2DSAdjust_R_RHW
{
    
private:
    
    TPZSandlerExtended m_Sandler;
    TPZFMatrix<STATE> m_I1_SqJ2;
    
    REAL m_Lval;
    REAL m_Rval;
    
    
public:
    
    /// Default constructor
    TF2DSAdjust_R_RHW();
    
    /// Copy constructor
    TF2DSAdjust_R_RHW(const TF2DSAdjust_R_RHW &other);
    
    /// Assignement constructor
    const TF2DSAdjust_R_RHW &operator=(const TF2DSAdjust_R_RHW &other);
    
    /// Destructor
    ~TF2DSAdjust_R_RHW();
    
    
    /// Get the L value
    REAL Lval(){return m_Lval;}
    
    /// Set the L value
    void SetLval(REAL lval){m_Lval = lval;}
    
    /// Get the R value
    REAL Rval(){return m_Rval;}
    
    /// Set the R value
    void SetRval(REAL rval){m_Rval = rval;}
    
    
    /// Method to initialize the class for R with possible data
    void PopulateR();
    
    /// Method to represent the error of cost function for NLopt
    STATE errorfunctionF2_R(const std::vector<STATE> &input);
    
    /// Method to represent the gradient of object function for NLopt
    void gradientfunctionF2_R(const std::vector<STATE> &input, std::vector<double> &grad);
    
    /// Method to adjust the parameters using NLopt optimization libarary
    void AdjustR();
    
    
    /// Second derivative (Hessian) of objective function for R
    void Hessian_R(TPZFMatrix<REAL> & Hessian, REAL I1val, REAL SqJ2Val);
    
    /// First derivative (Gradient) of objective function  for R
    void Residual_R(TPZFMatrix<REAL> &Residual, REAL I1val, REAL SqJ2Val);
    
    /// Method to assemble Hessian and Residual for R
    STATE AssembleR(TPZFMatrix<REAL> &I1_SqJ2, TPZFMatrix<REAL> &hessian, TPZFMatrix<REAL> &res);
    
    /// Method to adjust the parameter R using Quasi-Newton method
    void AdjustR2();
    
    /// Method to compute the X values (the inputs are A,B,C and R)
    STATE ComputeXval(STATE &a, STATE &b, STATE &c, STATE &l, STATE &r);
    
    /// Method to modifiy the parameter R
    void LoadCorrectionR(TPZFMatrix<REAL> &delx);
    
    /// Method to represent the Cap function
    STATE CapFunction(STATE &I1, STATE &SqJ2);
    
    
};


#endif /* TF2DSAdjust_R_RHW_h */
