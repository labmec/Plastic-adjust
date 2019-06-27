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
    TPZElasticResponse m_ER;
    
    REAL m_Lval;

    
    
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
    
    /// Method to initialize the class for R with possible data
    void PopulateR();
    
    /// Method to compute strian and stress
    void Compute(TPZVec<REAL> &eps, TPZVec<REAL> &eps_n, TPZVec<REAL> &sig, TPZVec<REAL> &sigtr);
    
    /// Method to compute trial stress
    void ComputedSigmaTr(TPZVec<TTestSection> &active, TPZFMatrix<REAL> &stresstr);
    
    /// Method to compute analytically the R value of DS
    REAL ComputeR_analytic(STATE &sigst1, STATE &sigst2);
    
    /// First derivative (Gradient) of objective function  for R (or jacobian)
    void grad_R(REAL & gradient, REAL sigTst1, REAL sigTst2, REAL sigst1, REAL sigst2);
    
    /// Objective function for R (or residual)
    void res_R(REAL & residual, REAL sigTst1, REAL sigTst2, REAL sigst1, REAL sigst2);
    
    /// Method to assemble gradient and Residual for R
    void AdjustL(TPZVec<TTestSection> &active, TPZFMatrix<REAL> &stresstr, TPZFMatrix<REAL> &Ldata);
    
    /// Method to modifiy the parameter L
    void LoadCorrectionL(REAL & delx);
    
    /// Method to adjust the parameter R using Gradient descent method
    void AdjustR(TPZVec<TTestSection> &active);
    
    
    
};


#endif /* TF2DSAdjust_R_RHW_h */
