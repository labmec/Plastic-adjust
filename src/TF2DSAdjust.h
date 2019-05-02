//
//  TF2DSAdjust.h
//  plastic_adjust
//
//  Created by Manouchehr Sanei on 4/11/19.
//
//

#ifndef TF2DSAdjust_h
#define TF2DSAdjust_h

#include <stdio.h>
#include "pzfmatrix.h"
#include "TPZSandlerExtended.h"
#include "StubFunctions.h"

class TF2DSAdjust
{
    
private:
    
    TPZElasticResponse m_ER;
    
    TPZSandlerExtended m_Sandler;
    TPZFMatrix<STATE> m_I1_SqJ2;
    
    
    REAL m_Wval;
    REAL m_Dval;
    
    REAL m_Lval;
    REAL m_Rval;
    
    
    
public:
    
    /// Default constructor
    TF2DSAdjust();
    
    /// Copy constructor
    TF2DSAdjust(const TF2DSAdjust &other);
    
    /// Assignement constructor
    const TF2DSAdjust &operator=(const TF2DSAdjust &other);
    
    /// Destructor
    ~TF2DSAdjust();
    
    
    /// Get the L value
    REAL Lval(){return m_Lval;}
    
    /// Set the L value
    void SetLval(REAL lval){m_Lval = lval;}
    
    /// Get the R value
    REAL Rval(){return m_Rval;}
    
    /// Set the R value
    void SetRval(REAL rval){m_Rval = rval;}
    
    
    /// Get the D value
    REAL Dval(){return m_Dval;}
    
    /// Set the D value
    void SetDval(REAL dval){m_Dval = dval;}
    
    
    /// Get the W value
    REAL Wval(){return m_Wval;}
    
    /// Set the R value
    void SetWval(REAL wval){m_Wval = wval;}
    
    
    /// Method to represent the error of cost function for L and R
    STATE errorfunctionF2(const std::vector<STATE> &input);
    
    
    /// Method to represent the Cap function
    STATE CapFunction(STATE &I1, STATE &SqJ2);
    
    /// Second derivative (Hessian) of objective function
    void HessianR(TPZFMatrix<REAL> & Hessian, REAL I1val, REAL SqJ2Val);
    
    /// First derivative (Gradient) of objective function
    void ResidualR(TPZFMatrix<REAL> &Residual, REAL I1val, REAL SqJ2Val);
    
    /// Method to assemble Hessian and Residual
    STATE AssembleR(TPZFMatrix<REAL> &I1_SqJ2, TPZFMatrix<REAL> &hessian, TPZFMatrix<REAL> &res);
    
    /// Method to adjust the parameters using Quasi-Newton method
    void AdjustR();
    
    STATE ComputeXval(STATE &a, STATE &b, STATE &c, STATE &l, STATE &r);
    
    /// Method to modifiy the parameters
    void LoadCorrectionR(TPZFMatrix<REAL> &delx);
    
    
    /// Second derivative (Hessian) of objective function
    void HessianL(TPZFMatrix<REAL> & Hessian, REAL I1val, REAL SqJ2Val);
    
    /// First derivative (Gradient) of objective function
    void ResidualL(TPZFMatrix<REAL> &Residual, REAL I1val, REAL SqJ2Val);
    
    
    void ComputeXval(TPZVec<TTestSection> &active, TPZFMatrix<REAL> &I1_SqJ2, std::vector<REAL> &X_data);
    
    void XvalvsEpsPlasticdata(TPZVec<TTestSection> &active, std::vector<REAL> &epsPV_data, TPZFMatrix<REAL> &I1_SqJ2, std::vector<REAL> &L_data, TPZFMatrix<REAL> &X_epsP);
    
    
    REAL ComputeDval_in();
    
    REAL ComputeWval_in();
    
    /// Method to modifiy the parameters
    void LoadCorrectionL(TPZFMatrix<REAL> &delx);
    
    /// Method to initialize the class with possible data
    void PopulateL();
    
    /// Method to initialize the class with possible data
    void PopulateR();
    
    
    /// Second derivative (Hessian) of objective function
    void HessianDW(TPZFMatrix<REAL> &Hessian, REAL &X, REAL &depsPv);
    
    /// First derivative (Gradient) of objective function
    void ResidualDW(TPZFMatrix<REAL> &Residual, REAL &X, REAL &depsPv);
    
    STATE AssembleDW(TPZFMatrix<REAL> &X_epsP, TPZFMatrix<REAL> &hessian, TPZFMatrix<REAL> &res);
    
    STATE DWcostFunction(REAL &X, REAL &depsPv);
    
    
    /// Method to adjust the parameters using Quasi-Newton method
    void AdjustDW();
    
    
    /// Method to modifiy the parameters
    void LoadCorrectionDW(TPZFMatrix<REAL> &delx);
    
    /// Method to adjust the parameters using Quasi-Newton method
    void AdjustX0();
    
    
    void FindMaxXvalue(std::vector<REAL> &X_data, REAL &max_X);
    
    /// Method to compute elastic strain
    void ComputedElasticStrain(TPZVec<TTestSection> &active, std::vector<REAL> &epsEV_data);
    
    /// Method to compute plastic strain
    void ComputedPlasticStrain(TPZVec<TTestSection> &active, std::vector<REAL> &epsPV_data);
    
    /// Sumation of two vector a and b
    std::vector<REAL> a_subtract_b(std::vector<REAL> & a, std::vector<REAL> & b);
    
    
    void ComputedInvariantStress(TPZVec<TTestSection> &active, TPZFMatrix<STATE> &I1SqJ2data);
};

#endif /* TF2DSAdjust_h */
