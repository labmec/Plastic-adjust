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
    
    
    /// Method to set elastic parameters
    void SetYoungPoisson(REAL E, REAL poisson)
    {
        m_ER.SetEngineeringData(E, poisson);
    }
    
    void SetLameDataElastic(REAL lambda, REAL mu)
    {
        m_ER.SetLameData(lambda, mu);
    }
    
    
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
    
    /// Second derivative (Hessian) of objective function for R
    void Hessian_R(TPZFMatrix<REAL> & Hessian, REAL I1val, REAL SqJ2Val);
    
    /// First derivative (Gradient) of objective function  for R
    void Residual_R(TPZFMatrix<REAL> &Residual, REAL I1val, REAL SqJ2Val);
    
    /// Method to assemble Hessian and Residual for R
    STATE AssembleR(TPZFMatrix<REAL> &I1_SqJ2, TPZFMatrix<REAL> &hessian, TPZFMatrix<REAL> &res);
    
    /// Method to adjust the parameter R using Quasi-Newton method
    void AdjustR();
    
    /// Method to compute the X values (the inputs are A,B,C and R)
    STATE ComputeXval(STATE &a, STATE &b, STATE &c, STATE &l, STATE &r);
    
    /// Method to modifiy the parameter R
    void LoadCorrectionR(TPZFMatrix<REAL> &delx);
    
    
    /// Second derivative (Hessian) of objective function for L (damage variable)
    void Hessian_LX(TPZFMatrix<REAL> & Hessian, REAL I1val, REAL SqJ2Val);
    
    /// First derivative (Gradient) of objective function for L (damage variable)
    void Residual_LX(TPZFMatrix<REAL> &Residual, REAL I1val, REAL SqJ2Val);
    
    
    /// Method to compute the X values from experimental data
    void ComputeXval(TPZVec<TTestSection> &active, std::vector<REAL> &I1data, std::vector<REAL> & SqJ2data, std::vector<REAL> &X_data);
    
    /// Method to combine the X values vs plastic strain deformation
    void XvalvsEpsPlasticdata(TPZVec<TTestSection> &active, std::vector<REAL> &epsPV_data, std::vector<REAL> &I1data, std::vector<REAL> & SqJ2data, std::vector<REAL> &X_data, TPZFMatrix<REAL> &X_epsP);
    
    /// Analytical function to compute the D value
    REAL ComputeDval_initial();
    
    /// Analytical function to compute the W value
    REAL ComputeWval_initial();
    
    /// Method to modifiy the parameter L
    void LoadCorrectionLX(TPZFMatrix<REAL> &delx);
    
    /// Method to initialize the class for L with possible data
    void PopulateDW();
    
    /// Method to initialize the class for R with possible data
    void PopulateR();
    
    
    /// Second derivative (Hessian) of objective function for parameters D and W
    void Hessian_DW(TPZFMatrix<REAL> &Hessian, REAL &X, REAL &depsPv);
    
    /// First derivative (Gradient) of objective function for parameters D and W
    void Residual_DW(TPZFMatrix<REAL> &Residual, REAL &X, REAL &depsPv);
    
    /// Method to assemble Hessian and Residual for parameters D and W
    STATE AssembleDW(TPZVec<TTestSection> &active, TPZFMatrix<REAL> &X_epsP, TPZFMatrix<REAL> &hessian, TPZFMatrix<REAL> &res);
    
    /// Method to make cost function for parameters D and W
    STATE DWcostFunction(REAL &X, REAL &depsPv);
    
    
    /// Method to adjust the parameters D and W using Quasi-Newton method
    void AdjustDW(TPZVec<TTestSection> &active);
    
    
    /// Method to modifiy the parameters D and W
    void LoadCorrectionDW(TPZFMatrix<REAL> &delx);
    
    /// Method to adjust the parameter X0 using Quasi-Newton method
    void AdjustX0();
    
    
    /// Method to find the maximum value of X
    void FindMaxXvalue(std::vector<REAL> &X_data, REAL &max_X);
    
    /// Method to compute elastic strain
    void ComputedElasticStrainLX(TPZVec<TTestSection> &active, std::vector<REAL> &epsEV_data);
    
    /// Method to compute plastic strain
    void ComputedPlasticStrainLX(TPZVec<TTestSection> &active, std::vector<REAL> &epsPV_data);
    
    /// Sumation of two vector a and b
    std::vector<REAL> a_subtract_b(std::vector<REAL> & a, std::vector<REAL> & b);
    
    
    void ComputedInvariantStressLX(TPZVec<TTestSection> &active, std::vector<REAL> &I1data, std::vector<REAL> & SqJ2data);
};

#endif /* TF2DSAdjust_h */
