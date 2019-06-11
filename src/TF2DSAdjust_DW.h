//
//  TF2DSAdjust_DW.h
//  plastic_adjust
//
//  Created by Manouchehr Sanei on 4/11/19.
//
//

#ifndef TF2DSAdjust_DW_h
#define TF2DSAdjust_DW_h

#include <stdio.h>
#include "pzfmatrix.h"
#include "StubFunctions.h"

class TF2DSAdjust_DW
{
    
private:
    
    TPZElasticResponse m_ER;
    
    TPZSandlerExtended m_Sandler;
    TPZFMatrix<STATE> m_X_epsPv;
    
    REAL m_Kbulk;
    REAL m_Dval;
    REAL m_Wval;
    REAL m_iWval;
   

public:
    
    /// Default constructor
    TF2DSAdjust_DW();
    
    /// Copy constructor
    TF2DSAdjust_DW(const TF2DSAdjust_DW &other);
    
    /// Assignement constructor
    const TF2DSAdjust_DW &operator=(const TF2DSAdjust_DW &other);
    
    /// Destructor
    ~TF2DSAdjust_DW();
    
    
    /// Set the Bulk modulus
    void SetBulkModulus(REAL K)
    {
        m_Kbulk = K;
    }
    
    
    /// Get the D value
    REAL Dval(){return m_Dval;}
    
    /// Set the D value
    void SetDval(REAL dval){m_Dval = dval;}
    
    
    /// Get the W value
    REAL Wval(){return m_Wval;}
    
    /// Set the W value
    void SetWval(REAL wval){m_Wval = wval;}
    
    /// Get the (iW = Log(DW)) value
    REAL iWval(){return m_iWval;}
    
    /// Set the (iW = Log(DW)) value
    void SetiWval(REAL iwval){m_iWval = iwval;}
    
    

    
    ///*********************************************************
    
    /// Method to compute elastic strain
    void ComputedElasticStrain(TPZVec<TTestSection> &active, std::vector<REAL> &epsEV_data);
    
    /// Method to compute plastic strain
    void ComputedPlasticStrain(TPZVec<TTestSection> &active, std::vector<REAL> &epsPV_data);
    
    
    /// Sumation of two vector a and b
    std::vector<REAL> a_subtract_b(std::vector<REAL> & a, std::vector<REAL> & b);
    
    /// Show data of X vs epsPv to be able to plot the results
    void Data_X_vs_EpsPv(TPZVec<TTestSection> &active, std::vector<REAL> &epsPV_data, std::vector<REAL> &X_data,TPZFMatrix<REAL> &X_epsPv);
    
    /// Method to initialize the elastic parameter
    void PopulateElastic();
    
    /// Plot  X vs epsPv
    void PlotXvsEpsPv(TPZVec<TTestSection> &active);
    
    
    ///*********************************************************

    /// Method to initialize the class for L with possible data
    void PopulateDW();
    
    /// Analytical function to compute the D value
    REAL ComputeDval_initial();
    
    /// Analytical function to compute the W value
    REAL ComputeWval_initial();
    
    
    /// Second derivative (Hessian) of objective function for parameters D and W
    void Hessian_DW(TPZFMatrix<REAL> &Hessian, REAL &X, REAL &depsPvX);
    
    /// First derivative (Gradient) of objective function for parameters D and W
    void Residual_DW(TPZFMatrix<REAL> &Residual, REAL &X, REAL &depsPvX);
    
    /// Method to assemble Hessian and Residual for parameters D and W
    STATE AssembleDW(TPZFMatrix<REAL> &X_epsPv, TPZFMatrix<REAL> &hessian, TPZFMatrix<REAL> &res);
    
    /// Method to make cost function for parameters D and W
    STATE costFunction_DW(REAL &X, REAL &depsPvX);
    
        
    /// Method to adjust the parameters D and W using Quasi-Newton method
    void AdjustDW();
    
    
    /// Method to modifiy the parameters D and W
    void LoadCorrection_DW(TPZFMatrix<REAL> &delx);
    


};

#endif /* TF2DSAdjust_DW_h */
