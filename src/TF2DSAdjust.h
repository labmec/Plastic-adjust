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
    
    REAL m_Wval;
    REAL m_Dval;
    
public:
    
    /// Default constructor
    TF2DSAdjust();
    
    /// Copy constructor
    TF2DSAdjust(const TF2DSAdjust &other);
    
    /// Assignement constructor
    const TF2DSAdjust &operator=(const TF2DSAdjust &other);
    
    /// Destructor
    ~TF2DSAdjust();
    
    /// Method to initialize the class with possible data
    void Populate();
    
    /// Method to compute elastic strain
    void ComputedElasticStrain(TPZVec<TTestSection> &active, std::vector<REAL> &epsEV_data);
    
    /// Method to compute plastic strain
    void ComputedPlasticStrain(TPZVec<TTestSection> &active, std::vector<REAL> &epsPV_data);
    
    /// Sumation of two vector a and b
    std::vector<REAL> a_subtract_b(std::vector<REAL> & a, std::vector<REAL> & b);
};

#endif /* TF2DSAdjust_h */
