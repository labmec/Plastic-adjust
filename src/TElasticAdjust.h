//
//  TElasticAdjust.h
//  plastic_adjust
//
//  Created by Philippe Devloo & Manouchehr Sanei on 03/03/19.
//
//

#ifndef TElasticAdjust_h
#define TElasticAdjust_h

#include <stdio.h>
#include "TPZElasticResponse.h"
#include "fadType.h"


class TTestSection;

class TElasticAdjust
{
public:
    enum EAdjustType {EYoungPoisson,ELambdaMu};
    
private:
    
    TPZElasticResponse fER;
    TPZVec<REAL> fSig0;
    
    REAL Assemble(TPZVec<TTestSection> &active, TPZFMatrix<REAL> &tangent, TPZFMatrix<REAL> &res);
    void LoadCorrection(TPZFMatrix<REAL> &delu);

    
public:
    
    /// Default constructor
    TElasticAdjust();
    
    /// Copy constructor
    TElasticAdjust(const TElasticAdjust &other);
    
    /// Assignement constructor
    const TElasticAdjust &operator=(const TElasticAdjust &other);
    
    /// Destructor
    ~TElasticAdjust();
    
    /// Method to present elastic response
    TPZElasticResponse ElasticResponse()
    {
        return fER;
    }
    
    /// Method to set elastic parameters
    void SetYoungPoisson(REAL E, REAL poisson)
    {
        fER.SetEngineeringData(E, poisson);
    }
    
    /// Method to compute sigma and dsigma
    void Compute(TPZVec<REAL> &eps, TPZVec<REAL> &sig, TPZVec<REAL> &dsigr, TPZVec<REAL> &dsiga);
    
    /// Method to adjust elastic parameters
    void Adjust(TPZVec<TTestSection> &active);

    
    /// Method to compute stress
    void ComputedSigma(TTestSection &sec, int i, TPZMatrix<REAL> &stress);
    
    
};



#endif /* TElasticAdjust_h */
