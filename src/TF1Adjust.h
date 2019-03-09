//
//  TF1Adjust.hpp
//  plastic_adjust
//
//  Created by Philippe Devloo on 03/03/19.
//

#ifndef TF1Adjust_hpp
#define TF1Adjust_hpp

#include <stdio.h>
#include "pzfmatrix.h"
#include "TPZSandlerExtended.h"
#include "StubFunctions.h"

class TF1Adjust
{
    
    TPZSandlerExtended fSandler;
    
    TPZFMatrix<STATE> fI1_SqJ2;
    
public:
    
    TF1Adjust(){}
    
    TF1Adjust(const TF1Adjust &copy) : fSandler(copy.fSandler), fI1_SqJ2(copy.fI1_SqJ2)
    {
        
    }
    
    TF1Adjust &operator=(const TF1Adjust &copy)
    {
        fSandler = copy.fSandler;
        fI1_SqJ2 = copy.fI1_SqJ2;
        return *this;
    }
    
    /// method to initialize the class with possible data
    void Populate();
    
    void Adjust(TPZVec<TTestSection> &active);
    
    void Adjust();
    
    void ComputedSigma(STATE I1Min, STATE I1max, STATE delI1, TPZMatrix<REAL> &I1_SqJ2);
    
    STATE errorfunction(const std::vector<STATE> &input);

    STATE errorfunction2(const std::vector<STATE> &input);
};

#endif /* TF1Adjust_hpp */
