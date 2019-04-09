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
    
    REAL m_A_val;
    REAL m_B_val;
    REAL m_C_val;

    
public:
    
    /// Default constructor
    TF1DSAdjust();
    
    /// Copy constructor
    TF1DSAdjust(const TF1DSAdjust &other);
    
    /// Assignement constructor
    const TF1DSAdjust &operator=(const TF1DSAdjust &other);
    
    /// Destructor
    ~TF1DSAdjust();
    
    
    /// Gets the A value
    REAL GetAval(){return m_A_val;}
    
    /// Gets the B value
    REAL GetBval(){return m_B_val;}
    
    /// Gets the C value
    REAL GetCval(){return m_C_val;}
    
    /// Function to represent mean value of B
    void B_F1_function(TPZFMatrix<REAL> &I1_SqJ2, REAL &bMean);
    
    /// Function to represent mean value of C
    void C_F1_function(TPZFMatrix<REAL> &I1_SqJ2, REAL &cMean);
    
    /// Function to represent mean value of A
    void A_F1_function(TPZFMatrix<REAL> &I1_SqJ2, REAL &aMean);
    
    
    /// Method to initialize the class with possible data
    void Populate();
    
    /// Method to adjust the parameters
    void Adjust();
    
    /// Method to represent the error of cost function
    STATE errorfunction(const std::vector<STATE> &input);
    

};


#endif /* TF1DSAdjust_h */
