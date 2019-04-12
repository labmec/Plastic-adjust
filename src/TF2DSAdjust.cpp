//
//  TF2DSAdjust.cpp
//  plastic_adjust
//
//  Created by Manouchehr Sanei on 4/11/19.
//
//

#include "TF2DSAdjust.h"
#include "StubFunctions.h"



TF2DSAdjust::TF2DSAdjust(){
    
    m_Wval= 0.;
    m_Dval= 0.;
    
}

TF2DSAdjust::TF2DSAdjust(const TF2DSAdjust &other) : m_Sandler(other.m_Sandler)
{
    m_Wval = other.m_Wval;
    m_Dval = other.m_Dval;

}

const TF2DSAdjust & TF2DSAdjust::operator=(const TF2DSAdjust &other)
{
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_Sandler = other.m_Sandler;
    
    m_Wval = other.m_Wval;
    m_Dval = other.m_Dval;
    
    return *this;
}

TF2DSAdjust::~TF2DSAdjust(){
    
}


void TF2DSAdjust::Populate()
{
    
    // Set elastic data
    TPZElasticResponse ER;
    REAL E  = 5237.29;
    REAL nu = 0.239828;
    ER.SetEngineeringData(E, nu);
    
     m_Sandler.MCormicRanchSand(m_Sandler);
    
}


void TF2DSAdjust::ComputedElasticStrain(TPZVec<TTestSection> &active, std::vector<REAL> &epsEV_data)
{
    REAL lamda = m_ER.Lambda();
    REAL mu    = m_ER.G();

    TPZFMatrix<REAL> deform,stress;
    for(int64_t i=0; i<active.size(); i++)
    {
        active[i].GetData(deform,stress);
        int64_t ndata = stress.Rows();
        TPZManVector<REAL,2> sig_measure(2);

        for(int64_t d = 0; d<ndata; d++)
        {
            sig_measure[0] = stress(d,0);
            sig_measure[1] = stress(d,1);
            
            REAL epsEV = (2*sig_measure[0]+sig_measure[1])/(2*mu+3*lamda);
            
            epsEV_data.push_back(epsEV);
   
        }
    }
}

void TF2DSAdjust::ComputedPlasticStrain(TPZVec<TTestSection> &active, std::vector<REAL> &epsPV_data)
{
    std::vector<REAL> epsElV_data;
    std::vector<REAL> epsTV_data;
    TPZFMatrix<REAL> deform,stress;
    ComputedElasticStrain(active, epsElV_data);
    
    for(int64_t i=0; i<active.size(); i++)
    {
        active[i].GetData(deform,stress);
        int64_t ndata = deform.Rows();
        TPZManVector<REAL,2> epst(2);
        
        for(int64_t d = 0; d<ndata; d++)
        {
            epst[0] = deform(d,0);
            epst[1] = deform(d,1);
            
            REAL epsTV = 2*epst[0]+epst[1];
            
            epsTV_data.push_back(epsTV);
            
        }
    }
    
    epsPV_data = a_subtract_b(epsTV_data,epsElV_data);
    
}

std::vector<REAL> TF2DSAdjust::a_subtract_b(std::vector<REAL> & a, std::vector<REAL> & b){
    if (a.size()!=b.size()) {
        DebugStop();
    }
    std::vector<REAL> y(a.size());
    for (int i = 0; i < a.size(); i++) {
        y[i] = a[i]-b[i];
    }
    return y;
}

