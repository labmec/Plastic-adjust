//
//  TF2DSAdjust_R_RHW.cpp
//  plastic_adjust
//
//  Created by Manouchehr Sanei on 6/25/19.
//
//

#include "TF2DSAdjust_R_RHW.h"
#include <math.h>
#include "StubFunctions.h"
#include "TF1DSAdjust.h"
#include "nlopt.hpp"



TF2DSAdjust_R_RHW::TF2DSAdjust_R_RHW(){
    
    m_Lval= 0.;
}

TF2DSAdjust_R_RHW::TF2DSAdjust_R_RHW(const TF2DSAdjust_R_RHW &other) : m_Sandler(other.m_Sandler)
{
    m_Lval = other.m_Lval;
    
}

const TF2DSAdjust_R_RHW & TF2DSAdjust_R_RHW::operator=(const TF2DSAdjust_R_RHW &other)
{
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_Sandler = other.m_Sandler;
    m_Lval = other.m_Lval;
    return *this;
}

TF2DSAdjust_R_RHW::~TF2DSAdjust_R_RHW(){
    
}

void TF2DSAdjust_R_RHW::PopulateR()
{
    /// Render the data for failure function of DiMaggio-Sandler model
    REAL a_val = 15.091;
    REAL b_val = 0.0284035;
    REAL c_val = 15.9158;
    
    /// Oedometric loading (6620,9618)
    REAL E  = 2214.77;
    REAL nu = 0.25141;
    
    
    /// Triaxial unloading (9618,9680)
//    REAL E  = 5237.29;
//    REAL nu = 0.239828;
    
    m_ER.SetEngineeringData(E, nu);
    
    m_Sandler.SetA(a_val);
    m_Sandler.SetB(b_val);
    m_Sandler.SetC(c_val);
    
    /// Initial guess for L parameter
    REAL l_val = -51.0092;
    SetLval(l_val);
    
}


void TF2DSAdjust_R_RHW::ComputedSigmaTr(TPZVec<TTestSection> &active, TPZFMatrix<REAL> &stresstr)
{
    TPZFMatrix<REAL> deform,stress_meas;
    stresstr.Zero();
    
    for(int64_t i=0; i<active.size(); i++)
    {
        active[i].GetData(deform,stress_meas);
    
        int64_t ndata = deform.Rows();
        stresstr.Resize(ndata-1, 2);
        TPZManVector<REAL,2> eps(2),eps_n(2),sig(2),sig_n(2),sigtr(2);
        
        for(int64_t d = 1; d<ndata; d++)
        {
            eps[0] = deform(d-1,0);
            eps[1] = deform(d-1,1);
            
            eps_n[0] = deform(d,0);
            eps_n[1] = deform(d,1);
            
            sig[0] = stress_meas(d-1,0);
            sig[1] = stress_meas(d-1,1);
            
            sig_n[0] = stress_meas(d,0);
            sig_n[1] = stress_meas(d,1);
            
            Compute(eps,eps_n,sig,sigtr);
            
            REAL difsig1 = sig_n[0]-sigtr[0];
            REAL difsig2 = sig_n[1]-sigtr[1];
            
            stresstr(d-1,0) = difsig1;
            stresstr(d-1,1) = difsig2;
        }

    }
}

void TF2DSAdjust_R_RHW::Compute(TPZVec<REAL> &eps, TPZVec<REAL> &eps_n, TPZVec<REAL> &sig, TPZVec<REAL> &sigtr)
{
    sigtr[0] = m_ER.Lambda()*(2.*(eps_n[0]-eps[0])+(eps_n[1]-eps[1]))+2.*m_ER.Mu()*(eps_n[0]-eps[0])+sig[0];
    sigtr[1] = m_ER.Lambda()*(2.*(eps_n[0]-eps[0])+(eps_n[1]-eps[1]))+2.*m_ER.Mu()*(eps_n[1]-eps[1])+sig[1];

}


REAL TF2DSAdjust_R_RHW::ComputeR_analytic(STATE &sigst1, STATE &sigst2){
    
    STATE A     = m_Sandler.A();
    STATE B     = m_Sandler.B();
    STATE C     = m_Sandler.C();
    STATE Lprev = Lval();
    
    STATE term1, term2;
    
    term1 = 2*(pow(Lprev,2) - 2*sqrt(3)*Lprev*sigst1 + 3*pow(sigst1,2));
    term2 = 2*pow(A,2) - 4*A*C*exp(B*Lprev) + 2*pow(C,2)*exp(2*B*Lprev) - pow(sigst2,2);
    
    STATE r = sqrt(term1/term2);
    
    return r;
}



void TF2DSAdjust_R_RHW::grad_R(REAL & gradient, REAL sigTst1, REAL sigTst2, REAL sigst1, REAL sigst2){

    
    STATE A     = m_Sandler.A();
    STATE B     = m_Sandler.B();
    STATE C     = m_Sandler.C();
    
    gradient = (2*sqrt(3)*B*C*exp(B*m_Lval)*(-m_Lval + sqrt(3)*sigTst1)*sigTst2*
                      (2*pow(A,2) - 4*A*C*exp(B*m_Lval) + 2*pow(C,2)*exp(2*B*m_Lval) - pow(sigTst2,2)))/
    (pow(A - C*exp(B*m_Lval),5)*(pow(m_Lval,2) - 2*sqrt(3)*m_Lval*sigTst1 + 3*pow(sigTst1,2))*
     (pow(sigTst2,2)/pow(A - C*exp(B*m_Lval),4) +
      (3*pow(-m_Lval + sqrt(3)*sigTst1,2)*pow(2*pow(A,2) - 4*A*C*exp(B*m_Lval) + 2*pow(C,2)*exp(2*B*m_Lval) -
                                             pow(sigTst2,2),2))/(pow(A - C*exp(B*m_Lval),4)*pow(pow(m_Lval,2) - 2*sqrt(3)*m_Lval*sigTst1 + 3*pow(sigTst1,2),2)))) -
    (sigTst2*((sqrt(3)*(-4*A*B*C*exp(B*m_Lval) + 4*B*pow(C,2)*exp(2*B*m_Lval))*(-m_Lval + sqrt(3)*sigTst1))/
              (pow(A - C*exp(B*m_Lval),2)*(pow(m_Lval,2) - 2*sqrt(3)*m_Lval*sigTst1 + 3*pow(sigTst1,2))) -
              (sqrt(3)*(2*m_Lval - 2*sqrt(3)*sigTst1)*(-m_Lval + sqrt(3)*sigTst1)*
               (2*pow(A,2) - 4*A*C*exp(B*m_Lval) + 2*pow(C,2)*exp(2*B*m_Lval) - pow(sigTst2,2)))/
              (pow(A - C*exp(B*m_Lval),2)*pow(pow(m_Lval,2) - 2*sqrt(3)*m_Lval*sigTst1 + 3*pow(sigTst1,2),2)) -
              (sqrt(3)*(2*pow(A,2) - 4*A*C*exp(B*m_Lval) + 2*pow(C,2)*exp(2*B*m_Lval) - pow(sigTst2,2)))/
              (pow(A - C*exp(B*m_Lval),2)*(pow(m_Lval,2) - 2*sqrt(3)*m_Lval*sigTst1 + 3*pow(sigTst1,2))) +
              (2*sqrt(3)*B*C*exp(B*m_Lval)*(-m_Lval + sqrt(3)*sigTst1)*
               (2*pow(A,2) - 4*A*C*exp(B*m_Lval) + 2*pow(C,2)*exp(2*B*m_Lval) - pow(sigTst2,2)))/
              (pow(A - C*exp(B*m_Lval),3)*(pow(m_Lval,2) - 2*sqrt(3)*m_Lval*sigTst1 + 3*pow(sigTst1,2)))))/
    (pow(A - C*exp(B*m_Lval),2)*(pow(sigTst2,2)/pow(A - C*exp(B*m_Lval),4) +
                                     (3*pow(-m_Lval + sqrt(3)*sigTst1,2)*pow(2*pow(A,2) - 4*A*C*exp(B*m_Lval) + 2*pow(C,2)*exp(2*B*m_Lval) -
                                                                            pow(sigTst2,2),2))/(pow(A - C*exp(B*m_Lval),4)*pow(pow(m_Lval,2) - 2*sqrt(3)*m_Lval*sigTst1 + 3*pow(sigTst1,2),2))));

}


void TF2DSAdjust_R_RHW::res_R(REAL & residual, REAL sigTst1, REAL sigTst2, REAL sigst1, REAL sigst2){

    STATE A     = m_Sandler.A();
    STATE B     = m_Sandler.B();
    STATE C     = m_Sandler.C();

    residual= - atan(sigst2/sigst1) + atan((sigTst2/pow(A - C*exp(B*m_Lval),2))/((sqrt(3)*(-m_Lval + sqrt(3)*sigTst1)*
                                                                                         (2*pow(A,2) - 4*A*C*exp(B*m_Lval) + 2*pow(C,2)*exp(2*B*m_Lval) - pow(sigTst2,2)))/
                                                                                        (pow(A - C*exp(B*m_Lval),2)*(pow(m_Lval,2) - 2*sqrt(3)*m_Lval*sigTst1 + 3*pow(sigTst1,2)))));
}



void TF2DSAdjust_R_RHW::AdjustL(TPZVec<TTestSection> &active, TPZFMatrix<REAL> &stresstr, TPZFMatrix<REAL> &Ldata)
{
    TPZFMatrix<REAL> deform,stress_meas;
    Ldata.Zero();
    
    for(int64_t i=0; i<active.size(); i++)
    {
        active[i].GetData(deform,stress_meas);
        ComputedSigmaTr(active,stresstr);
        
        int64_t ndata = stresstr.Rows();
        Ldata.Resize(ndata, 1);
        
        REAL sigst1,sigst2,sigTst1,sigTst2;
        TPZVec<REAL> dirRHW;
        REAL residual, gradient;
        
        
        for(int64_t d = 0; d<ndata; d++)
        {
            sigst1  = (2*stresstr(d,0)+stresstr(d,1))/sqrt(3);
            sigst2  = sqrt(2)*sqrt(((stresstr(d,1)-stresstr(d,0))*(stresstr(d,1)-stresstr(d,0)))/3);
            
            sigTst1 = (2*stress_meas(d,0)+stress_meas(d,1))/sqrt(3);
            sigTst2 = sqrt(2)*sqrt(((stress_meas(d,1)-stress_meas(d,0))*(stress_meas(d,1)-stress_meas(d,0)))/3);
            
            /// To compute L in RHW, Newton and Gradient descent are implemented.
            bool newton_method = true;
            
            int count = 0;
            REAL mdif = 0.01;
            REAL tol = mdif*1.e-6;
            REAL gamma = 0.1;
            
            while(mdif > tol && count <20)
            {
                if (newton_method) {
                    
                    res_R(residual,sigTst1,sigTst2,sigst1,sigst2);
                    grad_R(gradient,sigTst1,sigTst2,sigst1,sigst2);
                    REAL term1 = residual;
                    REAL term2 = gradient;
                    
                    if (abs(term2)<=1.0e-7)
                    {
                        if (term2<0.0){
                            term2=-1.0e-7;}else{
                            term2=1.0e-7;}
                    }
                    
                    REAL delx  = - gamma*(term1/term2);
                    
                    LoadCorrectionL(delx);
                    mdif = abs(delx);
                    count++;
                    
                }else{
                    
                    grad_R(gradient,sigTst1,sigTst2,sigst1,sigst2);
                    REAL term2 = gradient;
                    
                    REAL delx  = - gamma*(term2);
                    LoadCorrectionL(delx);
                    mdif = abs(delx);
                    count++;

                }
            }

            REAL lval  = Lval();
            Ldata(d,0) = lval;
        }
    }
    
    Ldata.Print("Ldata = ", std::cout,EMathematicaInput);
}


void TF2DSAdjust_R_RHW::LoadCorrectionL(REAL & delx)
{
    REAL l_val = Lval()+delx;
    SetLval(l_val);

}


void TF2DSAdjust_R_RHW::AdjustR(TPZVec<TTestSection> &active)
{
    TPZFMatrix<REAL> deform,stress_meas;
    TPZFMatrix<REAL> stresstr, Ldata;
    AdjustL(active,stresstr,Ldata);
    int64_t ndata = stresstr.Rows();
    TPZFMatrix<REAL> Rdata(ndata,2);
    
    for(int64_t i=0; i<active.size(); i++)
    {
        active[i].GetData(deform,stress_meas);
        REAL sigTst1,sigTst2;
        
        for(int64_t d = 0; d<ndata; d++)
        {
            sigTst1 = (2*stress_meas(d,0)+stress_meas(d,1))/sqrt(3);
            sigTst2 = sqrt(2)*sqrt(((stress_meas(d,1)-stress_meas(d,0))*(stress_meas(d,1)-stress_meas(d,0)))/3);
            
            REAL lval = Ldata(d,0);
            SetLval(lval);
            
            REAL rval = ComputeR_analytic(sigTst1, sigTst2);
            
            Rdata(d,0) = d+1;
            Rdata(d,1) = rval;
            
        }

    }
    
    Rdata.Print("Rdata = ",std::cout,EMathematicaInput);
}



