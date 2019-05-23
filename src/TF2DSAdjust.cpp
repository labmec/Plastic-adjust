//
//  TF2DSAdjust.cpp
//  plastic_adjust
//
//  Created by Manouchehr Sanei on 4/11/19.
//
//

#include <math.h>
#include "TF2DSAdjust.h"
#include "StubFunctions.h"
#include "TF1DSAdjust.h"




TF2DSAdjust::TF2DSAdjust(){
    
    m_Wval= 0.;
    m_Dval= 0.;
    m_Lval= 0.;
    m_Rval= 0.;
    
}

TF2DSAdjust::TF2DSAdjust(const TF2DSAdjust &other) : m_Sandler(other.m_Sandler), m_I1_SqJ2(other.m_I1_SqJ2)
{
    
    m_Wval = other.m_Wval;
    m_Dval = other.m_Dval;
    m_Lval = other.m_Lval;
    m_Rval = other.m_Rval;

}

const TF2DSAdjust & TF2DSAdjust::operator=(const TF2DSAdjust &other)
{
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_Sandler = other.m_Sandler;
    m_I1_SqJ2 = other.m_I1_SqJ2;
    
    m_Wval = other.m_Wval;
    m_Dval = other.m_Dval;
    m_Lval = other.m_Lval;
    m_Rval = other.m_Rval;
    
    return *this;
}

TF2DSAdjust::~TF2DSAdjust(){
    
}

/// First step: finding the R

void TF2DSAdjust::PopulateR()
{
    
    /// Render at least two points that touch the cap of DiMaggio-Sandler; temporary
    const int n_I1 = 2;
    
    REAL I1_data[n_I1] = {-56.9677, -73.3605};
    REAL F1_data[n_I1] = {0.0, 10.0713};
    
    m_I1_SqJ2.Redim(n_I1,2);
    
    for(int i=0; i<n_I1; i++)
    {
        REAL I1 = I1_data[i];
        REAL F1 = F1_data[i];
        m_I1_SqJ2(i,0) = I1;
        m_I1_SqJ2(i,1) = F1;
    }
    
    /// Render the data for failure function of DiMaggio-Sandler model
    REAL a_val = 15.091;
    REAL b_val = 0.0284035;
    REAL c_val = 15.9158;
   
    m_Sandler.SetA(a_val);
    m_Sandler.SetB(b_val);
    m_Sandler.SetC(c_val);
    
    
    /// Initial guess for L and R parameters
    REAL l_val = -60.9137409;
    SetLval(l_val);
    REAL r_val = 1.5;
    SetRval(r_val);

    
}

STATE TF2DSAdjust::CapFunction(STATE &I1, STATE &SqJ2){
    
    STATE A     = m_Sandler.A();
    STATE B     = m_Sandler.B();
    STATE C     = m_Sandler.C();
    STATE Lprev = Lval();
    STATE R     = Rval();

    /// The original DiMaggio is considered then the gamma = 1.0
    STATE gamma = 1.0;
    STATE term1, term2;
    
    term1 = (I1-Lprev) / (R * (A-(C * exp(B * Lprev))));
    term2 = (gamma * SqJ2) / (A-(C * exp(B * Lprev)));
    
    STATE f2 = (term1 * term1 + term2 * term2 - 1);
    
    return f2;
    
}

//STATE TF2DSAdjust::errorfunctionF2(const std::vector<STATE> &input)
//{
//    TPZFMatrix<REAL> I1_SqJ2 = m_I1_SqJ2;
//    int64_t n_data = m_I1_SqJ2.Rows();
//    
//    STATE errorF2=0;
//    
//    for(int i=0; i<n_data; i++)
//    {
//        STATE minF2 = CapFunction(I1_SqJ2(i,0), I1_SqJ2(i,1));
//        
//        errorF2 += minF2*minF2;
//        
//    }
//    
//    return errorF2;
//}


void TF2DSAdjust::Hessian_R(TPZFMatrix<REAL> & Hessian, REAL I1val, REAL SqJ2Val){
    
    Hessian.Resize(2,2);
    Hessian.Zero();
    
    STATE Aval = m_Sandler.A();
    STATE Bval = m_Sandler.B();
    STATE Cval = m_Sandler.C();
    
    
    Hessian(0,0) += 2*pow((-2*(I1val - m_Lval))/(pow(Aval - Cval*exp(Bval*m_Lval),2)*pow(m_Rval,2)) +
                            (2*Bval*Cval*exp(Bval*m_Lval)*pow(I1val - m_Lval,2))/(pow(Aval - Cval*exp(Bval*m_Lval),3)*pow(m_Rval,2)) +
                            (2*Bval*Cval*exp(Bval*m_Lval)*pow(SqJ2Val,2))/pow(Aval - Cval*exp(Bval*m_Lval),3),2) +
    2*(2/(pow(Aval - Cval*exp(Bval*m_Lval),2)*pow(m_Rval,2)) - (8*Bval*Cval*exp(Bval*m_Lval)*(I1val - m_Lval))/(pow(Aval - Cval*exp(Bval*m_Lval),3)*pow(m_Rval,2)) +
       (6*pow(Bval,2)*pow(Cval,2)*exp(2*Bval*m_Lval)*pow(I1val - m_Lval,2))/(pow(Aval - Cval*exp(Bval*m_Lval),4)*pow(m_Rval,2)) +
       (2*pow(Bval,2)*Cval*exp(Bval*m_Lval)*pow(I1val - m_Lval,2))/(pow(Aval - Cval*exp(Bval*m_Lval),3)*pow(m_Rval,2)) +
       (6*pow(Bval,2)*pow(Cval,2)*exp(2*Bval*m_Lval)*pow(SqJ2Val,2))/pow(Aval - Cval*exp(Bval*m_Lval),4) +
       (2*pow(Bval,2)*Cval*exp(Bval*m_Lval)*pow(SqJ2Val,2))/pow(Aval - Cval*exp(Bval*m_Lval),3))*
    (-1 + pow(I1val - m_Lval,2)/(pow(Aval - Cval*exp(Bval*m_Lval),2)*pow(m_Rval,2)) + pow(SqJ2Val,2)/pow(Aval - Cval*exp(Bval*m_Lval),2));
    
    Hessian(0,1) += (-4*pow(I1val - m_Lval,2)*((-2*(I1val - m_Lval))/(pow(Aval - Cval*exp(Bval*m_Lval),2)*pow(m_Rval,2)) +
                                             (2*Bval*Cval*exp(Bval*m_Lval)*pow(I1val - m_Lval,2))/(pow(Aval - Cval*exp(Bval*m_Lval),3)*pow(m_Rval,2)) +
                                             (2*Bval*Cval*exp(Bval*m_Lval)*pow(SqJ2Val,2))/pow(Aval - Cval*exp(Bval*m_Lval),3)))/(pow(Aval - Cval*exp(Bval*m_Lval),2)*pow(m_Rval,3)) +
    2*((4*(I1val - m_Lval))/(pow(Aval - Cval*exp(Bval*m_Lval),2)*pow(m_Rval,3)) -
       (4*Bval*Cval*exp(Bval*m_Lval)*pow(I1val - m_Lval,2))/(pow(Aval - Cval*exp(Bval*m_Lval),3)*pow(m_Rval,3)))*
    (-1 + pow(I1val - m_Lval,2)/(pow(Aval - Cval*exp(Bval*m_Lval),2)*pow(m_Rval,2)) + pow(SqJ2Val,2)/pow(Aval - Cval*exp(Bval*m_Lval),2));
    
    Hessian(1,0) += (-4*pow(I1val - m_Lval,2)*((-2*(I1val - m_Lval))/(pow(Aval - Cval*exp(Bval*m_Lval),2)*pow(m_Rval,2)) +
                                             (2*Bval*Cval*exp(Bval*m_Lval)*pow(I1val - m_Lval,2))/(pow(Aval - Cval*exp(Bval*m_Lval),3)*pow(m_Rval,2)) +
                                             (2*Bval*Cval*exp(Bval*m_Lval)*pow(SqJ2Val,2))/pow(Aval - Cval*exp(Bval*m_Lval),3)))/(pow(Aval - Cval*exp(Bval*m_Lval),2)*pow(m_Rval,3)) +
    (8*(I1val - m_Lval)*(-1 + pow(I1val - m_Lval,2)/(pow(Aval - Cval*exp(Bval*m_Lval),2)*pow(m_Rval,2)) + pow(SqJ2Val,2)/pow(Aval - Cval*exp(Bval*m_Lval),2)))/
    (pow(Aval - Cval*exp(Bval*m_Lval),2)*pow(m_Rval,3)) - (8*Bval*Cval*exp(Bval*m_Lval)*pow(I1val - m_Lval,2)*
                                                           (-1 + pow(I1val - m_Lval,2)/(pow(Aval - Cval*exp(Bval*m_Lval),2)*pow(m_Rval,2)) + pow(SqJ2Val,2)/pow(Aval - Cval*exp(Bval*m_Lval),2)))/
    (pow(Aval - Cval*exp(Bval*m_Lval),3)*pow(m_Rval,3));
    
    Hessian(1,1) += (8*pow(I1val - m_Lval,4))/(pow(Aval - Cval*exp(Bval*m_Lval),4)*pow(m_Rval,6)) +
    (12*pow(I1val - m_Lval,2)*(-1 + pow(I1val - m_Lval,2)/(pow(Aval - Cval*exp(Bval*m_Lval),2)*pow(m_Rval,2)) + pow(SqJ2Val,2)/pow(Aval - Cval*exp(Bval*m_Lval),2)))/
    (pow(Aval - Cval*exp(Bval*m_Lval),2)*pow(m_Rval,4));
    
}

void TF2DSAdjust::Residual_R(TPZFMatrix<REAL> &Residual, REAL I1val, REAL SqJ2Val){
    
    Residual.Resize(2,1);
    Residual.Zero();
    
    STATE Aval = m_Sandler.A();
    STATE Bval = m_Sandler.B();
    STATE Cval = m_Sandler.C();
    
    Residual(0,0) += 2*((-2*(I1val - m_Lval))/(pow(Aval - Cval*exp(Bval*m_Lval),2)*pow(m_Rval,2)) +
                        (2*Bval*Cval*exp(Bval*m_Lval)*pow(I1val - m_Lval,2))/(pow(Aval - Cval*exp(Bval*m_Lval),3)*pow(m_Rval,2)) +
                        (2*Bval*Cval*exp(Bval*m_Lval)*pow(SqJ2Val,2))/pow(Aval - Cval*exp(Bval*m_Lval),3))*
    (-1 + pow(I1val - m_Lval,2)/(pow(Aval - Cval*exp(Bval*m_Lval),2)*pow(m_Rval,2)) + pow(SqJ2Val,2)/pow(Aval - Cval*exp(Bval*m_Lval),2));
    
    Residual(1,0) += (-4*pow(I1val - m_Lval,2)*(-1 + pow(I1val - m_Lval,2)/(pow(Aval - Cval*exp(Bval*m_Lval),2)*pow(m_Rval,2)) + pow(SqJ2Val,2)/pow(Aval - Cval*exp(Bval*m_Lval),2)))/
    (pow(Aval - Cval*exp(Bval*m_Lval),2)*pow(m_Rval,3));

}

STATE TF2DSAdjust::AssembleR(TPZFMatrix<REAL> &I1_SqJ2, TPZFMatrix<REAL> &hessian, TPZFMatrix<REAL> &res)
{
    STATE errorR = 0 ;
    I1_SqJ2 = m_I1_SqJ2;
    int64_t n_data = m_I1_SqJ2.Rows();
    
    TPZFMatrix<REAL> hessianMatrix(2,2);
    TPZFMatrix<REAL> residual(2,1);
    
    hessian.Resize(2,2);
    res.Resize(2,1);
    
    hessian.Zero();
    res.Zero();
    
    for(int i=0; i<n_data; i++)
    {
        STATE minF2 = CapFunction(I1_SqJ2(i,0), I1_SqJ2(i,1));
        
        Hessian_R(hessianMatrix, I1_SqJ2(i,0), I1_SqJ2(i,1));
        Residual_R(residual, I1_SqJ2(i,0), I1_SqJ2(i,1));
        
        hessian += hessianMatrix;
        res     += residual;
        errorR += minF2*minF2;
        
    }
    return errorR;
}

void TF2DSAdjust::AdjustR()
{
    TPZFMatrix<REAL> I1_SqJ2 = m_I1_SqJ2;
    TPZFMatrix<REAL> res(2,1,0.);
    TPZFMatrix<REAL> hessian(2,2,0.);
    
    REAL error = AssembleR(I1_SqJ2, hessian, res);
    
    std::cout << "Incoming error " << error << std::endl;
    {
        std::ofstream out("hessian.nb");
        hessian.Print("hessian",out,EMathematicaInput);
        std::cout << std::endl;
    }
    int count = 0;
    REAL normres = Norm(res);
    std::cout << "Incoming normres = " << normres << std::endl;
    REAL tol = normres*1.e-8;
    
    while(normres > tol && count <10)
    {
        hessian.SolveDirect(res,ELU);
        res *= -1.;
        LoadCorrectionR(res);
        
        REAL error = AssembleR(I1_SqJ2, hessian, res);
        
        normres = Norm(res);
        std::cout << "error = " << error << " normres = " << normres << std::endl;
        count++;
    }
    
    std::cout << std::endl;
    std::cout << "found minimum of F2 parameters: ";
    std::cout << "R = (" << Rval() << ")"<< std::endl;
    
    STATE A  = m_Sandler.A();
    STATE B  = m_Sandler.B();
    STATE C  = m_Sandler.C();
    REAL r   = Rval();
    REAL l_0 = Lval();
    REAL X_0val = ComputeXval(A, B, C, r, l_0);
    
    /// The following X0 value can be used in AdjustX0 to compute modified X0
    std::cout << "Initial value of X0 = (" << X_0val << ")"<< std::endl;
    std::cout << std::endl;
    std::cout.flush();
}

STATE TF2DSAdjust::ComputeXval(STATE &a, STATE &b, STATE &c, STATE &r, STATE &l){
    
    STATE x = l-r*(a-c*exp(b*l));
    return x;
}

void TF2DSAdjust::LoadCorrectionR(TPZFMatrix<REAL> &delx)
{
    REAL l_val = Lval()+delx(0,0);
    SetLval(l_val);
    REAL r_val = Rval()+delx(1,0);
    SetRval(r_val);
    
}

/// Second step: finding the D and W
/// Reading the input data in order to compute the plastic strain
void TF2DSAdjust::ComputedElasticStrainLX(TPZVec<TTestSection> &active, std::vector<REAL> &epsEV_data)
{
//    TPZElasticResponse ER;
    REAL lamda = m_ER.Lambda();
    REAL mu    = m_ER.G();

    TPZFMatrix<REAL> deform,stress;
    for(int64_t i=0; i<active.size(); i++)
    {
        active[i].GetData(deform,stress);
        int64_t ndata = deform.Rows();
        TPZManVector<REAL,2> strain_measure(2);

        for(int64_t d = 0; d<ndata; d++)
        {
            strain_measure[0] = deform(d,0);
            strain_measure[1] = deform(d,1);
            
            REAL epsEV = (2*strain_measure[0]+strain_measure[1])/(2*mu+3*lamda);
            
            epsEV_data.push_back(epsEV);
   
        }
    }
}

void TF2DSAdjust::ComputedPlasticStrainLX(TPZVec<TTestSection> &active, std::vector<REAL> &epsPV_data)
{
    std::vector<REAL> epsElV_data;
    std::vector<REAL> epsTV_data;
    TPZFMatrix<REAL> deform,stress;
    ComputedElasticStrainLX(active, epsElV_data);
    
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

/// Reading the input data: invariant stress: I1 and SqJ2
void TF2DSAdjust::ComputedInvariantStressLX(TPZVec<TTestSection> &active, std::vector<REAL> &I1data, std::vector<REAL> & SqJ2data)
{
    TPZFMatrix<REAL> invstress;
    for(int64_t i=0; i<active.size(); i++)
    {
        active[i].GetInvStressData(invstress);
        int64_t ndata = invstress.Rows();
        
        for(int64_t j = 0; j<ndata; j++)
        {
            REAL I1   = invstress(j,0);
            REAL SqJ2 = invstress(j,0);
            
            I1data.push_back(I1);
            SqJ2data.push_back(SqJ2);
        }
    }
}

/// Methods to compute the L value in order to compute the X value from input data
void TF2DSAdjust::Hessian_LX(TPZFMatrix<REAL> & Hessian, REAL I1val, REAL SqJ2Val){
    
    Hessian.Resize(1,1);
    Hessian.Zero();
    
    STATE Aval = m_Sandler.A();
    STATE Bval = m_Sandler.B();
    STATE Cval = m_Sandler.C();
    REAL  rval = Rval();
    
    Hessian(0,0) += 2*pow((-2*(I1val - m_Lval))/(pow(Aval - Cval*exp(Bval*m_Lval),2)*pow(rval,2)) +
                            (2*Bval*Cval*exp(Bval*m_Lval)*pow(I1val - m_Lval,2))/
                            (pow(Aval - Cval*exp(Bval*m_Lval),3)*pow(rval,2)) +
                            (2*Bval*Cval*exp(Bval*m_Lval)*pow(SqJ2Val,2))/
                            pow(Aval - Cval*exp(Bval*m_Lval),3),2) +
    2*(2/(pow(Aval - Cval*exp(Bval*m_Lval),2)*pow(rval,2)) -
       (8*Bval*Cval*exp(Bval*m_Lval)*(I1val - m_Lval))/
       (pow(Aval - Cval*exp(Bval*m_Lval),3)*pow(rval,2)) +
       (6*pow(Bval,2)*pow(Cval,2)*exp(2*Bval*m_Lval)*pow(I1val - m_Lval,2))/
       (pow(Aval - Cval*exp(Bval*m_Lval),4)*pow(rval,2)) +
       (2*pow(Bval,2)*Cval*exp(Bval*m_Lval)*pow(I1val - m_Lval,2))/
       (pow(Aval - Cval*exp(Bval*m_Lval),3)*pow(rval,2)) +
       (6*pow(Bval,2)*pow(Cval,2)*exp(2*Bval*m_Lval)*pow(SqJ2Val,2))/
       pow(Aval - Cval*exp(Bval*m_Lval),4) +
       (2*pow(Bval,2)*Cval*exp(Bval*m_Lval)*pow(SqJ2Val,2))/
       pow(Aval - Cval*exp(Bval*m_Lval),3))*
    (-1 + pow(I1val - m_Lval,2)/(pow(Aval - Cval*exp(Bval*m_Lval),2)*pow(rval,2)) +
     pow(SqJ2Val,2)/pow(Aval - Cval*exp(Bval*m_Lval),2));
    
}

void TF2DSAdjust::Residual_LX(TPZFMatrix<REAL> &Residual, REAL I1val, REAL SqJ2Val){
    
    Residual.Resize(1,1);
    Residual.Zero();
    
    STATE Aval = m_Sandler.A();
    STATE Bval = m_Sandler.B();
    STATE Cval = m_Sandler.C();
    REAL  rval = Rval();
    
    Residual(0,0) += 2*((-2*(I1val - m_Lval))/(pow(Aval - Cval*exp(Bval*m_Lval),2)*pow(rval,2)) +
                        (2*Bval*Cval*exp(Bval*m_Lval)*pow(I1val - m_Lval,2))/
                        (pow(Aval - Cval*exp(Bval*m_Lval),3)*pow(rval,2)) +
                        (2*Bval*Cval*exp(Bval*m_Lval)*pow(SqJ2Val,2))/
                        pow(Aval - Cval*exp(Bval*m_Lval),3))*
    (-1 + pow(I1val - m_Lval,2)/(pow(Aval - Cval*exp(Bval*m_Lval),2)*pow(rval,2)) +
     pow(SqJ2Val,2)/pow(Aval - Cval*exp(Bval*m_Lval),2));
}


void TF2DSAdjust::ComputeXval(TPZVec<TTestSection> &active, std::vector<REAL> &I1data, std::vector<REAL> & SqJ2data, std::vector<REAL> &X_data)
{
    ComputedInvariantStressLX(active, I1data, SqJ2data);
    int n_data = I1data.size();
    TPZFMatrix<REAL> hessianMatrix(1,1);
    TPZFMatrix<REAL> residual(1,1);
    TPZFMatrix<REAL> hessian,res;
    hessian.Resize(1,1);
    res.Resize(1,1);
    hessian.Zero();
    res.Zero();
    
    STATE Aval = m_Sandler.A();
    STATE Bval = m_Sandler.B();
    STATE Cval = m_Sandler.C();
    REAL  rval = Rval();
    
    for(int i=0; i<n_data; i++)
    {
        Hessian_LX(hessianMatrix, I1data[i], SqJ2data[i]);
        Residual_LX(residual, I1data[i], SqJ2data[i]);
        
        hessian = hessianMatrix;
        res     = residual;
        
        int count = 0;
        REAL normres = Norm(res);
        REAL tol = normres*1.e-8;
        
        while(normres > tol && count <10)
        {
            hessian.SolveDirect(res,ELU);
            res *= -1.;
            LoadCorrectionLX(res);
            count++;
        }
        
        REAL lval = Lval();
        
        REAL xval = ComputeXval(Aval, Bval, Cval, rval, lval);
        
        X_data.push_back(xval);
    }
}

/// Method to modify L value in order to modify the X value from input data
void TF2DSAdjust::LoadCorrectionLX(TPZFMatrix<REAL> &delx)
{
    REAL l_val = Lval()+delx(0,0);
    SetLval(l_val);
}


void TF2DSAdjust::PopulateDW()
{
    // Set elastic data
//    TPZElasticResponse ER;
    REAL E  = 2214.77;
    REAL nu = 0.25141;
    SetYoungPoisson(E, nu);
    
    REAL lmbda = (E*nu)/((1+nu)*(1-2*nu));
    REAL mu    = E/(2*(1+nu));
    SetLameDataElastic(lmbda, mu);
    
    /// Render the data for failure function of DiMaggio-Sandler model
    REAL a_val = 15.091;
    REAL b_val = 0.0284035;
    REAL c_val = 15.9158;
    REAL r_val = 2.55954;
    
    m_Sandler.SetA(a_val);
    m_Sandler.SetB(b_val);
    m_Sandler.SetC(c_val);
    SetRval(r_val);
    
    /// Initial guess for L parameter
    REAL l_val = -60.9137409;
    SetLval(l_val);

}

/// Method to combine the X values vs plastic strain deformation
void TF2DSAdjust::XvalvsEpsPlasticdata(TPZVec<TTestSection> &active, std::vector<REAL> &epsPV_data, std::vector<REAL> &I1data, std::vector<REAL> & SqJ2data, std::vector<REAL> &X_data, TPZFMatrix<REAL> &X_epsP){
    
    ComputedPlasticStrainLX(active, epsPV_data);
    ComputeXval(active, I1data,SqJ2data, X_data);
    
    
    int64_t n_data = I1data.size();
    X_epsP.Resize(n_data, 2);

    for(int i=0; i<n_data; i++)
    {
        REAL xval = X_data[i];
        REAL epsplast = epsPV_data[i];
        
        X_epsP(i,0) = xval;
        X_epsP(i,1) = epsplast;
        
    }
        
}


REAL TF2DSAdjust::ComputeDval_initial(){
    
    TPZVec<TTestSection> active;
    std::vector<REAL> epsPV_data;
    TPZFMatrix<REAL> I1_SqJ2;
    std::vector<REAL> I1data, SqJ2data;
    std::vector<REAL> X_data;
    TPZFMatrix<REAL> X_epsP;
    
    XvalvsEpsPlasticdata(active, epsPV_data, I1data, SqJ2data, X_data, X_epsP);
    
    int64_t n_data = X_epsP.Rows();
    REAL sum_D_val = 0.0;
    
    for(int64_t i = 1; i < n_data-1; i++)
    {
        REAL dval1 = abs((X_epsP(i,1)-X_epsP(i-1,1))/(X_epsP(i,0)-X_epsP(i-1,0)));
        REAL dval2 = abs((X_epsP(i+1,1)-X_epsP(i,1))/(X_epsP(i+1,0)-X_epsP(i,0)));
        REAL num   = abs(log(dval1/dval2));
        REAL denum = abs(X_epsP(i,0)-X_epsP(i+1,0));

        REAL dval  = abs(num/denum);
        sum_D_val += dval;
    }
    
    return sum_D_val/(n_data-1);;
}


REAL TF2DSAdjust::ComputeWval_initial(){
    
    TPZVec<TTestSection> active;
    std::vector<REAL> epsPV_data;
    std::vector<REAL> I1data, SqJ2data;
    std::vector<REAL> X_data;
    TPZFMatrix<REAL> X_epsP;
    
    XvalvsEpsPlasticdata(active, epsPV_data, I1data, SqJ2data, X_data, X_epsP);
    
    REAL Dval = ComputeDval_initial();
    
    int64_t n_data = X_epsP.Rows();
    REAL sum_W_val = 0.0;
    
    for(int64_t i = 0; i < n_data; i++)
    {
        REAL num   = X_epsP(i,0);
        REAL denum = exp(Dval*X_epsP(i,1))-1;
        
        REAL dval  = abs(num/denum);
        sum_W_val += dval;
    }
    
    return sum_W_val/(n_data);;
}


void TF2DSAdjust::Hessian_DW(TPZFMatrix<REAL> &Hessian, REAL &X, REAL &depsPv){
    
    Hessian.Resize(2,2);
    Hessian.Zero();
    
    Hessian(0,0) += 2*pow(X,2);
    Hessian(0,1) += 2*X;
    Hessian(1,0) += 2*X;
    Hessian(1,1) += 2.;
    
}


void TF2DSAdjust::Residual_DW(TPZFMatrix<REAL> &Residual, REAL &X, REAL &depsPv){
    
    Residual.Resize(2,1);
    Residual.Zero();
    
    Residual(0,0) += -2*X*(-m_Wval - m_Dval*X + log(depsPv));
    Residual(1,0) += -2*(-m_Wval - m_Dval*X + log(depsPv));
    
}


STATE TF2DSAdjust::AssembleDW(TPZVec<TTestSection> &active, TPZFMatrix<REAL> &X_epsP, TPZFMatrix<REAL> &hessian, TPZFMatrix<REAL> &res)
{
    std::vector<REAL> epsPV_data;
    std::vector<REAL> I1data, SqJ2data;
    std::vector<REAL> X_data;
    
    XvalvsEpsPlasticdata(active, epsPV_data, I1data, SqJ2data, X_data, X_epsP);
    
    int64_t n_data = X_epsP.Rows();
    STATE errorDW = 0 ;
    
    TPZFMatrix<REAL> hessianMatrix(2,2);
    TPZFMatrix<REAL> residual(2,1);
    
    hessian.Resize(2,2);
    res.Resize(2,1);
    
    hessian.Zero();
    res.Zero();
    
    for(int i=1; i<n_data; i++)
    {
        REAL xval   = X_epsP(i,0);
        REAL depsPv = (X_epsP(i,1)-X_epsP(i-1,1))/(X_epsP(i,0)-X_epsP(i-1,0));
        
        STATE minDW = DWcostFunction(xval, depsPv);
        
        Hessian_DW(hessianMatrix, xval, depsPv);
        Residual_DW(residual, xval, depsPv);
        
        hessian += hessianMatrix;
        res     += residual;
        errorDW += minDW*minDW;
        
    }
    return errorDW;
}


void TF2DSAdjust::AdjustDW(TPZVec<TTestSection> &active)
{
    std::vector<REAL> epsPV_data;
    std::vector<REAL> I1data, SqJ2data;
    std::vector<REAL> X_data;
    TPZFMatrix<REAL>  X_epsP;
    
    XvalvsEpsPlasticdata(active, epsPV_data, I1data, SqJ2data, X_data, X_epsP);
    
    TPZFMatrix<REAL> res(2,1,0.);
    TPZFMatrix<REAL> hessian(2,2,0.);
    
    REAL error = AssembleDW(active, X_epsP, hessian, res);
    
    std::cout << "Incoming error " << error << std::endl;
    {
        std::ofstream out("hessian.nb");
        hessian.Print("hessian",out,EMathematicaInput);
        std::cout << std::endl;
    }
    int count = 0;
    REAL normres = Norm(res);
    std::cout << "Incoming normres = " << normres << std::endl;
    REAL tol = normres*1.e-8;
    
    while(normres > tol && count <10)
    {
        hessian.SolveDirect(res,ELU);
        res *= -1.;
        LoadCorrectionDW(res);
        
        REAL error = AssembleDW(active, X_epsP, hessian, res);
        
        normres = Norm(res);
        std::cout << "error = " << error << " normres = " << normres << std::endl;
        count++;
    }
    
    std::cout << std::endl;
    std::cout << "found minimum of F2 parameters: ";
    std::cout << "D = (" << Dval() << ")"<< std::endl;
    std::cout << std::endl;
    std::cout << "W = (" << Wval() << ")"<< std::endl;
    std::cout << std::endl;
    std::cout.flush();
}


STATE TF2DSAdjust::DWcostFunction(REAL &X, REAL &depsPv){
    
    REAL Dvalue = Dval();
    REAL Wvalue = Wval();
    
    STATE cost = -Wvalue - Dvalue*X + log(depsPv);
    return cost;
}

void TF2DSAdjust::LoadCorrectionDW(TPZFMatrix<REAL> &delx)
{
    REAL d_val = Dval()+delx(0,0);
    SetDval(d_val);
    REAL w_val = Wval()+delx(1,0);
    SetWval(w_val);
    
}

void TF2DSAdjust::AdjustX0()
{
    REAL a_val = 15.091;
    REAL b_val = 0.0284035;
    REAL c_val = 15.9158;
    REAL r_val = 2.55954;
    
    /// The quantity of L is obtained based on the AdjustR
    REAL l_val = -56.9677;
    
    REAL X0_1 = ComputeXval(a_val, b_val, c_val, r_val, l_val);
    
    TPZVec<TTestSection> active;
    TPZFMatrix<REAL>  I1_SqJ2;
    std::vector<REAL> X_data;
    std::vector<REAL> I1data, SqJ2data;
    REAL max_X;
    REAL deltaX;
    REAL X0;
    
    ComputeXval(active, I1data, SqJ2data, X_data);
    
    FindMaxXvalue(X_data, max_X);
    
    if (max_X==X0_1) {
        
        deltaX = 0.;
        
    } else if (max_X>X0_1){
        
        deltaX = X0_1-max_X;
        
    } else{
        
        deltaX = max_X-X0_1;
    };
    
    ComputedInvariantStressLX(active, I1data, SqJ2data);
   
    X0 = deltaX - 2*I1data[0];
    
    std::cout << std::endl;
    std::cout << "found minimum of F2 parameters: ";
    std::cout << "X0 = (" << X0 << ")"<< std::endl;
    std::cout << std::endl;
    std::cout.flush();
}

void TF2DSAdjust::FindMaxXvalue(std::vector<REAL> &X_data, REAL &max_X){
    
    int highNum = 0;
    int m;
    int64_t ndata = X_data.size();
    
    for (m = 0 ; m < ndata ; m++);
    {
        if (X_data[m] > highNum)
            highNum = X_data[m];
    }
    max_X = highNum;
}

