//
//  TF2DSAdjust_R.cpp
//  plastic_adjust
//
//  Created by Manouchehr Sanei on 6/10/19.
//
//

#include "TF2DSAdjust_R.h"
#include <math.h>
#include "StubFunctions.h"
#include "TF1DSAdjust.h"
#include "nlopt.hpp"



TF2DSAdjust_R::TF2DSAdjust_R(){
    
    m_Lval= 0.;
    m_Rval= 0.;
    
}

TF2DSAdjust_R::TF2DSAdjust_R(const TF2DSAdjust_R &other) : m_Sandler(other.m_Sandler), m_I1_SqJ2(other.m_I1_SqJ2)
{
    m_Lval = other.m_Lval;
    m_Rval = other.m_Rval;
    
}

const TF2DSAdjust_R & TF2DSAdjust_R::operator=(const TF2DSAdjust_R &other)
{
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_Sandler = other.m_Sandler;
    m_I1_SqJ2 = other.m_I1_SqJ2;
    
    m_Lval = other.m_Lval;
    m_Rval = other.m_Rval;
    
    return *this;
}

TF2DSAdjust_R::~TF2DSAdjust_R(){
    
}

void TF2DSAdjust_R::PopulateR()
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


STATE TF2DSAdjust_R::errorfunctionF2_R(const std::vector<STATE> &input)
{
    STATE L   = input[0];
    STATE R   = input[1];
    int64_t n_data = m_I1_SqJ2.Rows();
    
    STATE A   = m_Sandler.A();
    STATE B   = m_Sandler.B();
    STATE C   = m_Sandler.C();
    STATE gamma = 1.0;

    STATE errorF2=0;

    for(int i=0; i<n_data; i++)
    {
        REAL I1val    = m_I1_SqJ2(i,0);
        REAL SqJ2Val  = m_I1_SqJ2(i,1);
        
        STATE term1 = (I1val-L) / (R * (A-(C * exp(B * L))));
        STATE term2 = (gamma * SqJ2Val) / (A-(C * exp(B * L)));
        STATE minF2 = (term1 * term1 + term2 * term2 - 1);

        errorF2 += minF2*minF2;

    }

    return errorF2;
}



void TF2DSAdjust_R::gradientfunctionF2_R(const std::vector<STATE> &input, std::vector<double> &grad)
{
    STATE L   = input[0];
    STATE R   = input[1];
    int64_t n_data = m_I1_SqJ2.Rows();
    
    STATE Aval = m_Sandler.A();
    STATE Bval = m_Sandler.B();
    STATE Cval = m_Sandler.C();
    
    STATE grad0 = 0.;
    STATE grad1 = 0.;
    
    for(int i=0; i<n_data; i++)
    {
        REAL I1val    = m_I1_SqJ2(i,0);
        REAL SqJ2Val  = m_I1_SqJ2(i,1);
        
        
        grad0 += 2*((-2*(I1val - L))/(pow(Aval - Cval*exp(Bval*L),2)*pow(R,2)) +
                    (2*Bval*Cval*exp(Bval*L)*pow(I1val - L,2))/(pow(Aval - Cval*exp(Bval*L),3)*pow(R,2)) +
                    (2*Bval*Cval*exp(Bval*L)*pow(SqJ2Val,2))/pow(Aval - Cval*exp(Bval*L),3))*
        (-1 + pow(I1val - L,2)/(pow(Aval - Cval*exp(Bval*L),2)*pow(R,2)) + pow(SqJ2Val,2)/pow(Aval - Cval*exp(Bval*L),2));
        
        grad1 += (-4*pow(I1val - L,2)*(-1 + pow(I1val - L,2)/(pow(Aval - Cval*exp(Bval*L),2)*pow(R,2)) + pow(SqJ2Val,2)/pow(Aval - Cval*exp(Bval*L),2)))/
        (pow(Aval - Cval*exp(Bval*L),2)*pow(R,3));
    }
    
    grad[0] = grad0;
    grad[1] = grad1;
    
}


double myvfunctionR(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    TF2DSAdjust_R *loc = (TF2DSAdjust_R *) my_func_data;
    
    if (!grad.empty()) {
        
        loc->gradientfunctionF2_R(x, grad);
    }
    
    double err = loc->errorfunctionF2_R(x);
    for(int i=0; i<x.size(); i++) std::cout << "x[" << i << "]= " << x[i] << " ";
    std::cout << "error = " << err << std::endl;
    return err;
}

/// NLopt functions are: LD_SLSQP, LD_VAR1, LD_VAR2

void TF2DSAdjust_R::AdjustR()
{
    nlopt::opt opt(nlopt::LD_VAR2, 2);
    opt.set_min_objective(myvfunctionR, this);
    opt.set_xtol_rel(1e-6);
    std::vector<double> x(2,0.);
    std::vector<double> lb(2);
    lb[0] = 1.5*Lval();
    lb[1] = 0.;
    opt.set_lower_bounds(lb);
    std::vector<double> ub(2);
    ub[0] = 0.;
    ub[1] = 3.0*Rval();
    opt.set_upper_bounds(ub);
    
    /// Initialize the L and R value
    x[0] = Lval();
    x[1] = Rval();
    STATE Aval = m_Sandler.A();
    STATE Bval = m_Sandler.B();
    STATE Cval = m_Sandler.C();
    double minf;
    
    try{
        nlopt::result result = opt.optimize(x, minf);
        
        if (result < 0) {
            std::cerr << "nlopt failed: result = " << result << std::endl;
        } else {
            
            std::cout << std::endl;
            std::cout << "found minimum of F2 parameters: ";
            std::cout << std::endl;
            for(int i=1; i<x.size(); i++)
                std::cout << "R = " << x[i] << ", " << "X0 = " << ComputeXval(Aval,Bval,Cval,x[i],x[i-1]) << std::endl;
            
            std::cout << std::endl;
            std::cout<< std::setprecision(10) << "min_val = " << minf << std::endl;
        }
    }
    catch(std::exception &e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }
    
}



void TF2DSAdjust_R::Hessian_R(TPZFMatrix<REAL> & Hessian, REAL I1val, REAL SqJ2Val){
    
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

void TF2DSAdjust_R::Residual_R(TPZFMatrix<REAL> &Residual, REAL I1val, REAL SqJ2Val){
    
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

STATE TF2DSAdjust_R::AssembleR(TPZFMatrix<REAL> &I1_SqJ2, TPZFMatrix<REAL> &hessian, TPZFMatrix<REAL> &res)
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


void TF2DSAdjust_R::AdjustR2()
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
    std::cout << std::endl;
    std::cout << "R = " << Rval() << std::endl;
    
    STATE A  = m_Sandler.A();
    STATE B  = m_Sandler.B();
    STATE C  = m_Sandler.C();
    REAL r   = Rval();
    REAL l_0 = Lval();
    REAL X_0val = ComputeXval(A, B, C, r, l_0);
    
    /// Initial value of X0
    std::cout << "X0 = " << X_0val << std::endl;
    std::cout << std::endl;
    std::cout.flush();
}

STATE TF2DSAdjust_R::ComputeXval(STATE &a, STATE &b, STATE &c, STATE &r, STATE &l){
    
    STATE x = l-r*(a-c*exp(b*l));
    return x;
}


void TF2DSAdjust_R::LoadCorrectionR(TPZFMatrix<REAL> &delx)
{
    REAL l_val = Lval()+delx(0,0);
    SetLval(l_val);
    REAL r_val = Rval()+delx(1,0);
    SetRval(r_val);
    
}


STATE TF2DSAdjust_R::CapFunction(STATE &I1, STATE &SqJ2){
    
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
