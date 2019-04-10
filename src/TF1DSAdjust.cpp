//
//  TF1DSAdjust.cpp
//  plastic_adjust
//
//  Created by Manouchehr Sanei on 4/7/19.
//
//

#include "TF1DSAdjust.h"
#include "nlopt.hpp"


TF1DSAdjust::TF1DSAdjust(){
    
    m_Aval = 0.;
    m_Bval = 0.;
    m_Cval = 0.;

}

TF1DSAdjust::TF1DSAdjust(const TF1DSAdjust &other) : m_Sandler(other.m_Sandler), m_I1_SqJ2(other.m_I1_SqJ2)
{
    m_Aval = other.m_Aval;
    m_Bval = other.m_Bval;
    m_Cval = other.m_Cval;
}


const TF1DSAdjust & TF1DSAdjust::operator=(const TF1DSAdjust &other)
{
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_Sandler = other.m_Sandler;
    m_I1_SqJ2 = other.m_I1_SqJ2;
    
    m_Aval = other.m_Aval;
    m_Bval = other.m_Bval;
    m_Cval = other.m_Cval;
    
    return *this;
}

TF1DSAdjust::~TF1DSAdjust(){
    
}


void TF1DSAdjust::Populate()
{
    m_Sandler.MCormicRanchSand(m_Sandler);
    
    /// Render data; temporary
    const int n_I1 = 6;
    REAL I1_data[n_I1] = {-5.81226, -28.7936, -63.9478, -99.3083, -125.169, -129.247};
    REAL F1_data[n_I1] = {1.62448, 7.96679, 12.6667, 14.0312, 14.6975, 14.6444};
    
    m_I1_SqJ2.Redim(n_I1,2);
    
    for(int i=0; i<n_I1; i++)
    {
        REAL I1 = I1_data[i];
        REAL F1 = F1_data[i];
        m_I1_SqJ2(i,0) = I1;
        m_I1_SqJ2(i,1) = F1;
    }
    
    
    REAL b_val = computeB_F1();
    SetBval(b_val);
    REAL c_val = computeC_F1();
    SetCval(c_val);
    REAL a_val = computeA_F1();
    SetAval(a_val);
    
    m_Sandler.SetB(b_val);
    m_Sandler.SetC(c_val);
    m_Sandler.SetA(a_val);
    
    m_I1_SqJ2.Print("TestData",std::cout);
    std::cout << "Objective parameters A = " << m_Sandler.A() << " B = " << m_Sandler.B() << " C = " << m_Sandler.C() << std::endl;
}


REAL TF1DSAdjust::computeB_F1(){
    
    TPZFMatrix<REAL> I1_SqJ2 = m_I1_SqJ2;
    int64_t n_data = m_I1_SqJ2.Rows();
    REAL sum_B_val = 0.0;
    
    for(int64_t i = 1; i < n_data - 1; i++)
    {
        REAL SqJ2primeI1_n = abs((I1_SqJ2(i+1,1)-I1_SqJ2(i,1))/(I1_SqJ2(i+1,0)-I1_SqJ2(i,0)));
        REAL SqJ2primeI1   = abs((I1_SqJ2(i,1)-I1_SqJ2(i-1,1))/(I1_SqJ2(i,0)-I1_SqJ2(i-1,0)));
        REAL I1_val_denom  = abs(I1_SqJ2(i+1,0)-I1_SqJ2(i-1,0));
        REAL bvalue        = abs((log(SqJ2primeI1_n/SqJ2primeI1))/(I1_val_denom));
        sum_B_val += bvalue;
    }
    
    return sum_B_val/(n_data-2);;
}


REAL TF1DSAdjust::computeC_F1(){
    
    TPZFMatrix<REAL> I1_SqJ2 = m_I1_SqJ2;
    int64_t n_data = m_I1_SqJ2.Rows();
    REAL bval = Bval();
    REAL sum_C_val = 0.0;
    
    for(int64_t i = 1; i < n_data; i++)
    {
        REAL SqJ2_num  = (I1_SqJ2(i,1)-I1_SqJ2(i-1,1));
        REAL I1_denom  = exp(bval*I1_SqJ2(i,0))-exp(bval*I1_SqJ2(i-1,0));
        REAL cvalue    = abs((SqJ2_num)/(I1_denom));
        sum_C_val += cvalue;
    }
    
    return sum_C_val/(n_data-1);
}


REAL TF1DSAdjust::computeA_F1(){
    
    TPZFMatrix<REAL> &I1_SqJ2  = m_I1_SqJ2;
    int64_t n_data = m_I1_SqJ2.Rows();
    REAL bval = Bval();
    REAL cval = Cval();
    REAL sum_A_val = 0.0;
    
    for(int64_t i = 0; i < n_data; i++)
    {
        REAL avalue  = I1_SqJ2(i,1)+abs(cval)*exp(abs(bval)*I1_SqJ2(i,0));
        sum_A_val   += avalue;
    }
    
    return sum_A_val/n_data;
}


STATE TF1DSAdjust::errorfunction(const std::vector<STATE> &input)
{
    STATE B = input[0];
    STATE C = input[1];
    STATE A = input[2];
    STATE D = 0.;
    STATE W = 0.;
    STATE K(0.), G(0.), R(0.), phi(0.), N(0.), Psi(0.), kappa_0(0.);
    TPZSandlerExtended sandler(A,B,C,D,K,G,W,R,phi,N,Psi,kappa_0);
    
    STATE errorF1=0;
    int64_t n_data = m_I1_SqJ2.Rows();
    
    for(int i=0; i<n_data; i++)
    {
         STATE minF1 = sandler.F(this->m_I1_SqJ2(i,0))-this->m_I1_SqJ2(i,1);
        
        errorF1 += minF1*minF1;

    }

    return errorF1;
}


typedef struct {
    double a, b;
} my_constraint_data;

double myvfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    if (!grad.empty()) {
        DebugStop();
    }
    TF1DSAdjust *loc = (TF1DSAdjust *) my_func_data;
    double err = loc->errorfunction(x);
    for(int i=0; i<x.size(); i++) std::cout << "x[" << i << "]= " << x[i] << " ";
    std::cout << "error " << err << std::endl;
    return err;
}


void TF1DSAdjust::Adjust()
{
    nlopt::opt opt(nlopt::LN_NEWUOA_BOUND, 3);
    opt.set_min_objective(myvfunc, this);
    opt.set_xtol_rel(1e-6);
    std::vector<double> x(3,0.);
    std::vector<double> lb(3);
    lb[0] = 0.; lb[1] = 0.; lb[2] = 0.;
    opt.set_lower_bounds(lb);
    std::vector<double> ub(3);
    ub[0] = 2.*Bval();
    ub[1] = 2.*Cval();
    ub[2] = 2.*Aval();
    
    opt.set_upper_bounds(ub);
    
    /// Initialize the B and C value
    x[0] = Bval();
    x[1] = Cval();
    x[2] = Aval();
    
    double minf;
    
    try{
        nlopt::result result = opt.optimize(x, minf);
        
        if (result < 0) {
            std::cerr << "nlopt failed: result = " << result << std::endl;
        } else {
        
        std::cout << "found minimum at f(";
        for(int i=0; i<x.size(); i++) std::cout << x[i] << " ";
        std::cout << std::setprecision(10) << " min_val " << minf << std::endl;
        }
    }
    catch(std::exception &e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }
    
}


void TF1DSAdjust::Hessian(TPZFMatrix<REAL> & Hessian, REAL I1val, REAL SqJ2Val){
    
    Hessian.Zero();
    
    Hessian(0,0) +=  2.;
    Hessian(0,1) += -2*m_Cval*exp(m_Bval*I1val)*I1val;
    Hessian(0,2) += -2*exp(m_Bval*I1val);
    
    Hessian(1,0) += -2*m_Cval*exp(m_Bval*I1val)*I1val;
    Hessian(1,1) +=  2*pow(m_Cval,2)*exp(2*m_Bval*I1val)*pow(I1val,2)-2*m_Cval*exp(m_Bval*I1val)*pow(I1val,2)*(m_Aval-m_Cval*exp(m_Bval*I1val)-SqJ2Val);
    Hessian(1,2) +=  2*m_Cval*exp(2*m_Bval*I1val)*I1val-2*exp(m_Bval*I1val)*I1val*(m_Aval-m_Cval*exp(m_Bval*I1val)-SqJ2Val);
    
    Hessian(2,0) += -2*exp(m_Bval*I1val);
    Hessian(2,1) +=  2*m_Cval*exp(2*m_Bval*I1val)*I1val-2*exp(m_Bval*I1val)*I1val*(m_Aval-m_Cval*exp(m_Bval*I1val)-SqJ2Val);
    Hessian(2,2) +=  2*exp(2*m_Bval*I1val);

}

void TF1DSAdjust::Residual(TPZFMatrix<REAL> &Residual, REAL I1val, REAL SqJ2Val){
    
    Residual.Zero();
    
    Residual(0,0) +=  2*(m_Aval-m_Cval*exp(m_Bval*I1val)-SqJ2Val);
    Residual(1,0) += -2*m_Cval*exp(m_Bval*I1val)*I1val*(m_Aval-m_Cval*exp(m_Bval*I1val)-SqJ2Val);
    Residual(2,0) += -2*exp(m_Bval*I1val)*(m_Aval-m_Cval*exp(m_Bval*I1val)-SqJ2Val);

}


void TF1DSAdjust::Assemble(TPZFMatrix<REAL> &I1_SqJ2, TPZFMatrix<REAL> &hessian, TPZFMatrix<REAL> &res)
    {
    I1_SqJ2 = m_I1_SqJ2;
    int64_t n_data = m_I1_SqJ2.Rows();

    hessian.Zero();
    res.Zero();

    REAL I1val,SqJ2Val;

        for(int64_t j = 0; j<n_data; j++)
        {
            I1val = I1_SqJ2(j,0);
            SqJ2Val = I1_SqJ2(j,1);
            
            Residual(res, I1val, SqJ2Val);
            Hessian(hessian, I1val, SqJ2Val);
            
        }
}


void TF1DSAdjust::Adjust2()
{
    TPZFMatrix<REAL> I1_SqJ2 = m_I1_SqJ2;
    int64_t n_data = m_I1_SqJ2.Rows();
    TPZFMatrix<REAL> res(n_data,1,0.);
    TPZFMatrix<REAL> hessian(n_data,n_data,0.);
    
    std::vector<double> x(3,0.);
    REAL error = errorfunction(x);
    Assemble(I1_SqJ2, hessian, res);
    
    std::cout << "Incoming error " << error << std::endl;
    {
        std::ofstream out("tangent.nb");
        hessian.Print("tangent",out,EMathematicaInput);
        std::cout << "NOTHING COMING OUT?\n";
    }
    int count = 0;
    REAL normres = Norm(res);
    std::cout << "Incoming normres = " << normres << std::endl;
    REAL tol = normres*1.e-8;
    
    x[0] = Bval();
    x[1] = Cval();
    x[2] = Aval();
    
    while(normres > tol && count <10)
    {
        hessian.SolveDirect(res,ECholesky);
        res *= -1.;
        REAL error = errorfunction(x);
        Assemble(I1_SqJ2, hessian, res);
        normres = Norm(res);
        std::cout << "error = " << error << " normres = " << normres << std::endl;
        count++;
    }
    std::cout << "found minimum at f(";
    for(int i=0; i<x.size(); i++) std::cout << x[i] << " ";
    std::cout << std::setprecision(10) << " min_val " << error << std::endl;
    std::cout.flush();
}

