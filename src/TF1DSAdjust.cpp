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
    
    m_A_val = 0.;
    m_B_val = 0.;
    m_C_val = 0.;

}

TF1DSAdjust::TF1DSAdjust(const TF1DSAdjust &other) : m_Sandler(other.m_Sandler), m_I1_SqJ2(other.m_I1_SqJ2)
{
    m_A_val = other.m_A_val;
    m_B_val = other.m_B_val;
    m_C_val = other.m_C_val;
}


const TF1DSAdjust & TF1DSAdjust::operator=(const TF1DSAdjust &other)
{
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_Sandler = other.m_Sandler;
    m_I1_SqJ2 = other.m_I1_SqJ2;
    
    m_A_val = other.m_A_val;
    m_B_val = other.m_B_val;
    m_C_val = other.m_C_val;
    
    return *this;
}

TF1DSAdjust::~TF1DSAdjust(){
    
}


void TF1DSAdjust::B_F1_function(REAL &bMean, TPZFMatrix<REAL> &I1_SqJ2){
    
    I1_SqJ2        = m_I1_SqJ2;
    int64_t n_data = m_I1_SqJ2.Rows();
    REAL sum_B_val = 0.0;
    
    for(int64_t i = 1; i < n_data - 1; i++)
    {
        REAL SqJ2primeI1_n = abs((I1_SqJ2(i+1,1)-I1_SqJ2(i,1))/(I1_SqJ2(i+1,0)-I1_SqJ2(i,0)));
        REAL SqJ2primeI1   = abs((I1_SqJ2(i,1)-I1_SqJ2(i-1,1))/(I1_SqJ2(i,0)-I1_SqJ2(i-1,0)));
        REAL I1_val_denom  = abs(I1_SqJ2(i+1,0)-I1_SqJ2(i-1,0));
        REAL bvalue        = abs((log(SqJ2primeI1_n/SqJ2primeI1))/(I1_val_denom));
        
        sum_B_val += bvalue;
        m_B_val    = (sum_B_val)/(n_data-2);
    }
    
    bMean = GetBval();
    
}


void TF1DSAdjust::C_F1_function(REAL &cMean, TPZFMatrix<REAL> &I1_SqJ2){
    
    REAL bval = GetBval();
    
    I1_SqJ2        = m_I1_SqJ2;
    int64_t n_data = m_I1_SqJ2.Rows();
    REAL sum_C_val = 0.0;
    
    for(int64_t i = 1; i < n_data; i++)
    {
        REAL SqJ2_num  = (I1_SqJ2(i,1)-I1_SqJ2(i-1,1));
        REAL I1_denom  = exp(bval*I1_SqJ2(i,0))-exp(bval*I1_SqJ2(i-1,0));
        REAL cvalue    = abs((SqJ2_num)/(I1_denom));
        
        sum_C_val += cvalue;
        m_C_val    = (sum_C_val)/(n_data-1);
    }
    
    cMean = GetCval();
    
}


void TF1DSAdjust::A_F1_function(REAL &aMean, TPZFMatrix<REAL> &I1_SqJ2){
    
    REAL bval = GetBval();
    REAL cval = GetCval();
    
    I1_SqJ2        = m_I1_SqJ2;
    int64_t n_data = m_I1_SqJ2.Rows();
    REAL sum_A_val = 0.0;
    
    for(int64_t i = 0; i < n_data; i++)
    {
        REAL avalue  = I1_SqJ2(i,1)+abs(cval)*exp(abs(bval)*I1_SqJ2(i,0));
        
        sum_A_val   += avalue;
        m_A_val      = (sum_A_val)/(n_data);
    }
    
    aMean = GetAval();
    
}


STATE TF1DSAdjust::errorfunction(const std::vector<STATE> &input)
{
    STATE B = input[0];
    STATE C = exp(input[1]);
    STATE D = 0.;
    STATE W = 0.;
    STATE A = GetAval();
    STATE K(0.), G(0.), R(0.), phi(0.), N(0.), Psi(0.), kappa_0(0.);
    TPZSandlerExtended sandler(A,B,C,D,K,G,W,R,phi,N,Psi,kappa_0);
    int imin=0;
    STATE F1min = sandler.F(this->m_I1_SqJ2(imin,0))-log(abs(this->m_I1_SqJ2(imin,1)-GetAval()));
    for(int i=0; i<m_I1_SqJ2.Rows(); i++)
    {
        STATE FTrial = sandler.F(this->m_I1_SqJ2(imin,0))-log(abs(this->m_I1_SqJ2(imin,1)-GetAval()));
        if(FTrial < F1min)
        {
            imin = i;
            F1min = FTrial;
        }
    }
    A -= F1min;
    sandler.SetA(A);

    
    STATE error = 0.;
    for(int i=0; i<m_I1_SqJ2.Rows(); i++)
    {
        if(A<this->m_I1_SqJ2(i,1))
        {
            DebugStop();
        }
        STATE compare1 = log(A-this->m_I1_SqJ2(i,1));
        STATE compare2 = log(C)+this->m_I1_SqJ2(i,0)*B;
        error += (compare1-compare2)*(compare1-compare2);
    }
    return error;
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
    
        REAL a_val, b_val, c_val;
    
        B_F1_function(b_val, m_I1_SqJ2);
        C_F1_function(c_val, m_I1_SqJ2);
        A_F1_function(a_val, m_I1_SqJ2);

        m_Sandler.SetB(b_val);
        m_Sandler.SetC(c_val);
        m_Sandler.SetA(a_val);
    
        m_I1_SqJ2.Print("TestData",std::cout);
        std::cout << "Objective parameters A = " << m_Sandler.A() << " B = " << m_Sandler.B() << " C = " << m_Sandler.C() << std::endl;
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
    nlopt::opt opt(nlopt::LN_NEWUOA_BOUND, 2);
    opt.set_min_objective(myvfunc, this);
    opt.set_xtol_rel(1e-6);
    std::vector<double> x(2,0.);
    std::vector<double> lb(2);
    lb[0] = 0.; lb[1] = 0.;
    opt.set_lower_bounds(lb);
    std::vector<double> ub(2);
    ub[0] = 2.*GetBval();
    ub[1] = 2.*GetCval();
    opt.set_upper_bounds(ub);
    
    /// Initialize the B and C value
    x[0] = GetBval();
    x[1] = GetCval();
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


