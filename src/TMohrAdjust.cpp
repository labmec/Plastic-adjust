//
//  TMohrAdjust.cpp
//  plastic_adjust
//
//  Created by Manouchehr Sanei on 4/23/19.
//
//

#include "TMohrAdjust.h"
#include "nlopt.hpp"
#include<math.h>


TMohrAdjust::TMohrAdjust(){
    
    m_coh = 0.;
    m_phi = 0.;
    
}

TMohrAdjust::TMohrAdjust(const TMohrAdjust &other) : m_mohr(other.m_mohr), m_sign_tau(other.m_sign_tau)
{
    
    m_coh = other.m_coh;
    m_phi = other.m_phi;

}

const TMohrAdjust & TMohrAdjust::operator=(const TMohrAdjust &other)
{
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_coh = other.m_coh;
    m_phi = other.m_phi;
    
    return *this;
}

TMohrAdjust::~TMohrAdjust(){
    
}

void TMohrAdjust::Populate()
{
    /// Render data; temporary
    const int n_sign = 4;

    REAL sign_data[n_sign] = {1.8, 9.7, 21.5, 34.5};
    REAL tau_data[n_sign]  = {3.1, 7.4, 11.5, 14.0};
    
    m_sign_tau.Redim(n_sign, 2);
        
    for(int i=0; i<n_sign; i++)
    {
        REAL sig = sign_data[i];
        REAL tau = tau_data[i];
        m_sign_tau(i,0) = sig;
        m_sign_tau(i,1) = tau;
    }
    
    /// Initial guess
    REAL phi_val = compute_phi_initial();
    SetPhival(phi_val);
    REAL coh_val = compute_coh_initial();
    SetCohval(coh_val);
    
    m_mohr.SetPhi(phi_val);
    m_mohr.SetCohesion(coh_val);
    
    m_sign_tau.Print("TestData",std::cout);
    std::cout << "Objective parameters phi = " << m_mohr.Phi() << " coh = " << m_mohr.Cohesion() << std::endl;
}


REAL TMohrAdjust::compute_phi_initial(){
    
    TPZFMatrix<REAL> sign_tau = m_sign_tau;
    int64_t n_data = m_sign_tau.Rows();
    REAL sum_phi_val = 0.0;
    
    for(int64_t i = 1; i < n_data; i++)
    {
        REAL tauval  = (sign_tau(i,1)-sign_tau(i-1,1));
        REAL signval = (sign_tau(i,0)-sign_tau(i-1,0));
        REAL phival  = atan(abs(tauval/signval));
        sum_phi_val += (phival*180/3.1415);
    }
    
    return sum_phi_val/(n_data-1);
}


REAL TMohrAdjust::compute_coh_initial(){
    
    TPZFMatrix<REAL> sign_tau  = m_sign_tau;
    int64_t n_data = m_sign_tau.Rows();
    REAL phival = Phival();
    REAL sum_coh_val = 0.0;
    
    for(int64_t i = 0; i < n_data; i++)
    {
        REAL avalue  = sign_tau(i,1)-(sign_tau(i,0)*tan(phival*(3.1415/180)));
        sum_coh_val  += avalue;
    }
    
    return sum_coh_val/n_data;
}



STATE TMohrAdjust::errorfunctionMC(const std::vector<STATE> &input)
{
    
    REAL friction = input[0];
    REAL cohesion = input[1];

    STATE errorMC=0;
    int64_t n_data = m_sign_tau.Rows();
    
    for(int i=0; i<n_data; i++)
    {

        STATE minval = ComputeShearStress(this->m_sign_tau(i,0), friction, cohesion)-this->m_sign_tau(i,1);
        
        errorMC += minval*minval;
        
    }
    
    return errorMC;
}

void TMohrAdjust::gradientfunctionMC(const std::vector<STATE> &input, std::vector<double> &grad)
{
    STATE phi   = input[0];
    STATE coh   = input[1];
    int64_t n_data = m_sign_tau.Rows();
    
    STATE grad0 = 0.;
    STATE grad1 = 0.;
    
    for(int i=1; i<n_data; i++)
    {
        REAL sign = m_sign_tau(i,0);
        REAL tau  = m_sign_tau(i,1);
        
        grad0 += -2*sign*pow((1/cos(phi)),2)*(-m_coh+tau-sign*tan(phi));
        grad1 += -2*(-coh+tau-sign*tan(phi));
    }
    
    grad[0] = grad0;
    grad[1] = grad1;
    
}

double myvfunctionMC(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    TMohrAdjust *loc = (TMohrAdjust *) my_func_data;
    if (!grad.empty()) {
        
        loc->gradientfunctionMC(x, grad);
    }
    
    double err = loc->errorfunctionMC(x);
    for(int i=0; i<x.size(); i++) std::cout << "x[" << i << "]= " << x[i] << " ";
    std::cout << "error = " << err << std::endl;
    return err;
}


/// NLopt functions are: LN_NEWUOA_BOUND, LN_NEWUOA, LN_PRAXIS, LD_SLSQP

void TMohrAdjust::Adjust()
{
    nlopt::opt opt(nlopt::LN_NEWUOA_BOUND, 2);
    opt.set_min_objective(myvfunctionMC, this);
    opt.set_xtol_rel(1e-6);
    std::vector<double> x(2,0.);
    std::vector<double> lb(2);
    lb[0] = 0.; lb[1] = 0.;
    opt.set_lower_bounds(lb);
    std::vector<double> ub(2);
    ub[0] = 2.*Phival();
    ub[1] = 2.*Cohval();
    
    opt.set_upper_bounds(ub);
    
    /// Initialize the phi and cohesion value
    x[0] = Phival();
    x[1] = Cohval();

    
    double minf;
    
    try{
        nlopt::result result = opt.optimize(x, minf);
        
        if (result < 0) {
            std::cerr << "nlopt failed: result = " << result << std::endl;
        } else {
            
            std::cout << std::endl;
            std::cout << "found minimum of MC parameters : ";
            for(int i=1; i<x.size(); i++)
                std::cout << "phi = " << x[i-1] << ", " << "cohesion = " << x[i] << std::endl;
            std::cout << std::endl;
            std::cout<< std::setprecision(10) << "min_val = " << minf << std::endl;
        }
    }
    catch(std::exception &e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }
    
}


void TMohrAdjust::Hessian(TPZFMatrix<REAL> & Hessian, REAL signVal, REAL tauVal){
    
    Hessian.Resize(2,2);
    Hessian.Zero();
    
    Hessian(0,0) += 2*pow(signVal,2)*pow((1/cos(m_phi)),4)-4*signVal*pow((1/cos(m_phi)),2)*tan(m_phi)*(-m_coh +tauVal-signVal*tan(m_phi));
    Hessian(0,1) += 2*signVal*pow((1/cos(m_phi)),2);
    Hessian(1,0) += 2*signVal*pow((1/cos(m_phi)),2);
    Hessian(1,1) += 2;
    
}

void TMohrAdjust::Residual(TPZFMatrix<REAL> &Residual, REAL signVal, REAL tauVal){
    
    Residual.Resize(2,1);
    Residual.Zero();
    
    Residual(0,0) += -2*signVal*pow((1/cos(m_phi)),2)*(-m_coh+tauVal-signVal*tan(m_phi));
    Residual(1,0) += -2*(-m_coh+tauVal-signVal*tan(m_phi));
    
}

STATE TMohrAdjust::Assemble(TPZFMatrix<REAL> &sign_tau, TPZFMatrix<REAL> &hessian, TPZFMatrix<REAL> &res)
{
    STATE errorMC=0;
    sign_tau = m_sign_tau;
    int64_t n_data = m_sign_tau.Rows();
    
    TPZFMatrix<REAL> hessianMatrix(2,2);
    TPZFMatrix<REAL> residual(2,1);
    
    hessian.Resize(2,2);
    res.Resize(2,1);
    
    hessian.Zero();
    res.Zero();
    
    REAL friction = Phival();
    REAL cohesion = Cohval();
    
    for(int i=0; i<n_data; i++)
    {
        STATE minval = ComputeShearStress(sign_tau(i,0), friction, cohesion)-sign_tau(i,1);
    
        Hessian(hessianMatrix, sign_tau(i,0), sign_tau(i,1));
        Residual(residual, sign_tau(i,0), sign_tau(i,1));
        
        hessian += hessianMatrix;
        res += residual;
        errorMC += minval*minval;
    }
    
    return errorMC;
}


void TMohrAdjust::Adjust2()
{
    TPZFMatrix<REAL> sign_tau = m_sign_tau;
    
    TPZFMatrix<REAL> res(2,1,0.);
    TPZFMatrix<REAL> hessian(2,2,0.);
    
    REAL error = Assemble(sign_tau, hessian, res);
    
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
        hessian.SolveDirect(res,ECholesky);
        res *= -1.;
        LoadCorrection(res);
        
        REAL error = Assemble(sign_tau, hessian, res);
        
        normres = Norm(res);
        std::cout << "error = " << error << " normres = " << normres << std::endl;
        count++;
    }
    
    std::cout << std::endl;
    std::cout << "found minimum of MC parameters: ";
    std::cout << "phi = (" << Phival() << "), coh = (" << Cohval() << ")"<< std::endl;
    std::cout << std::endl;
    std::cout.flush();
    
}


void TMohrAdjust::LoadCorrection(TPZFMatrix<REAL> &delx)
{
    REAL phi = Phival()+delx(0,0);
    SetPhival(phi);
    REAL coh = Cohval()+delx(1,0);
    SetCohval(coh);

}

STATE TMohrAdjust::ComputeShearStress(REAL &sig_n, REAL &phi, REAL &coh){
    
    REAL tau = sig_n*tan(phi)+coh;
    
    return tau;
}
