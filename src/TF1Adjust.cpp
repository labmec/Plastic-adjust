//
//  TF1Adjust.cpp
//  plastic_adjust
//
//  Created by Philippe Devloo on 03/03/19.
//

#include "TF1Adjust.h"
#include "nlopt.hpp"

STATE TF1Adjust::errorfunction(const std::vector<STATE> &input)
{
    STATE B = input[0];
    STATE C = input[1]+1.;
    STATE D = 0.;
    STATE W = 0.;
    STATE A = 0.;
    STATE K(0.), G(0.), R(0.), phi(0.), N(0.), Psi(0.), kappa_0(0.);
    TPZSandlerExtended sandler(A,B,C,D,K,G,W,R,phi,N,Psi,kappa_0);
    int imin=0;
    STATE F1min = sandler.F(this->fI1_SqJ2(imin,0))-this->fI1_SqJ2(imin,1);
    for(int i=0; i<fI1_SqJ2.Rows(); i++)
    {
        STATE FTrial = sandler.F(this->fI1_SqJ2(i,0))-this->fI1_SqJ2(i,1);
        if(FTrial < F1min)
        {
            imin = i;
            F1min = FTrial;
        }
    }
    A -= F1min;
    sandler.SetA(A);
    STATE error = 0.;
    for(int i=0; i<fI1_SqJ2.Rows(); i++)
    {
        STATE FTrial = sandler.F(this->fI1_SqJ2(i,0))-this->fI1_SqJ2(i,1);
        if(FTrial < -1.e-10)
        {
            DebugStop();
        }
        error += fabs(FTrial);
    }
    return error;
}

STATE TF1Adjust::errorfunction2(const std::vector<STATE> &input)
{
    STATE B = input[0];
    STATE C = exp(input[1]);
    STATE D = 0.;
    STATE W = 0.;
    STATE A = 0.;
    STATE K(0.), G(0.), R(0.), phi(0.), N(0.), Psi(0.), kappa_0(0.);
    TPZSandlerExtended sandler(A,B,C,D,K,G,W,R,phi,N,Psi,kappa_0);
    int imin=0;
    STATE F1min = sandler.F(this->fI1_SqJ2(imin,0))-this->fI1_SqJ2(imin,1);
    for(int i=0; i<fI1_SqJ2.Rows(); i++)
    {
        STATE FTrial = sandler.F(this->fI1_SqJ2(i,0))-this->fI1_SqJ2(i,1);
        if(FTrial < F1min)
        {
            imin = i;
            F1min = FTrial;
        }
    }
    A -= F1min;
    sandler.SetA(A);
    STATE error = 0.;
    for(int i=0; i<fI1_SqJ2.Rows(); i++)
    {
        if(A<this->fI1_SqJ2(i,1))
        {
            DebugStop();
        }
        STATE compare1 = log(A-this->fI1_SqJ2(i,1));
        STATE compare2 = log(C)+this->fI1_SqJ2(i,0)*B;
        error += (compare1-compare2)*(compare1-compare2);
    }
    return error;
}

/// method to initialize the class with possible data
void TF1Adjust::Populate()
{
    fSandler.MCormicRanchSand(fSandler);
    const int imax = 10;
    const REAL IMIN = -10.;
    fI1_SqJ2.Redim(imax,2);
    for(int i=0; i<imax; i++)
    {
        REAL I1 = IMIN-IMIN*(i)/(imax-1);
        REAL F1 = fSandler.F(I1);
        fI1_SqJ2(i,0) = I1;
        fI1_SqJ2(i,1) = F1;
    }
    fI1_SqJ2.Print("TestData",std::cout);
    std::cout << "Objective parameters A = " << fSandler.A() << " B = " << fSandler.B() << " C = " << fSandler.C() <<
    " D = " << fSandler.D() <<
    " W = " << fSandler.W() << std::endl;
}



typedef struct {
    double a, b;
} my_constraint_data;

double myvfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    if (!grad.empty()) {
        DebugStop();
    }
    TF1Adjust *loc = (TF1Adjust *) my_func_data;
    double err = loc->errorfunction2(x);
    for(int i=0; i<x.size(); i++) std::cout << "x[" << i << "]= " << x[i] << " ";
    std::cout << "error " << err << std::endl;
    return err;
}


void TF1Adjust::Adjust()
{
    nlopt::opt opt(nlopt::LN_NEWUOA_BOUND, 2);
    opt.set_min_objective(myvfunc, this);
    opt.set_xtol_rel(1e-6);
    std::vector<double> x(2,0.);
    std::vector<double> lb(2);
    lb[0] = 0.; lb[1] = -5.;
    opt.set_lower_bounds(lb);
    std::vector<double> ub(2);
    ub[0] = 2.;
    ub[1] = 10.;
    opt.set_upper_bounds(ub);

    x[0] = 0.67;
    x[1] = log(0.18);
    x[0] = 0.0;
    x[1] = 0.0;
    double minf;
    
    try{
        nlopt::result result = opt.optimize(x, minf);
        if (result < 0) {
            std::cerr << "nlopt failed: result = " << result << std::endl;
        }
        std::cout << "found minimum at f(";
        for(int i=0; i<x.size(); i++) std::cout << x[i] << " ";
        std::cout << std::setprecision(10) << " min_val " << minf << std::endl;
    }
    catch(std::exception &e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }

}

