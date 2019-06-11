//
//  TF2DSAdjust_DW.cpp
//  plastic_adjust
//
//  Created by Manouchehr Sanei on 4/11/19.
//
//

#include <math.h>
#include "TF2DSAdjust_DW.h"
#include "StubFunctions.h"


TF2DSAdjust_DW::TF2DSAdjust_DW(){
    
    m_Kbulk = 0.;
    m_Wval  = 0.;
    m_Dval  = 0.;
    m_iWval = 0.;
    
}

TF2DSAdjust_DW::TF2DSAdjust_DW(const TF2DSAdjust_DW &other) : m_X_epsPv(other.m_X_epsPv)
{
    m_Kbulk = other.m_Kbulk;
    m_Wval  = other.m_Wval;
    m_Dval  = other.m_Dval;
    m_iWval = other.m_iWval;

}

const TF2DSAdjust_DW & TF2DSAdjust_DW::operator=(const TF2DSAdjust_DW &other)
{
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_X_epsPv = other.m_X_epsPv;
    m_Kbulk   = other.m_Kbulk;
    m_Wval    = other.m_Wval;
    m_Dval    = other.m_Dval;
    m_iWval   = other.m_iWval;
    
    return *this;
}

TF2DSAdjust_DW::~TF2DSAdjust_DW(){
    
}

/// Computing D and W are done in two steps:
/// First: Plot X vs epsPv (It is computed using the data of hydrostatic test)
/// Second: Select the right data of X vs epsPv from last step to adjust D and W


/// ************************ First Setp
void TF2DSAdjust_DW::ComputedElasticStrain(TPZVec<TTestSection> &active, std::vector<REAL> &epsEV_data)
{
    REAL K = m_Kbulk;

    TPZFMatrix<REAL> deform,stress;
    for(int64_t i=0; i<active.size(); i++)
    {
        active[i].GetData(deform,stress);
        int64_t ndata = deform.Rows();        
        REAL conf_stress;

        for(int64_t d = 0; d<ndata; d++)
        {
            conf_stress = stress(d,0);
            REAL epsEV = conf_stress/K;
            
            epsEV_data.push_back(epsEV);
   
        }
        
    }
}

void TF2DSAdjust_DW::ComputedPlasticStrain(TPZVec<TTestSection> &active, std::vector<REAL> &epsPv_data)
{
    std::vector<REAL> epsElV_data;
    std::vector<REAL> epsTV_data;
    TPZFMatrix<REAL> deform,stress;
    ComputedElasticStrain(active, epsElV_data);
    epsPv_data.assign(epsElV_data.size(), 0);
    
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
    
    std::vector<REAL> plastStrain = a_subtract_b(epsTV_data,epsElV_data);
    epsPv_data = plastStrain;
    
}

std::vector<REAL> TF2DSAdjust_DW::a_subtract_b(std::vector<REAL> & a, std::vector<REAL> & b){
    if (a.size()!=b.size()) {
        DebugStop();
    }
    std::vector<REAL> y(a.size());
    for (int i = 0; i < a.size(); i++) {
        y[i] = a[i]-b[i];
    }
    return y;
}


void TF2DSAdjust_DW::Data_X_vs_EpsPv(TPZVec<TTestSection> &active, std::vector<REAL> &epsPV_data, std::vector<REAL> &X_data, TPZFMatrix<REAL> &X_epsPv){
    
    TPZFMatrix<REAL> InvariantStress;
    for(int64_t i=0; i<active.size(); i++)
    {
        active[i].GetInvStressData(InvariantStress);
        int64_t ndata = InvariantStress.Rows();
        
        for(int64_t j = 0; j<ndata; j++)
        {
            REAL I1   = InvariantStress(j,0);
            X_data.push_back(I1);
        }
    }
    
    ComputedPlasticStrain(active, epsPV_data);
    
    REAL numdata = epsPV_data.size();
    
    X_epsPv.Resize(numdata, 2);
    
    for (int64_t k = 0; k < numdata; k++) {
        
        REAL epspvol = epsPV_data[k];
        REAL xval    = X_data[k];
        
        X_epsPv(k,0) = xval;
        X_epsPv(k,1) = epspvol;
        
    }
    
    X_epsPv.Print("X_epsPv = ",std::cout,EMathematicaInput);

}


void TF2DSAdjust_DW::PopulateElastic()
{
    // Set elastic data
    REAL K = 3191.68;
    SetBulkModulus(K);
    
}

void TF2DSAdjust_DW::PlotXvsEpsPv(TPZVec<TTestSection> &active)
{
    std::vector<REAL> epsPV_data;
    std::vector<REAL> X_data;
    TPZFMatrix<REAL> X_epsPv_data;
    
    
    Data_X_vs_EpsPv(active, epsPV_data, X_data, X_epsPv_data);
}



/// ************************ Second Setp


void TF2DSAdjust_DW::PopulateDW()
{

    /// Render the X vs epsPv from last step;
    const int n_I1 = 3;
    
    REAL X_data[n_I1]     = {-60.8142, -78.7041, -96.8918};
    REAL epsPv_data[n_I1] = {-0.0185137, -0.0210985, -0.0236779};
    
    m_X_epsPv.Redim(n_I1,2);
    
    for(int i=0; i<n_I1; i++)
    {
        REAL X = X_data[i];
        REAL epsPv = epsPv_data[i];
        m_X_epsPv(i,0) = X;
        m_X_epsPv(i,1) = epsPv;
    }
    
    /// Initial guess for parameters: D, W, and iW that is Log (D*W)
    REAL d_val = ComputeDval_initial();
    SetDval(d_val);
    
    REAL w_val = ComputeWval_initial();
    
    REAL iwval = log(d_val*w_val);
    SetiWval(iwval);
    
}


REAL TF2DSAdjust_DW::ComputeDval_initial(){
    
    TPZFMatrix<REAL> X_epsPv = m_X_epsPv;
    int64_t n_data = m_X_epsPv.Rows();
    REAL sum_D_val = 0.0;
    
    for(int64_t i = 1; i < n_data-1; i++)
    {
        REAL dval1 = (X_epsPv(i-1,1)-X_epsPv(i,1))/(X_epsPv(i-1,0)-X_epsPv(i,0));
        REAL dval2 = (X_epsPv(i+1,1)-X_epsPv(i,1))/(X_epsPv(i+1,0)-X_epsPv(i,0));
        REAL num   = log(dval1/dval2);
        REAL denum = ((X_epsPv(i-1,0)+X_epsPv(i,0))/2)-((X_epsPv(i+1,0)+X_epsPv(i,0))/2);
        
        REAL dval  = num/denum;
        sum_D_val += dval;
    }
    
    return sum_D_val/(n_data-1);;
}


REAL TF2DSAdjust_DW::ComputeWval_initial(){
    
    TPZFMatrix<REAL> X_epsPv = m_X_epsPv;
    int64_t n_data = m_X_epsPv.Rows();
    REAL dval = ComputeDval_initial();
    REAL sum_W_val = 0.0;
    
    for(int64_t i = 1; i < n_data; i++)
    {
        REAL num   = (X_epsPv(i,1)-X_epsPv(i-1,1))/(X_epsPv(i,0)-X_epsPv(i-1,0));
        REAL denum = dval*exp(dval*X_epsPv(i,0));
        
        REAL dval  = num/denum;
        sum_W_val += dval;
    }
    
    return sum_W_val/(n_data-1);;
}



void TF2DSAdjust_DW::Hessian_DW(TPZFMatrix<REAL> &Hessian, REAL &X, REAL &depsPvX){
    
    Hessian.Resize(2,2);
    Hessian.Zero();
    
    Hessian(0,0) += 2.;
    Hessian(0,1) += 2*X;
    Hessian(1,0) += 2*X;
    Hessian(1,1) += 2*pow(X,2);
    
}


void TF2DSAdjust_DW::Residual_DW(TPZFMatrix<REAL> &Residual, REAL &X, REAL &depsPvX){
    
    Residual.Resize(2,1);
    Residual.Zero();
    
    Residual(0,0) += -2*(-m_iWval - m_Dval*X + log(depsPvX));
    Residual(1,0) += -2*X*(-m_iWval - m_Dval*X + log(depsPvX));
    
}


STATE TF2DSAdjust_DW::AssembleDW(TPZFMatrix<REAL> &X_epsPv, TPZFMatrix<REAL> &hessian, TPZFMatrix<REAL> &res)
{
    
    X_epsPv = m_X_epsPv;
    int64_t n_data = m_X_epsPv.Rows();
    
    STATE errorDW = 0 ;
    
    TPZFMatrix<REAL> hessianMatrix(2,2);
    TPZFMatrix<REAL> residual(2,1);
    
    hessian.Resize(2,2);
    res.Resize(2,1);
    
    hessian.Zero();
    res.Zero();
    
    for(int i=1; i<n_data; i++)
    {
        REAL xval   = X_epsPv(i,0);
        REAL depsPvX = (X_epsPv(i,1)-X_epsPv(i-1,1))/(X_epsPv(i,0)-X_epsPv(i-1,0));
        
        STATE minDW = costFunction_DW(xval, depsPvX);
        
        Hessian_DW(hessianMatrix, xval, depsPvX);
        Residual_DW(residual, xval, depsPvX);
        
        hessian += hessianMatrix;
        res     += residual;
        errorDW += minDW*minDW;
        
        
    }
    return errorDW;
}



void TF2DSAdjust_DW::AdjustDW()
{
 
    TPZFMatrix<REAL> X_epsP = m_X_epsPv;
    TPZFMatrix<REAL> res(2,1,0.);
    TPZFMatrix<REAL> hessian(2,2,0.);

    REAL error = AssembleDW(X_epsP, hessian, res);
    
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
        LoadCorrection_DW(res);
        
        REAL error = AssembleDW(X_epsP, hessian, res);
        
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


STATE TF2DSAdjust_DW::costFunction_DW(REAL &X, REAL &depsPvX){
    
    REAL Dvalue  = Dval();
    REAL iWvalue = iWval();
    
    STATE cost = -iWvalue - Dvalue*X + log(depsPvX);
    return cost;
}

void TF2DSAdjust_DW::LoadCorrection_DW(TPZFMatrix<REAL> &delx)
{
    REAL d_val = Dval()+delx(0,0);
    SetDval(d_val);
    REAL iw_val = iWval()+delx(1,0);
    SetiWval(iw_val);
    
}


