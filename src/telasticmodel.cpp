#include "telasticmodel.h"
#include "StubFunctions.h"

TElasticModel::TElasticModel()
{

}


void TElasticModel::Compute(TPZVec<REAL> &eps, TPZVec<REAL> &sig, TPZVec<REAL> &dsigr, TPZVec<REAL> &dsiga)
{
    sig[0] = fER.fLambda*(2.*eps[0]+eps[1])+2.*fER.fMu*eps[0];
    sig[1] = fER.fLambda*(2.*eps[0]+eps[1])+2.*fER.fMu*eps[1];
    dsigr[0] = 2.*eps[0]+eps[1];
    dsigr[1] = 2.*eps[0];
    dsigr[2] = 1.;
    dsigr[3] = 0;
    dsiga[0] = 2.*eps[0]+eps[1];
    dsiga[1] = 2.*eps[1];
    dsiga[2] = 0.;
    dsiga[3] = 1.;
}

void TElasticModel::Adjust(TPZVec<TTestSection> &active)
{
    fSig0.resize(2*active.size());
    fSig0.Fill(0.);
    TPZFMatrix<REAL> res(2*(active.size()+1),1,0.);
    TPZFMatrix<REAL> tangent(res.Rows(),res.Rows(),0.);
    REAL error = Assemble(active,tangent,res);
    std::cout << "Incoming error " << error << std::endl;
    {
        std::ofstream out("tangent.nb");
        tangent.Print("tangent",out,EMathematicaInput);
        std::cout << "NOTHING COMING OUT?\n";
    }
    int count = 0;
    REAL normres = Norm(res);
    std::cout << "Incoming normres = " << normres << std::endl;
    REAL tol = normres*1.e-8;
    while(normres > tol && count <10)
    {
        tangent.SolveDirect(res,ECholesky);
        res *= -1.;
        LoadCorrection(res);
        error = Assemble(active,tangent,res);
        normres = Norm(res);
        std::cout << "error = " << error << " normres = " << normres << std::endl;
        count++;
    }
    fER.Print(std::cout);
    std::cout.flush();
}

REAL TElasticModel::Assemble(TPZVec<TTestSection> &active, TPZFMatrix<REAL> &tangent, TPZFMatrix<REAL> &res)
{
    REAL error = 0.;
    tangent.Zero();
    res.Zero();
    TPZFMatrix<REAL> deform,stress;
    for(int64_t i=0; i<active.size(); i++)
    {
        active[i].GetData(deform,stress);
        int64_t ndata = deform.Rows();
        TPZManVector<REAL,2> eps(2),sig_measure(2), sig_c(2);
        TPZManVector<REAL,4> dsigr(4), dsiga(4);
        for(int64_t d = 0; d<ndata; d++)
        {
            eps[0] = deform(d,0);
            eps[1] = deform(d,1);
            sig_measure[0] = stress(d,0);
            sig_measure[1] = stress(d,1);
            Compute(eps,sig_c,dsigr,dsiga);
            sig_c[0] += fSig0[i*2];
            sig_c[1] += fSig0[i*2+1];
            TPZManVector<FADREAL,2> sigfad(2), diffmeasure(2);
            {
                FADREAL a(4,sig_c[0]);
                for(int i=0; i<4; i++) a.fastAccessDx(i) = dsigr[i];
                sigfad[0] = a;
            }
            {
                FADREAL a(4,sig_c[1]);
                for(int i=0; i<4; i++) a.fastAccessDx(i) = dsiga[i];
                sigfad[1] = a;
            }
            diffmeasure[0] = (sig_measure[1]-sig_measure[0]-sigfad[1]+sigfad[0]);
            diffmeasure[1] = (sig_measure[1]+2.*sig_measure[0]-sigfad[1]-2.*sigfad[0]);
            FADREAL pointerror(4,0.);
            pointerror = diffmeasure[0]*diffmeasure[0];
            pointerror += diffmeasure[1]*diffmeasure[1];
            error += pointerror.val();
            res(0,0) += pointerror.fastAccessDx(0);
            res(1,0) += pointerror.fastAccessDx(1);
            res(2*i+2,0) += pointerror.fastAccessDx(2);
            res(2*i+3,0) += pointerror.fastAccessDx(3);
            TPZManVector<int64_t,4> dest(4);
            dest[0] = 0;
            dest[1] = 1;
            dest[2] = 2*i+2;
            dest[3] = 2*i+3;
            for(int k=0; k<4; k++) for(int l=0; l<4; l++)
            {
                tangent(dest[k],dest[l]) += 2.*diffmeasure[0].fastAccessDx(k)*diffmeasure[0].fastAccessDx(l)+
                        2.*diffmeasure[1].fastAccessDx(k)*diffmeasure[1].fastAccessDx(l);
            }
        }
    }
    return error;
}

REAL TElasticModel::Assemble2(TPZVec<TTestSection> &active, TPZFMatrix<REAL> &tangent, TPZFMatrix<REAL> &res)
{
    REAL error = 0.;
    tangent.Zero();
    res.Zero();
    TPZFMatrix<REAL> deform,stress;
    for(int64_t i=0; i<active.size(); i++)
    {
        active[i].GetData(deform,stress);
        int64_t ndata = deform.Rows();
        TPZManVector<REAL,2> eps(2),sig_measure(2), sig_c(2);
        TPZManVector<REAL,4> dsigr(4), dsiga(4);
        for(int64_t d = 0; d<ndata; d++)
        {
            eps[0] = deform(d,0);
            eps[1] = deform(d,1);
            sig_measure[0] = stress(d,0);
            sig_measure[1] = stress(d,1);
            Compute(eps,sig_c,dsigr,dsiga);
            sig_c[0] += fSig0[i*2];
            sig_c[1] += fSig0[i*2+1];
            error += 2*(sig_measure[0]-sig_c[0])*(sig_measure[0]-sig_c[0])+(sig_measure[1]-sig_c[1])*(sig_measure[1]-sig_c[1]);
            res(0,0) += 4.*(sig_c[0]-sig_measure[0])*dsigr[0]+2.*(sig_c[1]-sig_measure[0])*dsiga[0];
            res(1,0) += 4.*(sig_c[0]-sig_measure[0])*dsigr[1]+2.*(sig_c[1]-sig_measure[0])*dsiga[1];
            res(2*i+2,0) += 4.*(sig_c[0]-sig_measure[0])*dsigr[2]+2.*(sig_c[1]-sig_measure[0])*dsiga[2];
            res(2*i+3,0) += 4.*(sig_c[0]-sig_measure[0])*dsigr[3]+2.*(sig_c[1]-sig_measure[0])*dsiga[3];
            TPZManVector<int64_t,4> dest(4);
            dest[0] = 0;
            dest[1] = 1;
            dest[2] = 2*i+2;
            dest[3] = 2*i+3;
            for(int k=0; k<4; k++) for(int l=0; l<4; l++)
            {
                tangent(dest[k],dest[l]) += 4.*dsigr[k]*dsigr[l]+2.*dsiga[k]*dsiga[l];
            }
        }
    }
    return error;
}

void TElasticModel::LoadCorrection(TPZFMatrix<REAL> &delu)
{
    fER.fLambda += delu(0,0);
    fER.fMu += delu(1,0);
    for(int i=0; i<fSig0.size(); i++)
    {
        fSig0[i] += delu(i+2,0);
    }
}

void TElasticModel::ComputedSigma(TTestSection &sec, int i, TPZMatrix<REAL> &stress)
{
    TPZFMatrix<REAL> deform,stress_meas;
    sec.GetData(deform,stress_meas);
    stress.Redim(deform.Rows(),2);
    int64_t ndata = deform.Rows();
    TPZManVector<REAL,2> eps(2),sig_c(2);
    TPZManVector<REAL,4> dsigr(4),dsiga(4);
    for(int64_t d = 0; d<ndata; d++)
    {
        eps[0] = deform(d,0);
        eps[1] = deform(d,1);
        Compute(eps,sig_c,dsigr,dsiga);
        sig_c[0] += fSig0[i*2];
        sig_c[1] += fSig0[i*2+1];
        stress(d,0) = sig_c[0];
        stress(d,1) = sig_c[1];
    }
}
