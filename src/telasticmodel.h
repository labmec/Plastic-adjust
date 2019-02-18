#ifndef TELASTICMODEL_H
#define TELASTICMODEL_H

#include "TPZElasticResponse.h"
#include "fadType.h"

class TTestSection;

class TElasticModel
{
public:
    enum EAdjustType {EYoungPoisson,ELambdaMu};

private:

    TPZElasticResponse fER;
    TPZVec<REAL> fSig0;

    REAL Assemble(TPZVec<TTestSection> &active, TPZFMatrix<REAL> &tangent, TPZFMatrix<REAL> &res);

    REAL Assemble2(TPZVec<TTestSection> &active, TPZFMatrix<REAL> &tangent, TPZFMatrix<REAL> &res);

    void LoadCorrection(TPZFMatrix<REAL> &delu);

public:

    TElasticModel();

    TElasticModel(const TElasticModel &cp) : fER(cp.fER), fSig0(cp.fSig0)
    {

    }

    TPZElasticResponse ElasticResponse()
    {
        return fER;
    }

    TElasticModel &operator=(const TElasticModel &cp)
    {
        fER = cp.fER;
        fSig0 = cp.fSig0;
        return *this;
    }

    void SetYoungPoisson(REAL E, REAL poisson)
    {
        fER.SetUp(E, poisson);
    }

    void Compute(TPZVec<REAL> &eps, TPZVec<REAL> &sig, TPZVec<REAL> &dsigr, TPZVec<REAL> &dsiga);

    void Adjust(TPZVec<TTestSection> &active);

    void ComputedSigma(TTestSection &sec, int i, TPZMatrix<REAL> &stress);


};

#endif // TELASTICMODEL_H
