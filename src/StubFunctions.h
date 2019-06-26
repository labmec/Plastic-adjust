//
// Created by gus on 04/02/19.
//

#ifndef STUBFUNCTIONS_H
#define STUBFUNCTIONS_H

#include <iostream>
#include <string>
#include "pzfmatrix.h"
#include "pzstack.h"
#include "TPZElasticResponse.h"
#include "TPZSandlerExtended.h"



void open(std::string filePath, std::string nickname);

void model(std::string modelName, std::string test);

void select(std::string parameter, int initial, int final, std::string comment);

void adjust(std::string parameter);

void assign(std::string parameter, std::string fileNickname);

void clear();

class TTestData
{
public:

    enum EDataType {ETime, ESigr, ESiga, ESigaStar, EEpsr, EEpsa, EEpsaStar, ESigv, EEpsv, ESqJ2, EI1};

private:

    TPZFMatrix<double> fDeform;
    TPZFMatrix<double> fStress;
    TPZFMatrix<double> fInvariantStress;

    std::string fNickName;

    std::string fFileName;

    TPZVec<EDataType> fTypes;

    // 0 -> sigr
    // 1 -> siga
    // 2 -> epsr
    // 3 -> epsa
    // 4 -> I1
    // 5 -> SqJ2
    std::map<int,int> fRelevant;
    
    void ReadHeader(std::istream &input);

    void ReadLine(int64_t index, std::istream &input);

    // computes the number of lines in the file
    int64_t NumData();

public:

    TTestData(std::string &filename, std::string &nickname);

    TTestData();

    TTestData(const TTestData &cp);

    TTestData &operator=(const TTestData &cp);

    void GetData(int64_t first, int64_t last, TPZFMatrix<double> &deform, TPZFMatrix<double> &stress);
    void GetInvStressData(int64_t first, int64_t last, TPZFMatrix<double> &invariantStress);

};

class TTestSection
{
    TTestData *fOrigin;
    int64_t fFirst, fLast;
    std::string fComment;

public:

    TTestSection() : fOrigin(0), fFirst(-1), fLast(-1)
    {

    }

    TTestSection(TTestData *orig, int64_t first, int64_t last, const std::string &comment) : fOrigin(orig), fFirst(first),
        fLast(last), fComment(comment)
    {

    }

    TTestSection(const TTestSection &cp) : fOrigin(cp.fOrigin), fFirst(cp.fFirst), fLast(cp.fLast),fComment(cp.fComment)
    {

    }

    TTestSection &operator=(const TTestSection &cp)
    {
        fOrigin  = cp.fOrigin;
        fFirst   = cp.fFirst;
        fLast    = cp.fLast;
        fComment = cp.fComment;
        return *this;
    }

    void GetData(TPZFMatrix<double> &deform, TPZFMatrix<double> &stress)
    {
        fOrigin->GetData(fFirst,fLast,deform,stress);
    }
    
    void GetInvStressData(TPZFMatrix<double> &invariantStress)
    {
        fOrigin->GetInvStressData(fFirst,fLast,invariantStress);
    }
    
};

enum EAdjustModel {EElasticResponse, DiMaggioSandlerF2Response, DiMaggioSandlerF2findRinRHW};

struct TGlob
{
    TPZStack<TTestSection> fActive;

    std::map<std::string,TTestData> fMeasure;

    EAdjustModel fAdjust = EElasticResponse;

    TPZElasticResponse fER;
    
    TPZSandlerExtended m_Sandler;
};

extern TGlob glob;


#endif
