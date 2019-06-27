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

    TPZFMatrix<double> m_Deform;
    TPZFMatrix<double> m_Stress;
    TPZFMatrix<double> m_InvariantStress;

    std::string m_NickName;

    std::string m_FileName;

    TPZVec<EDataType> m_Types;

    // 0 -> sigr
    // 1 -> siga
    // 2 -> epsr
    // 3 -> epsa
    // 4 -> I1
    // 5 -> SqJ2
    std::map<int,int> m_Relevant;
    
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
    TTestData *m_Origin;
    int64_t fFirst, m_Last;
    std::string m_Comment;

public:

    TTestSection() : m_Origin(0), fFirst(-1), m_Last(-1)
    {

    }

    TTestSection(TTestData *orig, int64_t first, int64_t last, const std::string &comment) : m_Origin(orig), fFirst(first),
        m_Last(last), m_Comment(comment)
    {

    }

    TTestSection(const TTestSection &cp) : m_Origin(cp.m_Origin), fFirst(cp.fFirst), m_Last(cp.m_Last),m_Comment(cp.m_Comment)
    {

    }

    TTestSection &operator=(const TTestSection &cp)
    {
        m_Origin  = cp.m_Origin;
        fFirst   = cp.fFirst;
        m_Last    = cp.m_Last;
        m_Comment = cp.m_Comment;
        return *this;
    }

    void GetData(TPZFMatrix<double> &deform, TPZFMatrix<double> &stress)
    {
        m_Origin->GetData(fFirst,m_Last,deform,stress);
    }
    
    void GetInvStressData(TPZFMatrix<double> &invariantStress)
    {
        m_Origin->GetInvStressData(fFirst,m_Last,invariantStress);
    }
    
};

enum EAdjustModel {EElasticResponse, DiMaggioSandlerF2Response, DiMaggioSandlerF2findRinRHW};

struct TGlob
{
    TPZStack<TTestSection> m_Active;

    std::map<std::string,TTestData> fMeasure;

    EAdjustModel m_Adjust = EElasticResponse;

    TPZElasticResponse m_ER;
    
    TPZSandlerExtended m_Sandler;
};

extern TGlob glob;


#endif
