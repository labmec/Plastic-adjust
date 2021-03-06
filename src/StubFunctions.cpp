//
// Created by gus on 04/02/19.
//

#include "StubFunctions.h"
#include "pzstack.h"
#include "telasticmodel.h"

#include "TElasticAdjust.h"
#include "TF2DSAdjust_DW.h"
#include "TF2DSAdjust_R_RHW.h"




TGlob glob;

void open(std::string filePath, std::string nickname) {
    std::string complete = "../Inputdata/"+filePath;
    std::cout << "open(" << complete << ", " << nickname << ");" << std::endl;
    TTestData data(complete,nickname);
    glob.fMeasure[nickname] = data;
}

void model(std::string modelName, std::string test) {
    std::cout << "model(" << modelName << ", " << test << ");" << std::endl;
    if(modelName.compare("elasticity") == 0)
    {
        glob.m_Adjust = EElasticResponse;
    } else
    if(modelName.compare("CapDSplasticity") == 0)
    {
        glob.m_Adjust = DiMaggioSandlerF2Response;
    }
    if(modelName.compare("CapDSfindRinRHW") == 0)
    {
        glob.m_Adjust = DiMaggioSandlerF2findRinRHW;
    }
}

void select(std::string parameter, int initial, int final, std::string comment) {
    std::cout << "select(" << parameter << ", " << initial << ", " << final << ", " << comment << ");" << std::endl;
    if(glob.fMeasure.find(parameter) == glob.fMeasure.end())
    {
        std::cout << "SELECTION NOT FOUND ******\n";
    }
    else
    {
        TTestData *test = &glob.fMeasure[parameter];
        TTestSection selected(test,initial,final,comment);
        glob.m_Active.Push(selected);
    }
}

void adjust(std::string parameter) {
    std::cout << "adjust(" << parameter << ");" << std::endl;
    if(glob.m_Adjust == EElasticResponse)
    {
        TElasticAdjust model;
        model.Adjust(glob.m_Active);
        glob.m_ER = model.ElasticResponse();
        {
            int nactive = glob.m_Active.size();
            std::ofstream out("Adjust.nb");
            for(int i=0; i<nactive; i++)
            {
                TPZFMatrix<REAL> deform,stress_meas,stress_comp;
                glob.m_Active[i].GetData(deform,stress_meas);
                model.ComputedSigma(glob.m_Active[i],i,stress_comp);
                deform.Print("Deform = ",out,EMathematicaInput);
                stress_meas.Print("StressMeas = ",out,EMathematicaInput);
                stress_comp.Print("StressComp = ",out,EMathematicaInput);
            }
        }
    }
    if(glob.m_Adjust == DiMaggioSandlerF2Response)
    {
        
        TF2DSAdjust_DW model;
        model.PopulateElastic();
        model.PlotXvsEpsPv(glob.m_Active);
                
        {
            int nactive = glob.m_Active.size();
            std::ofstream out("Adjust.nb");
            for(int i=0; i<nactive; i++)
            {
                TPZFMatrix<REAL> deform,stress_meas,stress_comp,invariantStress;
                glob.m_Active[i].GetData(deform,stress_meas);
                glob.m_Active[i].GetInvStressData(invariantStress);
                
                deform.Print("Deform = ",out,EMathematicaInput);
                invariantStress.Print("InvariantStress = ",out,EMathematicaInput);
                
            }
        }
    }
    
    if(glob.m_Adjust == DiMaggioSandlerF2findRinRHW)
    {
        TF2DSAdjust_R_RHW model;
        model.PopulateR();
        model.AdjustR(glob.m_Active);
        
        {
            int nactive = glob.m_Active.size();
            std::ofstream out("Adjust.nb");
            for(int i=0; i<nactive; i++)
            {
                TPZFMatrix<REAL> deform,stress_meas;
                glob.m_Active[i].GetData(deform,stress_meas);
                
                deform.Print("Deform = ",out,EMathematicaInput);
                
            }
        }
    }
    
    
    else
    {
        std::cout << "I dont know what to adjust *****\n";
    }
}

void assign(std::string parameter, std::string fileNickname) {
    std::cout << "assign(" << parameter << ", " << fileNickname << ");" << std::endl;
}

void clear() {
    std::cout << "clear();" << std::endl;
}

TTestData::TTestData(std::string &filename, std::string &nickname) : m_NickName(nickname), m_FileName(filename)
{
    int64_t numdata = NumData();
    m_Deform.Redim(numdata,2);
    m_Stress.Redim(numdata,2);
    m_InvariantStress.Redim(numdata,2);
    std::ifstream input(filename);
    ReadHeader(input);
    for(int64_t i = 0; i<numdata; i++)
    {
        ReadLine(i,input);
    }

}
//TPZFMatrix<double> m_Deform;
//TPZFMatrix<double> m_Stress;

//std::string m_NickName;

//std::string m_FileName;

//TPZVec<EDataType> m_Types;

//// 0 -> sigc
//// 1 -> siga
//// 2 -> epsc
//// 3 -> epsa
//std::map<int,int> m_Relevant;



TTestData::TTestData() : m_Deform(), m_Stress(), m_InvariantStress(),m_NickName(), m_FileName(), m_Types(), m_Relevant()
{

}

TTestData::TTestData(const TTestData &cp) : m_Deform(cp.m_Deform), m_Stress(cp.m_Stress), m_InvariantStress(cp.m_InvariantStress), m_NickName(cp.m_NickName), m_FileName(cp.m_FileName), m_Types(cp.m_Types), m_Relevant(cp.m_Relevant)
{

}

TTestData &TTestData::operator=(const TTestData &cp)
{
    m_Deform   = cp.m_Deform;
    m_Stress   = cp.m_Stress;
    m_InvariantStress = cp.m_InvariantStress;
    m_NickName = cp.m_NickName;
    m_FileName = cp.m_FileName;
    m_Types    = cp.m_Types;
    m_Relevant = cp.m_Relevant;
    
    return *this;
}


void TTestData::GetData(int64_t first, int64_t last, TPZFMatrix<double> &deform, TPZFMatrix<double> &stress)
{
    deform.Redim(last-first+1,2);
    stress.Redim(last-first+1,2);
    for(int64_t i=first; i<=last; i++)
    {
        deform(i-first,0) = - m_Deform(i,0)*0.01;
        deform(i-first,1) = - m_Deform(i,1)*0.01;
        stress(i-first,0) = - m_Stress(i,0);
        stress(i-first,1) = - (this->m_Stress(i,1));
    }
    
    
}


void TTestData::GetInvStressData(int64_t first, int64_t last, TPZFMatrix<double> &invariantStress)
{
    invariantStress.Redim(last-first+1,2);
    for(int64_t i=first; i<=last; i++)
    {
        invariantStress(i-first,0) = - m_InvariantStress(i,0);
        invariantStress(i-first,1) = this->m_InvariantStress(i,1);
        
    }
}


void TTestData::ReadHeader(std::istream &input)
{
    std::string headerline;
    std::getline(input,headerline);
    TPZStack<std::string> titles;
    std::stringstream strstream(headerline);
    while(strstream)
    {
        std::string loc;
        strstream >> loc;
        if(loc.size()) titles.Push(loc);
    }
    this->m_Types.resize(titles.size());
    //Time	siga	epsc	epsa	epsv	sigc	siga*	SqJ2	I1

    for(int i=0; i<titles.size(); i++)
    {
        if(titles[i].compare("Time") == 0){
            m_Types[i] = ETime;
        }
        else if(titles[i].compare("siga") == 0)
        {
            m_Types[i] = ESiga;
        }
        else if(titles[i].compare("epsc") == 0)
        {
            m_Types[i] = EEpsr;
            m_Relevant[2] = i;
        }
        else if(titles[i].compare("epsa") == 0)
        {
            m_Types[i] = ESiga;
            m_Relevant[3] = i;
        }
        else if(titles[i].compare("epsv") == 0)
        {
            m_Types[i] = EEpsv;
        }
        else if(titles[i].compare("sigc") == 0)
        {
            m_Types[i] = ESigr;
            m_Relevant[0] = i;
        }
        else if(titles[i].compare("siga*") == 0)
        {
            m_Types[i] = ESigaStar;
            m_Relevant[1] = i;
        }
        else if(titles[i].compare("SqJ2") == 0)
        {
            m_Types[i] = ESqJ2;
            m_Relevant[5] = i;
        }
        else if(titles[i].compare("I1") == 0)
        {
            m_Types[i] = EI1;
            m_Relevant[4] = i;
        }
        else
        {
            std::cout << "I dont recognize " << titles[i] << std::endl;
        }
    }
    if(m_Relevant.size() != 6)
    {
        std::cout << "Didnt find all relevant titles\n";
    }
    else
    {
        for(auto it = m_Relevant.begin(); it != m_Relevant.end(); it++)
        {
            std::cout << "mapping index " << it->first << " to " << it->second << std::endl;
        }
    }
}

void TTestData::ReadLine(std::int64_t index, std::istream &input)
{
    std::string headerline;
    std::getline(input,headerline);
    TPZStack<REAL> values;
    std::stringstream strstream(headerline);
    for(int i=0; i< m_Types.size(); i++)
    {
        REAL val;
        strstream >> val;
        values.Push(val);
    }
    m_Stress(index,0)  = values[m_Relevant[0]];
    m_Stress(index,1)  = values[m_Relevant[1]];
    m_Deform(index,0)  = values[m_Relevant[2]];
    m_Deform(index,1)  = values[m_Relevant[3]];
    m_InvariantStress(index,0) = values[m_Relevant[4]];
    m_InvariantStress(index,1) = values[m_Relevant[5]];
}

/// computes the number of lines in the file
int64_t TTestData::NumData()
{
    std::ifstream input(m_FileName);
    std::string line;
    int64_t counter = 0;
    std::getline(input,line);
    if(!input)
    {
        std::cout << __PRETTY_FUNCTION__ << " Corrupted file\n";
        return -1;
    }
    while(input)
    {
        std::getline(input,line);
        if(input && line.size() > 0) counter++;
    }
    return counter;
}

