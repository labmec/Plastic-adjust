

// main.cpp : Defines the entry point for the console application.

#include <iostream>
#include <fstream>
#include <iomanip>

#include "json.hpp"
#include "StubFunctions.h"

/// Plasticity
#include "TPZPlasticStepPV.h"

/// Adjustment classes
#include "TF1DSAdjust.h"
#include "TMohrAdjust.h"
#include "TF2DSAdjust.h"


using json = nlohmann::json;
void printJSON(json commandFile, std::ostream& out);
void callMethods(json commandFile);
void translateToFunction(json singleCommand);


/// Method to present elastoplastic model
void LoadElastoPlasticModel();

// Read experimental data
TPZFMatrix<REAL> Read_Duplet(int n_data, std::string file);



int main() {
    
//    TF2DSAdjust F2;
//    F2.PopulateR();
//    F2.AdjustR();
//    return 0;
    
    LoadElastoPlasticModel();
    return 0;
    
//    TF1DSAdjust F1;
//    F1.Populate();
//    F1.Adjust2();
//    return 0;
    
    
//    TMohrAdjust MC;
//    MC.Populate();
//    MC.Adjust2();
//    return 0;
    
    
    std::string path;
    std::ifstream input;

    path = "../input_Oed8.json";
    input.open(path.c_str());

    while (!input.is_open()) {
        std::cout << ":: Oops! Let's try again with an actual file. Please write down its path:" << std::endl;
        std::cin >> path;
        input.open(path.c_str());
    }

    std::cout << ":: The JSON file was successfully read! Check it below: " << std::endl << std::endl;

    json commands;
    input >> commands;
//    printJSON(commands, std::cout);

    std::cout << std::endl << ":: The following methods would be called: " << std::endl << std::endl;
    callMethods(commands);
    

    return 0;
}

void printJSON(json commandFile, std::ostream& output) {
    output << std::setw(4) << commandFile << std::endl;
    output << std::flush;
}

void callMethods(json commandFile) {
    for (int i = 0; i < commandFile.size(); i++) {
        translateToFunction(commandFile[i]);
    }
}

void translateToFunction(json singleCommand) {
    std::string functionName = singleCommand[0].get<std::string>();

    if (functionName == "open") {
        std::string filePath = singleCommand[1].get<std::string>();
        std::string fileNickname = singleCommand[2].get<std::string>();

        open(filePath, fileNickname);
        return;
    }

    if (functionName == "model") {
        std::string modelName = singleCommand[1].get<std::string>();
        std::string test = singleCommand[2].get<std::string>();

        model(modelName, test);
        return;
    }

    if (functionName == "select") {
        std::string fileNickname = singleCommand[1].get<std::string>();
        int initialTime = singleCommand[2].get<int>();
        int finalTime = singleCommand[3].get<int>();
        std::string comment = singleCommand[4].get<std::string>();

        select(fileNickname, initialTime, finalTime, comment);
        return;
    }

    if (functionName == "adjust") {
        std::string parameter = singleCommand[1].get<std::string>();

        adjust(parameter);
        return;
    }

    if (functionName == "assign") {
        std::string parameter = singleCommand[1].get<std::string>();
        std::string fileNickname = singleCommand[2].get<std::string>();

        assign(parameter, fileNickname);
        return;
    }

    if (functionName == "clear") {
        clear();
        return;
    }

    std::cout << "Invalid function name, please check your JSON file." << std::endl;
    return;
}


void LoadElastoPlasticModel()
{
    
    // Experimental data
    std::string file_name;
    file_name = "/Users/manouchehr/Documents/GitHub/Plastic-adjust/exp_data/CP08.txt";
    int64_t n_data = 2000;
    TPZFMatrix<REAL> data = Read_Duplet(n_data, file_name);
    //    data.Print(std::cout);
    
    // DS Dimaggio Sandler PV
    TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse> LEDS;
    
    // LE Linear elastic response
    TPZElasticResponse ER;
    
    
        /// number CP08
        REAL E  = 2000.0; // MPa
        REAL nu = 0.2; // MPa
    
        STATE G = E / (2. * (1. + nu));
        STATE K = E / (3. * (1. - 2 * nu));
        REAL CA      = 14.0;
        REAL CB      = 0.020;
        REAL CC      = 13.0;
        REAL CD      = 0.013;
        REAL CR      = 2.0;
        REAL CW      = 0.04;
        REAL X_0     = -40.0;
        REAL phi = 0, psi = 1., N = 0;
    
        REAL Pc = X_0/3.0;
    
    
    ER.SetEngineeringData(E, nu);
    
    LEDS.SetElasticResponse(ER);
    LEDS.fYC.SetUp(CA, CB, CC, CD, K, G, CW, CR, phi, N, psi);
    
    TPZTensor<REAL> epsilon_t,sigma,sigma_target;
    sigma.Zero();
    epsilon_t.Zero();
    
    
    sigma.XX() = Pc;
    sigma.YY() = Pc;
    sigma.ZZ() = Pc;
    
    // Initial damage data
    REAL k_0;
    LEDS.InitialDamage(sigma, k_0);
    LEDS.fN.m_hardening = k_0;
    
    
    epsilon_t.Zero();
    
    TPZFNMatrix<2575,STATE> LEDS_epsilon_stress(n_data,4);
    for (int64_t id = 0; id < n_data; id++) {
        
        
        // For a given strain
        epsilon_t.XX() = -0.01*data(id,0);
        epsilon_t.YY() = -0.01*data(id,1);
        epsilon_t.ZZ() = -0.01*data(id,1);
        
        LEDS.ApplyStrainComputeSigma(epsilon_t, sigma_target);
        
        
        LEDS_epsilon_stress(id,0) = epsilon_t.XX();
        LEDS_epsilon_stress(id,1) = epsilon_t.YY();
        LEDS_epsilon_stress(id,2) = sigma_target.XX();
        LEDS_epsilon_stress(id,3) = sigma_target.YY();
    }
    
    LEDS_epsilon_stress.Print("data = ", std::cout,EMathematicaInput);
}

TPZFMatrix<REAL> Read_Duplet(int n_data, std::string file)
{
    TPZFMatrix<REAL> data(n_data,2);
    std::ifstream in(file.c_str());
    int count = 0;
    while(in)
    {
        in >> data(count,0);
        in >> data(count,1);
        
        count++;
        if (count == n_data)
        {
            break;
        }
    }
    return data;
}


