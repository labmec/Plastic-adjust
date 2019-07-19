

// main.cpp : Defines the entry point for the console application.

#include <iostream>
#include <fstream>
#include <iomanip>

#include "json.hpp"
#include "StubFunctions.h"

/// Plasticity
#include "TPZPlasticStepPV.h"
#include "TPZSandlerExtended.h"
#include "TPZYCMohrCoulombPV.h"

/// Adjustment classes
#include "TF1DSAdjust.h"
#include "TMohrAdjust.h"
#include "TF2DSAdjust_DW.h"
#include "TF2DSAdjust_R.h"


using json = nlohmann::json;
void printJSON(json commandFile, std::ostream& out);
void callMethods(json commandFile);
void translateToFunction(json singleCommand);


/// Method to present elastoplastic model
void LoadElastoPlasticModel();

/// Method to present Mohr-Coulomb model
void LoadMCModel();

// Read experimental data
TPZFMatrix<REAL> Read_Duplet(int n_data, std::string file);



int main() {
    
    
//    TF2DSAdjust_DW F2;
//    F2.PopulateDW();
//    F2.AdjustDW2();
//    
//    return 0;
//
    
//    TF2DSAdjust_R F2;
//    F2.PopulateR();
//    F2.AdjustR();
//    return 0;
    
    LoadElastoPlasticModel();
    return 0;
    
//    LoadMCModel();
//    return 0;
    
//    TF1DSAdjust F1;
//    F1.Populate();
//    F1.Adjust();
//    return 0;
    
    
//    TMohrAdjust MC;
//    MC.Populate();
//    MC.Adjust();
//    return 0;
    
    
    std::string path;
    std::ifstream input;

    path = "../input_CapF2DSfindRcp14.json";
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

void LoadMCModel()
{
    
    // Experimental data
    std::string file_name;
    file_name = "/Users/manouchehr/Documents/GitHub/Plastic-adjust/exp_data/CP10.txt";
    int64_t n_data = 701;
    TPZFMatrix<REAL> data = Read_Duplet(n_data, file_name);
    //     data.Print(std::cout);
    
    // MC Mohr Coloumb PV
    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
    
    // LE Linear elastic response
    TPZElasticResponse ER;
    
    /// number CP10
//    REAL E  = 4110.72; // MPa
//    REAL nu = 0.295099;
    
    REAL E  = 2693.74; // MPa
    REAL nu = 0.105264;
    
    
    ER.SetEngineeringData(E, nu);
    
    // Mohr Coulomb data; Nlopt
    REAL mc_cohesion    = 3.47113;
    REAL mc_phi         = 19.1662*M_PI/180;
    REAL mc_psi         = mc_phi; // because MS do not understand
    
    LEMC.SetElasticResponse(ER);
    LEMC.fYC.SetUp(mc_phi, mc_psi, mc_cohesion, ER);
    
    TPZTensor<REAL> epsilon_t,sigma,sigma_target;
    sigma.Zero();
    epsilon_t.Zero();
    
    epsilon_t.Zero();
    
    TPZFNMatrix<701,STATE> LEMC_epsilon_stress(n_data,4);
    
    for (int64_t id = 0; id < n_data; id++) {
        
        // For a given strain
        epsilon_t.XX() = -0.01*data(id,0);
        epsilon_t.YY() = -0.01*data(id,1);
        epsilon_t.ZZ() = -0.01*data(id,1);
        
        LEMC.ApplyStrainComputeSigma(epsilon_t, sigma_target);
        
        
        LEMC_epsilon_stress(id,0) = epsilon_t.XX();
        LEMC_epsilon_stress(id,1) = epsilon_t.YY();
        LEMC_epsilon_stress(id,2) = sigma_target.XX();
        LEMC_epsilon_stress(id,3) = sigma_target.YY();
    }
    
    LEMC_epsilon_stress.Print("data = ", std::cout,EMathematicaInput);
}

void LoadElastoPlasticModel()
{
    
    // Experimental data
    std::string file_name;
    file_name = "/Users/manouchehr/Documents/GitHub/Plastic-adjust/exp_data/CP14.txt";
    int64_t n_data = 4619;
    TPZFMatrix<REAL> data = Read_Duplet(n_data, file_name);
//     data.Print(std::cout);
    
    // DS Dimaggio Sandler PV
    TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse> LEDS;
    
    // LE Linear elastic response
    TPZElasticResponse ER;
    
    /// number 
    REAL E  = 2000.0; // MPa
    REAL nu = 0.25; // MPa
    
    STATE G = E / (2. * (1. + nu));
    STATE K = E / (3. * (1. - 2 * nu));
    REAL CA      = 16;
    REAL CB      = 0.04;
    REAL CC      = 14;
    REAL CD      = 0.02;
    REAL CR      = 3.0;
    REAL CW      = 0.2;
    REAL X_0     = -70.0;
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


