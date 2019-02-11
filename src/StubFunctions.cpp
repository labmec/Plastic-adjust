//
// Created by gus on 04/02/19.
//

#include "StubFunctions.h"

void open(std::string filePath, std::string nickname) {
    std::cout << "open(" << filePath << ", " << nickname << ");" << std::endl;
}

void model(std::string modelName, std::string test) {
    std::cout << "model(" << modelName << ", " << test << ");" << std::endl;
}

void select(std::string parameter, int initial, int final, std::string comment) {
    std::cout << "select(" << parameter << ", " << initial << ", " << final << ", " << comment << ");" << std::endl;
}

void adjust(std::string parameter) {
    std::cout << "select(" << parameter << ");" << std::endl;
}

void assign(std::string parameter, std::string fileNickname) {
    std::cout << "assign(" << parameter << ", " << fileNickname << ");" << std::endl;
}

void clear() {
    std::cout << "clear();" << std::endl;
}
