//
// Created by gus on 04/02/19.
//

#ifndef STUBFUNCTIONS_H
#define STUBFUNCTIONS_H

#include <iostream>
#include <string>


void open(std::string filePath, std::string nickname);

void model(std::string modelName, std::string test);

void select(std::string parameter, int initial, int final, std::string comment);

void adjust(std::string parameter);

void assign(std::string parameter, std::string fileNickname);

void clear();

#endif
