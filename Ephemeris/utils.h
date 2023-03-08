#pragma once
#include<string>
#include"ephemeris.h"

double randomreal();
vec randomdirection();
mat randommatrix();

double CalcTime();

std::string readline(FILE *fin);