#pragma once
#include<string>
#include"ephemeris.h"

uint64_t random64();
double randomreal();
double randomnormal();
vec randomdirection();
mat randommatrix();

double CalcTime();

std::string readline(FILE *fin);