#pragma once
#include "Normal.h"

class NormalInverseCdf :public Normal
{
public:
	NormalInverseCdf();
	NormalInverseCdf(double _mean, double _var, UniformGenerator* _gen);
	virtual double Generate();

};

double RationalApproximation(double t);
double NormalCDFInverse(double p);
