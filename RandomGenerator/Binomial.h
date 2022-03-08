#pragma once
#include "DiscreteGenerator.h"

class Binomial : public DiscreteGenerator
{
public:
	Binomial();
	Binomial(UniformGenerator* _gen, myLong _n,double _p);
	virtual double Generate();

private:
	double p;
	myLong n;
};

