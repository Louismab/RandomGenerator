#pragma once
#include "Exponential.h"

class ExponentialID : public Exponential
{

public:
	ExponentialID();
	ExponentialID(double _lambda, UniformGenerator* _gen);
	virtual double Generate();

};

