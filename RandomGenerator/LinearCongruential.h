#pragma once
#include "PseudoGenerator.h"

class LinearCongruential : public PseudoGenerator
{
public:
	LinearCongruential();
	LinearCongruential(myLong _seed, myLong _multiplier, myLong _increment, myLong _modulus);
	virtual double Generate();

	myLong get_Modulus();


protected:
	myLong Multiplier;
	myLong Increment;
	myLong Modulus;
};

