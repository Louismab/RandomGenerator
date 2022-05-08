#pragma once
#include "PseudoGenerator.h"
#include "LinearCongruential.h"

class EcuyerCombined : public PseudoGenerator
{
public:
	EcuyerCombined();
	virtual double Generate();


private:
	LinearCongruential Generator1;
	LinearCongruential Generator2;
};

