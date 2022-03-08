#pragma once
#include "DiscreteGenerator.h"

class HeadTail : public DiscreteGenerator
{
public:
	HeadTail();
	HeadTail(UniformGenerator* _gen);
	virtual double Generate();

};

