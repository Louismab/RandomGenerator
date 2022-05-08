#pragma once
#include "Normal.h"

class NormalBoxMuller_Quasi : public Normal
{
public:
	NormalBoxMuller_Quasi();
	NormalBoxMuller_Quasi(double _mean, double _var, UniformGenerator* _gen, UniformGenerator* _gen2);
	virtual double Generate();


private:
	bool NewSimulation = true;
	double SecondNormal;
	UniformGenerator* gen2;

};

