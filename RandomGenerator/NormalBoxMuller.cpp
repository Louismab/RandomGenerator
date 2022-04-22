#include "pch.h"
#include "NormalBoxMuller.h"
#define _USE_MATH_DEFINES
#include <math.h>



NormalBoxMuller::NormalBoxMuller()
{
}

NormalBoxMuller::NormalBoxMuller(double _mean, double _var, UniformGenerator* _gen)
	:Normal(_mean, _var, _gen)
{
}

double NormalBoxMuller::Generate()
{
	if (NewSimulation)
	{
		double firstUniform = generatorPtr->Generate();
		double secondUniform = generatorPtr->Generate();
		double R = pow(-2 * log(firstUniform), 0.5);
		double theta = 2 * M_PI * secondUniform;
		double firstNormal = R * cos(theta);
		SecondNormal = R * sin(theta);
		NewSimulation = false;
		return firstNormal * pow(var, 0.5) + mean;
	}
	else
	{
		NewSimulation = true;
		return SecondNormal * pow(var, 0.5) + mean;
	}
		
}
