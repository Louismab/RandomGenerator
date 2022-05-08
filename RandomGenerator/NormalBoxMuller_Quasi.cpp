#include "pch.h"
#include "NormalBoxMuller_Quasi.h"
#define _USE_MATH_DEFINES
#include <math.h>

NormalBoxMuller_Quasi::NormalBoxMuller_Quasi()
{
}

NormalBoxMuller_Quasi::NormalBoxMuller_Quasi(double _mean, double _var, UniformGenerator* _gen, UniformGenerator* _gen2)
	:Normal(_mean, _var, _gen), gen2(_gen2)
{
}

double NormalBoxMuller_Quasi::Generate()
{
	if (NewSimulation)
	{
		double firstUniform = generatorPtr->Generate();
		double secondUniform = gen2->Generate();
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