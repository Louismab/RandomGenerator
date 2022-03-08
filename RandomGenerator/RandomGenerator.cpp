#include "pch.h"
#include "RandomGenerator.h"
#include <math.h>

RandomGenerator::RandomGenerator()
{
}

double RandomGenerator::Mean(myLong nbSim)
{
	double result = 0. ;
	lastGeneratedNumbers = std::vector<double>(nbSim);

	for (myLong i = 0;i < nbSim;i++)
	{
		double currentNumber = Generate();
		result += currentNumber / nbSim;
		lastGeneratedNumbers[i]=currentNumber;   //we could use push_back with an empty vector at the beginning but not opti
	}

	return result;
}

double RandomGenerator::Variance(myLong nbSim)
{
	double result = 0.;
	double mean = Mean(nbSim);

	for (myLong i = 0;i < nbSim;i++)
	{
		double currentNumber = lastGeneratedNumbers[i];
		result += pow(currentNumber -mean,2)/nbSim;
	}

	return result;
}

std::vector<double> RandomGenerator::GenerateVector(int dimension)
{
	std::vector<double> result(dimension);
	for (int i = 0;i < dimension;i++)
	{
		result[i] = Generate();
	}
	return result;
}