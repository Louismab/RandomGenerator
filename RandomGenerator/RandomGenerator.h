#pragma once
#include <vector>

typedef unsigned long long myLong;

class RandomGenerator
{
public:
	RandomGenerator();
	virtual double Generate() = 0;
	std::vector<double> GenerateVector(int dimension);
	double Mean(myLong nbSim);
	double Variance(myLong nbSim);

	std::vector <double> lastGeneratedNumbers;
};

