#pragma once
#include "RandomGenerator.h"
#include "RandomProcess.h"

class Option
{
	public:
		Option();
		Option(RandomProcess* _process, double _K, std::vector<double> _r, double _T);
		virtual double ComputePrice(int NbSim, bool antithetic = false) = 0;

	protected:
		//double s;
		double K;
		std::vector<double> r;
		//double vol;
		double T;
		RandomProcess* process;



};

