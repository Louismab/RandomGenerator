#pragma once
#include "RandomGenerator.h"
#include "RandomProcess.h"

class Option
{
	public:
		Option();
		Option(RandomProcess* _process, double _K, std::vector<double> _r, double _T);
		virtual double ComputePrice(int NbSim, bool antithetic = false) = 0;
		virtual double ComputePrice_ControlVariate(int NbSim) = 0;
		virtual double ComputePrice_VDC(int NbSim) = 0;
		double calculate_variance();
		double calculate_mean();
		std::vector<double> calculate_ConfidenceInterval(double alpha = 0.99);

	protected:
		//double s;
		double K;
		std::vector<double> r;
		//double vol;
		double T;
		RandomProcess* process;

		std::vector<double> v; //stock the payoffs to calculte the variance



};

double fp(double x, double p);
