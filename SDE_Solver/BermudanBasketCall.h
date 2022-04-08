#pragma once
#include "Option.h"
#include <eigen-3.4.0/Eigen/Dense>

class BermudanBasketCall : public Option
{
public:
	BermudanBasketCall();
	BermudanBasketCall(RandomProcess* _process, double _K, std::vector<double> _r, double _T, std::vector<double> _weights, std::vector<double> _exeDates,  int _L = 5);

	double ComputePrice(int NbSim, bool antithetic = false);
	double ComputePrice_ControlVariate(int NbSim);

private:
	std::vector<double> exeDates;
	int L;
	std::vector<double> weights;
};



