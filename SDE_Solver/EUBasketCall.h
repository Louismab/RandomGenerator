#pragma once
#include "Option.h"

class EUBasketCall :public Option
{
public:
	EUBasketCall();
	EUBasketCall(RandomProcess* _process, double _K, std::vector<double> _r, double _T, std::vector<double> _weights);

	double ComputePrice(int NbSim, bool antithetic = false);
	//double ComputePrice_Antithetic(int NbSim);

private:
	std::vector<double> weights;
};

void print_vector(const std::vector<double>& v);