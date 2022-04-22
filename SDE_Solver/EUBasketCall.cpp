#include "pch.h"
#include "EUBasketCall.h"
#include <numeric>

EUBasketCall::EUBasketCall()
{
}

EUBasketCall::EUBasketCall(RandomProcess* _process, double _K, std::vector<double> _r, double _T, std::vector<double> _weights)
	: Option(_process, _K, _r, _T), weights(_weights)
{

}

double EUBasketCall::ComputePrice(int NbSim)
{
	double somme = 0.;
	double last_value;

	for (int n = 0; n < NbSim; ++n)
	{
		process->Simulate(0, T, T*365);
		std::vector<double> S_T = process->Get_ValueND(T);
		double WS_T = std::inner_product(std::begin(weights), std::end(weights), std::begin(S_T), 0.0);
		last_value = std::max(WS_T - K, 0.);
		somme = somme + last_value;
	}

	double price = std::exp(-r[0] * T) * (somme / NbSim); //mettre r en double
	return price;

}
