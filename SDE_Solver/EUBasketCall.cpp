#include "pch.h"
#include "EUBasketCall.h"
#include <numeric>
#include <iterator>

EUBasketCall::EUBasketCall()
{
}

EUBasketCall::EUBasketCall(RandomProcess* _process, double _K, std::vector<double> _r, double _T, std::vector<double> _weights)
	: Option(_process, _K, _r, _T), weights(_weights)
{

}

double EUBasketCall::ComputePrice(int NbSim,bool antithetic)
{
	double somme = 0.;
	double last_value;
	double WS_T;

	for (int n = 0; n < NbSim; ++n)
	{
		process->Simulate(0, T, T*365);
		std::vector<double> S_T = process->Get_ValueND(T,antithetic);
		
		if (antithetic)
		{
			process->Simulate_Antithetic(0, T, T * 365);
			std::vector<double> S_T = process->Get_ValueND(T, antithetic);
			//print_vector(S_T);
			int dim = S_T.size() / 2;
			std::vector<double> S_T_anti(dim);
			for (int i = 0;i < dim;i++)
			{
				S_T_anti[i] = (S_T[i] + S_T[i + dim]) / 2;
			}
			//print_vector(S_T_anti);
			WS_T = std::inner_product(std::begin(weights), std::end(weights), std::begin(S_T_anti), 0.0);
			std::cout << "WS_T: " << WS_T << std::endl;
		}
		else
		{
			process->Simulate(0, T, T * 365);
			std::vector<double> S_T = process->Get_ValueND(T, antithetic);
			//print_vector(S_T);
			WS_T = std::inner_product(std::begin(weights), std::end(weights), std::begin(S_T), 0.0);
			std::cout << "WS_T normal: " << WS_T << std::endl;
		}
		
		last_value = std::max(WS_T - K, 0.);
		somme = somme + last_value;
	}

	double price = std::exp(-r[0] * T) * (somme / NbSim); //mettre r en double
	return price;

}

void print_vector(const std::vector<double>& v)
{
	std::ostream_iterator<double> out_it(std::cout, ", ");
	std::copy(v.begin(), v.end(), out_it);
	std::cout << std::endl;
}


