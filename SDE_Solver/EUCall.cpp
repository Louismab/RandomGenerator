#include "pch.h"
#include "EUCall.h"
#include "BSEULER1D.h"
#include "Milstein1D.h"
#include "SinglePath.h"
#include <iostream>     
#include <algorithm> 
#include <cmath>

EUCall::EUCall()
{
}

EUCall::EUCall(RandomProcess* _process,  double _K, std::vector<double> _r,  double _T)
	: Option(_process, _K, _r, _T)
{

}

double EUCall::ComputePrice(int NbSim, bool antithetic)
{
	double somme = 0.;
	double last_value;

	for (int n = 0; n < NbSim; ++n)
	{
		process->Simulate(0, T, T*365);
		last_value = std::max(process->Get_Value(T) - K, 0.);
		somme = somme + last_value;
	}

	double price = std::exp(-r[0] * T) * (somme / NbSim);
	return price;

}

double ComputePrice_ControlVariate(int NbSim)
{
	return 0;
}