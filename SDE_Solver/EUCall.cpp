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

EUCall::EUCall(RandomProcess* _process, double _s, double _K, double _r, double _vol, double _T) : 
	Option(_process, _s, _K, _r, _vol, _T)
{

}

double EUCall::ComputePrice(int NbSim)
{
	double somme = 0.;
	double last_value;

	for (int n = 0; n < NbSim; ++n)
	{
		process -> Simulate(0, T, T);
		last_value = std::max(process->Get_Value(T) - K, 0.);
		somme = somme + last_value;
	}

	double price = std::exp(-r * T) * (somme / NbSim);
	return price;

}