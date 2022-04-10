#include "pch.h"
#include "Option.h"
#include < algorithm >
#include <numeric>

Option::Option()
{
}

Option::Option(RandomProcess* _process,  double _K, std::vector<double> _r,  double _T)
	: process(_process), K(_K), r(_r), T(_T),v(1)
{
	
}

double fp(double x, double p)
{
	return pow(x, p);
}

double Option::calculate_variance()
{
	double NbSim = v.size();
	double sum = std::accumulate(v.begin(), v.end(), 0.0);
	double mean = sum / v.size();

	//std::cout << mean << std::endl;

	std::vector<double> diff(v.size());
	std::transform(v.begin(), v.end(), diff.begin(), [mean](double x) { return x - mean; });
	double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
	double var = sq_sum / v.size();

	v.clear();
	return var/ NbSim;
}