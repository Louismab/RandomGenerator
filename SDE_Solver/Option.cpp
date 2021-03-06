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
	/*double sum = std::accumulate(v.begin(), v.end(), 0.0);
	double mean = sum / v.size();*/

	double mean = calculate_mean();

	/*std::vector<double> diff(v.size());
	std::transform(v.begin(), v.end(), diff.begin(), [mean](double x) { return x - mean; });
	double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
	double var = sq_sum / v.size();*/

	/*double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
	double stdev = std::sqrt(sq_sum / v.size() - mean * mean);
	double var = std::pow(stdev, 2);*/
	double somme=0;
	for (int i = 0;i < NbSim;i++)
	{
		somme += pow((v[i] - mean), 2);
	}

	double var = somme / NbSim;
	return var;
}

double Option::calculate_mean()
{
	double NbSim = v.size();
	double sum = std::accumulate(v.begin(), v.end(), 0.0);
	double mean = sum / v.size();

	return mean;
}

std::vector<double> Option::calculate_ConfidenceInterval(double alpha)
{
	std::vector<double> IC(2);
	double q;
	if (alpha == 0.99)
	{
		q = 2.576;
	}
	else if (alpha == 0.95)
	{
		q = 1.96;
	}
	else
	{
		std::cout << "Please choose alpha=0.99 or alpha=0.95" << std::endl;
	}

	double lb = calculate_mean() - q * (pow(calculate_variance() / v.size(), 0.5));
	double ub = calculate_mean() + q * (pow(calculate_variance() / v.size(), 0.5));

	IC[0] = lb;
	IC[1] = ub;

	return IC;
}