#include "pch.h"
#include "BS_ClosedForm.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream>


double norm_cdf(double x)
{
    return 0.5 * std::erfc(-x / std::sqrt(2));
}


double bs_price_call(double S, double strike, double volatility, double maturity, double r)
{
    double stddev = volatility * std::sqrt(maturity);
    double d1 = (std::log(S / strike) + (r + 0.5 * pow(volatility, 2)) * maturity) / stddev;
    double d2 = d1 - stddev;

    double price = S* norm_cdf(d1) - strike * std::exp(-r * maturity)* norm_cdf(d2);

    return price;

}

double bs_price_basket_call(double K, std::vector<double> r, double T, std::vector<double> S, Eigen::MatrixXd VCV)
{
    double vol_1 = VCV.diagonal().mean();
    double vol_2 = (1. / pow(S.size(), 2)) * VCV.sum();
    Eigen::VectorXd S_V = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(S.data(), S.size());
    Eigen::VectorXd F_V = S_V * exp(r[0] * T);
    double fwd = pow(F_V.prod(), 1. / S.size());
    double d_plus = (log(fwd / K) + 0.5 * vol_2 * T) / pow(vol_2 * T, 0.5);
    double d_minus = d_plus - pow(vol_2 * T, 0.5);
    double price = fwd * norm_cdf(d_plus) - K * norm_cdf(d_minus);

    return price;
}


