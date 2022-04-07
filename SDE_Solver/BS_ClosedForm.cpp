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

double norm_pdf(const double x) {
    return (1.0 / (pow(2 * M_PI, 0.5))) * exp(-0.5 * x * x);
}



double vanilla_payoff(double fwd, double strike, bool is_call)
{
    return std::max(is_call ? fwd - strike : strike - fwd, 0.);
}

double bs_time_value(double fwd, double strike, double volatility, double maturity)
{
    if (strike == 0.)
    {
        return 0.;
    }
    else
    {
        double stddev = volatility * std::sqrt(maturity);
        if (stddev == 0.)
        {
            return 0.;
        }
        double tmp = std::log(fwd / strike) / stddev;
        double d1 = tmp + 0.5 * stddev;
        double d2 = tmp - 0.5 * stddev;
        double res;
        if (fwd > strike)
        {
            res = strike * norm_cdf(-d2) - fwd * norm_cdf(-d1);
        }
        else
        {
            res = fwd * norm_cdf(d1) - strike * norm_cdf(d2);
        }
        if (res <= std::numeric_limits<double>::min())
        {
            res = 0.;
        }
        return res;
    }
}



double bs_price(double fwd, double strike, double volatility, double maturity, bool is_call)
{
    return vanilla_payoff(fwd, strike, is_call) + bs_time_value(fwd, strike, volatility, maturity);
}

double bs_price_call(double S, double strike, double volatility, double maturity, double r)
{
    double stddev = volatility * std::sqrt(maturity);
    double d1 = (std::log(S / strike) + (r + 0.5 * pow(volatility, 2)) * maturity) / stddev;
    double d2 = d1 - stddev;

    double price = S* norm_cdf(d1) - strike * std::exp(-r * maturity)* norm_cdf(d2);

    return price;

}


