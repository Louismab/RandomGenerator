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

double bs_time_value(double fwd, double strike, double volatility, double maturity, double r)
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

        double d1 = (std::log(fwd / strike) + (r + 0.5 * pow(volatility, 2)) * maturity) / stddev;
        double d2 = d1 - stddev;
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

double bs_price(double fwd, double strike, double volatility, double maturity, double r, bool is_call)
{
    return vanilla_payoff(fwd, strike, is_call) + bs_time_value(fwd, strike, volatility, maturity, r);
}

