#pragma once

//#ifndef CLOSED_FORM_H
//#define CLOSED_FORM_H

double vanilla_payoff(double fwd, double strike, bool is_call);
double bs_time_value(double fwd, double strike, double volatility, double maturity, double r);
double bs_price(double fwd, double strike, double volatility, double maturity, double r, bool is_call);

