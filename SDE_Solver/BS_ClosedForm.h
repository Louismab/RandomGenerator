#pragma once
#include <eigen-3.4.0/Eigen/Dense>
#include <vector>

//#ifndef CLOSED_FORM_H
//#define CLOSED_FORM_H

double bs_price_call(double S, double strike, double volatility, double maturity, double r);

double bs_price_basket_call(double K, std::vector<double> r, double T, std::vector<double> S, Eigen::MatrixXd VCV);
