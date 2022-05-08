// Testeur.cpp : Ce fichier contient la fonction 'main'. L'exécution du programme commence et se termine à cet endroit.
//

#include <iostream>
#include "LinearCongruential.h"
#include "EcuyerCombined.h"
#include "HeadTail.h"
#include "Bernoulli.h"
#include "Binomial.h"
#include "FiniteSet.h"
#include "NormalBoxMuller.h"
#include "NormalCLT.h"
#include <vector>
#include "PoissonAlgo1.h"
#include "PoissonAlgo2.h"
#include "ExponentialID.h"
#include "BSEULER1D.h"
#include "Milstein1D.h"
#include "EUCall.h"
#include <eigen-3.4.0/Eigen/Dense>
#include "BSEULERND.h"
#include <numeric>
#include "EUBasketCall.h"
#include "VanDerCorput.h"
#include <math.h>
#include <cmath>
#include "BS_ClosedForm.h"
#include "NormalInverseCdf.h"
#include "MilsteinND.h"
#include "BermudanCall.h"
#include "BermudanBasketCall.h"
#include <fstream>
#include <exception>
#include "NormalBoxMuller_Quasi.h"
#include <algorithm>
#include <random>

// Function in Order to Transfer our vectors to csv
void exportVectortoXl(std::string file, std::vector<double> V)
{
    std::ofstream myfile;
    myfile.open(file, std::ofstream::app);

    myLong vsize = V.size();
    for (myLong n = 0; n < vsize; n++)
    {
        myfile << V[n] << "\n";

    }
    myfile.close();
}
//*********************************************************************************************************
//************************************************MAIN*****************************************************

int main()
{
    RandomGenerator* Generator;
    LinearCongruential* Uniform = new LinearCongruential(27, 17, 43, 100);
    EcuyerCombined* Uniform2 = new EcuyerCombined();
    VanDerCorput* Vdc = new VanDerCorput();
    VanDerCorput* Vdc_3 = new VanDerCorput(3);
    NormalBoxMuller* Normal = new NormalBoxMuller(0, 1, Uniform2);
    NormalCLT* Normal_clt = new NormalCLT(0, 1, Uniform2);
    NormalInverseCdf* Normal_VDC = new NormalInverseCdf(0, 1, Vdc);
    NormalInverseCdf* Normal_Inv = new NormalInverseCdf(0, 1, Uniform2);
    NormalBoxMuller_Quasi* Normal_VDC_3 = new NormalBoxMuller_Quasi(0, 1, Vdc, Vdc_3);

    std::vector<double> v = { 1,2,3,4 };
    std::random_device rd;
    unsigned seed = 0;
    auto rng = std::default_random_engine{ seed };
    std::vector<double> x(v);
    std::shuffle(x.begin(), x.end(), rng);
    for (int i = 0;i < 4;i++)
    {
        std::cout << x[i] << std::endl;
    }
    for (int i = 0;i < 4;i++)
    {
        std::cout << v[i] << std::endl;
    }
    std::vector<double> y(v);
    rng = std::default_random_engine{ seed };
    std::shuffle(y.begin(), y.end(), rng);
    for (int i = 0;i < 4;i++)
    {
        std::cout << y[i] << std::endl;
    }

    /*std::cout << "Test Van der Corput" << std::endl;

    std::cout << "Mean= " << Vdc->Mean(10000) << std::endl;
    std::cout << "Variance= " << Vdc->Variance(10000) << std::endl;

    for (int i = 0;i < 10;++i)
    {
       std::cout << "Vdc base 67 : " << Vdc->Generate() << std::endl;
       //std::cout << "Vdc base 3 : " << Vdc_3->Generate() << std::endl;
    }*/


    std::cout << "Mean= " << Normal_VDC_3->Mean(10000) << std::endl;
    std::cout << "Variance= " << Normal_VDC_3->Variance(10000) << std::endl;

    /*std::vector < double> v1;
    for (int i = 0;i < 10000;++i)
    {
        v1.push_back(Normal_VDC_3->Generate());
    }
    exportVectortoXl("Normal_VDC_3.csv", v1);
    std::cout << "Mean= " << Normal_VDC_3->Mean(10000) << std::endl;
    std::cout << "Variance= " << Normal_VDC_3->Variance(10000) << std::endl;*/

    //1 dimension
    double s = 100.;
    double vol = 0.2;
    double maturity = 30./365.;
    double K = 100;
    double nb_simul = 10000;

    // Single Underlying
    double s = 100.;
    double vol = 0.2;
    std::vector<double> r = { 0.000 };

    // Basket Underlying
    std::vector<double> S = { 50.0, 150.0, 100.0 };
    std::vector<double> R = { 0.0, 0.0, 0.0 };
    std::vector<double> W = { 1. / 3, 1. / 3, 1. / 3 };
    Eigen::MatrixXd VCV{ {0.04,0.036,0.032},{0.036,0.04,0.036},{0.032,0.036,0.04} };
    //std::vector<double> c = { 1. / 4., 1. / 2.,3. / 4.,1. };
    std::vector<double> c = { 5. / 365., 15. / 365.,30. / 365. };

    //Display Parameters in the Command window
    std::cout << "##################  PARAMETERS: ################## \n" << std::endl;
    std::cout << " Maturity in year: T = " << maturity << "\n" << std::endl;
    std::cout << " Strike : K = " << K << "\n" << std::endl;
    std::cout << " Simulation number : nb_sim =" << nb_simul << "\n" << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << " ******* Single Underlying ******* \n" << std::endl;
    std::cout << " Single Underlying Spot at t=0 : S = " << s << "\n" << std::endl;
    std::cout << " Single Underlying Volatility at t=0 : sigma =" << vol << "\n" << std::endl;
    std::cout << " risk free rate : r =" << std::endl;
    print_vector(r);
    std::cout << "\n" << std::endl;
    std::cout << "\n" << std::endl;

    /*BermudanBasketCall* BermudanBasket_Euler = new BermudanBasketCall(EulerND, K, R, maturity, W, c, S, VCV);
    std::cout << "Price Bermudan Call Euler: \n" << std::endl;
    std::cout << "Price regular :" << BermudanBasket_Euler->ComputePrice(10000) << std::endl;
    std::cout << " (variance : " << BermudanBasket_Euler->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BermudanBasket_Euler->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BermudanBasket_Euler->calculate_ConfidenceInterval()[1] << std::endl;

    std::cout << "Price antithetic : " << BermudanBasket_Euler->ComputePrice(10000, true) << std::endl;
    std::cout << " (variance : " << BermudanBasket_Euler->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BermudanBasket_Euler->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BermudanBasket_Euler->calculate_ConfidenceInterval()[1] << std::endl;

    std::cout << "Price Control Variate: " << BermudanBasket_Euler->ComputePrice_ControlVariate(10000) << std::endl;
    std::cout << " (variance : " << BermudanBasket_Euler->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BermudanBasket_Euler->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BermudanBasket_Euler->calculate_ConfidenceInterval()[1] << std::endl;

    std::cout << "\n " << std::endl;

    /*BermudanBasketCall* BermudanBasket_Milstein = new BermudanBasketCall(MilsteinPD, K, R, maturity, W, c, S, VCV);
    std::cout << "Price Bermudan Call Milstein: \n" << std::endl;
    std::cout << "Price regular : " << BermudanBasket_Milstein->ComputePrice(10000) << std::endl;
    std::cout << " (variance : " << BermudanBasket_Milstein->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BermudanBasket_Milstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BermudanBasket_Milstein->calculate_ConfidenceInterval()[1] << std::endl;

    std::cout << "Price antithetic : " << BermudanBasket_Milstein->ComputePrice(10000, true) << std::endl;
    std::cout << " (variance : " << BermudanBasket_Milstein->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BermudanBasket_Milstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BermudanBasket_Milstein->calculate_ConfidenceInterval()[1] << std::endl;

    std::cout << "Price control variate : " << BermudanBasket_Milstein->ComputePrice_ControlVariate(10000) << std::endl;
    std::cout << " (variance : " << BermudanBasket_Milstein->calculate_variance() << ")" << std::endl << std::endl;
    std::cout << "lower bound : " << BermudanBasket_Milstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BermudanBasket_Milstein->calculate_ConfidenceInterval()[1] << std::endl;*/

    //graph Bermudan Basket call Euler

    /*std::vector < double> v_bermudan_basket_call_euler;
    std::vector < double> v_bermudan_basket_call_euler_anti;
    std::vector < double> v_bermudan_basket_call_euler_PCV;


    for (myLong i = 1; i <= 10000 / 10; i++)
    {
        std::cout << i << std::endl;
        v_bermudan_basket_call_euler.push_back(BermudanBasket_Euler->ComputePrice(i * 10));
        v_bermudan_basket_call_euler_anti.push_back(BermudanBasket_Euler->ComputePrice(i * 10, true));
        v_bermudan_basket_call_euler_PCV.push_back(BermudanBasket_Euler->ComputePrice_ControlVariate(i * 10));
    }
    exportVectortoXl("bermudan_euler_one_month_no_correl.csv", v_bermudan_basket_call_euler);
    exportVectortoXl("bermudan_euler_anti_one_month_no_correl.csv", v_bermudan_basket_call_euler_anti);
    exportVectortoXl("bermudan_euler_pcv_one_month_no_correl.csv", v_bermudan_basket_call_euler_PCV);

    //graph Bermudan call Milstein

    /*std::vector < double> v_bermudan_basket_call_milstein;
    std::vector < double> v_bermudan_basket_call_milstein_anti;
    std::vector < double> v_bermudan_basket_call_milstein_PCV;


    for (myLong i = 1; i <= 10000 / 10; i++)
    {
        std::cout << i << std::endl;
        v_bermudan_basket_call_milstein.push_back(BermudanBasket_Milstein->ComputePrice(i * 10));
        v_bermudan_basket_call_milstein_anti.push_back(BermudanBasket_Milstein->ComputePrice(i * 10, true));
        v_bermudan_basket_call_milstein_PCV.push_back(BermudanBasket_Milstein->ComputePrice_ControlVariate(i * 10));
    }
    exportVectortoXl("bermudan_milstein_one_month_no_correl.csv", v_bermudan_basket_call_milstein);
    exportVectortoXl("bermudan_milstein_anti_one_month_no_correl.csv", v_bermudan_basket_call_milstein_anti);
    exportVectortoXl("bermudan_milstein_pcv_one_month_no_correl.csv", v_bermudan_basket_call_milstein_PCV);*/

    /*std::vector < double> v1;
    for (int i = 0;i < 10000;++i)
    {
        v1.push_back(Vdc->Generate());
    }
    exportVectortoXl("Vdc.csv", v1);
    std::cout << "Mean= " << Vdc->Mean(10000) << std::endl;
    std::cout << "Variance= " << Vdc->Variance(10000) << std::endl;*/

    

    double price_cf = bs_price_call(s, K, vol, maturity, r[0]);
    std::cout << "Price BS closed form = " << price_cf << std::endl;

    std::cout << " Basket weigts by Underlying : " << "\n";
    print_vector(W);
    std::cout << "\n" << std::endl;

    std::cout << " Basket Underlying covariance matrix at t=0 : \n" << VCV << "\n" << std::endl;

    std::cout << "Price Call with Euler \n " << std::endl;
    EUCall* CallEuler = new EUCall(Euler, K, r, maturity);
    std::cout << "Price Call Euler: " << CallEuler->ComputePrice(10000) << std::endl;
    std::cout << " (variance : " << CallEuler->calculate_variance() << ")" << std::endl;
    std::cout << "Price Call Euler antithetic: " << CallEuler->ComputePrice(10000,true) << std::endl;
    std::cout << " (variance : " << CallEuler->calculate_variance() << ")" << std::endl;
    std::cout << "Price Call Euler Control Variate: " << CallEuler->ComputePrice_ControlVariate(1000) << std::endl;
    std::cout << " (variance : " << CallEuler->calculate_variance() << ")" << std::endl;
    EUCall* CallEuler_VDC = new EUCall(Euler_VDC, 100, r, maturity);
    std::cout << "Price Call Euler VDC: " << CallEuler_VDC->ComputePrice_VDC(5000) << std::endl;
    std::cout << " (variance : " << CallEuler_VDC->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << CallEuler_VDC->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << CallEuler_VDC->calculate_ConfidenceInterval()[1] << std::endl;

    /*std::vector<double> IC = CallEuler->calculate_ConfidenceInterval();
    std::cout << "lower bound : " << IC[0];
    std::cout << "upper bound : " << IC[1];*/

    std::cout << "Price Call with Milstein : \n " << std::endl;
    Milstein->Simulate(0, maturity, 30);
    EUCall* CallMilstein = new EUCall(Milstein, K, r, maturity);
    std::cout << "Price Call Miltsein : " << CallMilstein->ComputePrice(10000) << std::endl;
    std::cout << " (variance : " << CallMilstein->calculate_variance() << ")" << std::endl;
    std::cout << "Price Call Miltsein antithetic : " << CallMilstein->ComputePrice(10000,true) << std::endl;
    std::cout << " (variance : " << CallMilstein->calculate_variance() << ")" << std::endl;
    std::cout << "Price Call Miltsein Control Variate: " << CallMilstein->ComputePrice_ControlVariate(1000) << std::endl;
    std::cout << " (variance : " << CallMilstein->calculate_variance() << ")" << std::endl;
    EUCall* CallMilstein_VDC = new EUCall(Milstein_VDC, 100, r, maturity);
    std::cout << "Price Call Milstein VDC: " << CallMilstein_VDC->ComputePrice_VDC(5000) << std::endl;
    std::cout << " (variance : " << CallMilstein_VDC->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << CallMilstein_VDC->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << CallMilstein_VDC->calculate_ConfidenceInterval()[1] << std::endl;


    /*Euler_VDC->Simulate(0, 1, 365);
    double FinalValue= Euler_VDC->Get_Value(1);
    std::cout << "Final Value" << FinalValue << std::endl;
    for (double i = 0;i < 365;i++)
    {
        std::cout << i << ": " << Euler_VDC->Get_Value(1./365*i) << std::endl;
    }*/

    
    std::cout << "Price Call with Euler Van Der Corput \n " << std::endl;
    /*for (int i = 0;i < 10;i++)
    {
        Euler_VDC->Simulate(0, 30. / 365., 30);
        double last_value=Euler_VDC->Get_Value(30. / 365.);
        std::cout << "last value :" << last_value << std::endl;
    }*/

    

    //EUBasketCall* BasketCallMilstein = new EUBasketCall(MilsteinPD_VDC, K, R, maturity, W, S, VCV);


    std::cout << "\n " << std::endl;

    

    std::cout << "##################  COMPUTATION of Prices, Variance and IC  ################## \n" << std::endl;

    /*std::cout << "Price Call with Milstein : Van Der Corput \n " << std::endl;
    EUCall* CallMilstein_VDC = new EUCall(Milstein_VDC, K, r, maturity);
    std::cout << "Price Call Miltsein VDC : " << CallMilstein_VDC->ComputePrice(10000) << std::endl;
    std::cout << " (variance : " << CallMilstein->calculate_variance() << ")" << std::endl;
    std::cout << "Price Call Miltsein antithetic VDC : " << CallMilstein_VDC->ComputePrice(10000, true) << std::endl;
    std::cout << " (variance : " << CallMilstein->calculate_variance() << ")" << std::endl;
    std::cout << "Price Call Miltsein Control Variate VDC : " << CallMilstein_VDC->ComputePrice_ControlVariate(1000) << std::endl;
    std::cout << " (variance : " << CallMilstein->calculate_variance() << ")" << std::endl;*/

    // ---------------------------------------------------------------------*/

    EUBasketCall* BasketCallEuler = new EUBasketCall(EulerND, K, R, maturity, W, S, VCV);
    std::cout << "******** Price European Basket Call Euler:  ******** \n" << std::endl;
    std::cout << "Price regular :" << BasketCallEuler->ComputePrice(nb_simul) << std::endl;
    std::cout << " (variance : " << BasketCallEuler->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BasketCallEuler->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BasketCallEuler->calculate_ConfidenceInterval()[1] << std::endl;

    std::cout << "Price antithetic : " << BasketCallEuler->ComputePrice(nb_simul, true) << std::endl;
    std::cout << " (variance : " << BasketCallEuler->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BasketCallEuler->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BasketCallEuler->calculate_ConfidenceInterval()[1] << std::endl;

    std::cout << "Price Control Variate: " << BasketCallEuler->ComputePrice_ControlVariate(nb_simul) << std::endl;
    std::cout << " (variance : " << BasketCallEuler->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BasketCallEuler->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BasketCallEuler->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "Price Antithetic : " << BasketCallEuler->ComputePrice(10000, true) << std::endl;
    std::cout << " (variance : " << BasketCallEuler->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BasketCallEuler->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BasketCallEuler->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "Price Control Variate : " << BasketCallEuler->ComputePrice_ControlVariate(10000) << std::endl;
    std::cout << " (variance : " << BasketCallEuler->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BasketCallEuler->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BasketCallEuler->calculate_ConfidenceInterval()[1] << std::endl;
    EUBasketCall* BasketCallEuler_VDC = new EUBasketCall(EulerND_VDC, K, R, maturity, W, S, VCV);
    std::cout << "Price Basket Call Euler VDC: " << BasketCallEuler_VDC->ComputePrice_VDC(5000) << std::endl;
    std::cout << " (variance : " << BasketCallEuler_VDC->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BasketCallEuler_VDC->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BasketCallEuler_VDC->calculate_ConfidenceInterval()[1] << std::endl;

    std::cout << "\n " << std::endl;
    std::cout << "\n " << std::endl;

    EUBasketCall* BasketCallMilstein = new EUBasketCall(MilsteinPD, K, R, maturity, W, S, VCV);

    std::cout << "********   Price European Basket Call Milstein: ********  \n" << std::endl;
    std::cout << "Price regular : " << BasketCallMilstein->ComputePrice(nb_simul) << std::endl;
    std::cout << " (variance : " << BasketCallMilstein->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BasketCallMilstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BasketCallMilstein->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "Price Antithetic : " << BasketCallMilstein->ComputePrice(10000, true) << std::endl;
    std::cout << " (variance : " << BasketCallMilstein->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BasketCallMilstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BasketCallMilstein->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "Test ComputePrice_ControlVariate : " << BasketCallMilstein->ComputePrice_ControlVariate(10000) << std::endl;
    std::cout << " (variance : " << BasketCallMilstein->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BasketCallMilstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BasketCallMilstein->calculate_ConfidenceInterval()[1] << std::endl;
    EUBasketCall* BasketCallMilstein_VDC = new EUBasketCall(MilsteinPD_VDC, K, R, maturity, W, S, VCV);
    std::cout << "Price Basket Call Milstein VDC: " << BasketCallMilstein_VDC->ComputePrice_VDC(5000) << std::endl;
    std::cout << " (variance : " << BasketCallMilstein_VDC->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BasketCallMilstein_VDC->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BasketCallMilstein_VDC->calculate_ConfidenceInterval()[1] << std::endl;

    std::cout << "Price control variate : " << BasketCallMilstein->ComputePrice_ControlVariate(nb_simul) << std::endl;
    std::cout << " (variance : " << BasketCallMilstein->calculate_variance() << ")" << std::endl << std::endl;
    std::cout << "lower bound : " << BasketCallMilstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BasketCallMilstein->calculate_ConfidenceInterval()[1] << std::endl;

    //graph European Basket call Euler
    /*
    std::vector < double> v_euro_basket_call_euler;
    std::vector < double> v_euro_basket_call_euler_anti;
    std::vector < double> v_euro_basket_call_euler_PCV;


    std::cout << "\n " << std::endl;

    std::cout << "Bermudan Call" << std::endl;
    std::cout << "\n " << std::endl;
    
    
    BermudanCall* Bermudan_Euler = new BermudanCall(Euler, K, r, maturity, c);
    std::cout << "Price Bermudan Call Euler: \n" << std::endl;
    std::cout << "Price regular :" << Bermudan_Euler->ComputePrice(10000) << std::endl;
    std::cout << " (variance : " << Bermudan_Euler->calculate_variance() << ")" << std::endl;
    std::cout << "Price antithetic : " << Bermudan_Euler->ComputePrice(10000,true) << std::endl;
    std::cout << " (variance : " << Bermudan_Euler->calculate_variance() << ")" << std::endl;
    std::cout << "Price Control Variate: " << Bermudan_Euler->ComputePrice_ControlVariate(10000) << std::endl;
    std::cout << " (variance : " << Bermudan_Euler->calculate_variance() << ")" << std::endl;
    BermudanCall* Bermudan_EulerVDC = new BermudanCall(Euler_VDC, K, r, maturity, c);
    std::cout << "Price Bermudan Call Euler VDC: " << Bermudan_EulerVDC->ComputePrice_VDC(5000) << std::endl;
    std::cout << " (variance : " << Bermudan_EulerVDC->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << Bermudan_EulerVDC->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << Bermudan_EulerVDC->calculate_ConfidenceInterval()[1] << std::endl;

    std::cout << "\n " << std::endl;

    BermudanCall* Bermudan_Milstein = new BermudanCall(Milstein, K, r, maturity, c);
    std::cout << "Price Bermudan Call Milstein: \n" << std::endl;
    std::cout << "Price regular : " << Bermudan_Milstein->ComputePrice(10000) << std::endl;
    std::cout << " (variance : " << Bermudan_Milstein->calculate_variance() << ")" << std::endl;
    std::cout << "Price antithetic : " << Bermudan_Milstein->ComputePrice(10000,true) << std::endl;
    std::cout << " (variance : " << Bermudan_Milstein->calculate_variance() << ")" << std::endl;
    std::cout << "Price control variate : " << Bermudan_Milstein->ComputePrice_ControlVariate(10000) << std::endl;
    std::cout << " (variance : " << Bermudan_Milstein->calculate_variance() << ")" << std::endl;
    BermudanCall* Bermudan_MilsteinVDC = new BermudanCall(Milstein_VDC, K, r, maturity, c);
    std::cout << "Price Bermudan Call Milstein VDC: " << Bermudan_MilsteinVDC->ComputePrice_VDC(5000) << std::endl;
    std::cout << " (variance : " << Bermudan_MilsteinVDC->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << Bermudan_MilsteinVDC->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << Bermudan_MilsteinVDC->calculate_ConfidenceInterval()[1] << std::endl;

    std::cout << "\n " << std::endl;

    std::cout << "Bermudan Basket Call \n" << std::endl;

    BermudanBasketCall* BermudanBasket_Euler = new BermudanBasketCall(EulerND, K, R, maturity, W, c, S, VCV);
    std::cout << "Price Bermudan Basket Call Euler : \n " << std::endl;
    std::cout << "Price regular : " << BermudanBasket_Euler->ComputePrice(10000) << std::endl;
    std::cout << " (variance : " << BermudanBasket_Euler->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BermudanBasket_Euler->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BermudanBasket_Euler->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "Price antithetic : " << BermudanBasket_Euler->ComputePrice(10000,true) << std::endl;
    std::cout << " (variance : " << BermudanBasket_Euler->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BermudanBasket_Euler->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BermudanBasket_Euler->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "Price control variate : " << BermudanBasket_Euler->ComputePrice_ControlVariate(10000) << std::endl;
    std::cout << " (variance : " << BermudanBasket_Euler->calculate_variance() << ")" << std::endl;
    BermudanBasketCall* BermudanBasket_EulerVDC = new BermudanBasketCall(EulerND_VDC, K, R, maturity, W, c, S, VCV);
    std::cout << "Price Bermudan Basket Call Euler VDC: " << BermudanBasket_EulerVDC->ComputePrice_VDC(5000) << std::endl;
    std::cout << " (variance : " << BermudanBasket_EulerVDC->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BermudanBasket_EulerVDC->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BermudanBasket_EulerVDC->calculate_ConfidenceInterval()[1] << std::endl;
    
    std::cout << "\n " << std::endl;
    
    BermudanBasketCall* BermudanBasket_Milstein = new BermudanBasketCall(MilsteinPD, K, R, maturity, W, c, S, VCV);
    std::cout << "Price Bermudan Basket Call Milstein : \n" << std::endl;
    std::cout << "Price regular : " << BermudanBasket_Milstein->ComputePrice(10000) << std::endl;
    std::cout << " (variance : " << BermudanBasket_Milstein->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BermudanBasket_Milstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BermudanBasket_Milstein->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "Price antithetic : " << BermudanBasket_Milstein->ComputePrice(10000,true) << std::endl;
    std::cout << " (variance : " << BermudanBasket_Milstein->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BermudanBasket_Milstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BermudanBasket_Milstein->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "Price control variate : " << BermudanBasket_Milstein->ComputePrice_ControlVariate(10000) << std::endl;
    std::cout << " (variance : " << BermudanBasket_Milstein->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BermudanBasket_Milstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BermudanBasket_Milstein->calculate_ConfidenceInterval()[1] << std::endl;
    BermudanBasketCall* BermudanBasket_Milstein_VDC = new BermudanBasketCall(MilsteinPD_VDC, K, R, maturity, W, c, S, VCV);
    std::cout << "Price Bermudan Basket Call Milstein VDC: " << BermudanBasket_Milstein_VDC->ComputePrice_VDC(5000) << std::endl;
    std::cout << " (variance : " << BermudanBasket_Milstein_VDC->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BermudanBasket_Milstein_VDC->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BermudanBasket_Milstein_VDC->calculate_ConfidenceInterval()[1] << std::endl;

}




// Exécuter le programme : Ctrl+F5 ou menu Déboguer > Exécuter sans débogage
// Déboguer le programme : F5 ou menu Déboguer > Démarrer le débogage

// Astuces pour bien démarrer : 
//   1. Utilisez la fenêtre Explorateur de solutions pour ajouter des fichiers et les gérer.
//   2. Utilisez la fenêtre Team Explorer pour vous connecter au contrôle de code source.
//   3. Utilisez la fenêtre Sortie pour voir la sortie de la génération et d'autres messages.
//   4. Utilisez la fenêtre Liste d'erreurs pour voir les erreurs.
//   5. Accédez à Projet > Ajouter un nouvel élément pour créer des fichiers de code, ou à Projet > Ajouter un élément existant pour ajouter des fichiers de code existants au projet.
//   6. Pour rouvrir ce projet plus tard, accédez à Fichier > Ouvrir > Projet et sélectionnez le fichier .sln.

    /*for (int i = 0;i < 10;++i)
    {
        std::cout << Uniform->Generate() << std::endl;
    }

    std::cout << "Mean= " << Uniform->Mean(10000) << std::endl;*/


    /*for (int i = 0;i < 10;++i)
    {
        std::cout << i << std::endl;
        //v_euro_basket_call_euler.push_back(BasketCallEuler->ComputePrice(i * 20));
        //v_euro_basket_call_euler_anti.push_back(BasketCallEuler->ComputePrice(i * 20, true));
        v_euro_basket_call_euler_PCV.push_back(BasketCallEuler->ComputePrice_ControlVariate(i * 20));
    }
    //exportVectortoXl("Euro_BasketCallEuler_1m_simple.csv", v_euro_basket_call_euler);
    //exportVectortoXl("Euro_BasketCallEuler_1m_anti.csv", v_euro_basket_call_euler_anti);
    exportVectortoXl("Euro_BasketCallEuler_1m_pcv.csv", v_euro_basket_call_euler_PCV);
    */
    //graph European Basket call Milstein
    

    std::vector < double> v_euro_basket_call_milstein;
    std::vector < double> v_euro_basket_call_milstein_anti;
    std::vector < double> v_euro_basket_call_milstein_PCV;


    for (myLong i = 1; i <= 10000 / 20; i++)
    {
        std::cout << i << std::endl;
        v_euro_basket_call_milstein.push_back(BasketCallMilstein->ComputePrice(i * 20));
        v_euro_basket_call_milstein_anti.push_back(BasketCallMilstein->ComputePrice(i * 20, true));
        v_euro_basket_call_milstein_PCV.push_back(BasketCallMilstein->ComputePrice_ControlVariate(i * 20));
    }
    exportVectortoXl("Euro_BasketCallMilstein_1m_simple.csv", v_euro_basket_call_milstein);
    exportVectortoXl("Euro_BasketCallMilstein_1m_anti.csv", v_euro_basket_call_milstein_anti);
    exportVectortoXl("Euro_BasketCallMilstein_1m_pcv.csv", v_euro_basket_call_milstein_PCV);

    std::cout << "################## END ##################  " << std::endl;
    std::cout << "\n " << std::endl;
    
    
   
}