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

    std::cout << " ----- WELCOME TO OUR MC PRICER ----- \n" << std::endl;
    std::cout << " BELARBI Sofia, BERKOVICH Adam, MABILEAU Louis \n" << std::endl;


    // __________________  INPUT PARAMETERS  __________________

    // General Parameters
    double maturity = 30. / 365.;
    //double maturity = 1.;
    double K = 100;
    double nb_simul = 10000;
    double nb_simul_vdc = 5000; // Too long otherwise

    // Single Underlying
    double s = 100.;
    double vol = 0.2;
    std::vector<double> r = { 0.000 };

    // Basket Underlying
    std::vector<double> S = { 50.0, 150.0, 100.0 };
    std::vector<double> R = { 0.0, 0.0, 0.0 };
    std::vector<double> W = { 1. / 3, 1. / 3, 1. / 3 };
    Eigen::MatrixXd VCV{ {0.04,0.036,0.032},{0.036,0.04,0.036},{0.032,0.036,0.04} };
    double nb_stocks = S.size();

    // Bermudan Execution Dates
    std::vector<double> c = { 5. / 365., 15. / 365.,30. / 365. }; // For one month maturity 
    //std::vector<double> c = { 1. / 4., 1. / 2.,3. / 4.,1. }; // For one year maturity 

    //Display Parameters in the Command window
    std::cout << "__________________  PARAMETERS: __________________ \n" << std::endl;
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

    std::cout << " ******* Basket Underlying ******* \n" << std::endl;
    std::cout << " Basket Underlying Spot at t=0 : S1, S2, S3" << "\n";
    print_vector(S);
    std::cout << "\n" << std::endl;

    std::cout << " Basket weigts by Underlying : " << "\n";
    print_vector(W);
    std::cout << "\n" << std::endl;

    std::cout << " Basket Underlying covariance matrix at t=0 : \n" << VCV << "\n" << std::endl;

    std::cout << " Basket Underlying risk free rate at t=0 : " << "\n" << std::endl;
    print_vector(R);
    std::cout << "\n" << std::endl;

    std::cout << " ******* Bermudan ******* \n" << std::endl;
    std::cout << " Bermudan execution dates (in years) : " << "\n" << std::endl;
    print_vector(c);
    std::cout << "\n" << std::endl;


    // __________________  GENERATION  __________________
    std::cout << "---  Start of Generation \n" << std::endl;
    RandomGenerator* Generator;
    EcuyerCombined* Uniform2 = new EcuyerCombined();
    VanDerCorput* Vdc = new VanDerCorput();
    NormalBoxMuller* Normal = new NormalBoxMuller(0, 1, Uniform2);
    NormalCLT* Normal_clt = new NormalCLT(0, 1, Uniform2);
    NormalInverseCdf* Normal_VDC = new NormalInverseCdf(0, 1, Vdc);
    
    BSEULER1D* Euler = new BSEULER1D(Normal, s, r[0], vol);
    Milstein1D* Milstein = new Milstein1D(Normal, s, r[0], vol);
    BSEULER1D* Euler_VDC = new BSEULER1D(Normal_VDC, s, r[0], vol);
    Milstein1D* Milstein_VDC = new Milstein1D(Normal_VDC, s, r[0], vol);

    BSEULERND* EulerND = new BSEULERND(Normal_clt, S, R, VCV, 3);
    MilsteinND* MilsteinPD = new MilsteinND(Normal_clt, S, R, VCV, 3);
    BSEULERND* EulerND_VDC = new BSEULERND(Normal_VDC, S, R, VCV, 3);
    MilsteinND* MilsteinPD_VDC = new MilsteinND(Normal_VDC, S, R, VCV, 3);

    // Others: 
    /*
    LinearCongruential* Uniform = new LinearCongruential(27, 17, 43, 100);
    VanDerCorput* Vdc_3 = new VanDerCorput(3);
    NormalInverseCdf* Normal_Inv = new NormalInverseCdf(0, 1, Uniform2);
    NormalBoxMuller_Quasi* Normal_VDC_3 = new NormalBoxMuller_Quasi(0, 1, Vdc, Vdc_3);
    */

    std::cout << "---  End of Generation \n" << std::endl;

    // __________________  PRICE COMPUTATION  __________________
    std::cout << "__________________  COMPUTATION of Prices, Variance and IC  __________________ \n" << std::endl;

    // ********************************  EUROPEAN CALL  ********************************
    std::cout << "*******  EUROPEAN CALL  ******* \n" << std::endl;

    // CLOSED FORM
    double price_cf = bs_price_call(s, K, vol, maturity, r[0]);
    std::cout << "Price European Call BS Closed Form : " << price_cf << std::endl;
    std::cout << "\n" << std::endl;
    
    // EULER
    EUCall* CallEuler = new EUCall(Euler, K, r, maturity);
    EUCall* CallEuler_VDC = new EUCall(Euler_VDC, 100, r, maturity);

    std::cout << "Price European Call with Euler \n " << std::endl;

    std::cout << "Price European Call Euler simple: " << CallEuler->ComputePrice(nb_simul) << std::endl;
    std::cout << " - Variance : " << CallEuler->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << CallEuler->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << CallEuler->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price European Call Euler Antithetic: " << CallEuler->ComputePrice(nb_simul,true) << std::endl;
    std::cout << " - Variance : " << CallEuler->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << CallEuler->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << CallEuler->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price European Call Euler Control Variate: " << CallEuler->ComputePrice_ControlVariate(nb_simul) << std::endl;
    std::cout << " - Variance : " << CallEuler->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << CallEuler->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << CallEuler->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price European Call Euler VDC: " << CallEuler_VDC->ComputePrice_VDC(nb_simul_vdc) << std::endl;
    std::cout << " - Variance : " << CallEuler_VDC->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << CallEuler_VDC->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << CallEuler_VDC->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << "\n" << std::endl;
 
    //MILSTEIN
    EUCall* CallMilstein = new EUCall(Milstein, K, r, maturity);
    EUCall* CallMilstein_VDC = new EUCall(Milstein_VDC, 100, r, maturity);

    std::cout << "Price European Call with Milstein : \n " << std::endl;
    
    std::cout << "Price European Call Miltsein simple: " << CallMilstein->ComputePrice(nb_simul) << std::endl;
    std::cout << " - Variance : " << CallMilstein->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << CallMilstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << CallMilstein->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price European Call Miltsein Antithetic : " << CallMilstein->ComputePrice(nb_simul,true) << std::endl;
    std::cout << " - Variance : " << CallMilstein->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << CallMilstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << CallMilstein->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price European Call Miltsein Control Variate: " << CallMilstein->ComputePrice_ControlVariate(nb_simul) << std::endl;
    std::cout << " - Variance : " << CallMilstein->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << CallMilstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << CallMilstein->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << std::endl;
    
    std::cout << "Price European Call Milstein VDC: " << CallMilstein_VDC->ComputePrice_VDC(nb_simul_vdc) << std::endl;
    std::cout << " - Variance : " << CallMilstein_VDC->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << CallMilstein_VDC->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << CallMilstein_VDC->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << "\n" << std::endl;

    // ********************************  BERMUDAN CALL  ********************************
    std::cout << "*******  BERMUDAN CALL  ******* \n" << std::endl;

    // EULER
    BermudanCall* Bermudan_Euler = new BermudanCall(Euler, K, r, maturity, c);
    BermudanCall* Bermudan_EulerVDC = new BermudanCall(Euler_VDC, K, r, maturity, c);

    std::cout << "Price Bermudan Call with Euler: \n" << std::endl;

    std::cout << "Price Bermudan Call Euler simple:" << Bermudan_Euler->ComputePrice(nb_simul) << std::endl;
    std::cout << " - Variance : " << Bermudan_Euler->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << Bermudan_Euler->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << Bermudan_Euler->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price Bermudan Call Euler Antithetic : " << Bermudan_Euler->ComputePrice(nb_simul, true) << std::endl;
    std::cout << " - Variance : " << Bermudan_Euler->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << Bermudan_Euler->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << Bermudan_Euler->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price Bermudan Call Euler Control Variate: " << Bermudan_Euler->ComputePrice_ControlVariate(nb_simul) << std::endl;
    std::cout << " - Variance : " << Bermudan_Euler->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << Bermudan_Euler->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << Bermudan_Euler->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price Bermudan Call Euler VDC: " << Bermudan_EulerVDC->ComputePrice_VDC(nb_simul_vdc) << std::endl;
    std::cout << " - Variance : " << Bermudan_EulerVDC->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << Bermudan_EulerVDC->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << Bermudan_EulerVDC->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << "\n" << std::endl;

    //MILSTEIN
    BermudanCall* Bermudan_Milstein = new BermudanCall(Milstein, K, r, maturity, c);
    BermudanCall* Bermudan_MilsteinVDC = new BermudanCall(Milstein_VDC, K, r, maturity, c);

    std::cout << "Price Bermudan Call with Milstein: \n" << std::endl;

    std::cout << "Price Bermudan Call Milstein simple: : " << Bermudan_Milstein->ComputePrice(nb_simul) << std::endl;
    std::cout << " - Variance : " << Bermudan_Milstein->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << Bermudan_Milstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << Bermudan_Milstein->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price Bermudan Call Milstein Antithetic : " << Bermudan_Milstein->ComputePrice(nb_simul, true) << std::endl;
    std::cout << " - Variance : " << Bermudan_Milstein->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << Bermudan_Milstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << Bermudan_Milstein->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << std::endl;
    
    std::cout << "Price Bermudan Call Milstein Control Variate: : " << Bermudan_Milstein->ComputePrice_ControlVariate(nb_simul) << std::endl;
    std::cout << " - Variance : " << Bermudan_Milstein->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << Bermudan_Milstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << Bermudan_Milstein->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price Bermudan Call Milstein VDC: " << Bermudan_MilsteinVDC->ComputePrice_VDC(nb_simul_vdc) << std::endl;
    std::cout << " - Variance : " << Bermudan_MilsteinVDC->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << Bermudan_MilsteinVDC->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << Bermudan_MilsteinVDC->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << "\n" << std::endl;

    // ********************************  EUROPEAN BASKET CALL  ********************************
    std::cout << "*******  EUROPEAN BASKET CALL  ******* \n" << std::endl;

    // CLOSED FORM
    double price_basketcall_cf = bs_price_basket_call(K, R, maturity, S, VCV);
    std::cout << "Price European Basket Call BS Closed Form :  " << price_basketcall_cf << std::endl;
    std::cout << "\n " << std::endl;

    //EULER
    EUBasketCall* BasketCallEuler = new EUBasketCall(EulerND, K, R, maturity, W, S, VCV);
    EUBasketCall* BasketCallEuler_VDC = new EUBasketCall(EulerND_VDC, K, R, maturity, W, S, VCV);

    std::cout << "Price European Basket Call with Euler: \n" << std::endl;

    std::cout << "Price European Basket Call Euler simple:" << BasketCallEuler->ComputePrice(nb_simul) << std::endl;
    std::cout << " - Variance : " << BasketCallEuler->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << BasketCallEuler->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << BasketCallEuler->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price European Basket Call Euler Antithetic : " << BasketCallEuler->ComputePrice(nb_simul, true) << std::endl;
    std::cout << " - Variance : " << BasketCallEuler->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << BasketCallEuler->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << BasketCallEuler->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price European Basket Call Euler Control Variate: : " << BasketCallEuler->ComputePrice_ControlVariate(nb_simul) << std::endl;
    std::cout << " - Variance : " << BasketCallEuler->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << BasketCallEuler->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << BasketCallEuler->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price European Basket Call Euler VDC: " << BasketCallEuler_VDC->ComputePrice_VDC(nb_simul_vdc) << std::endl;
    std::cout << " - Variance : " << BasketCallEuler_VDC->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << BasketCallEuler_VDC->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << BasketCallEuler_VDC->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << "\n" << std::endl;

    //MILSTEIN
    EUBasketCall* BasketCallMilstein = new EUBasketCall(MilsteinPD, K, R, maturity, W, S, VCV);
    EUBasketCall* BasketCallMilstein_VDC = new EUBasketCall(MilsteinPD_VDC, K, R, maturity, W, S, VCV);

    std::cout << "Price European Basket Call with Milstein: \n" << std::endl;

    std::cout << "Price European Basket Call Milstein simple:" << BasketCallMilstein->ComputePrice(nb_simul) << std::endl;
    std::cout << " - Variance : " << BasketCallMilstein->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << BasketCallMilstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << BasketCallMilstein->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price European Basket Call Milstein Antithetic : " << BasketCallMilstein->ComputePrice(nb_simul, true) << std::endl;
    std::cout << " - Variance : " << BasketCallMilstein->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << BasketCallMilstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << BasketCallMilstein->calculate_ConfidenceInterval()[1] << std::endl; 
    std::cout << "\n" << std::endl;

    std::cout << "Price European Basket Call Milstein Control Variate: : " << BasketCallMilstein->ComputePrice_ControlVariate(nb_simul) << std::endl;
    std::cout << " - Variance : " << BasketCallMilstein->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << BasketCallMilstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << BasketCallMilstein->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price European Basket Call Milstein VDC: " << BasketCallMilstein_VDC->ComputePrice_VDC(nb_simul_vdc) << std::endl;
    std::cout << " - Variance : " << BasketCallMilstein_VDC->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << BasketCallMilstein_VDC->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << BasketCallMilstein_VDC->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << "\n" << std::endl;

    // ********************************  BERMUDAN BASKET CALL   ********************************
    std::cout << "*******  BERMUDAN BASKET CALL   ******* \n" << std::endl;

    //EULER
    BermudanBasketCall* BermudanBasket_Euler = new BermudanBasketCall(EulerND, K, R, maturity, W, c, S, VCV);
    BermudanBasketCall* BermudanBasket_EulerVDC = new BermudanBasketCall(EulerND_VDC, K, R, maturity, W, c, S, VCV);

    std::cout << "Price Bermudan Basket Call with Euler: \n" << std::endl;

    std::cout << "Price Bermudan Basket Call Euler simple:" << BermudanBasket_Euler->ComputePrice(nb_simul) << std::endl;
    std::cout << " - Variance : " << BermudanBasket_Euler->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << BermudanBasket_Euler->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << BermudanBasket_Euler->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price Bermudan Basket Call Euler Antithetic : " << BermudanBasket_Euler->ComputePrice(nb_simul, true) << std::endl;
    std::cout << " - Variance : " << BermudanBasket_Euler->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << BermudanBasket_Euler->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << BermudanBasket_Euler->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price Bermudan Basket Call Euler Control Variate: : " << BermudanBasket_Euler->ComputePrice_ControlVariate(nb_simul) << std::endl;
    std::cout << " - Variance : " << BermudanBasket_Euler->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << BermudanBasket_Euler->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << BermudanBasket_Euler->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price Bermudan Basket Call Euler VDC: " << BermudanBasket_EulerVDC->ComputePrice_VDC(nb_simul_vdc) << std::endl;
    std::cout << " - Variance : " << BermudanBasket_EulerVDC->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << BermudanBasket_EulerVDC->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << BermudanBasket_EulerVDC->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << "\n" << std::endl;

    //MILSTEIN
    BermudanBasketCall* BermudanBasket_Milstein = new BermudanBasketCall(MilsteinPD, K, R, maturity, W, c, S, VCV);
    BermudanBasketCall* BermudanBasket_Milstein_VDC = new BermudanBasketCall(MilsteinPD_VDC, K, R, maturity, W, c, S, VCV);

    std::cout << "Price Bermudan Basket Call with Milstein: \n" << std::endl;

    std::cout << "Price Bermudan Basket Call Milstein simple:" << BermudanBasket_Milstein->ComputePrice(nb_simul) << std::endl;
    std::cout << " - Variance : " << BermudanBasket_Milstein->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << BermudanBasket_Milstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << BermudanBasket_Milstein->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price Bermudan Basket Call Milstein Antithetic : " << BermudanBasket_Milstein->ComputePrice(nb_simul, true) << std::endl;
    std::cout << " - Variance : " << BermudanBasket_Milstein->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << BermudanBasket_Milstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << BermudanBasket_Milstein->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price Bermudan Basket Call Milstein Control Variate: : " << BermudanBasket_Milstein->ComputePrice_ControlVariate(nb_simul) << std::endl;
    std::cout << " - Variance : " << BermudanBasket_Milstein->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << BermudanBasket_Milstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << BermudanBasket_Milstein->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price Bermudan Basket Call Milstein VDC: " << BermudanBasket_Milstein_VDC->ComputePrice_VDC(nb_simul_vdc) << std::endl;
    std::cout << " - Variance : " << BermudanBasket_Milstein_VDC->calculate_variance() << std::endl;
    std::cout << " - IC Lower Bound : " << BermudanBasket_Milstein_VDC->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << " - IC Upper Bound : " << BermudanBasket_Milstein_VDC->calculate_ConfidenceInterval()[1] << std::endl;
    std::cout << "\n" << "\n" << std::endl;

    std::cout << "---  End of Computation \n" << std::endl;


}

// Exemple of transfer from C++ to CVS:
// Graph European call Euler
    /*

    std::vector < double> v_euro_call_euler;
    std::vector < double> v_euro_call_euler_anti;
    std::vector < double> v_euro_call_euler_PCV;


    for (myLong i = 1; i <= 10000 / 10; i++)
    {
        v_euro_call_euler.push_back(CallEuler->ComputePrice(i * 10));
        v_euro_call_euler_anti.push_back(CallEuler->ComputePrice(i * 10, true));
        v_euro_call_euler_PCV.push_back(CallEuler->ComputePrice_ControlVariate(i * 10));
    }
    exportVectortoXl("Euro_CallEuler_1Y_simple.csv", v_euro_call_euler);
    exportVectortoXl("Euro_CallEuler_1Y_anti.csv", v_euro_call_euler_anti);
    exportVectortoXl("Euro_CallEuler_1Y_pcv.csv", v_euro_call_euler_PCV);

    */
    