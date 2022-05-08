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
#include <string>
#include <typeinfo>

// Function in Order to Transfer our vectors to csv
void exportVectortoXl(std::string file, std::vector<double> V)
{
    std::ofstream myfile;

    myfile.open(file);
    
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


    // ***********************  INPUT PARAMETERS  ***********************

    // General Parameters
    double maturity = 30. / 365.;
    //double maturity = 1.;
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
    double nb_stocks = S.size();

    // Bermudan Execution Dates
    std::vector<double> c = { 5. / 365., 15. / 365.,30. / 365. };
    //std::vector<double> c = { 1. / 4., 1. / 2.,3. / 4.,1. };

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


    // ***********************  GENERATION  ***********************
    std::cout << "---  Start of Generation \n" << std::endl;
    RandomGenerator* Generator;
    LinearCongruential* Uniform = new LinearCongruential(27, 17, 43, 100);
    EcuyerCombined* Uniform2 = new EcuyerCombined();
    //VanDerCorput* Vdc = new VanDerCorput(3);
    NormalBoxMuller* Normal = new NormalBoxMuller(0, 1, Uniform2);
    NormalCLT* Normal_clt = new NormalCLT(0, 1, Uniform2);
    //NormalInverseCdf* Normal_VDC = new NormalInverseCdf(0, 1, Vdc);

    BSEULER1D* Euler = new BSEULER1D(Normal, s, r[0], vol);
    Milstein1D* Milstein = new Milstein1D(Normal, s, r[0], vol);

    BSEULERND* EulerND = new BSEULERND(Normal_clt, S, R, VCV, nb_stocks);
    MilsteinND* MilsteinPD = new MilsteinND(Normal_clt, S, R, VCV, nb_stocks);

    std::cout << "---  End of Generation \n" << std::endl;
    
    // ***********************  PRICE COMPUTATION  ***********************

    ///////////////////// EUROPEAN BASKET CALL /////////////////////////

    std::cout << "##################  COMPUTATION of Prices, Variance and IC  ################## \n" << std::endl;


    //EULER
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

    std::cout << "\n " << std::endl;
    std::cout << "\n " << std::endl;

    EUBasketCall* BasketCallMilstein = new EUBasketCall(MilsteinPD, K, R, maturity, W, S, VCV);

    std::cout << "********   Price European Basket Call Milstein: ********  \n" << std::endl;
    std::cout << "Price regular : " << BasketCallMilstein->ComputePrice(nb_simul) << std::endl;
    std::cout << " (variance : " << BasketCallMilstein->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BasketCallMilstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BasketCallMilstein->calculate_ConfidenceInterval()[1] << std::endl;

    std::cout << "Price antithetic : " << BasketCallMilstein->ComputePrice(nb_simul, true) << std::endl;
    std::cout << " (variance : " << BasketCallMilstein->calculate_variance() << ")" << std::endl;
    std::cout << "lower bound : " << BasketCallMilstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BasketCallMilstein->calculate_ConfidenceInterval()[1] << std::endl;

    std::cout << "Price control variate : " << BasketCallMilstein->ComputePrice_ControlVariate(nb_simul) << std::endl;
    std::cout << " (variance : " << BasketCallMilstein->calculate_variance() << ")" << std::endl << std::endl;
    std::cout << "lower bound : " << BasketCallMilstein->calculate_ConfidenceInterval()[0] << std::endl;
    std::cout << "upper bound : " << BasketCallMilstein->calculate_ConfidenceInterval()[1] << std::endl;

    //graph European Basket call Euler
    /*
    std::vector < double> v_euro_basket_call_euler;
    std::vector < double> v_euro_basket_call_euler_anti;
    std::vector < double> v_euro_basket_call_euler_PCV;


    for (myLong i = 1; i <= 10000 / 20; i++)
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