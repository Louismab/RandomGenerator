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

//Test adam
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

int main()
{
    RandomGenerator* Generator;
    LinearCongruential* Uniform = new LinearCongruential(27, 17, 43, 100);
    EcuyerCombined* Uniform2 = new EcuyerCombined();
    VanDerCorput* Vdc = new VanDerCorput();
    NormalBoxMuller* Normal = new NormalBoxMuller(0, 1, Uniform2);
    NormalCLT* Normal_clt = new NormalCLT(0, 1, Uniform2);
    NormalInverseCdf* Normal_VDC = new NormalInverseCdf(0, 1, Vdc);

    //1 dimension
    double s = 100.;
    double vol = 0.2;
    double maturity = 30. / 365.;
    double K = 100;
    std::vector<double> r = { 0.000 };


    //N dimensions
    std::vector<double> S = { 50.0, 150.0, 100.0 };
    std::vector<double> R = { 0.0, 0.0, 0.0 };
    std::vector<double> W = { 1./3, 1./3, 1./3 };
    Eigen::MatrixXd VCV{ {0.04,0,0},{0,0.04,0},{0,0,0.04} };
    //std::vector<double> c = { 1. / 4., 1. / 2.,3. / 4.,1. };
    std::vector<double> c = {5./365., 15. / 365.,30. / 365. };

    BSEULER1D* Euler = new BSEULER1D(Normal, s, r[0], vol);
    Milstein1D* Milstein = new Milstein1D(Normal, s, r[0], vol);
    BSEULER1D* Euler_VDC = new BSEULER1D(Normal_VDC, s, r[0], vol);
    Milstein1D* Milstein_VDC = new Milstein1D(Normal_VDC, s, r[0], vol);

    BSEULERND* EulerND = new BSEULERND(Normal_clt, S, R, VCV, 3);
    MilsteinND* MilsteinPD = new MilsteinND(Normal_clt, S, R, VCV, 3);
    BSEULERND* EulerND_VDC = new BSEULERND(Normal_VDC, S, R, VCV, 3);
    MilsteinND* MilsteinPD_VDC = new MilsteinND(Normal_VDC, S, R, VCV, 3);

    double price_cf = bs_price_call(s, K, vol, maturity, r[0]);
    std::cout << "Price BS closed form = " << price_cf << std::endl;

    double price_basketcall_cf = bs_price_basket_call(K, R, maturity, S, VCV);
    std::cout << "Price basket call closed form = " << price_basketcall_cf << std::endl;

    std::cout << "\n " << std::endl;

    std::cout << "Price Call with Euler \n " << std::endl;
    EUCall* CallEuler = new EUCall(Euler, K, r, maturity);
    std::cout << "Price Call Euler: " << CallEuler->ComputePrice(10000) << std::endl;
    std::cout << " (variance : " << CallEuler->calculate_variance() << ")" << std::endl;
    std::cout << "Price Call Euler antithetic: " << CallEuler->ComputePrice(10000,true) << std::endl;
    std::cout << " (variance : " << CallEuler->calculate_variance() << ")" << std::endl;
    std::cout << "Price Call Euler Control Variate: " << CallEuler->ComputePrice_ControlVariate(1000) << std::endl;
    std::cout << " (variance : " << CallEuler->calculate_variance() << ")" << std::endl;

    //graph call Euler
    std::vector < double> v1;
    v1.push_back(2);
    v1.push_back(4);
    exportVectortoXl("IC_graphs2.csv", v1);
    /*for (myLong i = 0;i <= 10000/10;i ++)
    {
        std::cout << i << std::endl;    
        v.push_back(CallEuler->ComputePrice(i*10));
    }
    std::vector<double> IC = CallEuler->calculate_ConfidenceInterval();
    std::cout << "lower bound : " << IC[0];
    std::cout << "upper bound : " << IC[1];
    exportVectortoXl("IC_graphs.xlsx", v);*/

    
    /*std::cout << "Price Call with Euler Van Der Corput \n " << std::endl;
    for (int i=0;i < 10;i++)
    {
        Euler_VDC->Simulate(0, 30. / 365., 30);
        double last_value=Euler_VDC->Get_Value(30. / 365.);
        std::cout << "last value :" << last_value << std::endl;
    }*/

    /*EUCall* CallEuler_VDC = new EUCall(Euler_VDC, 100, r, 30. / 365.);
    std::cout << "Mean Normal_VDC= " << Normal_VDC->Mean(10000) << std::endl;
    std::cout << "Variance Normal_VDC= " << Normal_VDC->Variance(10000) << std::endl;
    std::cout << "Price Call Euler VDC: " << CallEuler_VDC->ComputePrice(100000) << std::endl;
    std::cout << " (variance : " << CallEuler_VDC->calculate_variance() << ")" << std::endl;
    std::cout << "Price Call Euler antithetic VDC: " << CallEuler_VDC->ComputePrice(100000, true) << std::endl;
    std::cout << " (variance : " << CallEuler_VDC->calculate_variance() << ")" << std::endl;*/



    std::cout << "\n " << std::endl;

    /*std::cout << "Price Call with Milstein : \n " << std::endl;
    Milstein->Simulate(0, maturity, 30);
    EUCall* CallMilstein = new EUCall(Milstein, K, r, maturity);
    std::cout << "Price Call Miltsein : " << CallMilstein->ComputePrice(10000) << std::endl;
    std::cout << " (variance : " << CallMilstein->calculate_variance() << ")" << std::endl;
    std::cout << "Price Call Miltsein antithetic : " << CallMilstein->ComputePrice(10000,true) << std::endl;
    std::cout << " (variance : " << CallMilstein->calculate_variance() << ")" << std::endl;
    std::cout << "Price Call Miltsein Control Variate: " << CallMilstein->ComputePrice_ControlVariate(1000) << std::endl;
    std::cout << " (variance : " << CallMilstein->calculate_variance() << ")" << std::endl;

    std::cout << "\n " << std::endl;

    std::cout << "Price Call with Milstein : Van Der Corput \n " << std::endl;
    EUCall* CallMilstein_VDC = new EUCall(Milstein_VDC, K, r, maturity);
    std::cout << "Price Call Miltsein VDC : " << CallMilstein_VDC->ComputePrice(10000) << std::endl;
    std::cout << " (variance : " << CallMilstein->calculate_variance() << ")" << std::endl;
    std::cout << "Price Call Miltsein antithetic VDC : " << CallMilstein_VDC->ComputePrice(10000, true) << std::endl;
    std::cout << " (variance : " << CallMilstein->calculate_variance() << ")" << std::endl;
    std::cout << "Price Call Miltsein Control Variate VDC : " << CallMilstein_VDC->ComputePrice_ControlVariate(1000) << std::endl;
    std::cout << " (variance : " << CallMilstein->calculate_variance() << ")" << std::endl;

    // ---------------------------------------------------------------------

    EUBasketCall* BasketCallEuler = new EUBasketCall(EulerND, K, R, maturity, W, S, VCV);
    EUBasketCall* BasketCallMilstein = new EUBasketCall(MilsteinPD, K, R, maturity, W, S, VCV);

    std::cout << "\n " << std::endl;

    std::cout << "Price Basket Call with EulerND : \n" << std::endl;
    std::cout << "Price regular : " << BasketCallEuler->ComputePrice(10000) << std::endl;
    std::cout << " (variance : " << BasketCallEuler->calculate_variance() << ")" << std::endl;
    std::cout << "Price Antithetic : " << BasketCallEuler->ComputePrice(10000, true) << std::endl;
    std::cout << " (variance : " << BasketCallEuler->calculate_variance() << ")" << std::endl;
    std::cout << "Price Control Variate : " << BasketCallEuler->ComputePrice_ControlVariate(10000) << std::endl;
    std::cout << " (variance : " << BasketCallEuler->calculate_variance() << ")" << std::endl;

    std::cout << "\n " << std::endl;

    //----------------------------------------------------------

    std::cout << "Price Basket with MilsteinND : \n" << std::endl;
    std::cout << "Price regular : " << BasketCallMilstein->ComputePrice(10000) << std::endl;
    std::cout << " (variance : " << BasketCallMilstein->calculate_variance() << ")" << std::endl;
    std::cout << "Price Antithetic : " << BasketCallMilstein->ComputePrice(10000, true) << std::endl;
    std::cout << " (variance : " << BasketCallMilstein->calculate_variance() << ")" << std::endl;
    std::cout << "Test ComputePrice_ControlVariate : " << BasketCallMilstein->ComputePrice_ControlVariate(10000) << std::endl;
    std::cout << " (variance : " << BasketCallMilstein->calculate_variance() << ")" << std::endl;

    std::cout << "\n " << std::endl;

    

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

    std::cout << "\n " << std::endl;

    BermudanCall* Bermudan_Milstein = new BermudanCall(Milstein, K, r, maturity, c);
    std::cout << "Price Bermudan Call Milstein: \n" << std::endl;
    std::cout << "Price regular : " << Bermudan_Milstein->ComputePrice(10000) << std::endl;
    std::cout << " (variance : " << Bermudan_Milstein->calculate_variance() << ")" << std::endl;
    std::cout << "Price antithetic : " << Bermudan_Milstein->ComputePrice(10000,true) << std::endl;
    std::cout << " (variance : " << Bermudan_Milstein->calculate_variance() << ")" << std::endl;
    std::cout << "Price control variate : " << Bermudan_Milstein->ComputePrice_ControlVariate(10000) << std::endl;
    std::cout << " (variance : " << Bermudan_Milstein->calculate_variance() << ")" << std::endl;

    std::cout << "\n " << std::endl;*/

    std::cout << "Bermudan Basket Call \n" << std::endl;

    BermudanBasketCall* BermudanBasket_Euler = new BermudanBasketCall(EulerND, K, R, maturity, W, c, S, VCV);
    std::cout << "Price Bermudan Basket Call Euler : \n " << std::endl;
    std::cout << "Price regular : " << BermudanBasket_Euler->ComputePrice(10000) << std::endl;
    std::cout << " (variance : " << BermudanBasket_Euler->calculate_variance() << ")" << std::endl;
    std::cout << "Price antithetic : " << BermudanBasket_Euler->ComputePrice(10000,true) << std::endl;
    std::cout << " (variance : " << BermudanBasket_Euler->calculate_variance() << ")" << std::endl;
    std::cout << "Price control variate : " << BermudanBasket_Euler->ComputePrice_ControlVariate(10000) << std::endl;
    std::cout << " (variance : " << BermudanBasket_Euler->calculate_variance() << ")" << std::endl;
    
    std::cout << "\n " << std::endl;
    
    BermudanBasketCall* BermudanBasket_Milstein = new BermudanBasketCall(MilsteinPD, K, R, maturity, W, c, S, VCV);
    std::cout << "Price Bermudan Basket Call Milstein : \n" << std::endl;
    std::cout << "Price regular : " << BermudanBasket_Milstein->ComputePrice(10000) << std::endl;
    std::cout << " (variance : " << BermudanBasket_Milstein->calculate_variance() << ")" << std::endl;
    std::cout << "Price antithetic : " << BermudanBasket_Milstein->ComputePrice(10000,true) << std::endl;
    std::cout << " (variance : " << BermudanBasket_Milstein->calculate_variance() << ")" << std::endl;
    std::cout << "Price control variate : " << BermudanBasket_Milstein->ComputePrice_ControlVariate(10000) << std::endl;
    std::cout << " (variance : " << BermudanBasket_Milstein->calculate_variance() << ")" << std::endl;
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
        std::cout << Uniform2->Generate() << std::endl;
    }

    std::cout << "Mean2= " << Uniform2->Mean(10000) << std::endl;*/

    //ExponentialID* Expo = new ExponentialID(2, Uniform2);

    ////head tail
    //std::cout << "HeadTail" << std::endl;
    //Generator = new HeadTail (Uniform2);
    //for (int i = 0;i < 10;++i)
    //{
    //    std::cout << Generator->Generate() << std::endl;
    //}

    ////bernoulli
    //std::cout << "Bernoulli" << std::endl;
    //Generator = new Bernoulli(Uniform2,0.1);
    //for (int i = 0;i < 10;++i)
    //{
    //    std::cout << Generator->Generate() << std::endl;
    //}

    ////binomial
    //std::cout << "Binomial" << std::endl;
    //Generator = new Binomial(Uniform2, 100,0.1);
    //for (int i = 0;i < 10;++i)
    //{
    //    std::cout << Generator->Generate() << std::endl;
    //}

    ////finite set
    //std::cout << "Finite Set" << std::endl;
    //std::vector<double> _values = { 5, 10};
    //std::vector<double> _probas = { 0.9, 0.1};
    //Generator = new FiniteSet(Uniform2, _values, _probas);
    //for (int i = 0;i < 10;++i)
    //{
    //    std::cout << Generator->Generate() << std::endl;
    //}

    ////Poisson Algo 1
    //std::cout << "Poisson Algo 1" << std::endl;
    //Generator = new PoissonAlgo1(2,Uniform2);
    //std::cout << "Mean= " << Generator->Mean(10000) << std::endl;
    //std::cout << "Variance= " << Generator->Variance(10000) << std::endl;

    ////Poisson Algo 2
    //std::cout << "Poisson Algo 2" << std::endl;
    //Generator = new PoissonAlgo2(Expo);
    //std::cout << "Mean= " << Generator->Mean(10000) << std::endl;
    //std::cout << "Variance= " << Generator->Variance(10000) << std::endl;


    ////Exponential ID
    //std::cout << "Exponential ID" << std::endl;
    //Generator = new ExponentialID(2, Uniform2);
    //std::cout << "Mean= " << Generator->Mean(10000) << std::endl;
    //std::cout << "Variance= " << Generator->Variance(10000) << std::endl;

    //Normal BM
    //std::cout << "Normal BM" << std::endl;

    //for (int i = 0;i < 10;++i)
    //{
    //    std::cout << Normal->Generate() << std::endl;
    //}
    //std::cout << "Mean= " << Normal->Mean(10000) << std::endl;
    //std::cout << "Variance= " << Normal->Variance(10000) << std::endl;

    ////Normal CLT
    //std::cout << "Normal CLT" << std::endl;

    //for (int i = 0;i < 10;++i)
    //{
    //    std::cout << Generator->Generate() << std::endl;
    //}
    //std::cout << "Mean= " << Generator->Mean(10000) << std::endl;
    //std::cout << "Variance= " << Generator->Variance(10000) << std::endl;

    //std::cout << "BSEULER1D" << std::endl;

    //Euler->Simulate(0, 10, 10);
    //double FinalValue= Euler->Get_Value(10);
    //std::cout << "Final Value" << FinalValue << std::endl;
    //for (double i = 0;i < 11;i++)
    //{
    //    std::cout << i << ": " << Euler->Get_Value(i) << std::endl;
    //}

    //std::cout << "Milstein1D" << std::endl;

    //std::cout << "ok3" << std::endl;
    //FinalValue = Milstein->Get_Value(10);
    //std::cout << "Final Value" << FinalValue << std::endl;
    //for (double i = 0;i < 11;i++)
    //{
    //    std::cout << i << ": " << Milstein->Get_Value(i) << std::endl;
    //}

    //pricing call


    //Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2, 2);
    //std::cout << A << std::endl;




    //double x;
    //std::cin >> x;

    /*Milstein = new Milstein1D(Normal, 100, r[0], 0.05);
    EUCall* call = new EUCall(Milstein, 100, r, 30.0 / 365.0);
    //std::cout << "ok" << std::endl;
    std::cout << "Price Call 1 : " << call->ComputePrice(1000) << std::endl;

    Milstein = new Milstein1D(Normal, 50, r[0], 0.07);
    call = new EUCall(Milstein, 50, r, 30.0 / 365.0);
    std::cout << "Price Call 2 : " << call->ComputePrice(1000) << std::endl;

    Milstein = new Milstein1D(Normal, 60, r[0], 0.08);
    call = new EUCall(Milstein, 60, r, 30.0 / 365.0);
    std::cout << "Price Call 3 : " << call->ComputePrice(1000) << std::endl;

    std::cout << "Test Van der Corput" << std::endl;
    NormalInverseCdf* Normal2 = new NormalInverseCdf(0, 1, Vdc);
    std::cout << "Mean= " << Normal2->Mean(10000) << std::endl;
    std::cout << "Variance= " << Normal2->Variance(10000) << std::endl;
    */

    /*std::cout << "test antithetic" << std::endl;
    EulerND->Simulate_Antithetic(0, 30.0 / 365.0, 30);
    std::cout << "Final Value : " << std::endl;
    for (int i = 0; i < 6; i++)
    {
        std::cout << i << ": " << EulerND->Get_Value(30.0 / 365.0, i) << std::endl;
    }

    /*std::cout << "Test EulerND \n" << std::endl;
    EulerND->Simulate(0, 30.0 / 365.0, 30);
    std::cout << "Final Value Test EulerND : " << std::endl;
    for (int i = 0; i < 3; i++)
    {
        std::cout << i << ": " << EulerND->Get_Value(30.0 / 365.0, i) << std::endl;
    }

    std::cout << "\n " << std::endl;

    std::cout << "Test MilsteinND \n" << std::endl;
    MilsteinPD->Simulate(0, 30.0 / 365.0, 30);
    std::cout << "Final Value Test MilsteinND: " << std::endl;
    for (int i = 0; i < 3; i++)
    {
        std::cout << i << ": " << MilsteinPD->Get_Value(30.0 / 365.0, i) << std::endl;
    }

    std::cout << "\n " << std::endl; */

    /*//1 dimension
        std::vector<double> S1 = { 100 };
        std::vector<double> R1 = { 0.0 };
        std::vector<double> W1 = { 1. };
        Eigen::MatrixXd VCV1{ {0.04} };

        BSEULERND* Euler1stock = new BSEULERND(Normal_clt, S1, R1, VCV1, 1);
        EUBasketCall* BasketCallEuler1stock = new EUBasketCall(Euler1stock, 100, R1, maturity, W1, S1, VCV1);

        MilsteinND* Milstein1stock = new MilsteinND(Normal_clt, S1, R1, VCV1, 1);
        EUBasketCall* BasketCallMilstein1stock = new EUBasketCall(Milstein1stock, 100, R1, maturity, W1, S1, VCV1);

        /*std::cout << "Price with Euler Basket 1 stock : \n" << std::endl;
        std::cout << "Price regular : " << BasketCallEuler1stock->ComputePrice(10000) << std::endl;
        std::cout << "Price Antithetic : " << BasketCallEuler1stock->ComputePrice(10000, true) << std::endl;
        std::cout << "Price Control Variate : " << BasketCallEuler1stock->ComputePrice_ControlVariate(1000) << std::endl;

        std::cout << "Price with Milstein Basket 1 stock :\n" << std::endl;
        std::cout << "Price regular : " << BasketCallMilstein1stock->ComputePrice(10000) << std::endl;
        //std::cout << "Price Antithetic : " << BasketCallMilstein1stock->ComputePrice(1000, true) << std::endl;
        std::cout << "Price Control Variate : " << BasketCallMilstein1stock->ComputePrice_ControlVariate(10000) << std::endl;*/
    