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

///HELLO adam stp marche

// HELLO BIS TRES ZZZ

int main()
{
    RandomGenerator* Generator;
    LinearCongruential* Uniform = new LinearCongruential(27, 17, 43, 100);
    EcuyerCombined* Uniform2 = new EcuyerCombined();
    NormalBoxMuller* Normal = new NormalBoxMuller(0, 1, Uniform2);
    Generator = new NormalCLT(0, 1, Uniform2);

    //SDE Solver
    double s = 100.;
    double vol = 0.1;
    double r = 0.0001;

    BSEULER1D* Euler = new BSEULER1D(Normal, s, r, vol);

    Milstein1D* Milstein = new Milstein1D(Normal, s, r, vol);
    Milstein->Simulate(0, 10, 10);

    std::cout << "Pricing Call" << std::endl;
    EUCall* call = new EUCall(Milstein, s, 100, r, vol, 10);
    std::cout << "Price Call : " << call->ComputePrice(1000) << std::endl;

    //N dimensions
    std::cout << "3 dimensions: " << std::endl;
    std::vector<double> S = { 100.0, 50.0, 60.0 };
    std::vector<double> R = { 0.0, 0.0, 0.0 };
    Eigen::MatrixXd Vol{ {0.05,0.003,0.002},{0.003,0.07,0.001},{0.002,0.001,0.08} };

    BSEULERND* EulerND = new BSEULERND(Generator, S, R, Vol, 3);
    EulerND->Simulate(0, 30.0 / 365.0, 30);
    std::cout << "Final Value : " << std::endl;
    for (int i = 0; i < 3; i++)
    {
        std::cout << i << ": " << EulerND->Get_Value(30.0 / 365.0, i) << std::endl;
    }

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