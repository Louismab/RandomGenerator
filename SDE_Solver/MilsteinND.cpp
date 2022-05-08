#include "pch.h"
#include "MilsteinND.h"
#include <math.h>
#include <algorithm>
#include <random>

using namespace Eigen;


MilsteinND::MilsteinND()
{
}

MilsteinND::MilsteinND(RandomGenerator* _gen, std::vector<double> _s, std::vector<double> _r, Eigen::MatrixXd _VCV, int dim)
    : BSND(_gen, _s, _r, _VCV, dim)
{

    EigenSolver<MatrixXd> es(VCV);

    if (VCV.determinant() == 0) //if vol not definite positive
    {
        MatrixXd D = es.pseudoEigenvalueMatrix();
        MatrixXd P = es.pseudoEigenvectors();
        B = P * D;
    }

    else //vol is definite positive
    {
        B = VCV.llt().matrixL();
    }

}

void MilsteinND::Simulate(double start_time, double end_time, size_t nb_steps)
{

    double dt = (end_time - start_time) / nb_steps;
    Eigen::VectorXd last = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(s.data(), s.size());


    for (int i = 0;i < dim;i++)
    {
        paths[i] = new SinglePath(start_time, end_time, nb_steps);
        paths[i]->AddValue(last[i]);
    }


    for (int i = 0;i < nb_steps;i++)
    {
        //next = last + last * ( (r - 0.5 * pow(vol, 2.)) * dt + vol * dW + 0.5 * pow(vol, 2) * pow(dW, 2) );   
        std::vector<double> dW = gen->GenerateVector(dim);
        Eigen::VectorXd dW_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(dW.data(), dW.size()) * pow(dt, 0.5);

        Eigen::VectorXd R_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(r.data(), r.size());

        Eigen::MatrixXd B_square = B * B;

        Eigen::VectorXd T = 0.5 * VCV.diagonal();

        Eigen::VectorXd Z = B * dW_M;

        Eigen::VectorXd X = (R_M - T) * dt + Z + 0.5*(Z.cwiseProduct(Z));

        Eigen::VectorXd next = last + last.cwiseProduct(X);
        last = next;

        for (int j = 0;j < dim;j++)
        {
            paths[j]->AddValue(last[j]);
        }

    }
}

void MilsteinND::Simulate_Antithetic(double start_time, double end_time, size_t nb_steps)
{
    //For each dimension, create a single path and the associated antithetic path

    double dt = (end_time - start_time) / nb_steps;
    Eigen::VectorXd last = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(s.data(), s.size());
    Eigen::VectorXd last_anti = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(s.data(), s.size());

    for (int i = 0;i < dim;i++)
    {
        paths[i] = new SinglePath(start_time, end_time, nb_steps);
        paths[i]->AddValue(last[i]);

        paths_antithetic[i] = new SinglePath(start_time, end_time, nb_steps);
        paths_antithetic[i]->AddValue(last_anti[i]);
    }

    for (int i = 0;i < nb_steps;i++)
    {

        std::vector<double> dW = gen->GenerateVector(dim);
        Eigen::VectorXd dW_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(dW.data(), dW.size()) * pow(dt, 0.5);
        Eigen::VectorXd R_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(r.data(), r.size());
        Eigen::VectorXd dW_M_anti = -dW_M;

        Eigen::MatrixXd B_square = B * B;

        Eigen::VectorXd T = 0.5 * VCV.diagonal();

        Eigen::VectorXd Z = B * dW_M;
        Eigen::VectorXd Z_anti = B * dW_M_anti;

        Eigen::VectorXd X = (R_M - T) * dt + Z + 0.5 * (Z.cwiseProduct(Z));
        Eigen::VectorXd X_anti = (R_M - T) * dt + Z_anti + 0.5 * (Z_anti.cwiseProduct(Z_anti));

        Eigen::VectorXd next = last + last.cwiseProduct(X);
        Eigen::VectorXd next_anti = last_anti + last_anti.cwiseProduct(X_anti);

        last = next;
        last_anti = next_anti;

        for (int j = 0;j < dim;j++)
        {
            paths[j]->AddValue(last[j]);
            paths_antithetic[j]->AddValue(last_anti[j]);
        }

    }
}

void MilsteinND::Simulate_VDC(double start_time, double end_time, size_t nb_steps, myLong sim, myLong nbSim)
{
    double dt = (end_time - start_time) / nb_steps;
    Eigen::VectorXd last = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(s.data(), s.size());
    Eigen::MatrixXd x;
    std::vector<double> y;

    if (M_VDC.rows() != nbSim)
    {
        M_VDC.resize(nbSim, dim);
        for (size_t i = 0; i < nbSim; i++)
        {
            y = gen->GenerateVectorVDC(dim, i);
            for (size_t j = 0; j < dim; j++)
            {
                M_VDC(i, j) = y[j];
            }
        }

    }

    for (int i = 0;i < dim;i++)
    {
        paths[i] = new SinglePath(start_time, end_time, nb_steps);
        paths[i]->AddValue(last[i]);
    }

    Eigen::VectorXd dW_M = M_VDC.row(sim) * pow(dt, 0.5);

    Eigen::VectorXd R_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(r.data(), r.size());

    Eigen::MatrixXd B_square = B * B;

    Eigen::VectorXd T = 0.5 * VCV.diagonal();

    Eigen::VectorXd Z = B * dW_M;

    Eigen::VectorXd X = (R_M - T) * dt + Z + 0.5 * (Z.cwiseProduct(Z));

    Eigen::VectorXd next = last + last.cwiseProduct(X);
    last = next;

    for (int j = 0;j < dim;j++)
    {
        paths[j]->AddValue(last[j]);
    }

    for (int i = 0;i < nb_steps;i++)
    {
        x = M_VDC;

        for (unsigned int j = 0;j < dim;j++)
        {
            auto rng = std::default_random_engine{ i * (j + 1) };
            std::shuffle(x.col(j).begin(), x.col(j).end(), rng);
        }
        Eigen::VectorXd dW_M = x.row(sim) * pow(dt, 0.5);

        Eigen::VectorXd R_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(r.data(), r.size());

        Eigen::MatrixXd B_square = B * B;

        Eigen::VectorXd T = 0.5 * VCV.diagonal();

        Eigen::VectorXd Z = B * dW_M;

        Eigen::VectorXd X = (R_M - T) * dt + Z + 0.5 * (Z.cwiseProduct(Z));

        Eigen::VectorXd next = last + last.cwiseProduct(X);
        last = next;

        for (int j = 0;j < dim;j++)
        {
            paths[j]->AddValue(last[j]);
        }

    }
}
