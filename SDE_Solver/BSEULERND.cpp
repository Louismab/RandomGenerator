#include "pch.h"
#include "BSEULERND.h"
#include <math.h> 


using namespace Eigen;


BSEULERND::BSEULERND()
{
}

BSEULERND::BSEULERND(RandomGenerator* _gen, std::vector<double> _s, std::vector<double> _r, Eigen::MatrixXd _vol, int dim)
	: BSND(_gen, _s, _r, _vol, dim)
{
    
    EigenSolver<MatrixXd> es(vol);

    if (vol.determinant() == 0) //if vol not definite positive
    {
        MatrixXd D = es.pseudoEigenvalueMatrix();
        MatrixXd P = es.pseudoEigenvectors();
        B = P * D;
    }

    else //vol is definite positive
    {
        B = vol.llt().matrixL(); 
    }
    
}

void BSEULERND::Simulate(double start_time, double end_time, size_t nb_steps)
{
    double dt = (end_time - start_time) / nb_steps;
    Eigen::VectorXd last = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(s.data(), s.size());
    
    for (int i = 0;i < dim;i++)
    {
        paths[i]= new SinglePath(start_time, end_time, nb_steps);
        paths[i]->AddValue(last[i]);
    }


    for (int i = 0;i < nb_steps;i++)
    {

        std::vector<double> dW = gen->GenerateVector(dim)*pow(dt, 0.5);
        Eigen::VectorXd dW_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(dW.data(), dW.size());
        
        std::vector<double> M = r * dt;
       
        Eigen::VectorXd M_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(M.data(), M.size());
        Eigen::VectorXd Z = M_M + (B*dW_M);
        //std::cout << "dW " << dW_M << std::endl;
        //std::cout << "B " << B << std::endl;
        //std::cout << "Z: " << Z << std::endl;
        Eigen::VectorXd next = last + last.cwiseProduct(Z);
        last = next;
        //std::cout << "last: " << last << std::endl;
        
        for (int j = 0;j < dim;j++)
        {
            paths[j]->AddValue(last[j]);
        }

    } 
}

void BSEULERND::Simulate_Antithetic(double start_time, double end_time, size_t nb_steps)
{
    double dt = (end_time - start_time) / nb_steps;
    Eigen::VectorXd last = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(s.data(), s.size());
    Eigen::VectorXd last_anti = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(s.data(), s.size());

    paths.resize(dim * 2);

    for (int i = 0;i < dim*2;i++)
    {
        paths[i] = new SinglePath(start_time, end_time, nb_steps);
        if (i<dim)
        {
            paths[i]->AddValue(last[i]);
        }
        else
        {
            paths[i]->AddValue(last[i-dim]);
        }
    }

    for (int i = 0;i < nb_steps;i++)
    {

        std::vector<double> dW = gen->GenerateVector(dim) * pow(dt, 0.5);
        Eigen::VectorXd dW_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(dW.data(), dW.size());
        Eigen::VectorXd dW_M_anti = -dW_M;
        std::vector<double> M = r * dt;

        Eigen::VectorXd M_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(M.data(), M.size());
        Eigen::VectorXd Z = M_M + (B * dW_M);
        Eigen::VectorXd Z_anti = M_M + (B * dW_M_anti);

        Eigen::VectorXd next = last + last.cwiseProduct(Z);
        last = next;

        Eigen::VectorXd next_anti = last_anti + last_anti.cwiseProduct(Z_anti);
        last_anti = next_anti;

        for (int j = 0;j < dim;j++)
        {
            paths[j]->AddValue(last[j]);
            paths[j + dim]->AddValue(last_anti[j]);
        }

    }
}


/*int BSEULERND::get_dim()
{
    return dim;
}*/


std::vector<double> operator*(std::vector<double> lhs, double rhs)
{
    size_t dim = lhs.size();
    std::vector<double> result(dim);
    for (size_t i = 0; i < dim; i++)
    {
        result[i] = lhs[i] * rhs;
    }
    return result;
}