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
    Eigen::LLT<Eigen::MatrixXd> llt(vol); 
    EigenSolver<MatrixXd> es(vol);

    if (vol.determinant() == 0) //if vol not definite positive
    {
        MatrixXcd D = es.eigenvalues().asDiagonal();
        std::cout << D << std::endl;
        MatrixXcd P = es.eigenvectors();
        std::cout << P << std::endl;

        //B = P * D.sqrt();

        std::cout << P * D << std::endl;
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