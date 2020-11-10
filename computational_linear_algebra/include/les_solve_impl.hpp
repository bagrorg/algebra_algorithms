#ifndef _LES_SOLVER_IMPL_HPP
#define _LES_SOLVER_IMPL_HPP

#include "les_solve.hpp"

namespace lin_algebra{
    template<typename T>
    matrix<T> iteration_solve(matrix<T> &A, matrix<T> &b, long double eps){
        bool flag = !A.check_gersh_circles();
        matrix<T> x(0, A.get_rows(), 1);
        int bad_counts = 0;
        
        x.random_generate(A.max_elem());
    
        while ((x - A*x - b).norm() >= eps) {
            matrix<T> x_new(A * x + b);
            if(x_new.norm() >= x.norm() + 1) bad_counts++;
            else bad_counts = 0;

            if(bad_counts == 20 && flag) return matrix<T>(0, 1, 1);
            x = x_new;
        }
        return x;
    }
}

#endif