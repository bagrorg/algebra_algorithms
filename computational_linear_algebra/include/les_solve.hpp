#ifndef _LES_SOLVER_HPP
#define _LES_SOLVER_HPP

#include "matrix.hpp"


namespace lin_algebra {
    template<typename T>
    matrix<T> iteration_solve(matrix<T> &A, matrix<T> &b, long double eps);

}

#include "les_solve_impl.hpp"

#endif