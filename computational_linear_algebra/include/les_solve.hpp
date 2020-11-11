#ifndef _LES_SOLVER_HPP
#define _LES_SOLVER_HPP

#include "matrix.hpp"


namespace lin_algebra {

    template<typename T>
    class LES_solver : public arithmetic_restriction<T> {
    public:
        LES_solver(const matrix<T> &A, const matrix<T> &b);
        ~LES_solver() = default;

        matrix<T> gauss_seidel(T eps);
        matrix<T> iteration_solve(T eps);
        
    private:
        matrix<T> lower_matrix_les_solve(const matrix<T> &M, const matrix<T> &b);
        matrix<T> it_method_matrix_transform();
        std::pair<matrix<T>, matrix<T>> gauss_seidel_matrix_transform();

        bool is_null_diag();
        bool iteration_check(matrix<T> &x_new);


        int bad_count = 0, FAILURE_LIMIT = 20;
        bool gershgorin_check;
        matrix<T> LES_matrix, LES_vector;
        
        matrix<T> x;
    };
}

#include "les_solve_impl.hpp"

#endif