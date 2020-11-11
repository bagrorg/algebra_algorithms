#ifndef _LES_SOLVER_IMPL_HPP
#define _LES_SOLVER_IMPL_HPP

#include "les_solve.hpp"

namespace lin_algebra {
    template<typename T>
    LES_solver<T>::LES_solver(const matrix<T> &A, const matrix<T> &b) : LES_matrix(A), LES_vector(b),
                                                                                        x(matrix<T>(0, A.get_rows(), 1)) { }


    template<typename T>
    bool LES_solver<T>::iteration_check(matrix<T> &x_new) {
        if(x_new.norm() >= x.norm() + 1) bad_count++;
        else bad_count = 0;
        return bad_count == FAILURE_LIMIT;
    }

    template<typename T>
    bool LES_solver<T>::is_null_diag() {
        for(std::size_t i = 0; i < std::min(LES_matrix.get_rows(), LES_matrix.get_cols()); i++) {
            if(LES_matrix[i][i] != 0) return false;
        }
        return true;
    }

    template<typename T>
    matrix<T> LES_solver<T>::iteration_solve(T eps){
        matrix<T> tmp = it_method_matrix_transform();
        gershgorin_check = tmp.check_gersh_circles();
        x.random_generate(tmp.max_elem());

        while ((x - tmp * x - LES_vector).norm() >= eps) {
            matrix<T> x_new(tmp * x + LES_vector);
            if (!gershgorin_check && iteration_check(x_new)) return matrix<T>(0, 1, 1);
            x = x_new;
        }

        gershgorin_check = false;
        return x;
    }

    template<typename T>
    matrix<T> LES_solver<T>::gauss_seidel(T eps) {
        if (is_null_diag()) return matrix<T>(0, 1, 1);
        std::pair<matrix<T>, matrix<T>> LD_decomposition = gauss_seidel_matrix_transform();
        x.random_generate(LES_matrix.max_elem());

        while ((LES_matrix * x - LES_vector).norm() >= eps) {
            matrix<T> x_new(lower_matrix_les_solve(LD_decomposition.first, LD_decomposition.second * x + LES_vector));
            if (iteration_check(x_new)) return matrix<T>(0, 1, 1);
            x = x_new;
        }
        
        return x;
    }

    template<typename T>
    matrix<T> LES_solver<T>::it_method_matrix_transform() {
        matrix<T> tmp(LES_matrix);
        for(std::size_t i = 0; i < tmp.get_rows(); i++) {
            for(std::size_t j = 0; j < tmp.get_cols(); j++) {
                tmp[i][j] *= -1;
            }
            tmp[i][i] += 1;
        }

        return tmp;
    }

    template<typename T>
    std::pair<matrix<T>, matrix<T>> LES_solver<T>::gauss_seidel_matrix_transform() {
        matrix<T> lower_triang(0, LES_matrix.get_rows(), LES_matrix.get_cols()),
                    upper_triag(0, LES_matrix.get_rows(), LES_matrix.get_cols());
        
        for(std::size_t i = 0; i < LES_matrix.get_rows(); i++) {
            for(std::size_t j = 0; j < LES_matrix.get_cols(); j++) {
                if(i < j) upper_triag[i][j] = LES_matrix[i][j];
                else lower_triang[i][j] = LES_matrix[i][j];
            }
        }
        return {lower_triang, upper_triag.negate()};
    }

    template<typename T>
    matrix<T> LES_solver<T>::lower_matrix_les_solve(const matrix<T> &M, const matrix<T> &b) {
        matrix<T> solve(0, LES_matrix.get_rows(), 1), tmp(b);

        for(std::size_t i = 0; i < M.get_rows(); i++) {
            for(std::size_t j = 0; j <= i; j++) {
                if(j == i) solve[i][0] = tmp[i][0] / M[i][j];
                else tmp[i][0] -= solve[j][0] * M[i][j];
            }
        }

        return solve;
    }
}

#endif