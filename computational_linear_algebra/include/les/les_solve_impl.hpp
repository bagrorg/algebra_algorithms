#ifndef _LES_SOLVER_IMPL_HPP
#define _LES_SOLVER_IMPL_HPP

#include "les_solve.hpp"
#include <unistd.h>

unsigned int microseconds = 1000;



namespace lin_algebra {
    template<typename T>
    LES_solver<T>::LES_solver(const matrix <T> &A, const matrix <T> &b) : LES_matrix(A), LES_vector(b),
                                                                          x(matrix<T>(0, A.get_rows(), 1)) {}


    template<typename T>
    bool LES_solver<T>::iteration_check(matrix <T> &x_new) {
        //Если норма больше предыдущего то ++
        if (x_new.norm() >= x.norm() + 1) bad_count++;
            //Иначе обнуляем счетчик
        else bad_count = 0;
        return bad_count == FAILURE_LIMIT;
    }

    template<typename T>
    bool LES_solver<T>::is_null_diag() {
        for (std::size_t i = 0; i < std::min(LES_matrix.get_rows(), LES_matrix.get_cols()); i++) {
            if (LES_matrix[i][i] != 0) return false;
        }
        return true;
    }

    template<typename T>
    matrix <T> LES_solver<T>::iteration_solve(long double eps, matrix <T> x0, bool generate_random) {
        //Переделываем матрицу
        matrix<T> tmp = it_method_matrix_transform();
        //Проверяем, что все круги Гершгорина в единичной сфере
        gershgorin_check = tmp.gersh_in_unit_circle();
        //Если просят сгенерировать случайный, то генерируем
        if (generate_random) x = matrix<T>::random_generate(2, LES_matrix.get_rows(), 1);
            //Иначе ставим х0
        else x = x0;

        //Пока не удовлитворим условие на точность
        while ((x - tmp * x - LES_vector).norm() >= eps) {
            //Берем новый икс
            matrix<T> x_new(tmp * x + LES_vector);
            //Если круги Гершгорина не лежат в круге единичного радиуса и мы долго расходимся, то возвращаем 0
            if (!gershgorin_check && iteration_check(x_new)) return matrix<T>(0, 1, 1);
            //Если все хорошо, делаем переход
            x = x_new;
        }

        gershgorin_check = false;
        bad_count = 0;
        return x;
    }

    template<typename T>
    matrix <T> LES_solver<T>::gauss_seidel(long double eps, matrix <T> x0, bool generate_random) {
        //Если нулевая диагональ, то возвращаем 0
        if (is_null_diag()) return matrix<T>(0, 1, 1);
        //Делаем LD декомпозицию
        std::pair<matrix<T>, matrix<T>> LD_decomposition = gauss_seidel_matrix_transform();
        //Если надо - делаем рандомный
        if (generate_random) x = matrix<T>::random_generate(20, LES_matrix.get_rows(), 1);
        else x = x0;

        while ((LES_matrix * x - LES_vector).norm() >= eps) {
            //Итеративным переходом делаем новый x
            matrix<T> x_new(lower_matrix_les_solve(LD_decomposition.first, LD_decomposition.second * x + LES_vector));
            //Если мы долго расходимся по норме, то возвращаем 0
            if (iteration_check(x_new)) return matrix<T>(0, 1, 1);
            //Если все хорошо, то делаем переход
            x = x_new;
        }

        gershgorin_check = false;
        bad_count = 0;
        return x;
    }

    template<typename T>
    matrix <T> LES_solver<T>::it_method_matrix_transform() {
        //Преобразуем матрицу из A в такую чтобы x = Ax + b
        matrix<T> tmp, En(LES_matrix.get_rows());
        tmp = En - LES_matrix;
        return tmp;
    }

    template<typename T>
    std::pair<matrix < T>, matrix <T>> LES_solver<T>::gauss_seidel_matrix_transform() {
        //Раскладываем матрицу на L и D
        matrix<T> lower_triang(0, LES_matrix.get_rows(), LES_matrix.get_cols()),
                upper_triag(0, LES_matrix.get_rows(), LES_matrix.get_cols());

        for (std::size_t i = 0; i < LES_matrix.get_rows(); i++) {
            for (std::size_t j = 0; j < LES_matrix.get_cols(); j++) {
                if (i < j) upper_triag[i][j] = LES_matrix[i][j];
                else lower_triang[i][j] = LES_matrix[i][j];
            }
        }
        return {lower_triang, upper_triag.negate()};
    }

    template<typename T>
    matrix <T> LES_solver<T>::lower_matrix_les_solve(const matrix <T> &M, const matrix <T> &b) {
        //Заводим матрицу решений и решаем СЛУ за быстро сверху вниз
        matrix<T> solve(0, LES_matrix.get_rows(), 1), tmp(b);

        for (std::size_t i = 0; i < M.get_rows(); i++) {
            for (std::size_t j = 0; j <= i; j++) {
                if (j == i) solve[i][0] = tmp[i][0] / M[i][j];
                else tmp[i][0] -= solve[j][0] * M[i][j];
            }
        }

        return solve;
    }
}

#endif