#ifndef _LES_SOLVER_HPP
#define _LES_SOLVER_HPP

#include "../matrix/matrix.hpp"


namespace lin_algebra {

    template<typename T>
    class LES_solver : public real_and_complex_restriction<T> {
    public:
        //Матрица и вектор b : Ax = b
        LES_solver(const matrix<T> &A, const matrix<T> &b);

        ~LES_solver() = default;

        //Решение Гауссом-Зейдель (точность, начальный вектор, дан ли начальный вектор)
        matrix<T> gauss_seidel(long double eps, matrix<T> x0, bool generate_random);

        //Метод итераций (точность, начальный вектор, дан ли начальный вектор)
        matrix<T> iteration_solve(long double eps, matrix<T> x0, bool generate_random);

    private:
        //Решение нижнетреугольной матрицы за быстро
        matrix<T> lower_matrix_les_solve(const matrix<T> &M, const matrix<T> &b);

        //Перевод A из Ax = b в x = Ax - b
        matrix<T> it_method_matrix_transform();

        //Разбиваем на верхне и нижне треугольную
        std::pair<matrix<T>, matrix<T>> gauss_seidel_matrix_transform();

        //Является ли диагональ полностью нулевой
        bool is_null_diag();

        //Проверка расходимости для метода итераций
        bool iteration_check(matrix<T> &x_new);


        int bad_count = 0;
        const int FAILURE_LIMIT = 20;
        bool gershgorin_check;
        matrix<T> LES_matrix, LES_vector;

        matrix<T> x;
    };
}

#include "les_solve_impl.hpp"

#endif