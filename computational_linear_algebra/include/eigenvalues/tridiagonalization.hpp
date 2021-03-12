#ifndef TRIDIAGONALIZATION_HPP
#define TRIDIAGONALIZATION_HPP

#include "../matrix/matrix.hpp"

namespace lin_algebra {
    template<typename T>
    class triagonalize : public real_and_complex_restriction<T> {
    public:
        triagonalize(matrix<T> &A) : A_(A) {};
        ~triagonalize() = default;

        //Возвращает трездиагонализацию
        std::pair<matrix<T>, matrix<T>> start() {
            //Трехдиагональная матрица и матрица базиса, в котором A диагональна
            matrix<T> ret(A_), Q(A_.get_rows());

            for(int i = 0; i < A_.get_rows() - 1; i++) {
                //Вектор-столбец и вектор стандартного базиса
                matrix<T> v(0, A_.get_rows(), 1), ei(0, A_.get_rows(), 1);
                ei[i + 1][0] = 1;

                //Проставляем значения столбца
                for(int j = i + 1; j < v.get_rows(); j++) {
                    v[j][0] = ret[j][i];
                }
                //Вычитаем e_i+1 и нормализуем. Если вектор уже хороший (почти нулевой) то тогда это то что нам нужно
                if (v.norm() < 1e-12) continue;
                v = v.normilize() - ei;
                if (v.norm() < 1e-12) continue;
                v = v.normilize();
                //Создаем матрицу хаусхолдера и домножаем все на нее
                householder_matrix<T> H(A_.get_rows(), v);

                ret = H * ret;
                ret = H.mult_from_right(ret);
                Q = H * Q;
            }

            return {ret, Q};
        }
    private:
        matrix<T> A_;
    };
}

#endif //ALGEBRA_TRIDIAGONALIZATION_HPP
