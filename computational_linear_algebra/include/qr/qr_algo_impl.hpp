#ifndef _QR_ALGO_IMPL_HPP
#define _QR_ALGO_IMPL_HPP

#include "../matrix/matrix.hpp"
#include "qr_algo.hpp"

namespace lin_algebra {
    using std::pow;
	template<typename T>
	QR_algo<T>::QR_algo(matrix<T> const &A) : QR_matrix(A) {}

    template<typename T>
    std::pair<std::vector<givens_matrix < T>>, matrix <T>> QR_algo<T>::tridiag_algo() {
        //Создаем матрицу(которую будем делать верхнетреугольной)
        matrix<T> R(QR_matrix);
        //Список Матриц Гивенса
        std::vector<givens_matrix<T>> ret;

        for(std::size_t i = 0; i < R.get_cols() - 1; i++) {
            //Находим первый ненулевой элемент
            std::size_t j = i;
            std::size_t max_n = std::min(i + 2, QR_matrix.get_rows());
            while (j < max_n && R[j][i] == 0) j++;
            //Весь столбец нулевой
            if (j >= max_n) continue;

            //Операции аналогичные обычному разложению, но меньше итераций
            for(std::size_t k = j + 1; k < max_n; k++) {
                //Строим s и c и заводим матрицу гивенса
                T tmp_val = sqrt( ((R[j][i] * R[j][i]) + (R[k][i] * R[k][i])));
                givens_matrix<T> G(j, k, R[k][i] / tmp_val, R[j][i] / tmp_val, QR_matrix.get_rows());

                //Проводим домножение
                R = G * R;
                ret.push_back(G);
            }

            //Если элемент не диагональный
            if (j != i) {
                givens_matrix<T> G(i, j, 1, 0, QR_matrix.get_rows());
                R = G * R;
                ret.push_back(G);
            }
        }

        return {ret, R};
    }


	template<typename T>
	std::pair<matrix<T>, matrix<T>> QR_algo<T>::givens_algo() {
	    //Матрицы Q и R
		matrix<T> Q(QR_matrix.get_rows()), R(QR_matrix);

		for(std::size_t i = 0; i < R.get_cols(); i++) {
		    //Ищем первый ненулевой элемент
			std::size_t j = i;
			while (j < QR_matrix.get_rows() && R[j][i] == 0) j++;
			if (j >= QR_matrix.get_rows()) continue;

			//Строим матрицы гивенса, чтобы занулить столбец и домножаем на R и Q слева
			for(std::size_t k = j + 1; k < R.get_rows(); k++) {
				T tmp_val = std::sqrt(R[j][i] * R[j][i] + R[k][i] * R[k][i]);
				givens_matrix<T> G(j, k, R[k][i] / tmp_val, R[j][i] / tmp_val, QR_matrix.get_rows());

				R = G * R;
				Q = G * Q;
			}

			if (j != i) {
				givens_matrix<T> G(i, j, 1, 0, QR_matrix.get_rows());
				R = G * R;
				Q = G * Q;
			}
			
		}
		return {Q.transpose(), R};
	}


    template<typename T>
    std::pair<matrix<T>, matrix<T>> QR_algo<T>::housholder_algo() {
        //Будущие Q R
        matrix<T> Q(QR_matrix.get_rows()), R(QR_matrix);

        //Берем столбец в качестве вектора для отражения и проводим домножение на матрицы R и Q
        for(std::size_t i = 0; i < R.get_cols(); i++) {
            matrix<T> v(0, R.get_rows(), 1), ei(0, R.get_rows(), 1);
            ei[i][0] = 1;

            //Заполняем только нижнюю часть
            for(int j = i; j < v.get_rows(); j++) {
                v[j][0] = R[j][i];
            }
            //Если норма вектора маленькая, то он такой, как надо
            if (v.norm() < 1e-12) continue;
            v = v.normilize() - ei;
            if (v.norm() < 1e-12) continue;
            v = v.normilize();
            householder_matrix<T> H(R.get_rows(), v);

            R = H * R;
            Q = H * Q;
        }

        return {Q.transpose(), R};
    }
}

#endif