#ifndef _QR_ALGO_IMPL_HPP
#define _QR_ALGO_IMPL_HPP

#include "matrix.hpp"
#include "qr_algo.hpp"

namespace lin_algebra {
	template<typename T>
	QR_algo<T>::QR_algo(matrix<T> &A) : QR_matrix(A) {};

	template<typename T>
	std::pair<matrix<T>, matrix<T>> QR_algo<T>::givens_algo() {
		matrix<T> Q(QR_matrix.get_rows()), R(QR_matrix);

		for(std::size_t i = 0; i < R.get_cols(); i++) {
			std::size_t j = i;
			while (j < QR_matrix.get_rows() && R[j][i] == 0) j++;
			if (j >= QR_matrix.get_rows()) continue;

			for(std::size_t k = j + 1; k < R.get_rows(); k++) {
				T tmp_val = std::pow( (std::pow(R[j][i], 2) + std::pow(R[k][i], 2)), 0.5);
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
}

#endif