#ifndef _QR_ALGO_HPP
#define _QR_ALGO_HPP

#include "../matrix/matrix.hpp"
#include <vector>

namespace lin_algebra {
	template<typename T>
	class QR_algo : public real_and_complex_restriction<T> {
	public:
		QR_algo(matrix<T> const &A);

		//Алгоритм через Гивенса
		std::pair<matrix<T>, matrix<T>> givens_algo();
		//Алгоритм через Хаусхолдера
        std::pair<matrix<T>, matrix<T>> housholder_algo();
        //Быстрый алгоритм для трехдиагональных матриц
        std::pair<std::vector<givens_matrix<T>>, matrix<T>> tridiag_algo();
	private:
	    //Матрица, для которой проводится разложение
		matrix<T> QR_matrix;
	};
}

#include "qr_algo_impl.hpp"

#endif