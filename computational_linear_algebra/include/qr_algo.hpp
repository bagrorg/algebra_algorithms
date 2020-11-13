#ifndef _QR_ALGO_HPP
#define _QR_ALGO_HPP

#include "matrix.hpp"

namespace lin_algebra {
	template<typename T>
	class QR_algo : public arithmetic_restriction<T> {
	public:
		QR_algo(matrix<T> &A);

		std::pair<matrix<T>, matrix<T>> givens_algo();
	private:
		matrix<T> QR_matrix;
	};
}

#include "qr_algo_impl.hpp"

#endif