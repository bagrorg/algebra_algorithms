#ifndef _GIVENS_MATRIX_IMPL_HPP
#define _GIVENS_MATRIX_IMPL_HPP

#include "matrix.hpp"

namespace lin_algebra {
	template<typename T>
	matrix<T> givens_matrix<T>::operator*(matrix<T> const &m) {
		matrix<T> buff(m);
        for (std::size_t k = 0; k < buff.get_rows(); k++) {
            buff[_index_i][k] = m[_index_i][k] * _c + m[_index_j][k] * _s;
            buff[_index_j][k] = m[_index_i][k] * -_s + m[_index_j][k] * _c;
        }
        return buff;
	}    
}

#endif