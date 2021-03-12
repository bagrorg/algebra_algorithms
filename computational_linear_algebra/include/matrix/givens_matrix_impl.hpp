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

    template<typename T>
    givens_matrix <T> givens_matrix<T>::transpose() {
        return givens_matrix<T>(_index_i, _index_j, -_s, _c, _n);
    }

    template<typename T>
    matrix <T> givens_matrix<T>::mult_from_right(const matrix <T> &m) {
        matrix<T> buff(m);
        for (std::size_t k = 0; k < buff.get_rows(); k++) {
            buff[k][_index_i] = m[k][_index_i] * _c + m[k][_index_j] * _s;
            buff[k][_index_j] = m[k][_index_i] * -_s + m[k][_index_j] * _c;
        }
        return buff;
    }
}

#endif