#ifndef _HOUSEHOLDER_MATRIX_IMPL_HPP
#define _HOUSEHOLDER_MATRIX_IMPL_HPP

#include "matrix.hpp"

namespace lin_algebra {
	template<typename T>
	matrix<T> householder_matrix<T>::operator*(matrix<T> const &m) {
	    //Создаем временные матрицы
		matrix<T> buff(m), buff_2(0, 1, 1);

        buff_2 = v_.transpose() * buff; //v^T * M
        buff_2 = v_ * buff_2; // v * v^t * M
        buff_2 *= 2; // 2 * v * v^T * M

        return buff - buff_2; //M - 2 * v * v^T * M
	}

    template<typename T>
    matrix <T> householder_matrix<T>::mult_from_right(const matrix <T> &m) {
        //Операции аналогичные домножению слева, но с другой стороны
        matrix<T> buff(m), buff_2(0, 1, 1);

        buff_2 = (buff * 2 * v_) * (v_.transpose());

        return buff - buff_2;
    }
}

#endif