#ifndef _MATRIX_IMPL_HPP
#define _MATRIX_IMPL_HPP

#include "matrix.hpp"
#include <algorithm>

namespace lin_algebra{
    template<typename T>
    matrix<T>::matrix(T val, std::size_t r, std::size_t c) {
        _rows = r;
        _cols = c;
        _data = new T *[_rows];
        _data[0] = new T[_rows * _cols];
        for (std::size_t i = 1; i != _rows; i++) {
            _data[i] = _data[i - 1] + _cols;
        }

        for(std::size_t i = 0; i < r; i++) {
            for(std::size_t j = 0; j < c; j++) {
                _data[i][j] = val;
            }
        }
    }

    template<typename T>
    matrix<T>::matrix(const matrix &m) {
        _data = new T *[m._rows];
        _data[0] = new T[m._rows * m._cols];
        for (std::size_t i = 1; i != m._rows; i++) {
            _data[i] = _data[i - 1] + m._cols;
        }
        _rows = m._rows;
        _cols = m._cols;

        for (std::size_t i = 0; i < _rows; i++) {
            for (std::size_t j = 0; j < _cols; j++) {
                _data[i][j] = m._data[i][j];
            }
        }
    }

    template<typename T>
    matrix<T>::matrix(std::size_t n) : matrix(0, n, n) {
        for(std::size_t i = 0; i < _rows; i++) {
            _data[i][i] = 1;
        }
    }

    template<typename T>
    std::size_t matrix<T>::get_rows() const { return _rows; }

    template<typename T>
    std::size_t matrix<T>::get_cols() const { return _cols; }

    template<typename T>
    void matrix<T>::set(std::size_t i, std::size_t j, T val) {
        _data[i][j] = val;
    }

    template<typename T>
    T* matrix<T>::operator[](std::size_t index) {
        return _data[index];
    }

    template<typename T>
    const T* matrix<T>::operator[](std::size_t index) const {
        return _data[index];
    }

    template<typename T>
    T matrix<T>::get(std::size_t i, std::size_t j) const {
        return _data[i][j];
    }

    template<typename T>
    void matrix<T>::swap(matrix &m) {
        std::swap(_rows, m._rows);
        std::swap(_cols, m._cols);
        std::swap(_data, m._data);
    }

    template<typename T>
    matrix<T> &matrix<T>::operator=(matrix m) {
        swap(m);
        return *this;
    }

    template<typename T>
    bool matrix<T>::operator==(matrix const &m) {
        if (_cols != m._cols || _rows != m._rows) {
            return false;
        }
        for (std::size_t i = 0; i < _cols; i++) {
            for (std::size_t j = 0; j < _rows; j++) {
                if (_data[i][j] != m._data[i][j]) {
                    return false;
                }
            }
        }
        return true;
    }

    template<typename T>
    bool matrix<T>::operator!=(matrix const &m) {
        return !(*this == m);
    }

    template<typename T>
    matrix<T> &matrix<T>::operator+=(matrix const &m) {
        for (std::size_t i = 0; i < _cols; i++) {
            for (std::size_t j = 0; j < _rows; j++) {
                _data[i][j] += m._data[i][j];
            }
        }
        return *this;
    }

    template<typename T>
    matrix<T> &matrix<T>::operator-=(matrix const &m) {
        for (std::size_t i = 0; i < _cols; i++) {
            for (std::size_t j = 0; j < _rows; j++) {
                _data[i][j] -= m._data[i][j];
            }
        }
        return *this;
    }

    template<typename T>
    matrix<T> &matrix<T>::operator*=(matrix const &m) {
        return *this = *this * m;
    }

    template<typename T>
    matrix<T> matrix<T>::operator+(matrix const &m) {
        matrix<T> ret(*this);
        ret += m;
        return ret;
    }

    template<typename T>
    matrix<T> matrix<T>::operator-(matrix const &m) {
        matrix<T> ret(*this);
        ret -= m;
        return ret;
    }

    template<typename T>
    matrix<T> matrix<T>::operator*(matrix const &m) {
        matrix<T> buff(0, _rows, m._cols);
        for (std::size_t i = 0; i < _rows; i++) {
            for (std::size_t j = 0; j < m._cols; j++) {
                for (std::size_t k = 0; k < _cols; k++) {
                    buff._data[i][j] += _data[i][k] * m._data[k][j];
                }
            }
        }
        return buff;
    }

    template<typename T>
    T matrix<T>::max_elem() {
        T ans = _data[0][0];
        for(std::size_t i = 0; i < _rows; i++) {
            for(std::size_t j = 0; j < _cols; j++) {
                ans = std::max(ans, _data[i][j]);
            }
        }

        return ans;
    }

    template<typename T>
    T matrix<T>::norm() {
        T ans = 0;
        for(std::size_t i = 0; i < _rows; i++) {
            for(std::size_t j = 0; j < _cols; j++) {
                ans += _data[i][j] * _data[i][j];
            }
        }

        return std::pow(ans, 0.5);
    }

    template<typename T>
    matrix<T>::~matrix() {
        delete [] _data[0];
        delete [] _data;
    }

    template<typename T>
    std::ostream &operator<<(std::ostream &out, const matrix<T> &m) {
        for (std::size_t i = 0; i < m.get_rows(); i++) {
            for (std::size_t j = 0; j < m.get_cols(); j++) {
                out << m[i][j] << ' ';
            }
            out << std::endl;
        }

        return out;
    }

    template<typename T>
    std::istream &operator>>(std::istream &in, matrix<T> &m) {
        for (std::size_t i = 0; i < m.get_rows(); i++) {
            for (std::size_t j = 0; j < m.get_cols(); j++) {
                in >> m[i][j];
            }
        }

        return in;
    }

    template <typename T>
    void matrix<T>::random_generate(T abs_limit) {
        srand(time(0));

        for(std::size_t i = 0; i < _rows; i++){
            for(std::size_t j = 0; j < _cols; j++) {    
                _data[i][j] = static_cast<T> (rand()) / (static_cast<T> (RAND_MAX / (2 * abs_limit))) - abs_limit;
            }
        }
    }

    template<typename T>
    matrix<T> matrix<T>::negate() {
        matrix<T> tmp(*this);
        for(std::size_t i = 0; i < _rows; i++) {
            for(std::size_t j = 0; j < _cols; j++) {
                tmp[i][j] *= -1;
            }
        }

        return tmp;
    }

    template<typename T>
    bool matrix<T>::check_gersh_circles() {
        T cent, rad;
        for(std::size_t i = 0; i < _rows; i++) {
            cent = _data[i][i];
            rad = 0;
            for(std::size_t j = 0; j < _cols; j++) {
                rad += std::abs(_data[i][j]);
            }
            rad -= cent;
            
            if(std::abs(cent + rad) >= 1 || std::abs(cent - rad) >= 1){
                return false;
            }
        }

        return true;
    }

    template<typename T>
    matrix<T> matrix<T>::transpose() {
        matrix<T> tmp(0, _cols, _rows);
        for(std::size_t i = 0; i < _rows; i++) {
            for(std::size_t j = 0; j < _cols; j++) {
                tmp[j][i] = _data[i][j];
            }
        }
        return tmp;
    }
}

#endif