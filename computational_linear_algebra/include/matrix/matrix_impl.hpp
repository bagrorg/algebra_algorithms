#ifndef _MATRIX_IMPL_HPP
#define _MATRIX_IMPL_HPP

#include <iomanip>
#include "matrix.hpp"
#include <algorithm>

namespace lin_algebra {
    using std::abs, std::sqrt, std::pow;
    template<typename T>
    matrix<T>::matrix(T val, std::size_t r, std::size_t c) {
        //Выставляем значения для rows, cols и выделяем память
        _rows = r;
        _cols = c;
        _data = new T *[_rows];
        _data[0] = new T[_rows * _cols];
        for (std::size_t i = 1; i != _rows; i++) {
            _data[i] = _data[i - 1] + _cols;
        }


        //Заполняем все значением val
        for (std::size_t i = 0; i < r; i++) {
            for (std::size_t j = 0; j < c; j++) {
                _data[i][j] = val;
            }
        }
    }

    template<typename T>
    matrix<T>::matrix(const matrix &m) {
        //Выделяем память
        _data = new T *[m._rows];
        _data[0] = new T[m._rows * m._cols];
        for (std::size_t i = 1; i != m._rows; i++) {
            _data[i] = _data[i - 1] + m._cols;
        }
        //Копируем rows, cols
        _rows = m._rows;
        _cols = m._cols;

        //Копируем элементы матрицы
        for (std::size_t i = 0; i < _rows; i++) {
            for (std::size_t j = 0; j < _cols; j++) {
                _data[i][j] = m._data[i][j];
            }
        }
    }

    template<typename T>
    matrix<T>::matrix(std::size_t n) : matrix(0, n, n) {
        //Создаем матрицу с 0 везде и ставим 1 на диагонали
        for (std::size_t i = 0; i < _rows; i++) {
            _data[i][i] = 1;
        }
    }

    template<typename T>
    std::size_t matrix<T>::get_rows() const { return _rows; }

    template<typename T>
    matrix <T> matrix<T>::normilize() {
        //Подсчитываем норму матрицы/вектора
        double tmp = this->norm();
        //В случае, если норма 0 возвращаем сам вектор/матрицу
        if (static_cast<T> (tmp) == static_cast<T> (0)) return *this;

        //Делаем копию матрицы
        matrix<T> buff(*this);
        //Делим все коэффициенты матрицы на норму
        buff *= static_cast<T> (1.0 / tmp);

        return buff;
    }

    template<typename T>
    std::size_t matrix<T>::get_cols() const { return _cols; }

    template<typename T>
    void matrix<T>::set(std::size_t i, std::size_t j, T val) {
        _data[i][j] = val;
    }

    template<typename T>
    T *matrix<T>::operator[](std::size_t index) {
        if (index >= _rows) throw matrix_exception(BAD_INDEX);
        return _data[index];
    }

    template<typename T>
    const T *matrix<T>::operator[](std::size_t index) const {
        if (index >= _rows) throw matrix_exception(BAD_INDEX);
        return _data[index];
    }

    template<typename T>
    matrix<T>::~matrix() {
        delete[] _data[0];
        delete[] _data;
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
    matrix <T> &matrix<T>::operator=(matrix m) {
        swap(m);
        return *this;
    }


    template<typename T>
    bool matrix<T>::operator==(matrix const &m) {
        //Если есть неравенство по столбцам/строкам то сразу говорим нет
        if (_cols != m._cols || _rows != m._rows) {
            return false;
        }

        //Проверяем, что все элементы матриц совпадают с погрешностью
        for (std::size_t i = 0; i < _rows; i++) {
            for (std::size_t j = 0; j < _cols; j++) {
                if (abs(_data[i][j] - m._data[i][j]) > 1e-7) {
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
    matrix <T> &matrix<T>::operator+=(matrix const &m) {
        if (m._cols != _cols || m._rows != _rows) throw matrix_exception(ADD_OR_MINUS_WRONG_SIZE);
        //Прибавляем ко всем элементам матрицы элементы матрицы m
        for (std::size_t i = 0; i < _rows; i++) {
            for (std::size_t j = 0; j < _cols; j++) {
                _data[i][j] += m._data[i][j];
            }
        }
        return *this;
    }

    template<typename T>
    matrix <T> &matrix<T>::operator-=(matrix const &m) {
        if (m._cols != _cols || m._rows != _rows) throw matrix_exception(ADD_OR_MINUS_WRONG_SIZE);
        //Отнимаем ко всем элементам матрицы элементы матрицы m
        for (std::size_t i = 0; i < _rows; i++) {
            for (std::size_t j = 0; j < _cols; j++) {
                _data[i][j] -= m._data[i][j];
            }
        }
        return *this;
    }

    template<typename T>
    matrix <T> &matrix<T>::operator*=(matrix const &m) {
        if (_cols != m._rows) throw matrix_exception(MULT_WRONG_SIZE);
        return *this = *this * m;
    }

    template<typename T>
    matrix <T> &matrix<T>::operator*=(T m) {
        return *this = *this * m;
    }

    template<typename T>
    matrix <T> matrix<T>::operator+(matrix const &m) {
        if (m._cols != _cols || m._rows != _rows) throw matrix_exception(ADD_OR_MINUS_WRONG_SIZE);
        matrix<T> ret(*this);
        ret += m;
        return ret;
    }

    template<typename T>
    matrix <T> matrix<T>::operator-(matrix const &m) {
        if (m._cols != _cols || m._rows != _rows) throw matrix_exception(ADD_OR_MINUS_WRONG_SIZE);
        matrix<T> ret(*this);
        ret -= m;
        return ret;
    }

    template<typename T>
    matrix <T> matrix<T>::operator*(matrix const &m) {
        if (_cols != m._rows) throw matrix_exception(MULT_WRONG_SIZE);
        //Перемножаем матрицы
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
        //Ищем максимум
        for (std::size_t i = 0; i < _rows; i++) {
            for (std::size_t j = 0; j < _cols; j++) {
                ans = std::max(ans, _data[i][j]);
            }
        }

        return ans;
    }

    template<typename T>
    long double matrix<T>::norm() {
        long double ans = 0;
        //Прибавляем к ans квадраты модулей
        for (std::size_t i = 0; i < _rows; i++) {
            for (std::size_t j = 0; j < _cols; j++) {
                ans += abs(_data[i][j]) * abs(_data[i][j]);
            }
        }
        //Извлекаем квадрат
        return sqrt(ans);
    }

    template<typename T>
    std::ostream &operator<<(std::ostream &out, const matrix <T> &m) {
        for (std::size_t i = 0; i < m.get_rows(); i++) {
            for (std::size_t j = 0; j < m.get_cols(); j++) {
                out << m[i][j] << ' ';
            }
            out << std::endl;
        }

        return out;
    }

    template<typename T>
    std::istream &operator>>(std::istream &in, matrix <T> &m) {
        for (std::size_t i = 0; i < m.get_rows(); i++) {
            for (std::size_t j = 0; j < m.get_cols(); j++) {
                in >> m[i][j];
            }
        }

        return in;
    }

    template<typename T>
    matrix <T> matrix<T>::random_generate(double abs_limit, std::size_t n, std::size_t m) {
        //Выставляем семечко
        srand(time(0));
        //Создаем матрицу
        matrix<T> tmp(0, n, m);

        for (std::size_t i = 0; i < n; i++) {
            for (std::size_t j = 0; j < m; j++) {
                //Генерируем случайное число по модулю меньшее abs_limit
                tmp[i][j] = static_cast<T> (
                        static_cast<double> (rand()) / (static_cast<double> (RAND_MAX / (abs_limit * 2))) - abs_limit);
            }
        }

        return tmp;
    }

    template<typename T>
    matrix <T> matrix<T>::negate() {
        matrix<T> tmp(*this);
        return tmp * -1;
    }

    template<typename T>
    bool matrix<T>::gersh_in_unit_circle() {
        //Центр и радиус
        long double cent, rad;

        for (std::size_t i = 0; i < _rows; i++) {
            //Центр - диагональный элемент
            cent = _data[i][i];
            rad = 0;
            //Суммируем модули в строке
            for (std::size_t j = 0; j < _cols; j++) {
                rad += abs(_data[i][j]);
            }
            //Вычитаем модуль центра(посчитали лишним когда суммировали строку)
            rad -= abs(cent);
            //Если вышли из шарика, то возвращаем false
            if (abs(cent + rad) >= 1 || abs(cent - rad) >= 1) {
                return false;
            }
        }

        return true;
    }

    template<typename T>
    matrix <T> matrix<T>::transpose() {
        matrix<T> tmp(0, _cols, _rows);
        for (std::size_t i = 0; i < _rows; i++) {
            for (std::size_t j = 0; j < _cols; j++) {
                tmp[j][i] = _data[i][j];
            }
        }
        return tmp;
    }

    template<typename T>
    matrix <T> matrix<T>::operator*(T m) {
        matrix<T> buff(*this);
        for (std::size_t i = 0; i < _rows; i++) {
            for (std::size_t j = 0; j < _cols; j++) {
                buff[i][j] *= m;
            }
        }

        return buff;
    }

    template<typename T>
    bool matrix<T>::gersh_small_circle(double eps) {
        //Центр и радиус
        long double cent, rad;

        for (std::size_t i = 0; i < _rows; i++) {
            cent = _data[i][i];
            rad = 0;
            for (std::size_t j = 0; j < _cols; j++) {
                rad += abs(_data[i][j]);
            }
            rad -= abs(cent);
            //Если радиус большой, то возвращаем false
            if (rad >= eps) return false;
        }

        return true;
    }

    template<typename T>
    matrix <T> matrix<T>::delete_last_row_col() {
        //Создаем матрицу и заполняем всем, кроме краев
        matrix<T> ret(0, _rows - 1, _cols - 1);
        for(int i = 0; i < ret.get_rows(); i++) {
            for(int j = 0; j < ret.get_cols(); j++) {
                ret[i][j] = _data[i][j];
            }
        }

        return ret;
    }
}

#endif