#ifndef _MATRIX_HPP
#define _MATRIX_HPP

#include <cstdio>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include <complex>


namespace lin_algebra{
    template<class T> struct is_complex : std::false_type {};
    template<class T> struct is_complex<std::complex<T>> : std::true_type {};

    template <typename T, typename std::enable_if<std::is_arithmetic<T>::value || is_complex<T>::value>::type* = nullptr>
    class real_and_complex_restriction {};

    template<typename T>
    class matrix : public real_and_complex_restriction<T> {
    public:
        //Cоздается пустая матрица
        matrix() : matrix(1) {};
        //Создается единичная матрица размером n x n
        matrix(std::size_t n);
        //Создается матрица размером r x c с val на всех позициях
        matrix(T val, std::size_t r, std::size_t c);
        //Копирование
        matrix(const matrix &m);
        ~matrix();

        //Позволяет поменять матрицы местами (нужно для копирования)
        void swap(matrix &m);

        //Возвращают размеры матрицы
        std::size_t get_rows() const;
        std::size_t get_cols() const;

        //Сеттеры и геттеры (по сути не нужны, так как есть оператор [][])
        void set(std::size_t i, std::size_t j, T val);
        T get(std::size_t i, std::size_t j) const;

        //+ - * над матрицами
        matrix<T> operator+(matrix const &m);
        matrix<T> operator-(matrix const &m);
        matrix<T> operator*(matrix const &m);

        //Операторы доступа к элементам матрицы
        T* operator[](std::size_t index);
        const T* operator[](std::size_t index) const;


        //Операции над матрицами
        matrix<T> &operator+=(matrix const &m);
        matrix<T> &operator-=(matrix const &m);
        matrix<T> &operator*=(matrix const &m);
        matrix<T> &operator*=(T m);
        matrix<T> operator*(T m);
        matrix<T> &operator=(matrix m);

        //Сравнение матриц
        bool operator==(matrix const &m);
        bool operator!=(matrix const &m);

        //Функция возвращает A / ||A||
        matrix<T> normilize();

        //Максимальный элемент матрицы
        T max_elem();

        //Норма (евклидова если вектор, фробениуса иначе)
        long double norm();

        //Проверка, что все круги Гершгорина лежат в единичной сфере
        bool gersh_in_unit_circle();

        //Проверка, что круги гершгорина имеют радиус миньший eps
        bool gersh_small_circle(double eps);

        //Возвращает матрицу с bij = - aij
        matrix<T> negate();

        //Возвращает транспонированную матрицу
        matrix<T> transpose();

        template <typename U>
        friend std::istream &operator>>(std::istream &in, matrix<U> &m);
        template <typename U>
        friend std::ostream &operator<<(std::ostream &out, const matrix<U> &m);

        matrix<T> delete_last_row_col();

        //Возвращает матрицу с случайными элементами, по модулю меньшими abs_limit
        static matrix<T> random_generate(double abs_limit, std::size_t n, std::size_t m);
    protected:
        std::size_t _rows;
        std::size_t _cols;
        T** _data;
    };

    template<typename T>
    class givens_matrix {
    public:
        //Принимает на вход позиции i,j, s,c и размер матрицы
        givens_matrix(std::size_t i, std::size_t j, T s, T c, std::size_t n) : _index_i(i), _index_j(j),
                                                                                _n(n), _s(s), _c(c) {};
        ~givens_matrix() = default;

        //Быстрое умножение
        matrix<T> operator*(matrix<T> const &m);
        //Транспонирование
        givens_matrix<T> transpose();
        //Домножение справа (в плюсах через оператор * только слева можно сделать, справа можно бы было через матрицу обычную реализовать, но там бы был код пострашнее)
        matrix<T> mult_from_right(matrix<T> const &m);
    private:
        std::size_t _index_i, _index_j;
        std::size_t _n;
        T _s, _c;
    };

    template<typename T>
    class householder_matrix {
    public:
        //Принимаем на вход нормированный вектор v и размер матрицы
        householder_matrix(std::size_t n, matrix<T> v) : size_(n), v_(v) {};
        ~householder_matrix() = default;

        //быстрое умножение
        matrix<T> operator*(matrix<T> const &m);
        //Аналогично домножение справа
        matrix<T> mult_from_right(matrix<T> const &m);
    private:
        std::size_t size_;
        matrix<T> v_;
    };


    //Эксепшны
    enum ERROR_TYPE {
        MULT_WRONG_SIZE = 0,
        ADD_OR_MINUS_WRONG_SIZE,
        BAD_INDEX
    };

    class matrix_exception : public std::exception {
    public:
        explicit matrix_exception(ERROR_TYPE type) : flag(type) {};

        const char *what() const noexcept override {
            switch (flag) {
                case MULT_WRONG_SIZE :
                    return "Size error while multiplication";
                case ADD_OR_MINUS_WRONG_SIZE :
                    return "Size error while plus or minus";
                case BAD_INDEX :
                    return "Index error";
                default:
                    return "Unknown error";
            }
        }

    private:
        ERROR_TYPE flag;
    };
}

#include "matrix_impl.hpp"
#include "givens_matrix_impl.hpp"
#include "householder_matrix_impl.hpp"
#endif