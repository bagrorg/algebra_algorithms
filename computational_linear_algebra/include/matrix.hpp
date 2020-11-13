#ifndef _MATRIX_HPP
#define _MATRIX_HPP

#include <cstdio>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>

namespace lin_algebra{
    template <typename T, typename std::enable_if<std::is_arithmetic<T>::value>::type* = nullptr>
    class arithmetic_restriction {};

    template<typename T>
    class matrix : public arithmetic_restriction<T> {
    public:
        matrix() = default;
        matrix(std::size_t n);
        matrix(T val, std::size_t r, std::size_t c);
        matrix(const matrix &m);
        ~matrix();

        std::size_t get_rows() const;
        std::size_t get_cols() const;

        void set(std::size_t i, std::size_t j, T val);
        T get(std::size_t i, std::size_t j) const;
        void swap(matrix &m);
        
        matrix<T> operator+(matrix const &m);
        matrix<T> operator-(matrix const &m);
        matrix<T> operator*(matrix const &m);

        T* operator[](std::size_t index);
        const T* operator[](std::size_t index) const;

        matrix<T> &operator+=(matrix const &m);
        matrix<T> &operator-=(matrix const &m);
        matrix<T> &operator*=(matrix const &m);
        matrix<T> &operator=(matrix m);

        bool operator==(matrix const &m);
        bool operator!=(matrix const &m);

        T max_elem();
        T norm();
        bool check_gersh_circles();
        matrix<T> negate();
        matrix<T> transpose();

        template <typename U>
        friend std::istream &operator>>(std::istream &in, matrix<U> &m);
        
        template <typename U>
        friend std::ostream &operator<<(std::ostream &out, const matrix<U> &m);

        void random_generate(T abs_limit);
    protected:
        std::size_t _rows;
        std::size_t _cols;
        T **_data;  //кукуха поехала
    };

    template<typename T>
    class givens_matrix {
    public:
        givens_matrix(std::size_t i, std::size_t j, T s, T c, std::size_t n) : _index_i(i), _index_j(j),
                                                                                _n(n), _s(s), _c(c) {};
        ~givens_matrix() = default;

        matrix<T> operator*(matrix<T> const &m);
    private:
        std::size_t _index_i, _index_j;
        std::size_t _n;
        T _s, _c;
    };
    
}

#include "matrix_impl.hpp"
#include "givens_matrix_impl.hpp"

#endif